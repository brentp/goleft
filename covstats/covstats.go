package covstats

import (
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"

	arg "github.com/alexflint/go-arg"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/brentp/goleft/samplename"
	"github.com/brentp/xopen"
)

var cli = struct {
	N       int      `arg:"-n,help:number of reads to sample for length"`
	Regions string   `arg:"-r,help:optional bed file to specify target regions"`
	Bams    []string `arg:"positional,required,help:bam for which to estimate coverage"`
}{N: 1000000}

func pcheck(e error) {
	if e != nil {
		panic(e)
	}
}

func readCoverage(path string) int {
	fh, err := xopen.Ropen(path)
	pcheck(err)
	cov := 0
	for {
		line, err := fh.ReadString('\n')

		if err == io.EOF {
			break
		}
		line = strings.TrimSuffix(line, "\n")
		toks := strings.SplitN(line, "\t", 5)
		s, err := strconv.Atoi(toks[1])
		pcheck(err)
		e, err := strconv.Atoi(toks[2])
		pcheck(err)
		cov += e - s
	}
	return cov
}

func meanStd(arr []int) (mean, std float64) {
	l := float64(len(arr))
	for _, a := range arr {
		mean += float64(a) / l
	}
	for _, a := range arr {
		std += math.Pow(float64(a)-mean, 2) / l
	}
	return mean, math.Sqrt(std)
}

// Stats hold info about a bam returned from `BamStats`
type Stats struct {
	InsertMean float64
	InsertSD   float64
	// 5th percentile of insert size
	InsertPct5 int
	// 95th percentile of insert size
	InsertPct95      int
	TemplateMean     float64
	TemplateSD       float64
	ReadLengthMean   float64
	ReadLengthMedian float64
	// ProportionBad is the proportion of reads that were Dup|QCFail
	ProportionBad            float64
	ProportionUnmapped       float64
	ProportionProperlyPaired float64
	ProportionDuplicate      float64

	MaxReadLength int
}

func (s Stats) String() string {
	return fmt.Sprintf("%.2f\t%.2f\t%d\t%d\t%.2f\t%.2f", s.InsertMean, s.InsertSD, s.InsertPct5, s.InsertPct95, s.TemplateMean, s.TemplateSD)
}

// number of reads to skip to avoid crap at start of chrom
const skipReads = 100000

// BamStats takes bam reader sample N well-behaved sites and return the coverage and insert-size info
func BamStats(br *bam.Reader, n int) Stats {
	br.Omit(bam.AllVariableLengthData)
	sizes := make([]int, 0, 2*n)
	insertSizes := make([]int, 0, n)
	templateLengths := make([]int, 0, n)
	var nBad, nUnmapped int
	for i := 0; i < skipReads; i++ {
		_, err := br.Read()
		if err == io.EOF {
			log.Fatal("covmed: not enough reads to sample for bam stats")
		}
	}
	s := Stats{}
	var k int

	for len(insertSizes) < n {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		pcheck(err)
		if rec.Flags&sam.Unmapped != 0 {
			nUnmapped++
			continue
		}
		k++
		if rec.Flags&(sam.Duplicate|sam.QCFail) != 0 {
			if rec.Flags&sam.Duplicate != 0 {
				s.ProportionDuplicate++
			}
			nBad++
			continue
		}
		if rec.Flags&sam.ProperPair != 0 {
			s.ProportionProperlyPaired++
		}
		if len(sizes) < 2*n {
			_, read := rec.Cigar.Lengths()
			sizes = append(sizes, read)
		} else {
			// single end reads have no pairs so we have to skip.
			if len(insertSizes) == 0 {
				break
			}
		}

		if rec.Pos < rec.MatePos && rec.Flags&sam.ProperPair == sam.ProperPair && len(rec.Cigar) == 1 && rec.Cigar[0].Type() == sam.CigarMatch {
			insertSizes = append(insertSizes, rec.MatePos-rec.End())
			templateLengths = append(templateLengths, rec.TempLen)
		}
	}

	sort.Ints(sizes)

	if len(sizes) > 0 {
		s.ProportionBad = float64(nBad) / float64(k)
		s.ProportionDuplicate = s.ProportionDuplicate / float64(k)
		s.ProportionProperlyPaired = s.ProportionProperlyPaired / float64(k)
		s.ProportionUnmapped = float64(nUnmapped) / float64(k)
		s.ReadLengthMedian = float64(sizes[(len(sizes)-1)/2]) - 1
		s.ReadLengthMean, _ = meanStd(sizes)
		sort.Ints(sizes)
		s.MaxReadLength = sizes[len(sizes)-1]
	}

	if len(insertSizes) > 0 {
		s.InsertMean, s.InsertSD = meanStd(insertSizes)
		s.TemplateMean, s.TemplateSD = meanStd(templateLengths)
		sort.Ints(insertSizes)
		l := float64(len(insertSizes) - 1)
		s.InsertPct5 = insertSizes[int(0.05*l+0.5)]
		s.InsertPct95 = insertSizes[int(0.95*l+0.5)]

	}
	return s
}

// Main is called from the dispatcher
func Main() {
	fmt.Fprintln(os.Stdout, "coverage\tinsert_mean\tinsert_sd\tinsert_5th\tinsert_95th\ttemplate_mean\ttemplate_sd\tpct_unmapped\tpct_bad_reads\tpct_duplicate\tpct_proper_pair\tread_length\tbam\tsample")

	arg.MustParse(&cli)
	for _, bamPath := range cli.Bams {

		fh, err := os.Open(bamPath)
		pcheck(err)

		brdr, err := bam.NewReader(fh, 2)
		pcheck(err)

		names := strings.Join(samplename.Names(brdr.Header()), ",")
		if names == "" {
			names = "<no-read-groups>"
		}

		ifh, ierr := os.Open(bamPath + ".bai")
		if ierr != nil {
			// if .bam.bai didn't exist, check .bai
			ifh, err = os.Open(bamPath[:len(bamPath)-4] + ".bai")
		}
		pcheck(err)

		idx, err := bam.ReadIndex(ifh)
		pcheck(err)

		genomeBases := 0
		mapped := uint64(0)
		sizes := BamStats(brdr, cli.N)
		var notFound []string
		for _, ref := range brdr.Header().Refs() {
			genomeBases += ref.Len()
			stats, ok := idx.ReferenceStats(ref.ID())
			if !ok {
				if !strings.Contains(ref.Name(), "random") && ref.Len() > 10000 {
					notFound = append(notFound, ref.Name())
				}
				continue
			}
			mapped += stats.Mapped
		}
		if len(notFound) > 0 {
			fmt.Fprintf(os.Stderr, "chromosomes: %s not found in %s\n", strings.Join(notFound, ","), bamPath)
		}
		if cli.Regions != "" {
			genomeBases = readCoverage(cli.Regions)
		}

		// TODO: check that reads are from coverage regions.
		coverage := (1 - sizes.ProportionBad) * float64(mapped) * sizes.ReadLengthMean / float64(genomeBases)

		fmt.Fprintf(os.Stdout, "%.2f\t%s\t%.2f\t%.1f\t%.1f\t%.1f\t%d\t%s\t%s\n", coverage, sizes.String(), 100*sizes.ProportionUnmapped,
			100*sizes.ProportionBad,
			100*sizes.ProportionDuplicate,
			100*sizes.ProportionProperlyPaired,
			sizes.MaxReadLength, bamPath, names)
	}
}
