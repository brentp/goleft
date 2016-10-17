package covmed

import (
	"fmt"
	"io"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"

	arg "github.com/alexflint/go-arg"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/brentp/xopen"
)

var cli = struct {
	N       int    `arg:"-n,help:number of reads to sample for length"`
	Bam     string `arg:"positional,required,help:bam for which to estimate coverage"`
	Regions string `arg:"positional,help:optional bed file to specify target regions"`
}{N: 100000}

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

// Main is called from the dispatcher
func Main() {

	arg.MustParse(&cli)
	log.Println(cli.Bam)

	fh, err := os.Open(cli.Bam)
	pcheck(err)

	brdr, err := bam.NewReader(fh, 2)
	pcheck(err)

	ifh, ierr := os.Open(cli.Bam + ".bai")
	if ierr != nil {
		// if .bam.bai didn't exist, check .bai
		ifh, err = os.Open(cli.Bam[:len(cli.Bam)-4] + ".bai")
	}
	pcheck(err)

	idx, err := bam.ReadIndex(ifh)
	pcheck(err)

	genomeBases := 0
	mapped := uint64(0)
	sizes := make([]int, 0, cli.N)
	for _, ref := range brdr.Header().Refs() {
		stats, ok := idx.ReferenceStats(ref.ID())
		if !ok {
			fmt.Fprintf(os.Stderr, "chromosome: %s not found in %s\n", ref.Name(), cli.Bam)
			continue
		}
		genomeBases += ref.Len()
		mapped += stats.Mapped

	}
	if cli.Regions != "" {
		genomeBases = readCoverage(cli.Regions)
	}
	// TODO: check that reads are from coverage regions.
	for len(sizes) < cli.N {
		rec, err := brdr.Read()
		if rec.Flags&(sam.Secondary|sam.Supplementary|sam.Unmapped) != 0 {
			continue
		}
		if err == io.EOF {
			break
		}
		pcheck(err)
		_, read := rec.Cigar.Lengths()
		sizes = append(sizes, read)
	}
	sort.Ints(sizes)

	readSize := float64(sizes[(len(sizes)-1)/2]) - 1
	coverage := float64(mapped) * readSize / float64(genomeBases)

	fmt.Fprintf(os.Stdout, "%.2f\n", coverage)
}
