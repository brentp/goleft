package indexcov

import (
	"fmt"
	"log"
	"os"
	"reflect"
	"strings"
	"unsafe"

	arg "github.com/alexflint/go-arg"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/sam"
)

var cli = struct {
	Chrom string   `arg:"-c,help:optional chromosome to extract depth. default is entire genome."`
	Bam   []string `arg:"positional,required,help:bam(s) for which to estimate coverage"`
}{}

func chunksize(c bgzf.Chunk) int64 {
	return c.End.File - c.Begin.File
}

func getRefs(idx *bam.Index) []RefIndex {
	refs := reflect.ValueOf(*idx).FieldByName("idx").FieldByName("Refs")
	ptr := unsafe.Pointer(refs.Pointer())
	return (*(*[1 << 28]RefIndex)(ptr))[:refs.Len()]
}

// MaxCN is the maximum normalized value.
var MaxCN = float32(6)

// NormalizedDepth returns a list of numbers for the normalized depth of the given region.
// Values are scaled to have a mean of 1. If end is 0, the full chromosome is returned.
func NormalizedDepth(idx *bam.Index, refID int, start int, end int) []float32 {
	refs := getRefs(idx)

	last := refs[0]
	// get the last chromosome with any data.
	totalIntervals := 0
	for k := 0; k < len(refs)-1; k++ {
		totalIntervals += len(refs[k].Intervals)
		if len(refs[k].Intervals) > 0 {
			last = refs[k]
		}
	}
	// this gives the total file size.
	size := float64(last.Intervals[len(last.Intervals)-1].File + TileWidth)
	ref := refs[refID]

	meanSizePerTile := size / float64(totalIntervals)

	si, ei := start/TileWidth, end/TileWidth
	if end == 0 || ei >= len(ref.Intervals) {
		ei = len(ref.Intervals) - 1
	}
	if ei <= si {
		return nil
	}
	depths := make([]float32, 0, ei-si)
	for i, o := range ref.Intervals[si:ei] {
		depths = append(depths, float32(float64(ref.Intervals[si+i+1].File-o.File)/meanSizePerTile))
		if depths[i] > MaxCN {
			depths[i] = MaxCN
		}
	}
	return depths
}

func getRef(b *bam.Reader, chrom string) *sam.Reference {
	refs := b.Header().Refs()
	if strings.HasPrefix(chrom, "chr") {
		chrom = chrom[3:]
	}
	for _, ref := range refs {
		if chrom == ref.Name() {
			return ref
		}
		if strings.HasPrefix(ref.Name(), "chr") {
			if chrom == ref.Name()[3:] {
				return ref
			}
		}
	}
	return nil
}

func getShortName(b string) string {

	fh, err := os.Open(b)
	if err != nil {
		log.Fatal(err)
	}
	defer fh.Close()
	br, err := bam.NewReader(fh, 2)
	if err != nil {
		log.Fatal(err)
	}
	defer br.Close()
	m := make(map[string]bool)
	for _, rg := range br.Header().RGs() {
		m[rg.Get(sam.Tag([2]byte{'S', 'M'}))] = true
	}
	if len(m) > 1 {
		log.Println("warning: more than one tag for %s", b)
	}
	for sm := range m {
		return sm
	}
	vs := strings.Split(b, "/")
	v := vs[len(vs)-1]
	vs = strings.SplitN(v, ".", 1)
	return vs[len(vs)-1]
}

// Main is called from the goleft dispatcher
func Main() {

	arg.MustParse(&cli)

	rdr, err := os.Open(cli.Bam[0])
	if err != nil {
		log.Println(cli.Bam[0])
		panic(err)
	}
	brdr, err := bam.NewReader(rdr, 2)
	if err != nil {
		panic(err)
	}

	var refs []*sam.Reference

	if cli.Chrom != "" {
		refs = append(refs, getRef(brdr, cli.Chrom))
	} else {
		refs = brdr.Header().Refs()
	}
	rdr.Close()
	brdr.Close()
	if refs == nil {
		panic(fmt.Sprintf("unable to find chromosome: %s", cli.Chrom))
	}

	var idxs []*bam.Index
	names := make([]string, 0, len(cli.Bam))

	for _, b := range cli.Bam {

		rdr, err = os.Open(b + ".bai")
		if err != nil {
			panic(err)
		}

		idx, err := bam.ReadIndex(rdr)
		if err != nil {
			panic(err)
		}
		idxs = append(idxs, idx)
		names = append(names, getShortName(b))
	}

	depths := make([][]float32, len(idxs))
	fmt.Printf("#chrom\tstart\tend\t%s\n", strings.Join(names, "\t"))
	for _, ref := range refs {
		chrom := ref.Name()

		for k, idx := range idxs {
			depths[k] = NormalizedDepth(idx, ref.ID(), 0, ref.Len())
		}
		for i := 0; i < len(depths[0]); i++ {
			fmt.Printf("%s\t%d\t%d\t%s\n", chrom, i*16384, (i+1)*16384, depthsFor(depths, i))
		}
	}
}

func depthsFor(depths [][]float32, i int) string {
	s := make([]string, len(depths))
	for j := 0; j < len(depths); j++ {
		s[j] = fmt.Sprintf("%.3g", depths[j][i])
	}
	return strings.Join(s, "\t")
}
