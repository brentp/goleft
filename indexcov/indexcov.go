package indexcov

import (
	"fmt"
	"log"
	"os"
	"reflect"
	"sort"
	"strings"
	"unsafe"

	arg "github.com/alexflint/go-arg"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/sam"
)

var cli = struct {
	Prefix string   `arg:"-p,required,help:prefix for output files"`
	Chrom  string   `arg:"-c,help:optional chromosome to extract depth. default is entire genome."`
	Bam    []string `arg:"positional,required,help:bam(s) for which to estimate coverage"`
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

	// sizes is used to get the median.
	sizes := make([]uint32, 0, 16384)
	// get the last chromosome with any data.
	for k := 0; k < len(refs)-1; k++ {
		for i, iv := range refs[k].Intervals[1:] {
			sizes = append(sizes, uint32(iv.File-refs[k].Intervals[i].File))
		}
	}
	// this gives the total file size.
	ref := refs[refID]

	// we get the median as it's more stable than mean.
	sort.Slice(sizes, func(i, j int) bool { return sizes[i] < sizes[j] })
	medianSizePerTile := float64(sizes[len(sizes)/2])

	si, ei := start/TileWidth, end/TileWidth
	if end == 0 || ei >= len(ref.Intervals) {
		ei = len(ref.Intervals) - 1
	}
	if ei <= si {
		return nil
	}
	depths := make([]float32, 0, ei-si)
	for i, o := range ref.Intervals[si:ei] {
		depths = append(depths, float32(float64(ref.Intervals[si+i+1].File-o.File)/medianSizePerTile))
		if depths[i] > MaxCN {
			depths[i] = MaxCN
		}
	}
	return depths
}

const slots = 200

func tint(f float32) int {
	if v := int(f); v < slots {
		return v
	}
	return slots - 1
}

// DepthsAt calculates the count of items in depths that are at 100 * d
func DepthsAt(depths []float32, counts []int) {
	if len(counts) != slots {
		panic(fmt.Sprintf("indexcov: expecting counts to be length %d", slots))
	}
	for _, d := range depths {
		counts[tint(d*100+0.5)]++
	}
}

func DepthsROC(counts []int) []float32 {
	totals := make([]int, len(counts))
	totals[len(totals)-1] = counts[len(totals)-1]
	for i := len(totals) - 2; i >= 0; i-- {
		totals[i] = totals[i+1] + counts[i]
	}
	max := float32(totals[0])
	roc := make([]float32, len(counts))
	for i := 0; i < len(roc); i++ {
		roc[i] = float32(totals[i]) / max
	}
	return roc
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
	br, err := bam.NewReader(fh, 1)
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

func getWriter(prefix string) (*bgzf.Writer, error) {
	fh, err := os.Create(fmt.Sprintf("%s.indexcov.bed.gz", prefix))
	if err != nil {
		return nil, err
	}
	return bgzf.NewWriter(fh, 1), nil
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

	counts := make([][]int, len(idxs))
	depths := make([][]float32, len(idxs))

	bgz, err := getWriter(cli.Prefix)
	if err != nil {
		panic(err)
	}

	fmt.Fprintf(bgz, "#chrom\tstart\tend\t%s\n", strings.Join(names, "\t"))
	for i, ref := range refs {
		chrom := ref.Name()

		for k, idx := range idxs {
			depths[k] = NormalizedDepth(idx, ref.ID(), 0, ref.Len())
			if i == 0 {
				counts[k] = make([]int, slots)
			}
			DepthsAt(depths[k], counts[k])
		}
		for i := 0; i < len(depths[0]); i++ {
			fmt.Fprintf(bgz, "%s\t%d\t%d\t%s\n", chrom, i*16384, (i+1)*16384, depthsFor(depths, i))
		}
	}
	bgz.Close()

	writeROCs(counts, names, cli.Prefix)

}

func writeROCs(counts [][]int, names []string, prefix string) {
	fh, err := os.Create(fmt.Sprintf("%s.indexcov.roc", prefix))
	if err != nil {
		panic(err)
	}
	fmt.Fprintf(fh, "cov\t%s\n", strings.Join(names, "\t"))
	rocs := make([][]float32, len(counts))

	for i, scount := range counts {
		rocs[i] = DepthsROC(scount)
	}
	nSamples := len(names)

	vals := make([]string, nSamples)

	for i := 0; i < slots; i++ {
		for k := 0; k < nSamples; k++ {
			vals[k] = fmt.Sprintf("%.2f", rocs[k][i])
		}
		fmt.Fprintf(fh, "%.2f\t%s\n", float64(i)/100, strings.Join(vals, "\t"))
	}
}

func depthsFor(depths [][]float32, i int) string {
	s := make([]string, len(depths))
	for j := 0; j < len(depths); j++ {
		if i >= len(depths[j]) {
			s[j] = "0"
		} else {
			s[j] = fmt.Sprintf("%.3g", depths[j][i])
		}
	}
	return strings.Join(s, "\t")
}
