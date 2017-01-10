package indexcov

import (
	"fmt"
	"html/template"
	"io"
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
	chartjs "github.com/brentp/go-chartjs"
	"github.com/brentp/go-chartjs/types"
)

// Ploidy indicates the expected ploidy of the samples.
var Ploidy = 2

var cli = struct {
	Prefix    string   `arg:"-p,required,help:prefix for output files"`
	IncludeGL bool     `arg:"-e,help:plot GL chromosomes like: GL000201.1 which are not plotted by default"`
	Sex       []string `arg:"-X,help:name of the sex chromosome(s) used to infer sex; The first will be used to populate the sex column in a ped file."`
	Chrom     string   `arg:"-c,help:optional chromosome to extract depth. default is entire genome."`
	Bam       []string `arg:"positional,required,help:bam(s) for which to estimate coverage"`
}{Sex: []string{"X", "Y"}}

func getRefs(idx *bam.Index) []RefIndex {
	refs := reflect.ValueOf(*idx).FieldByName("idx").FieldByName("Refs")
	ptr := unsafe.Pointer(refs.Pointer())
	return (*(*[1 << 28]RefIndex)(ptr))[:refs.Len()]
}

// MaxCN is the maximum normalized value.
var MaxCN = float32(6)

// Index wraps a bam.Index to cache calculated values.
type Index struct {
	*bam.Index

	//mu                *sync.RWMutex
	medianSizePerTile float64
	refs              []RefIndex
}

func vOffset(o bgzf.Offset) int64 {
	return o.File<<16 | int64(o.Block)
}

// NormalizedDepth returns a list of numbers for the normalized depth of the given region.
// Values are scaled to have a mean of 1. If end is 0, the full chromosome is returned.
func (x *Index) NormalizedDepth(refID int, start int, end int) []float32 {

	if x.medianSizePerTile == 0.0 {
		x.refs = getRefs(x.Index)

		// sizes is used to get the median.
		sizes := make([]uint64, 0, 16384)
		// get the last chromosome with any data.
		for k := 0; k < len(x.refs)-1; k++ {
			if len(x.refs[k].Intervals) < 2 {
				continue
			}
			for i, iv := range x.refs[k].Intervals[1:] {
				sizes = append(sizes, uint64(vOffset(iv)-vOffset(x.refs[k].Intervals[i])))
			}
		}

		// we get the median as it's more stable than mean.
		sort.Slice(sizes, func(i, j int) bool { return sizes[i] < sizes[j] })
		x.medianSizePerTile = float64(sizes[len(sizes)/2])
		// if we have a single chunk of a chrom, then we get a lot of zeros so we address that here.
		if x.medianSizePerTile == 0 {
			i := len(sizes) / 2
			for ; i < len(sizes) && sizes[i] == 0; i++ {
			}
			sizes = sizes[i:]
			x.medianSizePerTile = float64(sizes[len(sizes)/2])
		}
	}
	ref := x.refs[refID]

	si, ei := start/TileWidth, end/TileWidth
	if end == 0 || ei >= len(ref.Intervals) {
		ei = len(ref.Intervals) - 1
	}
	if ei <= si {
		return nil
	}
	depths := make([]float32, 0, ei-si)
	for i, o := range ref.Intervals[si:ei] {
		depths = append(depths, float32(float64(vOffset(ref.Intervals[si+i+1])-vOffset(o))/x.medianSizePerTile))
		if depths[i] > MaxCN {
			depths[i] = MaxCN
		}
	}
	return depths
}

const slots = 70

// with 0.5, we'll get centered at 1 and max of 2.
// so the max is 1/slotsMid
const slotsMid = float64(2) / float64(3)

func tint(f float32) int {
	if v := int(f); v < slots {
		return v
	}
	return slots - 1
}

// CountsAtDepth calculates the count of items in depths that are at 100 * d
func CountsAtDepth(depths []float32, counts []int) {
	if len(counts) != slots {
		panic(fmt.Sprintf("indexcov: expecting counts to be length %d", slots))
	}
	for _, d := range depths {
		counts[tint(d*(slots*float32(slotsMid))+0.5)]++
	}
}

// CountsROC returns a slice that indicates the cumulative proportion of
// 16KB chunks that were at least (normalized) depth given by their index.
func CountsROC(counts []int) []float32 {
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
	fh, err := os.Create(fmt.Sprintf("%s-indexcov.bed.gz", prefix))
	if err != nil {
		return nil, err
	}
	return bgzf.NewWriter(fh, 1), nil
}

func zero(ints []int) {
	for i := range ints {
		ints[i] = 0
	}
}

// Main is called from the goleft dispatcher
func Main() {

	chartjs.XFloatFormat = "%.0f"
	arg.MustParse(&cli)
	if strings.HasSuffix(cli.Prefix, "/") {
		cli.Prefix = cli.Prefix + "qc"
	}

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

	var idxs []*Index
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
		idxs = append(idxs, &Index{Index: idx})
		names = append(names, getShortName(b))
	}

	counts := make([][]int, len(idxs))
	depths := make([][]float32, len(idxs))

	bgz, err := getWriter(cli.Prefix)
	if err != nil {
		panic(err)
	}

	rfh, err := os.Create(fmt.Sprintf("%s-indexcov.roc", cli.Prefix))
	if err != nil {
		panic(err)
	}

	// we plot all the coverage roc charts in a single html file.
	charts := make([]chartjs.Chart, 0, len(refs))
	sexes := make(map[string][]float64)

	fmt.Fprintf(bgz, "#chrom\tstart\tend\t%s\n", strings.Join(names, "\t"))
	for ir, ref := range refs {
		chrom := ref.Name()
		// Some samples may not have all the data, so we always take the longest sample for printing.
		longest, longesti := 0, 0

		for k, idx := range idxs {
			depths[k] = idx.NormalizedDepth(ref.ID(), 0, ref.Len())
			if len(depths[k]) > longest {
				longesti = k
			}
			if ir == 0 {
				counts[k] = make([]int, slots)
			} else {
				zero(counts[k])
			}

			CountsAtDepth(depths[k], counts[k])
		}

		for _, x := range cli.Sex {
			if x == chrom && len(depths[longesti]) > 0 {
				sexes[chrom] = GetCN(depths)
			}
		}

		for i := 0; i < len(depths[longesti]); i++ {
			fmt.Fprintf(bgz, "%s\t%d\t%d\t%s\n", chrom, i*16384, (i+1)*16384, depthsFor(depths, i))
		}
		if len(depths[longesti]) > 0 {
			c := writeROCs(counts, names, chrom, cli.Prefix, rfh)
			if cli.IncludeGL || !strings.HasPrefix(chrom, "GL") {
				charts = append(charts, c)
				if err := plotDepths(depths, names, chrom, cli.Prefix); err != nil {
					panic(err)
				}
			}
			if len(charts) > 1 {
				charts[len(charts)-1].Options.Legend = &chartjs.Legend{Display: types.False}
			}
		}
	}
	chartjs.XFloatFormat = "%.2f"
	saveCharts(fmt.Sprintf("%s-indexcov-roc.html", cli.Prefix), "", charts...)
	bgz.Close()
	writePed(sexes, cli.Sex, names, cli.Prefix)
}

func writePed(sexes map[string][]float64, keys []string, samples []string, prefix string) {
	if len(sexes) == 0 {
		return
	}
	for _, k := range keys {
		if _, ok := sexes[k]; !ok {
			fmt.Printf("chromosome %s not found. not writing ped\n")
			os.Exit(1)
		}
	}
	sexes["_inferred"] = make([]float64, len(sexes[keys[0]]))
	f, err := os.Create(fmt.Sprintf("%s-indexcov.ped", prefix))
	if err != nil {
		panic(err)
	}
	hdr := make([]string, len(keys))
	for i, k := range keys {
		hdr[i] = "CN" + k
	}
	fmt.Fprintf(f, "#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\t%s\n", strings.Join(hdr, "\t"))

	tmpl := "unknown\t%s\t-9\t-9\t%d\t-9\t"
	for i, s := range samples {
		inferred := int(0.5 + sexes[keys[0]][i])
		fmt.Fprintf(f, tmpl, s, inferred)
		sexes["_inferred"][i] = float64(inferred)
		s := make([]string, 0, len(keys))
		for _, k := range keys {
			s = append(s, fmt.Sprintf("%.2f", sexes[k][i]))
		}
		fmt.Fprintln(f, strings.Join(s, "\t"))
	}
	if len(keys) > 1 {
		chart, customjs, err := plotSex(sexes, keys[:2], samples)
		if err != nil {
			panic(err)
		}
		saveCharts(fmt.Sprintf("%s-indexcov-sex.html", cli.Prefix), customjs, chart)
	}
}

// GetCN returns an float per sample estimating the number of copies of that chromosome.
// It is a very crude estimate, but that's what indexcov is and it tends to work well.
func GetCN(depths [][]float32) []float64 {
	if depths == nil {
		return nil
	}
	meds := make([]float64, 0, len(depths))
	for _, d := range depths {
		tmp := make([]float32, 0, len(d))
		for _, dp := range d {
			// exclude sites that are exactly 0 as these are the centromere.
			if dp != 0 {
				tmp = append(tmp, dp)
			}
		}
		sort.Slice(tmp, func(i, j int) bool { return tmp[i] < tmp[j] })
		med := float64(float32(Ploidy) * tmp[int(float64(len(tmp))*0.5)])
		meds = append(meds, med)
	}
	return meds
}

func saveCharts(path string, customjs string, charts ...chartjs.Chart) {
	wtr, err := os.Create(path)
	if err != nil {
		panic(err)
	}
	defer wtr.Close()
	if err := chartjs.SaveCharts(wtr, map[string]interface{}{"height": 800, "width": 800, "custom": template.JS(customjs)}, charts...); err != nil {
		panic(err)
	}
}

func getROCs(counts [][]int) [][]float32 {
	rocs := make([][]float32, len(counts))

	for i, scount := range counts {
		rocs[i] = CountsROC(scount)
	}
	return rocs

}

func writeROCs(counts [][]int, names []string, chrom string, prefix string, fh io.Writer) chartjs.Chart {
	rocs := getROCs(counts)
	chart, err := plotROCs(rocs, names, chrom)
	if err != nil {
		panic(err)
	}
	fmt.Fprintf(fh, "#chrom\tcov\t%s\n", strings.Join(names, "\t"))
	nSamples := len(names)

	vals := make([]string, nSamples)

	for i := 0; i < slots; i++ {
		for k := 0; k < nSamples; k++ {
			vals[k] = fmt.Sprintf("%.2f", rocs[k][i])
		}
		fmt.Fprintf(fh, "%s\t%.2f\t%s\n", chrom, float64(i)/(slots*slotsMid), strings.Join(vals, "\t"))
	}
	return chart
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
