package main

import (
	"fmt"
	"html/template"
	"io"
	"math"
	"math/rand"
	"os"
	"strconv"
	"strings"

	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/mat"

	"github.com/brentp/faidx"
	chartjs "github.com/brentp/go-chartjs"
	"github.com/brentp/go-chartjs/types"
	"github.com/brentp/goleft/dcnv/debiaser"
	"github.com/brentp/xopen"
	"go4.org/sort"
)

// Intervals is the wrapper for a slice of intervals.
type Intervals struct {
	Chrom string
	// 0-based starts
	Starts []uint32
	// 1-based ends
	GCs           []float64
	SeqComplexity []float64
	Ends          []uint32

	// shape is n-sites * n-samples
	Depths  *mat.Dense
	_depths []float64

	GCB *GcBounds

	sampleMedians []float64
	sampleScalars []float64
	Samples       []string
}

func (ivs Intervals) NSamples() int {
	_, c := ivs.Depths.Dims()
	return c
}

func mustAtoi(a string) uint32 {
	v, err := strconv.Atoi(a)
	if err != nil {
		panic(err)
	}
	return uint32(v)
}

func mustAtof(a string) float64 {
	v, err := strconv.ParseFloat(a, 64)
	if err != nil {
		panic(err)
	}
	return v
}

func imin(a, b int) int {
	if a < b {
		return a
	}
	return b
}

type GcBounds struct {
	Min float64
	Max float64
}

func (ivs *Intervals) addFromLine(l string, fa *faidx.Faidx, fp *faidx.FaPos) {
	toks := strings.Split(l, "\t")
	toks[len(toks)-1] = strings.TrimSpace(toks[len(toks)-1])
	// subtract $n bases since GC before will afffect reads here.
	s, e := mustAtoi(toks[1]), mustAtoi(toks[2])
	if s < 250 {
		s = 250
	}
	st, err := fa.Stats(toks[0], int(s-250), int(e))
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
	}

	if ivs.GCB != nil && (st.GC > ivs.GCB.Max || st.GC < ivs.GCB.Min) {
		return
	}

	ivs.Starts = append(ivs.Starts, s)
	ivs.Ends = append(ivs.Ends, e)

	for c := 3; c < len(toks); c++ {
		d := mustAtof(toks[c])
		//d /= float64(iv.End - iv.Start)
		ivs._depths = append(ivs._depths, d)
	}

	ivs.GCs = append(ivs.GCs, st.GC)
}

// SampleMedians gets the Median log2 values for each sample.
func (ivs *Intervals) SampleMedians() []float64 {
	r, _ := ivs.Depths.Dims()
	depths := make([]float64, r)

	ivs.sampleMedians = make([]float64, ivs.NSamples())
	for sampleI := 0; sampleI < ivs.NSamples(); sampleI++ {
		// sorting the extracted array is much faster.
		mat.Col(depths, sampleI, ivs.Depths)
		// lop off the lower depths (for exome).
		// and then normalized on the median above that lower bound.
		sort.Slice(depths, func(i, j int) bool { return depths[i] < depths[j] })
		var k int
		for k = 0; k < len(depths) && depths[k] == 0; k++ {
		}
		ivs.sampleMedians[sampleI] = depths[k:][int(0.65*float64(len(depths)-k))]
	}
	return ivs.sampleMedians
}

// NormalizeBySampleMedian divides each depth by the sample median.
func (ivs *Intervals) NormalizeBySampleMedian() {
	meds := ivs.SampleMedians()
	r, _ := ivs.Depths.Dims()
	mat := ivs.Depths
	for i := 0; i < r; i++ {
		row := mat.RawRowView(i)
		floats.Div(row, meds)
	}
}

func median(b []float64) float64 {
	a := make([]float64, len(b))
	copy(a, b)

	sort.Slice(a, func(i, j int) bool { return a[i] < a[j] })
	if len(a)%2 == 0 {
		am := a[len(a)/2]
		bm := a[len(a)/2+1]
		return (am + bm) / 2
	}
	return a[len(a)/2]
}

type MatFn func(*mat.Dense)

func (ivs *Intervals) ReadRegions(path string, fasta string) {
	fai, err := faidx.New(fasta)
	if err != nil {
		panic(err)
	}
	fp := &faidx.FaPos{}

	rdr, err := xopen.Ropen(path)
	if err != nil {
		panic(err)
	}
	i := 0
	for {
		line, err := rdr.ReadString('\n')
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		if i == 0 && (line[0] == '#' || strings.HasPrefix(line, "chrom")) {
			ivs.Samples = strings.Split(strings.TrimSpace(line), "\t")[3:]
			continue
		}
		if i == 0 || i == 1 {
			ivs.Chrom = string(line[:strings.Index(line, "\t")])
			fp.Chrom = ivs.Chrom
		}
		i++
		ivs.addFromLine(line, fai, fp)
	}
	ns := len(ivs.Samples)
	ivs.Depths = mat.NewDense(len(ivs._depths)/ns, ns, ivs._depths)
}

func (ivs *Intervals) Write(n int) {
	meds := ivs.SampleMedians()
	r, c := ivs.Depths.Dims()
	lmeds := make([]float64, c)
	for i, v := range meds {
		lmeds[i] = math.Log2(v)
	}
	_s := make([]string, ivs.NSamples())
	formatIV := func(depths []float64) string {
		for k := 0; k < len(_s); k++ {
			_s[k] = fmt.Sprintf("%.1f", depths[k])
		}
		return strings.Join(_s, "\t")
	}
	for i := 0; i < r; i++ {
		if i == n {
			break
		}
		fmt.Printf("%s:%d-%d\t%s\n", ivs.Chrom, ivs.Starts[i], ivs.Ends[i], formatIV(ivs.Depths.RawRowView(i)))
	}
}

func Pipeliner(mat *mat.Dense, fns ...MatFn) {
	for _, fn := range fns {
		fn(mat)
	}
}

// truncate depth values above this to cnMax
const cnMax = 2.5

type vs struct {
	xs []float64
	ys []float64
}

func (v *vs) Xs() []float64 {
	return v.xs
}

func (v *vs) Ys() []float64 {
	return v.ys
}

func (v *vs) Rs() []float64 {
	return nil
}

func (v *vs) Len() int {
	return len(v.xs)
}

// make it meet gonum/plot plotter.XYer

func (v *vs) XY(i int) (x, y float64) {
	return v.xs[i], v.ys[i]
}

func asValues(vals []float64, multiplier float64) chartjs.Values {

	// skip until we find non-zero.
	v := vs{xs: make([]float64, 0, len(vals)), ys: make([]float64, 0, len(vals))}
	seenNonZero := false
	for i, r := range vals {
		if r == 0 && !seenNonZero {
			continue
		}
		seenNonZero = true
		v.xs = append(v.xs, float64(i)*multiplier)
		if r > cnMax {
			r = cnMax
		}
		v.ys = append(v.ys, float64(r))
	}
	return &v
}

func randomColor(s int) *types.RGBA {
	rand.Seed(int64(s))
	return &types.RGBA{
		R: uint8(26 + rand.Intn(230)),
		G: uint8(26 + rand.Intn(230)),
		B: uint8(26 + rand.Intn(230)),
		A: 240}
}

func plotDepths(depths *mat.Dense, samples []string, chrom string, base string) error {
	chart := chartjs.Chart{Label: chrom}
	xa, err := chart.AddXAxis(chartjs.Axis{Type: chartjs.Linear, Position: chartjs.Bottom, ScaleLabel: &chartjs.ScaleLabel{FontSize: 16, LabelString: "position on " + chrom, Display: chartjs.True}})
	if err != nil {
		return err
	}
	ya, err := chart.AddYAxis(chartjs.Axis{Type: chartjs.Linear, Position: chartjs.Left,
		Tick:       &chartjs.Tick{Min: 0, Max: 2.5},
		ScaleLabel: &chartjs.ScaleLabel{FontSize: 16, LabelString: "scaled coverage", Display: chartjs.True}})
	if err != nil {
		return err
	}

	w := 0.4
	r, nsamples := depths.Dims()
	if nsamples > 30 {
		w = 0.3
	}
	if nsamples > 50 {
		w = 0.2
	}
	depth := make([]float64, r)
	for i := 0; i < nsamples; i++ {
		mat.Col(depth, i, depths)
		xys := asValues(depth, 16384)
		//log.Println(chrom, samples[i], len(xys.Xs()))
		c := randomColor(i)
		dataset := chartjs.Dataset{Data: xys, Label: samples[i], Fill: chartjs.False, PointRadius: 0, BorderWidth: w,
			BorderColor: c, BackgroundColor: c, SteppedLine: chartjs.True, PointHitRadius: 6}
		dataset.XAxisID = xa
		dataset.YAxisID = ya
		chart.AddDataset(dataset)
	}
	chart.Options.Responsive = chartjs.False
	chart.Options.Tooltip = &chartjs.Tooltip{Mode: "nearest"}
	wtr, err := os.Create(fmt.Sprintf("%s-depth-%s.html", base, chrom))
	if err != nil {
		return err
	}
	link := template.HTML(`<a href="index.html">back to index</a>`)
	if err := chart.SaveHTML(wtr, map[string]interface{}{"width": 850, "height": 550, "customHTML": link}); err != nil {
		return err
	}
	if err := wtr.Close(); err != nil {
		return err
	}
	return nil
}

func main() {

	bed := os.Args[1]
	fasta := os.Args[2]
	ivs := &Intervals{}

	ivs.ReadRegions(bed, fasta)

	db := debiaser.GeneralDebiaser{}
	db.Window = 9
	db.Vals = make([]float64, len(ivs.GCs))
	copy(db.Vals, ivs.GCs)
	Pipeliner(ivs.Depths, db.Sort, db.Debias, db.Unsort)

	ivs.NormalizeBySampleMedian()

	nsites, nsamples := ivs.Depths.Dims()
	dps := make([]string, nsamples)
	fdp, err := os.Create("depth.bed")
	if err != nil {
		panic(err)
	}
	plotDepths(ivs.Depths, ivs.Samples, ivs.Chrom, "dd")
	fmt.Fprintf(fdp, "#chrom\tstart\tend\t%s\n", strings.Join(ivs.Samples, "\t"))
	for i := 0; i < nsites; i++ {
		iv := ivs.Depths.RawRowView(i)
		for si := range dps {
			dps[si] = fmt.Sprintf("%.2f", iv[si])
		}
		fmt.Fprintf(fdp, "%s\t%d\t%d\t%s\n", ivs.Chrom, ivs.Starts[i], ivs.Ends[i], strings.Join(dps, "\t"))
	}
}
