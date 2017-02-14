package main

import (
	"fmt"
	"image/color"
	"io"
	"log"
	"math"
	"os"
	"strconv"
	"strings"

	"github.com/brentp/faidx"
	"github.com/brentp/goleft/dcnv/debiaser"
	"github.com/brentp/goleft/dcnv/scalers"
	"github.com/brentp/goleft/emdepth"
	"github.com/brentp/xopen"
	"github.com/gonum/matrix/mat64"
	"github.com/gonum/plot"
	"github.com/gonum/plot/plotter"
	"github.com/gonum/plot/vg"
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
	Depths  *mat64.Dense
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
	st, err := fa.Stats(toks[0], int(s-250), int(e+250))
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
	}

	if ivs.GCB != nil && st.GC > ivs.GCB.Max || st.GC < ivs.GCB.Min {
		return
	}

	ivs.Starts = append(ivs.Starts, s)
	ivs.Ends = append(ivs.Ends, e)

	/*
		last := len(ivs.Starts) - 1
		fp.Start = int(ivs.Starts[last] - 100)
		fp.End = int(ivs.Ends[last] + 100)
		fa.Q(fp)
		//ivs.SeqComplexity = append(ivs.SeqComplexity, 1-float64(fp.Duplicity()))
	*/

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
		mat64.Col(depths, sampleI, ivs.Depths)
		// lop off the lower depths (for exome).
		// and then normalized on the median above that lower bound.
		sort.Slice(depths, func(i, j int) bool { return depths[i] < depths[j] })
		var k int
		for k = 0; k < len(depths) && depths[k] == 0; k++ {
		}
		ivs.sampleMedians[sampleI] = depths[k:][int(0.65*float64(len(depths)-k))]
	}
	med := median(ivs.sampleMedians)
	ivs.sampleScalars = make([]float64, 0, len(ivs.sampleMedians))
	for _, sm := range ivs.sampleMedians {
		ivs.sampleScalars = append(ivs.sampleScalars, med/sm)
	}
	return ivs.sampleMedians
}

func (ivs *Intervals) SampleScalars() []float64 {
	ivs.SampleMedians()
	return ivs.sampleScalars
}

// CorrectBySampleMedian subtracts the sample median from each sample.
func (ivs *Intervals) CorrectBySampleMedian() {
	scalars := ivs.SampleMedians()
	r, c := ivs.Depths.Dims()
	mat := ivs.Depths
	for i := 0; i < r; i++ {
		row := mat.RawRowView(i)
		for j := 0; j < c; j++ {
			row[j] /= scalars[j]
		}
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

type MatFn func(*mat64.Dense)

func Pipeliner(mat *mat64.Dense, fns ...MatFn) {
	for _, fn := range fns {
		fn(mat)
	}
}

// A plotter plots the before and after a MatFn is Applied
type Plotter struct {
	Idxs []int
}

// Wrap will plot a line of Values from a mat64.Dense before and after calling fn.
func (p *Plotter) Wrap(fn MatFn, Xs []float64, xlabel, ylabel, path string) MatFn {
	if p.Idxs == nil {
		p.Idxs = append(p.Idxs, 0)
	}
	pl, err := plot.New()
	if err != nil {
		panic(err)
	}
	pl.X.Label.Text = xlabel
	pl.Y.Label.Text = ylabel

	plotFn := func(mat *mat64.Dense, name string, color color.RGBA) {
		r, _ := mat.Dims()
		if Xs == nil {
			log.Fatal("must specify xs")
			Xs = make([]float64, r)
			for i := 0; i < r; i++ {
				Xs[i] = float64(i)
			}
		}

		col := make([]float64, r)
		for k, idx := range p.Idxs {
			mat64.Col(col, idx, mat)
			px := make(plotter.XYs, len(Xs))
			for j, x := range Xs {
				px[j].X = x
				px[j].Y = col[j]
			}
			l, err := plotter.NewScatter(px)
			if err != nil {
				panic(err)
			}

			l.GlyphStyle.Color = color
			l.GlyphStyle.Radius = vg.Length(vg.Points(1))
			l.Color = color
			pl.Add(l)
			if k == 0 {
				pl.Legend.Add(name, l)
			}
		}

	}
	wfn := func(mat *mat64.Dense) {
		plotFn(mat, "before", color.RGBA{225, 0, 4, 255})
		fn(mat)
		plotFn(mat, "after", color.RGBA{4, 0, 225, 255})
		if err := pl.Save(vg.Length(12)*vg.Inch, vg.Length(8)*vg.Inch, path); err != nil {
			panic(err)
		}

	}
	return wfn
}

// mean without the lowest and highest values.
func mean(in []float64) float64 {

	min := math.MaxFloat32
	max := float64(0)
	s := float64(0)
	for _, v := range in {
		if v > max {
			max = v
		}
		if v < min {
			min = v
		}

		s += v
	}
	//s -= (min + max)
	return s / float64(len(in))
}

func pLess(Depths []float64, val float64) float64 {
	var c float64
	for _, d := range Depths {
		if d < val {
			c++
		}
	}
	return c / float64(len(Depths))
}

// CallCopyNumbers returns Intervals for which any sample has non-zero copy-number
func (ivs *Intervals) CallCopyNumbers() {
	samples := ivs.Samples

	cache := &emdepth.Cache{}
	nskip := 0
	r, c := ivs.Depths.Dims()
	row32 := make([]float32, c)
	for i := 0; i < r; i++ {
		row := ivs.Depths.RawRowView(i)
		//fmt.Fprintf(os.Stdout, "%s\t%d\t%d\t%s\n", ivs.Chrom, iv.Start, iv.End, formatIV(iv))
		/*
			if pLess(row, 7) > 0.5 {
				nskip++
				continue
			}
			if mean(row) < 15 {
				nskip++
				continue
			}
		*/
		for k, d := range row {
			row32[k] = float32(d)
		}

		em := emdepth.EMDepth(row32, emdepth.Position{Start: ivs.Starts[i], End: ivs.Ends[i]})
		cnvs := cache.Add(em)
		ivs.printCNVs(cnvs, samples)
	}

	ivs.printCNVs(cache.Clear(nil), samples)
	log.Println("skipped:", nskip)

	//fmt.Fprintf(os.Stdout, "%s\t%d\t%d\t%s\t%s\n", ivs.Chrom, last.Start, last.End, formatCns(emdepth.EMDepth(last.AdjustedDepths)), formatFloats(last.AdjustedDepths))
}

const MinSize = 500

func (ivs *Intervals) printCNVs(cnvs []*emdepth.CNV, samples []string) {
	if len(cnvs) == 0 {
		return
	}
	fs := make([]string, 0, len(samples))
	fjoin := func(sl []float32) string {
		fs = fs[:0]
		for _, v := range sl {
			fs = append(fs, fmt.Sprintf("%.2f", v))
		}
		return strings.Join(fs, ",")
	}
	ijoin := func(sl []int) string {
		fs = fs[:0]
		for _, v := range sl {
			fs = append(fs, strconv.Itoa(v))
		}
		return strings.Join(fs, ",")
	}
	sort.Slice(cnvs, func(i, j int) bool { return cnvs[i].Position[0].Start < cnvs[j].Position[0].Start })
	for _, cnv := range cnvs {
		l := len(cnv.Position) - 1
		if cnv.Position[l].End-cnv.Position[0].Start < MinSize {
			continue
		}
		sample := samples[cnv.SampleI]
		cn := bestCN(cnv.CN)
		if cn == 2 {
			continue
		}
		fmt.Fprintf(os.Stdout, "%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%d\n", ivs.Chrom, cnv.Position[0].Start, cnv.Position[l].End,
			cn, sample, ijoin(cnv.CN), fjoin(cnv.Depth), fjoin(cnv.Log2FC), cnv.PSize)
	}
}

func bestCN(cns []int) int {
	var m float64
	for _, v := range cns {
		m += float64(v)
	}
	m /= float64(len(cns))
	return int(m + 0.5)
}

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
	ivs.Depths = mat64.NewDense(len(ivs._depths)/ns, ns, ivs._depths)
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

func main() {

	/*
		f, err := os.Create("dcnv.cpu.pprof")
		if err != nil {
			panic(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	*/

	bed := os.Args[1]
	fasta := os.Args[2]
	ivs := &Intervals{GCB: &GcBounds{Min: 0.25, Max: 0.75}}

	ivs.ReadRegions(bed, fasta)

	//ivs.CorrectBySampleMedian()

	//db := debiaser.GeneralDebiaser{}
	db := debiaser.ChunkDebiaser{
		ScoreWindow: 0.01}
	db.Window = 1
	db.Vals = make([]float64, len(ivs.GCs))
	copy(db.Vals, ivs.GCs)
	zsc := &scalers.ZScore{}
	l2 := &scalers.Log2{}
	sc := l2
	_ = sc
	//Pipeliner(ivs.Depths, sc.Scale)
	pl := Plotter{Idxs: []int{5}}
	_ = pl

	_, _ = zsc, l2
	//sc.Scale(ivs.Depths)

	// Correct by GC
	db.Window = 1
	Pipeliner(ivs.Depths, db.Sort, pl.Wrap(db.Debias, db.Vals, "GC", "normalized depth", "gc.png"), db.Unsort)
	//copy(db.Vals, ivs.SeqComplexity)
	//Pipeliner(ivs.Depths, db.Sort, pl.Wrap(db.Debias, db.Vals, "complexity", "normalized depth", "cpx.png"), db.Unsort)
	//sc.UnScale(ivs.Depths)

	fmt.Println("#chrom\tstart\tend\tcn\tsample\tcns\tdepths\tl2s\tn-regions")
	ivs.CallCopyNumbers()
}
