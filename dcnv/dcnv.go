package main

import (
	"fmt"
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
	"github.com/gonum/matrix"
	"github.com/gonum/matrix/mat64"
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

func (ivs *Intervals) addFromLine(l string, fa *faidx.Faidx) {
	toks := strings.Split(l, "\t")
	toks[len(toks)-1] = strings.TrimSpace(toks[len(toks)-1])

	ivs.Starts = append(ivs.Starts, mustAtoi(toks[1]))
	ivs.Ends = append(ivs.Ends, mustAtoi(toks[2]))

	for c := 3; c < len(toks); c++ {
		d := mustAtof(toks[c])
		//d /= float64(iv.End - iv.Start)
		ivs._depths = append(ivs._depths, d)
	}
	// subtract $n bases since GC before will afffect reads here.
	last := len(ivs.Starts) - 1
	st, err := fa.Stats(toks[0], int(ivs.Starts[last])-100, int(ivs.Ends[last])-100)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
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
	scalars := ivs.SampleScalars()
	r, c := ivs.Depths.Dims()
	mat := ivs.Depths
	for i := 0; i < r; i++ {
		row := mat.RawRowView(i)
		for j := 0; j < c; j++ {
			row[j] *= scalars[j]
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

func (ivs *Intervals) CorrectByGC(window int) {
	// TODO: put the debiaser on the Intervals object append
	// figure out how to make this Scale->Sort->Debias->Unsort->Unscale more sane.
	db := debiaser.GeneralDebiaser{
		Posns:  ivs.Starts,
		Vals:   ivs.GCs,
		Window: window}
	zsc := &scalers.ZScore{}
	Pipeliner(ivs.Depths, zsc.Scale, db.Sort, db.Debias, db.Unsort, zsc.UnScale)
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
		if pLess(row, 7) > 0.5 {
			nskip++
			continue
		}
		if mean(row) < 15 {
			nskip++
			continue
		}
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
		fmt.Fprintf(os.Stdout, "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%d\n", ivs.Chrom, cnv.Position[0].Start, cnv.Position[l].End,
			sample, ijoin(cnv.CN), fjoin(cnv.Depth), fjoin(cnv.Log2FC), cnv.PSize)
	}
}

func allEqual(a, b []int) bool {
	for i, v := range a {
		if v != b[i] {
			return false
		}
	}
	return true
}

func (ivs *Intervals) SVD(n int) {
	zscore := &scalers.ZScore{}
	zscore.Scale(ivs.Depths)
	var svd mat64.SVD
	if ok := svd.Factorize(ivs.Depths, matrix.SVDThin); !ok {
		panic("error with SVD")
	}

	// get svd and zero out first n components.
	s, u, v := extractSVD(&svd)
	for i := 0; i < n; i++ {
		s[i] = 0
	}
	sigma := mat64.NewDense(len(s), len(s), nil)
	for i := 0; i < len(s); i++ {
		sigma.Set(i, i, s[i])
	}

	ivs.Depths.Product(u, sigma, v)
	zscore.UnScale(ivs.Depths)
}

func max32(a, b float64) float64 {
	if a > b {
		return a
	}
	return b
}

func extractSVD(svd *mat64.SVD) (s []float64, u, v *mat64.Dense) {
	var um, vm mat64.Dense
	um.UFromSVD(svd)
	vm.VFromSVD(svd)
	s = svd.Values(nil)
	return s, &um, &vm
}

func (ivs *Intervals) ReadRegions(path string, fasta string) {
	fai, err := faidx.New(fasta)
	if err != nil {
		panic(err)
	}
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
		}
		i++
		ivs.addFromLine(line, fai)
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

	window := 19
	bed := os.Args[1]
	fasta := os.Args[2]
	ivs := &Intervals{}
	ivs.ReadRegions(bed, fasta)
	fmt.Fprintln(os.Stderr, ivs.Samples)
	_ = window

	ivs.CorrectBySampleMedian()
	ivs.CorrectByGC(window)
	ivs.SVD(7)
	log.Println(ivs.SampleMedians())

	nsites, nsamples := ivs.Depths.Dims()
	dps := make([]string, nsamples)
	fmt.Printf("#chrom\tstart\tend\t%s\n", strings.Join(ivs.Samples, "\t"))
	for i := 0; i < nsites; i++ {
		iv := ivs.Depths.RawRowView(i)
		for si := range dps {
			dps[si] = fmt.Sprintf("%.2f", iv[si])
		}
		fmt.Printf("%s\t%d\t%d\t%s\n", ivs.Chrom, ivs.Starts[i], ivs.Ends[i], strings.Join(dps, "\t"))
	}
}
