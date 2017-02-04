package main

import (
	"fmt"
	"io"
	"log"
	"math"
	"math/rand"
	"os"
	"strconv"
	"strings"
	"time"

	"github.com/JaderDias/movingmedian"
	"github.com/brentp/faidx"
	"github.com/brentp/goleft/dcnv/scalers"
	"github.com/brentp/goleft/emdepth"
	"github.com/brentp/xopen"
	"github.com/gonum/matrix"
	"github.com/gonum/matrix/mat64"
	"go4.org/sort"
)

// Interval is the struct used by dcnv
type Interval struct {
	Start          uint32
	End            uint32
	Depths         []float32
	GC             float32
	AdjustedDepths []float32
}

// Intervals is the wrapper for a slice of intervals.
type Intervals struct {
	Chrom         string
	Intervals     []*Interval
	sampleMedians []float32
	sampleScalars []float32
	Samples       []string
}

func (ivs Intervals) NSamples() int {
	return len(ivs.Intervals[0].Depths)
}

func sortByGC(ivs *Intervals) {
	sort.Slice(ivs.Intervals, func(i, j int) bool { return ivs.Intervals[i].GC < ivs.Intervals[j].GC })
}

func sortRandom(ivs *Intervals) {
	t := time.Now()
	regions := ivs.Intervals
	rand.Seed(int64(t.Nanosecond()))
	for i := range regions {
		j := rand.Intn(i + 1)
		regions[i], regions[j] = regions[j], regions[i]
	}
}

func mustAtoi(a string) uint32 {
	v, err := strconv.Atoi(a)
	if err != nil {
		panic(err)
	}
	return uint32(v)
}

func mustAtof(a string) float32 {
	v, err := strconv.ParseFloat(a, 32)
	if err != nil {
		panic(err)
	}
	return float32(v)
}

func intervalFromLine(l string, fa *faidx.Faidx) *Interval {
	toks := strings.Split(l, "\t")
	toks[len(toks)-1] = strings.TrimSpace(toks[len(toks)-1])
	iv := &Interval{Start: mustAtoi(toks[1]), End: mustAtoi(toks[2]),
		Depths:         make([]float32, 0, len(toks)-3),
		AdjustedDepths: make([]float32, 0, len(toks)-3),
	}
	for c := 3; c < len(toks); c++ {
		d := mustAtof(toks[c])
		//d /= float32(iv.End - iv.Start)
		iv.Depths = append(iv.Depths, d)
		iv.AdjustedDepths = append(iv.AdjustedDepths, d)
	}
	// subtract $n bases since GC before will afffect reads here.
	st, err := fa.Stats(toks[0], int(iv.Start)-100, int(iv.End))
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
	}
	iv.GC = float32(st.GC)
	return iv
}

// SampleMedians gets the Median log2 values for each sample.
func (ivs *Intervals) SampleMedians() []float32 {
	regions := ivs.Intervals
	depths := make([]float32, 0, len(regions))

	ivs.sampleMedians = make([]float32, ivs.NSamples())
	for sampleI := 0; sampleI < ivs.NSamples(); sampleI++ {
		// sorting the extracted array is much faster.
		depths = depths[:0]
		for _, d := range regions {
			if d.AdjustedDepths[sampleI] > 0 {
				depths = append(depths, d.AdjustedDepths[sampleI])
			}
		}
		// lop off the lower depths (for exome).
		// and then normalized on the median above that lower bound.
		sort.Slice(depths, func(i, j int) bool { return depths[i] < depths[j] })
		ivs.sampleMedians[sampleI] = depths[int(0.65*float64(len(depths)))]
	}
	med := median(ivs.sampleMedians)
	ivs.sampleScalars = make([]float32, 0, len(ivs.sampleMedians))
	for _, sm := range ivs.sampleMedians {
		ivs.sampleScalars = append(ivs.sampleScalars, med/sm)
	}
	return ivs.sampleMedians
}

func (ivs *Intervals) SampleScalars() []float32 {
	ivs.SampleMedians()
	return ivs.sampleScalars
}

// CorrectBySampleMedian subtracts the sample median from each sample.
func (ivs *Intervals) CorrectBySampleMedian() {
	scalars := ivs.SampleScalars()
	for _, r := range ivs.Intervals {
		for i, s := range scalars {
			r.AdjustedDepths[i] *= s
		}
	}
}

func median(b []float32) float32 {
	a := make([]float32, len(b))
	for i, bb := range b {
		a[i] = bb
	}

	sort.Slice(a, func(i, j int) bool { return a[i] < a[j] })
	if len(a)%2 == 0 {
		am := a[len(a)/2]
		bm := a[len(a)/2+1]
		return (am + bm) / 2
	}
	return a[len(a)/2]
}

// CorrectByGC sorts so that Intervals with similar GC are grouped together
// and then docs a moving median correction on the log2 of the coverage.
func (ivs *Intervals) CorrectByGC(window int) {
	// sort random to make sure adjacent true sites are randomized away from each other.
	sortRandom(ivs)
	sortByGC(ivs)
	zscore := &scalers.ZScore{}
	r, _ := len(ivs.Intervals), len(ivs.Samples)
	dps := make([]float64, r)
	z := ivs.IntoMatrix(nil)
	zscore.Scale(z)
	for sampleI := 0; sampleI < ivs.NSamples(); sampleI++ {
		mat64.Col(dps, sampleI, z)
		corrected := correctByMovingMedian(window, dps)
		z.SetCol(sampleI, corrected)
	}
	zscore.UnScale(z)
	ivs.FromMatrix(z)
}

// after sorting by GC, this is used to adjust zs to subtract bias (subtract the median).
// it's assumed that this is sorted by whatever effect is to be removed.
func correctByMovingMedian(window int, dps []float64) []float64 {

	mm := movingmedian.NewMovingMedian(window)
	mid := (window-1)/2 + 1
	for i := 0; i < mid; i++ {
		mm.Push(dps[i])
	}

	var i int
	for i = 0; i < len(dps)-mid; i++ {
		mm.Push(dps[i+mid])
		dps[i] -= mm.Median()
	}
	for ; i < len(dps); i++ {
		dps[i] -= mm.Median()
	}
	return dps
}

// mean without the lowest and highest values.
func mean(in []float32) float32 {

	min := float32(math.MaxFloat32)
	max := float32(0)
	s := float32(0)
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
	return s / float32(len(in))
}

// all2 returns true if all values in the slice are == 2.
func all2(cns []int) bool {
	for _, c := range cns {
		if c != 2 {
			return false
		}
	}
	return true
}

// all2 returns true if all values in the slice are == 2.
func pLess(Depths []float32, val float32) float32 {
	var c float32
	for _, d := range Depths {
		if d < val {
			c++
		}
	}
	return c / float32(len(Depths))
}

// CallCopyNumbers returns Intervals for which any sample has non-zero copy-number
func (ivs *Intervals) CallCopyNumbers() {
	ivs.SortByPosition()
	samples := ivs.Samples

	cache := &emdepth.Cache{}
	nskip := 0
	for _, iv := range ivs.Intervals {
		//fmt.Fprintf(os.Stdout, "%s\t%d\t%d\t%s\n", ivs.Chrom, iv.Start, iv.End, formatIV(iv))
		if pLess(iv.Depths, 7) > 0.5 {
			nskip++
			continue
		}
		if mean(iv.Depths) < 15 {
			nskip++
			continue
		}

		em := emdepth.EMDepth(iv.AdjustedDepths, emdepth.Position{Start: iv.Start, End: iv.End})
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

func (ivs *Intervals) SortByPosition() {
	sort.Slice(ivs.Intervals, func(i, j int) bool { return ivs.Intervals[i].Start < ivs.Intervals[j].Start })
}

func svdmean(vs []float32) float64 {
	s := float64(0)
	for _, v := range vs {
		s += float64(v)
	}
	return s / float64(len(vs))
}

func svdsd(vs []float32, mean float64) float64 {
	sd := float64(0)
	for _, v := range vs {
		sd += math.Pow(float64(v)-mean, 2)
	}
	return math.Sqrt(sd / float64(len(vs)-1))
}

func (ivs *Intervals) FromMatrix(mat *mat64.Dense) {
	for i, iv := range ivs.Intervals {
		row := mat.RawRowView(i)
		for k, v := range row {
			iv.AdjustedDepths[k] = max32(0, float32(v))
		}
	}
}

func (ivs *Intervals) IntoMatrix(mat *mat64.Dense) *mat64.Dense {
	if mat == nil {
		mat = mat64.NewDense(len(ivs.Intervals), len(ivs.Samples), nil)
	}
	for i, iv := range ivs.Intervals {
		row := mat.RawRowView(i)
		for i, d := range iv.AdjustedDepths {
			row[i] = float64(d)
		}
	}
	return mat
}

func (ivs *Intervals) SVD(n int) {
	zscore := &scalers.ZScore{}
	log.Println(ivs.Intervals[0].AdjustedDepths)
	z := ivs.IntoMatrix(nil)
	zscore.Scale(z)
	var svd mat64.SVD
	log.Println(z.RowView(0))
	if ok := svd.Factorize(z, matrix.SVDThin); !ok {
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

	z.Product(u, sigma, v)
	zscore.UnScale(z)
	ivs.FromMatrix(z)
}

func max32(a, b float32) float32 {
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
	m := make([]*Interval, 0, 100000)
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
		iv := intervalFromLine(line, fai)
		if medLE1(iv.Depths) {
			continue
		}
		m = append(m, iv)
	}
	ivs.Intervals = m
}

func medLE1(dps []float32) bool {
	var n float64
	for _, v := range dps {
		if v <= 20 {
			n++
		}
	}
	return n/float64(len(dps)) >= 0.5
}

func (ivs *Intervals) Write(n int) {
	meds := ivs.SampleMedians()
	lmeds := make([]float32, len(meds))
	for i, v := range meds {
		lmeds[i] = float32(math.Log2(float64(v)))
	}
	_s := make([]string, len(ivs.Intervals[0].Depths))
	formatIV := func(i *Interval) string {
		for k := 0; k < len(_s); k++ {
			_s[k] = fmt.Sprintf("%.0f:%.1f", i.Depths[k], i.AdjustedDepths[k])
		}
		return strings.Join(_s, "\t")
	}
	ivs.SortByPosition()
	for i, iv := range ivs.Intervals {
		if i == n {
			break
		}
		fmt.Printf("%s:%d-%d\t%s\n", ivs.Chrom, iv.Start, iv.End, formatIV(iv))
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
	//ivs.CorrectByGC(window)
	ivs.SortByPosition()
	ivs.SVD(7)
	log.Println(ivs.SampleMedians())

	nsites := len(ivs.Intervals)
	dps := make([]string, len(ivs.Samples))
	fmt.Printf("#chrom\tstart\tend\t%s\n", strings.Join(ivs.Samples, "\t"))
	for i := 0; i < nsites; i++ {
		iv := ivs.Intervals[i]
		for si := range dps {
			dps[si] = fmt.Sprintf("%.2f", iv.AdjustedDepths[si])
		}
		fmt.Printf("%s\t%d\t%d\t%s\n", ivs.Chrom, iv.Start, iv.End, strings.Join(dps, "\t"))
	}

	//	ivs.CallCopyNumbers()
}
