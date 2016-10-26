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
	"github.com/brentp/goleft/emdepth"
	"github.com/brentp/xopen"
	"go4.org/sort"
)

// Interval is the struct used by dcnv
type Interval struct {
	Chrom          string
	Start          uint32
	End            uint32
	Depths         []float32
	Log2s          []float32
	GC             float32
	AdjustedDepths []float32
	cns            []int
}

func (i *Interval) copy(cns []int) *Interval {
	n := len(i.Depths)
	c := &Interval{Chrom: i.Chrom, Start: i.Start, End: i.End, Depths: make([]float32, n),
		Log2s: make([]float32, n), GC: i.GC, AdjustedDepths: make([]float32, n)}
	for k := 0; k < n; k++ {
		c.Log2s[k] = i.Log2s[k]
		c.Depths[k] = i.Depths[k]
		c.AdjustedDepths[k] = i.AdjustedDepths[k]
	}
	c.cns = cns
	return c
}

// combine the data from 2 intervals.
func (i *Interval) update(b *Interval) {
	n := len(i.Depths)
	ilen := float32(i.End - i.Start)
	blen := float32(b.End - b.Start)
	tot := ilen + blen
	for k := 0; k < n; k++ {
		i.Log2s[k] = (ilen*i.Log2s[k] + blen*b.Log2s[k]) / tot
		i.Depths[k] = (ilen*i.Depths[k] + blen*b.Depths[k]) / tot
		i.AdjustedDepths[k] = (ilen*i.AdjustedDepths[k] + blen*b.AdjustedDepths[k]) / tot
	}

}

// Intervals is the wrapper for a slice of intervals.
type Intervals struct {
	Intervals     []*Interval
	sampleMedians []float32
}

func (ivs Intervals) NSamples() int {
	return len(ivs.Intervals[0].Depths)
}

func sortByGC(regions []*Interval) {
	sort.Slice(regions, func(i, j int) bool { return regions[i].GC < regions[j].GC })
}

func sortRandom(regions []*Interval) {
	t := time.Now()
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
	iv := &Interval{Chrom: toks[0], Start: mustAtoi(toks[1]), End: mustAtoi(toks[2]),
		Depths: make([]float32, 0, len(toks)-3),
		Log2s:  make([]float32, 0, len(toks)-3),
	}
	for c := 3; c < len(toks); c++ {
		d := mustAtof(toks[c])
		//d /= float32(iv.End - iv.Start)
		//iv.Depths = uint32(d)
		iv.Depths = append(iv.Depths, d)
		if d == 0 {
			iv.Log2s = append(iv.Log2s, -math.MaxFloat32)
		} else {
			iv.Log2s = append(iv.Log2s, float32(math.Log2(float64(d))))
		}
	}
	st, err := fa.Stats(iv.Chrom, int(iv.Start), int(iv.End))
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
	}
	iv.GC = float32(st.GC)
	return iv
}

// CalcSampleMedians sets the Median log2 values for each sample.
func (ivs *Intervals) SampleMedians() []float32 {
	regions := ivs.Intervals
	if len(ivs.sampleMedians) == 0 {
		ivs.sampleMedians = make([]float32, ivs.NSamples())
		for sampleI := 0; sampleI < ivs.NSamples(); sampleI++ {
			sort.Slice(regions, func(i, j int) bool { return regions[i].Log2s[sampleI] < regions[j].Log2s[sampleI] })
			ivs.sampleMedians[sampleI] = regions[len(regions)/2].Log2s[sampleI]
		}
	}
	return ivs.sampleMedians
}

// CorrectBySampleMedian subtracts the sample median from each sample.
func (ivs *Intervals) CorrectBySampleMedian() {
	meds := ivs.SampleMedians()
	for i := 0; i < ivs.NSamples(); i++ {
		m := meds[i]
		for _, r := range ivs.Intervals {
			r.Log2s[i] -= m
		}
	}
}

// SetAdjustedDepths sets the AdjustedDepths for each sample based on the log2.
func (ivs Intervals) SetAdjustedDepths() {
	meds := ivs.SampleMedians()
	for _, i := range ivs.Intervals {
		if len(i.AdjustedDepths) == 0 {
			i.AdjustedDepths = make([]float32, len(i.Depths))
		}
		for k, l2 := range i.Log2s {
			i.AdjustedDepths[k] = meds[k] + float32(math.Pow(2, float64(l2)))
		}
	}
}

// CorrectByGC sorts so that Intervals with similar GC are grouped together
// and then docs a moving median correction on the log2 of the coverage.
func (ivs *Intervals) CorrectByGC(window int) {
	// sort random to make sure adjacent true sites are randomized away from each other.
	sortRandom(ivs.Intervals)
	sortByGC(ivs.Intervals)
	for sampleI := 0; sampleI < len(ivs.Intervals[0].Depths); sampleI++ {
		correctByMovingMedian(ivs, window, sampleI)
	}
}

// after sorting by GC, this is used to adjust log2s to subtract bias (subtract the median).
func correctByMovingMedian(ivs *Intervals, window int, sampleI int) {
	regions := ivs.Intervals
	mm := movingmedian.NewMovingMedian(window)
	mid := (window-1)/2 + 1
	for i := 0; i < mid; i++ {
		mm.Push(float64(regions[i].Log2s[sampleI]))
	}

	var i int
	for i = 0; i < len(regions)-mid; i++ {
		mm.Push(float64(regions[i+mid].Log2s[sampleI]))
		regions[i].Log2s[sampleI] -= float32(mm.Median())
	}
	for ; i < len(regions); i++ {
		regions[i].Log2s[sampleI] -= float32(mm.Median())
	}
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
func all0(Depths []float32) bool {
	for _, d := range Depths {
		if d != 0 {
			return false
		}
	}
	return true
}

// CallCopyNumbers returns Intervals for which any sample has non-zero copy-number
func (ivs *Intervals) CallCopyNumbers() {
	ivs.SortByPosition()

	formatCns := func(cns []int) string {
		s := make([]string, 0, len(cns))
		for _, c := range cns {
			s = append(s, strconv.Itoa(c))
		}
		return strings.Join(s, "\t")
	}
	formatFloats := func(vs []float32) string {
		s := make([]string, 0, len(vs))
		for _, v := range vs {
			s = append(s, fmt.Sprintf("%.1f", v))
		}
		return strings.Join(s, "\t")
	}

	r := ivs.Intervals
	//cache := make([]*Intervals, 0, 8)
	var last *Interval
	for _, iv := range r {
		if all0(iv.Depths) {
			continue
		}

		cns := emdepth.EMDepth(iv.AdjustedDepths)
		if all2(cns) {
			continue
		}
		if last == nil {
			last = iv.copy(cns)
		} else if last.End == iv.Start && allEqual(cns, last.cns) {
			last.update(iv)
		} else {
			fmt.Fprintf(os.Stdout, "%s\t%d\t%d\t%s\t%s\n", last.Chrom, last.Start, last.End, formatCns(emdepth.EMDepth(last.AdjustedDepths)), formatFloats(last.AdjustedDepths))
			last = iv.copy(cns)
		}
	}
	fmt.Fprintf(os.Stdout, "%s\t%d\t%d\t%s\t%s\n", last.Chrom, last.Start, last.End, formatCns(emdepth.EMDepth(last.AdjustedDepths)), formatFloats(last.AdjustedDepths))
	log.Println(len(r))
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

func readRegions(path string, fasta string) *Intervals {
	fai, err := faidx.New(fasta)
	if err != nil {
		log.Fatal(err)
	}
	m := make([]*Interval, 0, 100000)
	rdr, err := xopen.Ropen(path)
	if err != nil {
		log.Fatal(err)
	}
	i := 0
	for {
		line, err := rdr.ReadString('\n')
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatal(err)
		}
		if i == 0 && line[0] == '#' || strings.HasPrefix(line, "chrom") {
			continue
		}
		i++
		m = append(m, intervalFromLine(line, fai))
	}
	return &Intervals{Intervals: m}
}

func main() {

	window := 5
	bed := os.Args[1]
	fasta := os.Args[2]
	ivs := readRegions(bed, fasta)
	_ = window
	ivs.CorrectBySampleMedian()
	// TODO: correct by Median across all samples
	ivs.CorrectByGC(window)
	ivs.SetAdjustedDepths()
	ivs.CallCopyNumbers()

}
