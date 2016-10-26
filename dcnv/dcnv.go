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
		d /= float32(iv.End - iv.Start)
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

// after sorting be GC, this is used to adjust log2s to subtract bias (subtract the median).
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
	s -= (min + max)
	return s / float32(len(in)-2)
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

	r := ivs.Intervals
	for _, iv := range r {
		if all0(iv.Depths) {
			continue
		}

		cns := emdepth.EMDepth(iv.AdjustedDepths)
		if all2(cns) {
			continue
		}
		fmt.Fprintf(os.Stdout, "%s:%d-%d %v %v %v\n", iv.Chrom, iv.Start, iv.End, cns, iv.AdjustedDepths, iv.Depths)
	}
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

	window := 11
	bed := os.Args[1]
	fasta := os.Args[2]
	ivs := readRegions(bed, fasta)
	ivs.CorrectBySampleMedian()
	ivs.CorrectByGC(window)
	ivs.SetAdjustedDepths()
	ivs.CallCopyNumbers()

}
