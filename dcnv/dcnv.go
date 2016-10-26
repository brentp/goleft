package main

import (
	"fmt"
	"io"
	"log"
	"math"
	"math/rand"
	"os"
	"sort"
	"strconv"
	"strings"
	"time"

	"github.com/JaderDias/movingmedian"
	"github.com/brentp/faidx"
	"github.com/brentp/goleft/emdepth"
	"github.com/brentp/xopen"
)

type interval struct {
	chrom         string
	start         uint32
	end           uint32
	depth         []float32
	log2          []float32
	GC            float32
	adjustedDepth []float32
	out           int
}

func (i *interval) setAdjusted() {
	if len(i.adjustedDepth) == 0 {
		i.adjustedDepth = make([]float32, len(i.depth))
	}
	for k, l2 := range i.log2 {
		i.adjustedDepth[k] = float32(math.Pow(2, float64(l2)))
	}
}

func (i *interval) String() string {
	s := fmt.Sprintf("%.1f", i.log2)
	if len(s) < 4 {
		s += " "
	}
	return s
}

func sortByGC(regions []*interval) {
	sort.Slice(regions, func(i, j int) bool { return regions[i].GC < regions[j].GC })
}
func sortRandom(regions []*interval) {
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

func intervalFromLine(l string, fa *faidx.Faidx) *interval {
	toks := strings.Split(l, "\t")
	iv := &interval{chrom: toks[0], start: mustAtoi(toks[1]), end: mustAtoi(toks[2]),
		depth: make([]float32, 0, len(toks)-3),
		log2:  make([]float32, 0, len(toks)-3),
	}
	for c := 3; c < len(toks); c++ {
		d := mustAtof(toks[3])
		//d *= float32(iv.end - iv.start)
		//iv.depth = uint32(d)
		iv.depth = append(iv.depth, d)
		iv.log2 = append(iv.log2, float32(math.Log2(float64(d))))
	}
	st, err := fa.Stats(iv.chrom, int(iv.start), int(iv.end))
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
	}
	iv.GC = float32(st.GC)
	return iv
}

func medianLog2BySample(regions []*interval, sampleI int) float32 {
	sort.Slice(regions, func(i, j int) bool { return regions[i].log2[sampleI] < regions[j].log2[sampleI] })
	return regions[len(regions)/2].log2[sampleI]
}

// CorrectBySampleMedian subtracts the sample median from each sample.
func CorrectBySampleMedian(regions []*interval) {
	n := len(regions[0].depth)
	for i := 0; i < n; i++ {
		m := medianLog2BySample(regions, i)
		for _, r := range regions {
			r.log2[i] -= m
		}
	}
}

// SetAdjustedDepths sets the adjustedDepth for each sample based on the log2.
func SetAdjustedDepths(regions []*interval) {
	for _, r := range regions {
		r.setAdjusted()
	}
}

// CorrectByGC sorts so that intervals with similar GC are grouped together
// and then docs a moving median correction on the log2 of the coverage.
func CorrectByGC(regions []*interval, window int) {
	// sort random to make sure adjacent true sites are randomized away from each other.
	sortRandom(regions)
	sortByGC(regions)
	for sampleI := 0; sampleI < len(regions[0].depth); sampleI++ {
		correctByMovingMedian(regions, window, sampleI)
	}
}

// after sorting be GC, this is used to adjust log2s to subtract bias (subtract the median).
func correctByMovingMedian(regions []*interval, window int, sampleI int) {
	mm := movingmedian.NewMovingMedian(window)
	mid := (window-1)/2 + 1
	for i := 0; i < mid; i++ {
		mm.Push(float64(regions[i].log2[sampleI]))
	}

	var i int
	for i = 0; i < len(regions)-mid; i++ {
		mm.Push(float64(regions[i+mid].log2[sampleI]))
		regions[i].log2[sampleI] -= float32(mm.Median())
	}
	for ; i < len(regions); i++ {
		regions[i].log2[sampleI] -= float32(mm.Median())
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

// CallCopyNumbers returns intervals for which any sample has non-zero copy-number
func CallCopyNumbers(r []*interval) {
	for i, iv := range r {
		cns := emdepth.EMDepth(iv.adjustedDepth)
		if all2(cns) {
			continue
		}
		fmt.Fprintln(os.Stderr, i, iv, cns)
	}
}

func readRegions(path string, fasta string) []*interval {
	fai, err := faidx.New(fasta)
	if err != nil {
		log.Fatal(err)
	}
	m := make([]*interval, 0, 100000)
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
	return m
}

func main() {

	window := 11
	bed := os.Args[1]
	fasta := os.Args[2]
	ivs := readRegions(bed, fasta)
	CorrectBySampleMedian(ivs)
	CorrectByGC(ivs, window)
	SetAdjustedDepths(ivs)
	CallCopyNumbers(ivs)

}
