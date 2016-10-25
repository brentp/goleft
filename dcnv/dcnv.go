package main

import (
	"fmt"
	"math"
	"math/rand"
	"sort"
	"strconv"
	"strings"
	"time"

	"github.com/JaderDias/movingmedian"
)

type interval struct {
	chrom string
	start uint32
	end   uint32
	depth float32
	log2  float32
	GC    float32
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

func intervalFromLine(l string) *interval {
	toks := strings.SplitN(l, "\t", 5)

	iv := &interval{chrom: toks[0], start: mustAtoi(toks[1]), end: mustAtoi(toks[2])}
	d := mustAtof(toks[3])
	//d *= float32(iv.end - iv.start)
	//iv.depth = uint32(d)
	iv.depth = d
	iv.log2 = float32(math.Log2(float64(iv.depth)))
	iv.GC = mustAtof(toks[4])
	return iv
}

func medianLog2(regions []*interval) float32 {
	sort.Slice(regions, func(i, j int) bool { return regions[i].log2 < regions[j].log2 })
	return regions[len(regions)/2].log2
}

func correctBySampleMedian(regions []*interval) {
	m := medianLog2(regions)
	for _, i := range regions {
		i.log2 -= m
	}
}

func correctByGC(regions []*interval, window int) {
	// sort random to make sure adjacent true sites are randomized away from each other.
	sortRandom(regions)
	sortByGC(regions)
	correctByMovingMedian(regions, window)
}

func correctByMovingMedian(regions []*interval, window int) {
	mm := movingmedian.NewMovingMedian(window)
	mid := (window-1)/2 + 1
	for i := 0; i < mid; i++ {
		mm.Push(float64(regions[i].log2))
	}

	var i int
	for i = 0; i < len(regions)-mid; i++ {
		mm.Push(float64(regions[i+mid].log2))
		regions[i].log2 -= float32(mm.Median())
	}
	for ; i < len(regions); i++ {
		regions[i].log2 -= float32(mm.Median())
	}

	// read in sample
	// log2
	// normalize each sample's log2 to its median
	// sort random
	// sort by gc
	// normalize my movingmedian
}

func main() {

}
