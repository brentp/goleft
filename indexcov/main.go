package main

import (
	"fmt"
	"os"

	"github.com/biogo/hts/bam"
)

const TileWidth = 0x4000

// NormalizedDepth returns a list of numbers for the normalized depth of the given region.
// Values are scaled to have a mean of 1. If end is 0, the full chromosome is returned.
func NormalizedDepth(idx *bam.Index, refID int, start int, end int) []float32 {
	in := idx.Index()

	last := in.Refs[0]
	// get the last chromosome with any data.
	for k := len(in.Refs) - 1; k > 0; k-- {
		if len(in.Refs[k].Intervals) > 0 {
			last = in.Refs[k]
			break
		}
	}
	// this gives the total file size.
	size := float64(last.Intervals[len(last.Intervals)-1].File + TileWidth)
	meanSizePerTile := size / float64(TileWidth)

	si, ei := start/TileWidth, end/TileWidth
	ref := in.Refs[refID]
	if end == 0 {
		ei = len(ref.Intervals) - 1
	}
	if ei >= len(ref.Intervals) {
		ei = len(ref.Intervals) - 1
	}
	depths := make([]float32, 0, ei-si)
	for i, o := range ref.Intervals[si:ei] {
		depths = append(depths, float32(float64(ref.Intervals[si+i+1].File-o.File)/meanSizePerTile))
	}
	return depths
}

func mean(a []float32) float32 {
	s := float32(0)
	for _, v := range a {
		s += v / float32(len(a))
	}
	return s
}

func main() {
	rdr, err := os.Open("/data/human/15-0022949.sub.bam.bai")
	if err != nil {
		panic(err)
	}

	idx, err := bam.ReadIndex(rdr)
	if err != nil {
		panic(err)
	}

	depths := NormalizedDepth(idx, 19, 0, 0)
	fmt.Println(mean(depths))
	depths = NormalizedDepth(idx, 20, 0, 0)
	fmt.Println(mean(depths))

	depths = NormalizedDepth(idx, 21, 0, 0)
	fmt.Println(mean(depths))

	for _, r := range idx.Index().Refs {
		if r.Stats != nil {
			fmt.Println(r.Stats)
		}
	}
}
