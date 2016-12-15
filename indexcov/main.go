package main

import (
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"

	"github.com/biogo/hts/bam"
)

const TileWidth = 0x4000

// NormalizedDepth returns a list of numbers for the normalized depth of the given region.
// Values are scaled to have a mean of 1. If end is 0, the full chromosome is returned.
func NormalizedDepth(idx *bam.Index, refID int, start int, end int) []float32 {
	in := idx.Index()

	last := in.Refs[0]
	// get the last chromosome with any data.
	totalIntervals := 0
	for k := 0; k < len(in.Refs)-1; k++ {
		totalIntervals += len(in.Refs[k].Intervals)
		if len(in.Refs[k].Intervals) > 0 {
			last = in.Refs[k]
		}
	}
	// this gives the total file size.
	size := float64(last.Intervals[len(last.Intervals)-1].File + TileWidth)
	ref := in.Refs[refID]

	meanSizePerTile := size / float64(totalIntervals)
	log.Println(meanSizePerTile, len(last.Intervals))

	si, ei := start/TileWidth, end/TileWidth
	if end == 0 || ei >= len(ref.Intervals) {
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

func getTid(b *bam.Reader, chrom string) int {
	refs := b.Header().Refs()
	if strings.HasPrefix(chrom, "chr") {
		chrom = chrom[3:]
	}
	for i, ref := range refs {
		if chrom == ref.Name() {
			return i
		}
		if strings.HasPrefix(ref.Name(), "chr") {
			if chrom == ref.Name()[3:] {
				return i
			}
		}
	}
	return -1
}

func main() {
	chrom := os.Args[1]
	rdr, err := os.Open(os.Args[2])
	if err != nil {
		panic(err)
	}
	brdr, err := bam.NewReader(rdr, 1)
	if err != nil {
		panic(err)
	}

	tid := getTid(brdr, chrom)
	rdr.Close()
	brdr.Close()
	if tid == -1 {
		panic(fmt.Sprintf("unable to find chromosome: %s", chrom))
	}

	rdr, err = os.Open(os.Args[2] + ".bai")
	if err != nil {
		panic(err)
	}

	idx, err := bam.ReadIndex(rdr)
	if err != nil {
		panic(err)
	}

	s, err := strconv.Atoi(chrom)
	if err != nil {
		panic(err)
	}

	depths := NormalizedDepth(idx, s-1, 0, 0)
	//fmt.Println(mean(depths))
	for i := 0; i < len(depths); i++ {
		fmt.Printf("%s\t%d\t%d\t%.3f\n", chrom, i*16384, (i+1)*16384, depths[i])
	}
}
