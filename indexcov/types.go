package indexcov

import (
	"reflect"
	"unsafe"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
)

const (
	// TileWidth is the length of the interval tiling used
	// in BAI and tabix indexes.
	TileWidth = 0x4000

	// StatsDummyBin is the bin number of the reference
	// statistics bin used in BAI and tabix indexes.
	StatsDummyBin = 0x924a
)

type oRefIndex struct {
	Bins      []bin
	Stats     *referenceStats
	Intervals []bgzf.Offset
}

type bin struct {
	Bin    uint32
	Chunks []bgzf.Chunk
}

type referenceStats struct {
	// Chunk is the span of the indexed BGZF
	// holding alignments to the reference.
	Chunk bgzf.Chunk

	// Mapped is the count of mapped reads.
	Mapped uint64

	// Unmapped is the count of unmapped reads.
	Unmapped uint64
}

func getRefs(idx *bam.Index) [][]int64 {
	refs := reflect.ValueOf(*idx).FieldByName("idx").FieldByName("Refs")
	ptr := unsafe.Pointer(refs.Pointer())

	ret := (*(*[1 << 28]oRefIndex)(ptr))[:refs.Len()]
	// save some memory.
	m := make([][]int64, len(ret))
	for i, r := range ret {
		m[i] = make([]int64, len(r.Intervals))
		for k, iv := range r.Intervals {
			m[i][k] = vOffset(iv)
		}
		r.Bins, r.Intervals = nil, nil
	}
	return m
}
