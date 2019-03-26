package indexcov

import (
	"log"
	"reflect"
	"unsafe"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/bgzf/index"
	"github.com/biogo/hts/csi"
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

type csiRefIndex struct {
	bins  []csiBin
	stats *index.ReferenceStats
}

type csiBin struct {
	bin     uint32
	left    bgzf.Offset
	records uint64
	chunks  []bgzf.Chunk
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

func getCSISizes(idx *csi.Index) ([][]int64, uint64, uint64) {
	var mapped, unmapped uint64
	refs := reflect.ValueOf(*idx).FieldByName("refs")
	ptr := unsafe.Pointer(refs.Pointer())
	ret := (*(*[1 << 28]csiRefIndex)(ptr))[:refs.Len()]

	m := make([][]int64, len(ret))
	for i, r := range ret {
		st, ok := idx.ReferenceStats(i)
		if ok {
			mapped += st.Mapped
			unmapped += st.Unmapped
		} else {
			log.Printf("no reference stats found for %dth reference", i)
		}
		if len(r.bins) < 2 {
			m[i] = make([]int64, 0)
			continue
		}
		m[i] = make([]int64, len(r.bins)-1)
		for k, iv := range r.bins[1:] {
			m[i][k] = vOffset(iv.left) - vOffset(r.bins[k].left)
			if m[i][k] < 0 {
				panic("expected positive change in vOffset")
			}
		}
	}
	return m, mapped, unmapped
}

func getSizes(idx *bam.Index) ([][]int64, uint64, uint64) {
	var mapped, unmapped uint64
	refs := reflect.ValueOf(*idx).FieldByName("idx").FieldByName("Refs")
	ptr := unsafe.Pointer(refs.Pointer())

	ret := (*(*[1 << 28]oRefIndex)(ptr))[:refs.Len()]
	// save some memory.
	m := make([][]int64, len(ret))
	for i, r := range ret {
		st, ok := idx.ReferenceStats(i)
		if ok {
			mapped += st.Mapped
			unmapped += st.Unmapped
		} else {
			log.Printf("no reference stats found for %dth reference", i)
		}
		if len(r.Intervals) < 2 {
			m[i] = make([]int64, 0)
			continue
		}
		m[i] = make([]int64, len(r.Intervals)-1)
		for k, iv := range r.Intervals[1:] {
			m[i][k] = vOffset(iv) - vOffset(r.Intervals[k])
			if m[i][k] < 0 {
				panic("expected positive change in vOffset")
			}
		}
		r.Bins, r.Intervals = nil, nil
	}
	return m, mapped, unmapped
}
