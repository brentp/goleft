package indexcov

import (
	"log"
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

func getSizes(idx *bam.Index) ([][]int64, uint64, uint64) {
	var mapped, unmapped uint64
	refs := reflect.ValueOf(*idx).FieldByName("idx").FieldByName("Refs")
	ptr := unsafe.Pointer(refs.Pointer())

	ret := (*(*[1 << 28]oRefIndex)(ptr))[:refs.Len()]
	// save some memory.
	m := make([][]int64, len(ret))
	n_messages := 0
	for i, r := range ret {
		st, ok := idx.ReferenceStats(i)
		if ok {
			mapped += st.Mapped
			unmapped += st.Unmapped
		} else {
			if n_messages <= 10 {
				log.Printf("no reference stats found for %dth reference chromosome", i)
			}
			if n_messages == 10 {
				log.Printf("not reporting further chromosomes without stats. %d", i)
			}
			n_messages += 1
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
