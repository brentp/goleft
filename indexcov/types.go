package indexcov

import "github.com/biogo/hts/bgzf"

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
