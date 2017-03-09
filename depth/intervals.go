package depth

import (
	"bufio"
	"io"
	"log"

	"github.com/biogo/store/interval"
	"github.com/brentp/xopen"
)

// Integer-specific intervals
type irange struct {
	Start, End int
	UID        uintptr
}

func (i irange) Overlap(b interval.IntRange) bool {
	// Half-open interval indexing.
	return i.End > b.Start && i.Start < b.End
}
func (i irange) ID() uintptr              { return i.UID }
func (i irange) Range() interval.IntRange { return interval.IntRange{i.Start, i.End} }

// Overlaps checks for overlaps without pulling intervals from the tree.
func Overlaps(tree *interval.IntTree, start, end int) bool {
	if tree == nil {
		return false
	}

	q := irange{Start: start, End: end, UID: uintptr(tree.Len())}

	overlaps := false
	tree.DoMatching(func(iv interval.IntInterface) bool {
		overlaps = true
		return true
	}, q)
	return overlaps

}

// ReadTree takes a bed file and returns map of trees.
func ReadTree(p string) map[string]*interval.IntTree {
	if p == "" {
		return nil
	}
	r, err := xopen.Ropen(p)
	if err != nil {
		panic(err)
	}
	tree := make(map[string]*interval.IntTree, 10)
	br := bufio.NewReader(r)
	k := 0

	for {
		line, err := br.ReadBytes('\n')
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}

		chrom, start, end := chromStartEndFromLine(line)
		if _, ok := tree[chrom]; !ok {
			tree[chrom] = &interval.IntTree{}
		}
		tree[chrom].Insert(irange{start, end, uintptr(k)}, false)
		k++

	}
	log.Printf("read %d intervals into interval tree", k)
	return tree
}
