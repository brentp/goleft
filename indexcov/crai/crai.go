package crai

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// note: for index cov, just need alnStart, alnSpan and sliceLen
// will need to scale sliceLen by 16384/alnSpan and then artifically partition
// into 16KB chunks?

// Slice holds the index information for a particular cram slice
type Slice struct {
	alnStart int64
	alnSpan  int64
	// Container start byte offset in the file
	containerStart int64
	// Slice start byte offset in the container data (‘blocks’)
	sliceStart int64
	sliceLen   int32
}

func (s Slice) Start() int64 {
	return s.alnStart
}

func (s Slice) SliceBytes() int32 {
	return s.sliceLen
}

func (s Slice) Span() int64 {
	return s.alnSpan
}

type Index struct {
	Slices [][]Slice
}

func ReadIndex(r io.Reader) (*Index, error) {
	b := bufio.NewReader(r)

	idx := &Index{Slices: make([][]Slice, 0, 2)}
	iline := 1

	for line, err := b.ReadString('\n'); err == nil; line, err = b.ReadString('\n') {
		parts := strings.Split(strings.TrimSpace(line), "\t")
		if len(parts) != 6 {
			return nil, fmt.Errorf("crai: expected 6 fields in index, got %d", len(parts))
		}

		si, err := strconv.Atoi(parts[0])
		if err != nil {
			return nil, fmt.Errorf("crai: unable to parse seqID (%s) at line %d", parts[0], line)
		}
		for i := len(idx.Slices); i <= si; i++ {
			idx.Slices = append(idx.Slices, make([]Slice, 0, 16))
		}

		sl := Slice{}
		if alnStart, err := strconv.Atoi(parts[1]); err != nil {
			return nil, fmt.Errorf("crai: unable to parse alignment start (%s) at line %d", parts[1], iline)
		} else {
			sl.alnStart = int64(alnStart)
		}

		if alnSpan, err := strconv.Atoi(parts[2]); err != nil {
			return nil, fmt.Errorf("crai: unable to parse alignment span (%s) at line %d", parts[2], iline)
		} else {
			sl.alnSpan = int64(alnSpan)
		}

		if containerStart, err := strconv.Atoi(parts[3]); err != nil {
			return nil, fmt.Errorf("crai: unable to parse alignment container start (%s) at line %d", parts[3], iline)
		} else {
			sl.containerStart = int64(containerStart)
		}

		if sliceStart, err := strconv.Atoi(parts[4]); err != nil {
			return nil, fmt.Errorf("crai: unable to parse alignment slice start (%s) at line %d", parts[4], iline)
		} else {
			sl.sliceStart = int64(sliceStart)
		}

		if sliceLen, err := strconv.Atoi(parts[5]); err != nil {
			return nil, fmt.Errorf("crai: unable to parse alignment slice length (%s) at line %d", parts[5], iline)
		} else {
			sl.sliceLen = int32(sliceLen)
		}
		idx.Slices[si] = append(idx.Slices[si], sl)

		iline++
	}
	return idx, nil
}
