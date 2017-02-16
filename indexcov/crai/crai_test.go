package crai_test

import (
	"strings"
	"testing"

	"github.com/brentp/goleft/indexcov/crai"
)

func TestCrai(t *testing.T) {

	s := strings.NewReader(`0	1	2	3	4	5
0	10	20	30	40	50
`)
	cr, err := crai.ReadIndex(s)
	if err != nil {
		t.Fatal(err)
	}

	if len(cr.Slices) != 1 {
		t.Fatalf("expected 1 chromosome, got %d", len(cr.Slices))
	}

	sl := cr.Slices[0]

	if len(sl) != 2 {
		t.Fatalf("expected 2 slices, got %d", len(sl))
	}

	for _, t := range sl {
		_, _, _ = t.Start(), t.Span(), t.SliceBytes()
	}

}
