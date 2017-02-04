package debiaser_test

import (
	"reflect"
	"testing"

	"github.com/brentp/goleft/dcnv/debiaser"
	"github.com/gonum/matrix/mat64"
)

func TestGeneralDebias(t *testing.T) {

	g := debiaser.GeneralDebiaser{Vals: []float64{0, 1, 2, 3, 10, 9, 8, 6},
		Posns: []uint32{0, 1, 2, 3, 4, 5, 6, 7},
	}

	mat := mat64.NewDense(100, len(g.Posns), nil)
	r, c := mat.Dims()
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			mat.Set(i, j, float64(i)*100+float64(j)/10)
		}
	}

	g.Sort(mat)

	if !reflect.DeepEqual(g.Posns, []uint32{0, 1, 2, 3, 7, 6, 5, 4}) {
		t.Fatalf("got %s", g.Posns)
	}
	if !reflect.DeepEqual(g.Vals, []float64{0, 1, 2, 3, 6, 8, 9, 10}) {
		t.Fatalf("got %s", g.Vals)
	}

	g.Unsort(mat)
	if !reflect.DeepEqual(g.Posns, []uint32{0, 1, 2, 3, 4, 5, 6, 7}) {
		t.Fatalf("got %s", g.Posns)
	}
	if !reflect.DeepEqual(g.Vals, []float64{0, 1, 2, 3, 10, 9, 8, 6}) {
		t.Fatalf("got %s", g.Vals)
	}

}
