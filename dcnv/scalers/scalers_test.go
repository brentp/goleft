package scalers_test

import (
	"log"
	"math"
	"testing"

	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/stat"

	"github.com/brentp/goleft/dcnv/scalers"
)

const eps = 0.001

func TestZScoreRoundTrip(t *testing.T) {

	zsc := &scalers.ZScore{}
	imat := mat.NewDense(10, 10, nil)
	r, c := imat.Dims()
	if r != 10 || c != 10 {
		t.Fatal("unexpected size")
	}
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			imat.Set(i, j, float64((i+1)*(j+1)))
		}
	}

	zsc.Scale(imat)
	for i := 0; i < r; i++ {
		if math.Abs(stat.Mean(mat.RawRowView(i), nil)) > eps {
			t.Fatalf("expected 0, got %f", stat.Mean(mat.RawRowView(i), nil))
		}
	}
	zsc.UnScale(mat)
	for i := 0; i < r; i++ {
		row := mat.RawRowView(i)
		for j := 0; j < c; j++ {
			if math.Abs(row[j]-float64((i+1)*(j+1))) > eps {
				log.Fatalf("expected: %f, got %f", float64((i+1)*(j+1)), row[j])
			}
		}
	}

}

func TestLog2RoundTrip(t *testing.T) {

	zsc := &scalers.Log2{}
	imat := mat.NewDense(10, 10, nil)
	r, c := imat.Dims()
	if r != 10 || c != 10 {
		t.Fatal("unexpected size")
	}
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			imat.Set(i, j, float64((i+1)*(j+1)))
		}
	}

	zsc.Scale(imat)
	for i := 0; i < r; i++ {
		if floats.Min(imat.RawRowView(i)) < 0 {
			t.Fatalf("log2:expected >0, got %f", stat.Mean(imat.RawRowView(i), nil))
		}
	}
	zsc.UnScale(imat)
	for i := 0; i < r; i++ {
		row := imat.RawRowView(i)
		for j := 0; j < c; j++ {
			if math.Abs(row[j]-float64((i+1)*(j+1))) > eps {
				log.Fatalf("expected: %f, got %f", float64((i+1)*(j+1)), row[j])
			}
		}
	}

}
