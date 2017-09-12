package debiaser_test

import (
	"fmt"
	"math/rand"
	"reflect"
	"testing"

	"gonum.org/v1/gonum/mat"

	"github.com/brentp/goleft/dcnv/debiaser"
	"github.com/brentp/goleft/dcnv/scalers"
)

func fillMatrix(mat *mat.Dense) {
	r, c := imat.Dims()
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			imat.Set(i, j, 100+float64(i)*100+10*float64(j)+float64(j)/10+10*rand.Float64())
		}
	}
}

func printMatrix(imat *mat.Dense) {
	r, c := imat.Dims()
	for i := 0; i < r; i++ {
		fmt.Printf("[")
		row := imat.RawRowView(i)
		for j := 0; j < c; j++ {
			if j < c-1 {
				fmt.Printf("%.2f, ", row[j])
			} else {
				fmt.Printf("%.2f", row[j])
			}
		}
		fmt.Printf("]\n")
	}
}

func TestGeneralDebiasSort(t *testing.T) {

	g := debiaser.GeneralDebiaser{Vals: []float64{0, 1, 2, 3, 10, 9, 8, 6}}

	imat := mat.NewDense(len(g.Vals), 4, nil)
	cpy := mat.NewDense(len(g.Vals), 4, nil)
	fillMatrix(imat)
	cpy.Copy(imat)

	g.Sort(imat)
	/*
		for i := 0; i < r; i++ {
			fmt.Printf("%v\n", mat.RawRowView(i))
		}
	*/

	/*
		if !reflect.DeepEqual(g.Posns, []uint32{0, 1, 2, 3, 7, 6, 5, 4}) {
			t.Fatalf("got %s", g.Posns)
		}
	*/
	if !reflect.DeepEqual(g.Vals, []float64{0, 1, 2, 3, 6, 8, 9, 10}) {
		t.Fatalf("got %s", g.Vals)
	}
	if reflect.DeepEqual(cpy, mat) {
		t.Fatalf("expected different matrix after sort")
	}

	g.Unsort(imat)
	/*
		if !reflect.DeepEqual(g.Posns, []uint32{0, 1, 2, 3, 4, 5, 6, 7}) {
			t.Fatalf("got %s", g.Posns)
		}
	*/
	if !reflect.DeepEqual(g.Vals, []float64{0, 1, 2, 3, 10, 9, 8, 6}) {
		t.Fatalf("got %s", g.Vals)
	}
	/*
		fmt.Println("\nafter unsort")
		for i := 0; i < r; i++ {
			fmt.Printf("%v -- %v\n", mat.RawRowView(i), cpy.RawRowView(i))
		}
	*/
	if !reflect.DeepEqual(cpy, imat) {
		t.Fatalf("expected indentical matrix after unsort")
	}

}

func TestGeneralDebias(t *testing.T) {
	g := debiaser.GeneralDebiaser{Vals: []float64{0, 1, 2, 3, 10, 9, 8, 6, 5, 4, 3, 2, 1, 0.5},
		Window: 1,
		//Posns:  []uint32{0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14},
	}
	imat := mat.NewDense(len(g.Vals), 7, nil)
	fillMatrix(imat)
	var zscore scalers.ZScore
	zscore.Scale(imat)
	zscore.UnScale(imat)
	// check z-scaling
	//fmt.Println("after")
	//printMatrix(mat)
	zscore.Scale(imat)
	g.Sort(imat)
	g.Debias(imat)
	zscore.UnScale(imat)
	//fmt.Println("after")
	//printMatrix(mat)
	/*
		if reflect.DeepEqual(g.Posns, []uint32{0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14}) {
			t.Fatal("expected unsorted posns")
		}
	*/

	g.Unsort(imat)
	/*
		if !reflect.DeepEqual(g.Posns, []uint32{0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14}) {
			t.Fatal("expected sorted posns")
		}
	*/

}
