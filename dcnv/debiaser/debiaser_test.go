package debiaser_test

import (
	"fmt"
	"math/rand"
	"reflect"
	"testing"

	"github.com/brentp/goleft/dcnv/debiaser"
	"github.com/brentp/goleft/dcnv/scalers"
	"github.com/gonum/matrix/mat64"
)

func fillMatrix(mat *mat64.Dense) {
	r, c := mat.Dims()
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			mat.Set(i, j, 100+float64(i)*100+10*float64(j)+float64(j)/10+10*rand.Float64())
		}
	}
}

func printMatrix(mat *mat64.Dense) {
	r, c := mat.Dims()
	for i := 0; i < r; i++ {
		fmt.Printf("[")
		row := mat.RawRowView(i)
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

	mat := mat64.NewDense(len(g.Vals), 4, nil)
	cpy := mat64.NewDense(len(g.Vals), 4, nil)
	fillMatrix(mat)
	cpy.Copy(mat)

	g.Sort(mat)
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

	g.Unsort(mat)
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
	if !reflect.DeepEqual(cpy, mat) {
		t.Fatalf("expected indentical matrix after unsort")
	}

}

func TestGeneralDebias(t *testing.T) {
	g := debiaser.GeneralDebiaser{Vals: []float64{0, 1, 2, 3, 10, 9, 8, 6, 5, 4, 3, 2, 1, 0.5},
		Window: 1,
		//Posns:  []uint32{0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14},
	}
	mat := mat64.NewDense(len(g.Vals), 7, nil)
	fillMatrix(mat)
	var zscore scalers.ZScore
	zscore.Scale(mat)
	zscore.UnScale(mat)
	// check z-scaling
	//fmt.Println("after")
	//printMatrix(mat)
	zscore.Scale(mat)
	g.Sort(mat)
	g.Debias(mat)
	zscore.UnScale(mat)
	//fmt.Println("after")
	//printMatrix(mat)
	/*
		if reflect.DeepEqual(g.Posns, []uint32{0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14}) {
			t.Fatal("expected unsorted posns")
		}
	*/

	g.Unsort(mat)
	/*
		if !reflect.DeepEqual(g.Posns, []uint32{0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14}) {
			t.Fatal("expected sorted posns")
		}
	*/

}
