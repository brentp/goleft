package debiaser

import (
	"fmt"
	"log"
	"math"
	"sort"

	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/mat"

	"github.com/JaderDias/movingmedian"
)

// Debiaser implements inplace removal of bias from a mat64 of (scaled) values.
type Debiaser interface {
	Debias(*mat.Dense)
}

// Sorter provides method to sort and then unsort a mat64
type Sorter interface {
	Sort(*mat.Dense)
	Unsort(*mat.Dense)
}

// SortedDebiaser is useful when we need to: sort, then unbias based on that sort order, then unsort.
// An example of this is for GC bias.
type SortedDebiaser interface {
	Sorter
	Debiaser
}

// GeneralDebiaser is an implementation of a SortedDebiaser that can be used simply by setting attributes.int
// Window is the size of the moving window for correction.
// Usage is to Call GeneralDebiaser.Sort() then Debias(), then Unsort(). Presumbly, Those calls will be flanked
// by a scaler call such as to scaler.ZScore.
type GeneralDebiaser struct {
	Vals   []float64
	Window int
	inds   []int
	tmp    *mat.Dense
}

func (g *GeneralDebiaser) setTmp(r, c int) {
	if g.tmp == nil {
		g.tmp = mat.NewDense(r, c, nil)
	} else {
		gr, gc := g.tmp.Dims()
		if gr != r || gc != c {
			g.tmp = mat.NewDense(r, c, nil)
		}
	}
}

// Sort sorts the rows in mat according the order in g.Vals.
func (g *GeneralDebiaser) Sort(mat *mat.Dense) {
	if g.inds == nil {
		g.inds = make([]int, len(g.Vals))
	}
	floats.Argsort(g.Vals, g.inds)
	r, c := mat.Dims()
	g.setTmp(r, c)

	changed := false
	for ai, bi := range g.inds {
		if ai != bi {
			changed = true
		}
		g.tmp.SetRow(ai, mat.RawRowView(bi))
	}
	if !changed {
		log.Println("WARNING: no change after sorting. This usually means .Vals is unset or same as previous run")
	}
	// copy g.tmp into mat
	mat.Copy(g.tmp)
}

// Unsort reverts the values to be position sorted.
func (g *GeneralDebiaser) Unsort(mat *mat.Dense) {
	if g.inds == nil {
		panic("unsort: must call sort first")
	}
	r, c := mat.Dims()
	g.setTmp(r, c)
	// copy mat into g.tmp
	g.tmp.Copy(mat)
	tmp := make([]float64, len(g.Vals))
	for ai, bi := range g.inds {
		mat.SetRow(bi, g.tmp.RawRowView(ai))
		tmp[bi] = g.Vals[ai]
	}
	g.Vals = tmp
}

// Debias by subtracting moving median in each sample.
// It's assumed that g.Sort() has been called before this and that g.Unsort() will be called after.
// It's also assumed that the values in mat have been scaled, for example by a `scaler.ZScore`.
func (g *GeneralDebiaser) Debias(imat *mat.Dense) {
	r, c := imat.Dims()
	col := make([]float64, r)
	for sampleI := 0; sampleI < c; sampleI++ {
		mat.Col(col, sampleI, imat)

		mm := movingmedian.NewMovingMedian(g.Window)
		mid := (g.Window-1)/2 + 1
		for i := 0; i < mid; i++ {
			mm.Push(col[i])
		}
		for i := 0; i < mid; i++ {
			col[i] /= math.Max(mm.Median(), 1)
		}

		var i int
		for i = mid; i < len(col)-mid; i++ {
			mm.Push(col[i+mid])
			col[i] /= math.Max(mm.Median(), 1)
		}
		for ; i < len(col); i++ {
			col[i] /= math.Max(mm.Median(), 1)
		}
		imat.SetCol(sampleI, col)
	}
}

type ChunkDebiaser struct {
	GeneralDebiaser
	// ScoreWindow defines the range of Vals used per window.
	// E.g. if this is 0.1 then all values from 0.25-0.35 will be normalized to the median of
	// Depths occuring in that range.
	ScoreWindow float64
}

func (cd *ChunkDebiaser) Debias(imat *mat.Dense) {
	if cd.ScoreWindow == 0 {
		panic("must set ChunkDebiaser.ScoreWindow")
	}
	r, c := imat.Dims()
	col := make([]float64, r)

	slices := make([]int, 1, 100)
	v0 := cd.Vals[0]
	for i := 0; i < len(cd.Vals); i++ {
		if cd.Vals[i]-v0 > cd.ScoreWindow {
			v0 = cd.Vals[i]
			slices = append(slices, i)
		}
	}
	slices = append(slices, len(cd.Vals))
	dpSubset := make([]float64, 0, len(cd.Vals))

	for sampleI := 0; sampleI < c; sampleI++ {
		mat.Col(col, sampleI, imat)
		for i := 1; i < len(slices); i++ {
			si, ei := slices[i-1], slices[i]
			dpSubset = dpSubset[:(ei - si)]
			copy(dpSubset, col[si:ei])
			sort.Float64s(dpSubset)
			var k int
			for ; k < len(dpSubset) && dpSubset[k] == 0; k++ {
			}
			median := dpSubset[(ei-si-k)/2]

			if median > 0 {
				for j := si; j < ei; j++ {
					col[j] /= median
				}
			}
		}
		imat.SetCol(sampleI, col)
	}
}

type SVD struct {
	MinVariancePct float64
}

func (isvd *SVD) Debias(imat *mat.Dense) {
	var svd mat.SVD
	if ok := svd.Factorize(imat, mat.SVDThin); !ok {
		panic("error with SVD")
	}

	// get svd and zero out first n components.
	s, u, v := extractSVD(&svd)
	sum := floats.Sum(s)
	str := "variance:"

	var n int
	for n = 0; n < 15 && 100*s[n]/sum > isvd.MinVariancePct; n++ {
		str += fmt.Sprintf(" %.2f", 100*s[n]/sum)
	}
	log.Println(str)

	sigma := mat.NewDense(len(s), len(s), nil)
	// leave the first n as 0 and set the rest.
	for i := n; i < len(s); i++ {
		sigma.Set(i, i, s[i])
	}
	imat.Product(u, sigma, v.T())
}

func extractSVD(svd *mat.SVD) (s []float64, u, v *mat.Dense) {
	var um *mat.Dense
	var vm *mat.Dense
	svd.UTo(um)
	svd.VTo(vm)
	s = svd.Values(nil)
	return s, um, vm
}
