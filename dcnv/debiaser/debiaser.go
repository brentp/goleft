package debiaser

import (
	"github.com/gonum/floats"
	"github.com/gonum/matrix/mat64"
)

// Debiaser implements inplace removal of bias from a mat64 of (scaled) values.
type Debiaser interface {
	Debias(*mat64.Dense)
}

// Sorter provides method to sort and then unsort a mat64
type Sorter interface {
	Sort(*mat64.Dense)
	Unsort(*mat64.Dense)
}

// SortedDebiaser is useful when we need to: sort, then unbias based on that sort order, then unsort.
// An example of this is for GC bias.
type SortedDebiaser interface {
	Sorter
	Debiaser
}

type GeneralDebiaser struct {
	Vals   []float64
	Posns  []uint32
	Window int
	inds   []int
	tmp    *mat64.Dense
}

func (g *GeneralDebiaser) setTmp(r, c int) {
	if g.tmp == nil {
		g.tmp = mat64.NewDense(r, c, nil)
	} else {
		gr, gc := g.tmp.Dims()
		if gr != r || gc != c {
			g.tmp = mat64.NewDense(r, c, nil)
		}
	}
}

func (g *GeneralDebiaser) Sort(mat *mat64.Dense) {
	if g.inds == nil {
		g.inds = make([]int, len(g.Vals))
	}
	floats.Argsort(g.Vals, g.inds)
	r, c := mat.Dims()
	g.setTmp(r, c)
	posns := make([]uint32, len(g.Posns))
	copy(posns, g.Posns)

	for ai, bi := range g.inds {
		g.tmp.SetRow(ai, mat.RawRowView(bi))
		g.Posns[ai] = posns[bi]
	}
	// copy g.tmp into mat
	mat.Copy(g.tmp)
}

func (g *GeneralDebiaser) Unsort(mat *mat64.Dense) {
	if g.inds == nil {
		panic("unsort: must call sort first")
	}
	r, c := mat.Dims()
	g.setTmp(r, c)
	// copy mat into g.tmp
	g.tmp.Copy(mat)
	tmp := make([]float64, len(g.Vals))
	posns := make([]uint32, len(g.Posns))
	copy(posns, g.Posns)
	for ai, bi := range g.inds {
		mat.SetRow(ai, g.tmp.RawRowView(bi))
		tmp[ai] = g.Vals[bi]
		g.Posns[ai] = posns[bi]
	}
	g.Vals = tmp
}
