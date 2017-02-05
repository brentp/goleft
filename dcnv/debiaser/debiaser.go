package debiaser

import (
	"github.com/JaderDias/movingmedian"
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

// GeneralDebiaser is an implementation of a SortedDebiaser that can be used simply by setting attributes.int
// Vals and Posns should have equal length.
// Window is the size of the moving window for correction.
// Usage is to Call GeneralDebiaser.Sort() then Debias(), then Unsort(). Presumbly, Those calls will be flanked
// by a scaler call such as to scaler.ZScore.
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

// Sort sorts the rows in mat according the order in g.Vals.
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

// Unsort reverts the values to be position sorted.
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
		mat.SetRow(bi, g.tmp.RawRowView(ai))
		tmp[bi] = g.Vals[ai]
		g.Posns[bi] = posns[ai]
	}
	g.Vals = tmp
}

// Debias by subtracting moving median in each sample.
// It's assumed that g.Sort() has been called before this and that g.Unsort() will be called after.
// It's also assumed that the values in mat have been scaled, for example by a `scaler.ZScore`.
func (g *GeneralDebiaser) Debias(mat *mat64.Dense) {
	r, c := mat.Dims()
	col := make([]float64, r)
	for sampleI := 0; sampleI < c; sampleI++ {
		mat64.Col(col, sampleI, mat)

		mm := movingmedian.NewMovingMedian(g.Window)
		mid := (g.Window-1)/2 + 1
		for i := 0; i < mid; i++ {
			mm.Push(col[i])
		}

		var i int
		for i = 0; i < len(col)-mid; i++ {
			mm.Push(col[i+mid])
			col[i] -= mm.Median()
		}
		for ; i < len(col); i++ {
			col[i] -= mm.Median()
		}
		mat.SetCol(sampleI, col)
	}
}
