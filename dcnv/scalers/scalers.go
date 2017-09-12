// Package scalers holds the interface for scaling depths to standardized scores.
package scalers

import (
	"math"
	"sort"

	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/stat"
)

// Scaler allows transformation and back of the depths.
// As an example, see the `ZScore` struct. Usually, these
// will be 0-centered
type Scaler interface {
	// Scale Converts from AdjustedDepth to a scaled value
	Scale(*mat.Dense)
	UnScale(*mat.Dense)
}

var _ Scaler = &ZScore{}
var _ Scaler = &Log2{}

// ZScore implements the Scaler interface for StdScore (z-score)
type ZScore struct {
	means []float64
	sds   []float64
}

// UnScale converts back to depths.
func (z *ZScore) UnScale(a *mat.Dense) {
	r, _ := a.Dims()
	for i := 0; i < r; i++ {
		row := a.RawRowView(i)
		for c, v := range row {
			row[c] = math.Max(0, v*z.sds[i]+z.means[i])
		}
	}
}

// Scale converts from depths to z-scores.
func (z *ZScore) Scale(a *mat.Dense) {
	r, _ := a.Dims()
	z.means = make([]float64, r)
	z.sds = make([]float64, r)
	// convert to z-score
	for i := 0; i < r; i++ {
		row := a.RawRowView(i)
		m, sd := stat.MeanStdDev(row, nil)
		for c, d := range row {
			row[c] = (d - m) / sd
		}
		z.means[i] = m
		z.sds[i] = sd
	}
}

type RowCentered struct {
	Centerer func([]float64) float64
	centers  []float64
}

type ColCentered struct {
	Centerer func([]float64) float64
	centers  []float64
}

func (rc *RowCentered) Scale(a *mat.Dense) {
	r, _ := a.Dims()
	if rc.centers == nil {
		rc.centers = make([]float64, 0, r)
	}
	rc.centers = rc.centers[:0]
	for i := 0; i < r; i++ {
		row := a.RawRowView(i)
		rc.centers = append(rc.centers, rc.Centerer(row))
		for c := range row {
			row[c] -= rc.centers[i]
		}
	}

}

func (cc *ColCentered) Scale(a *mat.Dense) {
	r, c := a.Dims()
	if cc.centers == nil {
		cc.centers = make([]float64, 0, c)
	}
	cc.centers = cc.centers[:0]
	col := make([]float64, r)
	for i := 0; i < c; i++ {
		mat.Col(col, i, a)
		cc.centers = append(cc.centers, cc.Centerer(col))
		for c := range col {
			col[c] -= cc.centers[i]
		}
		a.SetCol(i, col)
	}

}

func (rc *RowCentered) UnScale(a *mat.Dense) {
	r, _ := a.Dims()
	for i := 0; i < r; i++ {
		row := a.RawRowView(i)
		cnt := rc.centers[i]
		for j := range row {
			row[j] += cnt
		}
	}
}

func (cc *ColCentered) UnScale(a *mat.Dense) {
	r, c := a.Dims()
	col := make([]float64, r)
	for i := 0; i < c; i++ {
		mat.Col(col, i, a)
		cnt := cc.centers[i]
		for j := range col {
			col[j] += cnt
		}
		a.SetCol(i, col)
	}
}

func gmean(vs []float64) float64 {
	os := make([]float64, len(vs))
	copy(os, vs)
	sort.Float64s(os)
	return os[len(os)/2]
	return stat.Mean(vs, nil)
}

// Log2 implements Scaler interface to perform log2 transformation on depths.
type Log2 struct {
	CC *ColCentered
}

// Scale converts from depths to log2s
func (l *Log2) Scale(a *mat.Dense) {
	r, _ := a.Dims()
	if l.CC == nil {
		l.CC = &ColCentered{Centerer: gmean}
	}
	for i := 0; i < r; i++ {
		row := a.RawRowView(i)
		for c, d := range row {
			row[c] = math.Log2(1 + d)
		}
	}
	l.CC.Scale(a)
}

// UnScale converts from log2s to depths
func (l *Log2) UnScale(a *mat.Dense) {
	r, _ := a.Dims()
	l.CC.UnScale(a)
	for i := 0; i < r; i++ {
		row := a.RawRowView(i)
		for c, d := range row {
			row[c] = math.Pow(2, d)
		}
	}
}
