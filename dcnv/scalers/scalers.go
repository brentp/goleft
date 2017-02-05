// Package scalers holds the interface for scaling depths to standardized scores.
package scalers

import (
	"math"
	"sort"

	"github.com/gonum/matrix/mat64"
	"github.com/gonum/stat"
)

// Scaler allows transformation and back of the depths.
// As an example, see the `ZScore` struct. Usually, these
// will be 0-centered
type Scaler interface {
	// Scale Converts from AdjustedDepth to a scaled value
	Scale(*mat64.Dense)
	UnScale(*mat64.Dense)
}

var _ Scaler = &ZScore{}
var _ Scaler = &Log2{}

// ZScore implements the Scaler interface for StdScore (z-score)
type ZScore struct {
	means []float64
	sds   []float64
}

// UnScale converts back to depths.
func (z *ZScore) UnScale(a *mat64.Dense) {
	r, _ := a.Dims()
	for i := 0; i < r; i++ {
		row := a.RawRowView(i)
		for c, v := range row {
			row[c] = math.Max(0, v*z.sds[i]+z.means[i])
		}
	}
}

// Scale converts from depths to z-scores.
func (z *ZScore) Scale(a *mat64.Dense) {
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

func (rc *RowCentered) Scale(a *mat64.Dense) {
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

func (rc *RowCentered) UnScale(a *mat64.Dense) {
	r, _ := a.Dims()
	for i := 0; i < r; i++ {
		row := a.RawRowView(i)
		cnt := rc.centers[i]
		for j := range row {
			row[j] += cnt
		}
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
	RC *RowCentered
}

// Scale converts from depths to log2s
func (l *Log2) Scale(a *mat64.Dense) {
	r, _ := a.Dims()
	if l.RC == nil {
		l.RC = &RowCentered{Centerer: gmean}
	}
	for i := 0; i < r; i++ {
		row := a.RawRowView(i)
		for c, d := range row {
			if d > 0 {
				row[c] = math.Log2(d)
			} else {
				// use -6 because much lower numbers screw the SVD scaling.
				row[c] = -6
			}
		}
	}
	l.RC.Scale(a)
}

// UnScale converts from log2s to depths
func (l *Log2) UnScale(a *mat64.Dense) {
	r, _ := a.Dims()
	l.RC.UnScale(a)
	for i := 0; i < r; i++ {
		row := a.RawRowView(i)
		for c, d := range row {
			row[c] = math.Pow(2, d)
		}
	}
}
