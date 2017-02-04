// Package scalers holds the interface for scaling depths to standardized scores.
package scalers

import (
	"math"

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

// Log2 implements Scaler interface to perform log2 transformation on depths.
type Log2 struct {
}

// Scale converts from depths to log2s
func (l *Log2) Scale(a *mat64.Dense) {
	r, _ := a.Dims()
	for i := 0; i < r; i++ {
		row := a.RawRowView(i)
		for c, d := range row {
			row[c] = math.Log2(d)
		}
	}
}

// UnScale converts from log2s to depths
func (l *Log2) UnScale(a *mat64.Dense) {
	r, _ := a.Dims()
	for i := 0; i < r; i++ {
		row := a.RawRowView(i)
		for c, d := range row {
			row[c] = math.Pow(2, d)
		}
	}
}
