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
	Scale()
	UnScale()
	// Values should return the existing array when possible; otherwise, it must create one of size r, c.
	Values(r, c int) *mat64.Dense
}

// ZScore implements the Scaler interface for StdScore (z-score)
type ZScore struct {
	means []float64
	sds   []float64
	mat   *mat64.Dense
}

// Values returns existing matrix or a new one of the given size.
func (z *ZScore) Values(r, c int) *mat64.Dense {
	if z.mat == nil {
		z.mat = mat64.NewDense(r, c, nil)
	}
	return z.mat
}

// UnScale converts back to depths.
func (z *ZScore) UnScale() {
	a := z.mat
	r, _ := a.Dims()
	for i := 0; i < r; i++ {
		row := a.RawRowView(i)
		for c, v := range row {
			// TODO: make sure this checks out.
			row[c] = math.Max(0, v*z.sds[i]+z.means[i])
		}
	}
}

// Scale converts from depths to z-scores.
func (z *ZScore) Scale() {
	a := z.mat
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
