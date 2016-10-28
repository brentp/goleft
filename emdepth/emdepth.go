// Package emdepth uses a simplified EM algorithm to assign copy-number given a set of depths.
// it is based off the implementation in cn.mops.
// Like cn.mops, it works best with more samples and it assumes that most samples will
// have copy-number 2.
// emdepth consists of a single function EMDepth that iteratively assigns depths to copy-numbers, adjusts
// the center of each copy-number bin. and re-assigns...
// This package does no normalization and therefore expects incoming data to be normalized.
package emdepth

import (
	"sort"
	"sync"
)

// mean of all except highest and lowest values.
func mean32(a []float32) float32 {
	var sum float32
	for _, v := range a {
		sum += v
	}
	return sum / float32(len(a))
}

func mean64(a []float64) float64 {
	if len(a) == 0 {
		return 0
	}
	var sum float64
	for _, v := range a {
		sum += v
	}
	return sum / float64(len(a))
}

func abs(a float64) float64 {
	if a < 0 {
		return -a
	}
	return a
}

func abs32(a float32) float32 {
	if a < 0 {
		return -a
	}
	return a
}

const maxCN = 8
const maxiter = 10
const eps = 0.001

// used to test for convergence return the sum of abs differences
// and the largest absolute difference.
func summaxdiff(a, b []float64) (float64, float64) {
	var m, sum float64
	for i, av := range a {
		d := abs(av - b[i])
		sum += d
		if d > m {
			m = d
		}
	}
	return sum, m
}

func getBins(nSamples int) [][]float64 {
	binned := make([][]float64, maxCN)
	for i := 0; i < maxCN; i++ {
		if i == 2 {
			binned[i] = make([]float64, 0, nSamples)
		} else {
			binned[i] = make([]float64, 0)
		}
	}
	return binned
}

var binPool *sync.Pool

func init() {
	binPool = &sync.Pool{New: func() interface{} {
		return getBins(32)
	}}
}

// EMDepth returns a slice of integer copy-numbers (CNs) corresponding to the given normalized
// depths. It uses a simple EM to assign depths to copy-numbers bins with a preference for CN 2.
// And to adjust the mean depth of each bin after each iteration.
func EMDepth(depths []float32) []int {

	m := float64(mean32(depths))
	// lambda are the centers of each CN
	lambda := make([]float64, maxCN)
	// put each sample into the bin they belong to.
	// we re-use these since they can take a lot of mem.
	binned := binPool.Get().([][]float64)

	// keep lastCenters to check for convergence.
	lastCenters := make([]float64, len(lambda))

	lambda[0], lambda[2] = eps*m, m
	lastCenters[0] = lambda[0]

	// EXPECTATION:
	// initialize centers of each bin relative to copy-number 2.
	for i := 1; i < len(lambda); i++ {
		lambda[i] = lambda[2] * (float64(i) / 2)
	}

	// iterate at most maxiter times. stop early if the largest depth change
	// from 1 iter to the next is < 0.5 or if the total of all changes is < delta.
	sumd, maxd := float64(100), float64(100)
	for iter := 0; iter < maxiter && sumd > eps && maxd > 0.5; iter++ {
		for i := 1; i < len(lambda); i++ {
			lastCenters[i] = lambda[i]
			binned[i] = binned[i][:0]
		}
		binned[0] = binned[0][:0]
		//fmt.Println(getCNs(depths, centers), centers)
		// MAXIMIZATION
		// put samples in the bin they are closest 2.
		for _, df := range depths {
			d := float64(df)
			// most common case of copy-number 2.
			if d > lambda[1] && d < lambda[3] {
				if abs(d-lambda[2]) < abs(d-lambda[1]) && abs(d-lambda[2]) < abs(d-lambda[3]) {
					binned[2] = append(binned[2], d)
					continue
				}
			}
			idx := sort.SearchFloat64s(lambda, d)
			if idx == 0 {
				binned[0] = append(binned[0], d)
				continue
			}
			if idx == len(lambda) {
				binned[idx-1] = append(binned[idx-1], d)
				continue
			}
			// idx will always be the index of the value >= d so we check if idx or idx-1 is better.
			if abs(d-lambda[idx]) < abs(d-lambda[idx-1]) {
				binned[idx] = append(binned[idx], d)
			} else {
				binned[idx-1] = append(binned[idx-1], d)
			}
		}

		// EXPECTATION:
		// adjust copy-number 2 depth.
		lambda[2] = mean64(binned[2])
		if lambda[2] == 0 {
			n := float64(len(depths))
			// we exclude the top bin to avoid over-adjusting lambda[2] for extreme outliers.
			for i := 1; i < len(lambda)-1; i++ {
				bin := binned[i]
				pdepth := float64(len(bin)) / n
				if lambda[i] < eps {
					lambda[i] = eps
				}
				lambda[2] += mean64(bin) * (2 / float64(i)) * pdepth
			}
		}

		for i := 1; i < len(lambda); i++ {
			// adjust the expected depths of other copy-numbers based on that from CN2.
			lambda[i] = lambda[2] * (float64(i) / 2)
		}
		// make CN 2 more likely by expanding the range between CN1 and CN3.
		span := lambda[2] - lambda[1]
		lambda[1] -= (span / 2)
		lambda[3] += (span / 2)
		sumd, maxd = summaxdiff(lambda, lastCenters)
	}
	binPool.Put(binned)
	return getCNs(depths, lambda)
}

func distp5(v float64) float64 {
	if v < 0.5 {
		return 0.5 - v
	}
	return v - 0.5
}

// for copy-number 1 and 3, we check the CDF to make sure they are much better than CN2.
func adjustCNs(cns []int, depths []float32, centers []float64) []int {

	for i, cn := range cns {
		if cn == 1 || cn == 3 {
			dk := int(0.5 + depths[i])
			o, o2 := pmf(dk, centers[cn]), pmf(dk, centers[2])
			if o*0.95 < o2 {
				cns[i] = 2
			}
		}
	}
	return cns
}

func getCNs(depths []float32, centers []float64) []int {
	cns := make([]int, len(depths))
	for i, df := range depths {
		d := float64(df)
		idx := sort.SearchFloat64s(centers, d)
		if idx == 0 {
			cns[i] = 0
			continue
		}
		if idx == len(centers) {
			cns[i] = len(centers)
			continue
		}
		// idx will always be the index of the value >= d so we check if idx or idx-1 is better.
		if abs(d-centers[idx]) < abs(d-centers[idx-1]) {
			cns[i] = idx
		} else {
			cns[i] = idx - 1
		}

	}
	cns = adjustCNs(cns, depths, centers)
	return cns
}
