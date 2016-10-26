// Package emdepth uses a simplified EM algorithm to assign copy-number given a set of depths.
// it is based off the implementation in cn.mops.
// Like cn.mops, it works best with more samples and it assumes that most samples will
// have copy-number 2.
// emdepth consists of a single function EMDepth that iteratively assigns depths to copy-numbers, adjusts
// the center of each copy-number bin. and re-assigns...
// This package does no normalization and therefore expects incoming data to be normalized.
package emdepth

import (
	"math"
	"sort"

	stats "github.com/r0fls/gostats"
)

// mean of all except highest and lowest values.
func mean32(a []float32) float32 {
	var sum float32
	min := float32(math.MaxFloat32)
	max := float32(0)
	for _, v := range a {
		sum += v
		if v < min {
			min = v
		}
		if v > max {
			max = v
		}
	}
	sum -= max
	sum -= min
	return sum / float32(len(a)-2)
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

const maxCN = 8
const delta = 0.1
const maxiter = 20

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

// EMDepth returns a slice of integer copy-numbers (CNs) corresponding to the given normalized
// depths. It uses a simple EM to assign depths to copy-numbers bins with a preference for CN 2.
// And to adjust the mean depth of each bin after each iteration.
func EMDepth(depths []float32) []int {

	m := float64(mean32(depths))
	centers := make([]float64, maxCN)
	// put each sample into the bin they belong to.
	binned := make([][]float64, maxCN)
	for i := 0; i < len(centers); i++ {
		if i == 2 {
			binned[i] = make([]float64, 0, len(depths))
		} else {
			binned[i] = make([]float64, 0)
		}
	}

	// keep lastCenters to check for convergence.
	lastCenters := make([]float64, len(centers))

	centers[0] = 0.005 * m
	centers[2] = m
	lastCenters[0] = centers[0]
	// EXPECTATION:
	// initialize centers of each bin relative to copy-number 2.
	for i := 1; i < len(centers); i++ {
		centers[i] = centers[2] * (float64(i) / 2)
	}
	// iterate at most maxiter times. stop early if the largest depth change
	// from 1 iter to the next is < 0.5 or if the total of all changes is < delta.
	sumd, maxd := float64(100), float64(100)
	for iter := 0; iter < maxiter && sumd > delta && maxd > 0.5; iter++ {
		for i := 1; i < len(centers); i++ {
			lastCenters[i] = centers[i]
		}
		//fmt.Println(getCNs(depths, centers), centers)
		// MAXIMIZATION
		// put samples in the bin they are closest 2.
		for _, df := range depths {
			d := float64(df)
			// most common case of copy-number 2.
			if d > centers[1] && d < centers[3] {
				if abs(d-centers[2]) < abs(d-centers[1]) && abs(d-centers[2]) < abs(d-centers[3]) {
					binned[2] = append(binned[2], d)
					continue
				}
			}
			idx := sort.SearchFloat64s(centers, d)
			if idx == 0 {
				binned[0] = append(binned[0], d)
				continue
			}
			if idx == len(centers) {
				binned[idx-1] = append(binned[idx-1], d)
				continue
			}
			// idx will always be the index of the value >= d so we check if idx or idx-1 is better.
			if abs(d-centers[idx]) < abs(d-centers[idx-1]) {
				binned[idx] = append(binned[idx], d)
			} else {
				binned[idx-1] = append(binned[idx-1], d)
			}
		}

		// EXPECTATION:
		// adjust copy-number 2 depth.
		centers[2] = mean64(binned[2])
		centers[0] = 0.005 * centers[2]
		for i := 1; i < len(centers); i++ {
			// adjust the expected depths of other copy-numbers based on that from CN2.
			centers[i] = centers[2] * (float64(i) / 2)
		}
		// make CN 2 more likely by expanding the range between CN1 and CN3.
		centers[1] *= 1 / 1.3
		centers[3] *= 1.3
		sumd, maxd = summaxdiff(centers, lastCenters)

	}
	return getCNs(depths, centers)
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
			p := stats.Poisson(float64(centers[cn])).Cdf(int(0.5 + depths[i]))
			p2 := stats.Poisson(float64(centers[2])).Cdf(int(0.5 + depths[i]))
			d := distp5(p)
			d2 := distp5(p2)
			// TODO: make this a parameter.
			if d2-0.05 < d {
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
	return adjustCNs(cns, depths, centers)
}
