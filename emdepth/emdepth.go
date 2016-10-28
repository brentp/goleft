// Package emdepth uses a simplified EM algorithm to assign copy-number given a set of depths.
// it is based off the implementation in cn.mops.
// Like cn.mops, it works best with more samples and it assumes that most samples will
// have copy-number 2.
// emdepth consists of a single function EMDepth that iteratively assigns depths to copy-numbers, adjusts
// the center of each copy-number bin. and re-assigns...
// This package does no normalization and therefore expects incoming data to be normalized.
package emdepth

import (
	"fmt"
	"math"
	"sort"

	stats "github.com/r0fls/gostats"
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

const maxCN = 5
const delta = 0.001
const maxiter = 10

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

const eps = 0.001

func pmf(k int, pmean float64) float64 {
	// TODO: pre-calc exp(-pmean) if math.Exp show up in bench.
	return math.Pow(pmean, float64(k)) * math.Exp(-pmean) / math.Gamma(float64(k+1))
}

// top right of eqn 5.
func pdepth(depths []float32, cn int, lambda float32) []float32 {
	// i indexes copy-number, k indexes sample
	beta := float64(cn) / 2 * float64(lambda)
	if cn == 0 {
		beta = 0.001 / 2 * float64(lambda)
	}
	//em := math.Exp(-beta)
	ps := make([]float32, len(depths))
	for k, d := range depths {
		ps[k] = float32(pmf(int(d+0.5), beta))
	}
	return ps
}

func sum(a []float32) float32 {
	var s float32
	for _, v := range a {
		s += v
	}
	return s
}

// returns a square alpha ik.
func estep(alpha []float32, depths []float32, lambda float32, aik [][]float32) [][]float32 {
	// i index copy-numbers
	//	k indexes sample
	// N is number of samples.
	if len(aik) == 0 {
		// best-to-reuse as this can take a lot of memory.
		aik = make([][]float32, len(alpha))
	}
	//expm := math.Exp(float64(-lambda))

	denom := make([]float32, len(depths))
	n := len(alpha)
	// calculte prob for cn 2.
	// eqn(1), eqn(5) denom.
	for k, d := range depths {
		for i := 0; i < n; i++ {
			denom[k] += alpha[i] * float32(pmf(int(d+0.5), float64(float32(i)/2*lambda)))
		}
	}

	// calculate prob for other cns.
	for cn := range alpha {
		// TODO: send aik[cn] to pdeth to avoid allocation.
		aik[cn] = pdepth(depths, cn, lambda)
		//fmt.Println(lambda, cn, aik[cn])
		// eqn 5 from cn.mops.
		for k, ad := range aik[cn] {
			aik[cn][k] = alpha[cn] * ad / denom[k]
		}
	}
	return aik
}

func mstep(depths []float32, adepths [][]float32) (alpha []float32, lambda float32) {
	n := len(adepths)
	alpha = make([]float32, n)
	N := float32(len(depths))
	// eqn 6, 7
	G := float32(11)
	var lambdaDenom float32
	ys := float32(len(alpha)) + G

	alphaDenom := 1 + 1/N*float32(ys-float32(n))

	for cn, ai := range adepths {
		amean := mean32(ai)
		yi := float32(1)
		if cn == 2 {
			yi = 1 + G
		}
		top := amean + 1/N*(yi-1)
		alpha[cn] = top / alphaDenom

		// eqn 7
		if cn == 0 {
			lambdaDenom += amean * eps / 2
		} else {
			lambdaDenom += amean * float32(cn) / 2
		}

	}
	return alpha, mean32(depths) / lambdaDenom
}

func Mops(depths []float32) []int {
	// alpha[i] is percentage of samples with CNi
	// lambda is mean read-count for CN2
	// x is depth
	// p(x|i) is likely that read-count x is from CNi == 1/x! * e^-(i/2 * lambda) * i/2*lambda
	alpha := make([]float32, maxCN)
	for i := 0; i < len(alpha); i++ {
		alpha[i] = eps
	}
	alpha[2] = 1.0 - 5*eps*float32(len(alpha)-1)

	nlambda := mean32(depths)
	lambda := float32(math.MaxFloat32) / 10.0
	var aik [][]float32

	for n := 0; abs32(lambda-nlambda) > 0.01 && n < 10; n++ {
		lambda = nlambda
		aik = estep(alpha, depths, lambda, aik)
		alpha, nlambda = mstep(depths, aik)
	}
	for _, row := range aik {
		fmt.Println(row)
	}
	return []int{}
}

// EMDepth returns a slice of integer copy-numbers (CNs) corresponding to the given normalized
// depths. It uses a simple EM to assign depths to copy-numbers bins with a preference for CN 2.
// And to adjust the mean depth of each bin after each iteration.
func EMDepth(depths []float32) []int {

	m := float64(mean32(depths))
	// lambda are the centers of each CN
	lambda := make([]float64, maxCN)
	// put each sample into the bin they belong to.
	binned := make([][]float64, maxCN)
	for i := 0; i < len(lambda); i++ {
		if i == 2 {
			binned[i] = make([]float64, 0, len(depths))
		} else {
			binned[i] = make([]float64, 0)
		}
	}

	// keep lastCenters to check for convergence.
	lastCenters := make([]float64, len(lambda))

	lambda[0] = eps * m
	lambda[2] = m
	lastCenters[0] = lambda[0]
	// EXPECTATION:
	// initialize centers of each bin relative to copy-number 2.
	for i := 1; i < len(lambda); i++ {
		lambda[i] = lambda[2] * (float64(i) / 2)
	}
	// iterate at most maxiter times. stop early if the largest depth change
	// from 1 iter to the next is < 0.5 or if the total of all changes is < delta.
	sumd, maxd := float64(100), float64(100)
	for iter := 0; iter < maxiter && sumd > delta && maxd > 0.5; iter++ {
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
			for i := 1; i < len(lambda); i++ {
				bin := binned[i]
				pdepth := float64(len(bin)) / n
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
	return getCNs(depths, lambda)
}

func CNINI(cns []int, depths []float32, centers []float64) float32 {
	var in float32
	for i, cni := range cns {
		if cni == 0 {
			continue
		}
		v := float32(i)
		if v == 0 {
			v = eps
		}
		in += float32(cni) / float32(abs(math.Log(float64(v/2))))
	}
	return in
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
	cns = adjustCNs(cns, depths, centers)
	return cns
}
