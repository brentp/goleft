package emdepth

import "math"

var gammas []float64

func init() {
	gammas = make([]float64, 1000)
	for i := 0; i < 1000; i++ {
		gammas[i], _ = math.Lgamma(float64(i + 1))
	}
}

func pmf(k int, mu float64) float64 {
	sign := 1
	var lg float64
	if k < len(gammas) {
		lg = gammas[k]
	} else {
		lg, sign = math.Lgamma(float64(k + 1))
	}
	return math.Exp(float64(k)*math.Log(mu) - float64(sign)*lg - mu)
}

func cdf(mu float64, k int) float64 {
	var tot float64
	for ; k >= 0; k-- {
		tot += pmf(k, mu)
	}
	return tot
}
