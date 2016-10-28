// Package mops implements the EM algorithm described in the cn.mops paper.
package mops

import "math"

// mean of all except highest and lowest values.
func mean32(a []float32) float32 {
	var sum float32
	for _, v := range a {
		sum += v
	}
	return sum / float32(len(a))
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
const delta = 0.001
const maxiter = 10

const eps = 0.001

func pmf(k int, pmean float64) float64 {
	// TODO: pre-calc exp(-pmean) if math.Exp show up in bench.
	return math.Pow(pmean, float64(k)) * math.Exp(-pmean) / math.Gamma(float64(k+1))
}

// top right of eqn 5.
func pdepth(depths []float32, cn int, lambda float32, out []float32) {
	// i indexes copy-number, k indexes sample
	beta := float64(cn) / 2 * float64(lambda)
	if cn == 0 {
		beta = 0.001 / 2 * float64(lambda)
	}
	//em := math.Exp(-beta)
	for k, d := range depths {
		out[k] = float32(pmf(int(d+0.5), beta))
	}
}

// returns a square alpha ik.
func estep(alpha []float32, depths []float32, lambda float32, aik [][]float32) [][]float32 {
	// i index copy-numbers
	//	k indexes sample
	// N is number of samples.
	N := len(depths)
	if len(aik) == 0 {
		// best-to-reuse as this can take a lot of memory.
		aik = make([][]float32, len(alpha))
		for i := range aik {
			aik[i] = make([]float32, N)
		}
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
		// aik[cn] is set by this call to avoid allocation by re-use
		pdepth(depths, cn, lambda, aik[cn])
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

type Mopped struct {
	aik   [][]float32
	alpha []float32
}

func (m *Mopped) Gain() float32 {
	// cn.mops eqn(8)
	var ig float32
	for i, ai := range m.aik {
		v := float64(i)
		if i == 0 {
			v = float64(eps)
		}
		ig += mean32(ai) * float32(abs(math.Log(float64(v/2))))
	}
	return ig
}

func Mops(depths []float32) *Mopped {
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

	// em iterations.
	for n := 0; abs32(lambda-nlambda) > 0.01 && n < 10; n++ {
		lambda = nlambda
		aik = estep(alpha, depths, lambda, aik)
		alpha, nlambda = mstep(depths, aik)
	}
	return &Mopped{aik: aik, alpha: alpha}
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
