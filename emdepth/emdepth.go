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
	"sync"
)

// mean of all except highest and lowest values.
func median32(a []float32) float64 {
	// make a copy so we dont modify original order
	b := make([]float64, len(a))
	for i, v := range a {
		b[i] = float64(v)
	}
	sort.Float64s(b)
	if len(b)%2 == 1 {
		return b[len(b)/2]
	}
	return (b[len(b)/2] + b[len(b)/2+1]) / 2
}

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

type Position struct {
	Start uint32
	End   uint32
}

// EMDepth returns a slice of integer copy-numbers (CNs) corresponding to the given normalized
// depths. It uses a simple EM to assign depths to copy-numbers bins with a preference for CN 2.
// And to adjust the mean depth of each bin after each iteration.
// I is a unique identifier for the depths that can be use to assiciate the returned struct with
// position info.
func EMDepth(depths []float32, p Position) *EMD {

	m := median32(depths)
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
			for i := 1; i < len(lambda); i++ {
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
	return &EMD{Lambda: lambda, Depths: depths, Position: p}
}

// EMD holds the posterior maximum depth for each copy-number
type EMD struct {
	Lambda   []float64
	Depths   []float32
	Position Position
}

// Same returns
// 1) slice of indexes for samples that have the same (non 2 CN),
// 2) a slice indicating the sample-set that has changed state in o in a way that
// ends a CNV, either by going back to state2 or by switching from a state > 2
// less <2.
// 3)a float indicating the proportion of samples that in the same copy-number
// state in both.
func (e *EMD) Same(o *EMD) (non2 []int, changed []int, pct float64) {
	ofc := o.Log2FC()
	var nSame float64
	non2 = make([]int, 0, 2)
	changed = make([]int, 0, 1)
	for i, ee := range e.Log2FC() {
		oo := ofc[i]
		if ee > -0.99 && ee < 0.59 && oo > -0.99 && oo < 0.59 {
			nSame++
			continue
		}
		// same direction.
		if (oo >= 0.59 && ee >= 0.59) || (oo <= -0.99 && ee <= -0.99) {
			non2 = append(non2, i)
			nSame++
		} else {
			changed = append(changed, i)
		}
	}
	return non2, changed, nSame / float64(len(e.Depths))
}

// Log2FC is the fold-change relative to copy-number 2.
func (e *EMD) Log2FC() []float64 {
	m := make([]float64, 0, len(e.Depths))
	for _, d := range e.Depths {
		m = append(m, math.Log2(float64(d)/e.Lambda[2]))
	}
	return m
}

// CN finds the posterior maximum CN for each sample.
func (e *EMD) CN() []int {
	cns := make([]int, len(e.Depths))
	for i, df := range e.Depths {
		cns[i] = e.Type(df)
	}
	return cns
}

// Type returns the copy-number for the given depth
func (e *EMD) Type(d float32) int {
	df := float64(d)
	idx := sort.SearchFloat64s(e.Lambda, df)
	var cn int
	if idx == 0 {
		cn = 0
	} else if idx == len(e.Lambda) {
		cn = len(e.Lambda)
	} else {
		// idx will always be the index of the value >= d so we check if idx or idx-1 is better.
		if abs(df-e.Lambda[idx]) < abs(df-e.Lambda[idx-1]) {
			cn = idx
		} else {
			cn = idx - 1
		}
	}
	return e.adjustCN(cn, df)
}

// adjustCN for copy-number 1 and 3. Use poisson PMF to make sure they
// are better than CN2. If not, set to CN2.
func (e *EMD) adjustCN(cn int, depth float64) int {

	// TODO: do this for all cn +/- 1.
	if cn == 1 || cn == 3 {
		dk := int(0.5 + depth)
		o, o2 := pmf(dk, e.Lambda[cn]), pmf(dk, e.Lambda[2])
		if o*0.95 < o2 {
			cn = 2
		}
	}
	return cn
}

const ProportionSame = 0.6

// Cache is a way to track depth states.eps
// As new items are added to the cache, the value from EMD.Same
// is compared.
// Each time an EMD is added, a slice of Regions that have ended is returned.
// The cache is cleared if the proportion of samples that changed state is > ProportionSame
type Cache struct {
	last *EMD
	// for each sample (key), where is it non-CN2?
	cnvs map[int][]*EMD
}

// CNV holds the information for a CNV for a single sample.
type CNV struct {
	SampleI  int
	Depth    []float32
	Position []Position
	Log2FC   []float32
	CN       []int
	PSize    int
}

func (c *Cache) Add(e *EMD) []*CNV {
	// TODO: here we check e.Position and eject anything that doesnt have and End with 3 * len(position).
	if c.last == nil {
		if c.cnvs == nil {
			c.cnvs = make(map[int][]*EMD, 200)
		}
		c.last = e
	}
	ret := c.Clear(&e.Position)
	noncn2, _, _ := c.last.Same(e)
	for _, si := range noncn2 {
		c.cnvs[si] = append(c.cnvs[si], e)
	}
	c.last = e
	return ret

}

func (c *Cache) Clear(p *Position) []*CNV {
	if p == nil {
		p = &Position{End: c.last.Position.End + 10000, Start: c.last.Position.Start + 10000}
	}
	cnvs := make([]*CNV, 0, 5)
	keys := make([]int, len(c.cnvs))
	L := p.End - p.Start
	for si, emd := range c.cnvs {
		// not a big enough gap yet.
		if p.Start-emd[len(emd)-1].Position.End < 3*L {
			continue
		}
		cnvs = append(cnvs, makecnvs(emd, si))
		cnvs[len(cnvs)-1].PSize = len(c.cnvs)
		keys = append(keys, si)
	}
	for _, key := range keys {
		delete(c.cnvs, key)
	}
	return cnvs
}

// merge individal aberrant depth calls into CNVs.
func makecnvs(es []*EMD, sampleI int) *CNV {
	var cnv *CNV
	// in here, we know we have adjacent calls from  same sample.
	for i, es := range es {
		fc := es.Log2FC()[sampleI]
		if fc > -0.5 && fc < 0.3 {
			continue
		}
		cn := es.Type(es.Depths[sampleI])
		if i == 0 {
			cnv = &CNV{SampleI: sampleI, CN: []int{cn}, Depth: []float32{es.Depths[sampleI]},
				Position: []Position{es.Position}, Log2FC: []float32{float32(fc)}}
			continue
		}
		cnv.CN = append(cnv.CN, cn)
		cnv.Depth = append(cnv.Depth, es.Depths[sampleI])
		cnv.Position = append(cnv.Position, es.Position)
		cnv.Log2FC = append(cnv.Log2FC, float32(fc))
	}
	return cnv
}

func weightedmean(a Position, b Position, af float32, bf float32) float32 {
	la := float32(a.End - a.Start)
	lb := float32(b.End - b.Start)
	return (af*float32(la) + bf*float32(lb)) / (la + lb)
}
