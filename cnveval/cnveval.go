// Package cnveval provides a way to evalute CNVs based on a truth-set.
package cnveval

import (
	"fmt"
	"math"
	"sort"
)

// CNV indicates the region, sample and copy-number of a region.
type CNV struct {
	Chrom  string
	Start  int
	End    int
	Sample string
	CN     int
}

// Truth indicates the region, samples and copy-number for a truth-set.
type Truth struct {
	Chrom   string
	Start   int
	End     int
	Samples []string
	CN      int
	used    map[string]bool
}

func (t Truth) String() string {
	return fmt.Sprintf("%s:%d-%d[%d]", t.Chrom, t.Start, t.End, t.CN)
}

func (c CNV) String() string {
	return fmt.Sprintf("%s:%d-%d[%d] in %s", c.Chrom, c.Start, c.End, c.CN, c.Sample)
}

func (t Truth) equals(o Truth) bool {
	return t.Start == o.Start && t.End == o.End
}

// SC is a size-class
type SC int

var (
	// Any is used to avoid any size-class designation when reporting results.
	Any    SC
	Small  SC = 20000
	Medium SC = 100000
	Large  SC = math.MaxInt64
)

func (s SC) Order() int {
	switch s {
	case Small:
		return 0
	case Medium:
		return 1
	case Large:
		return 2
	case Any:
		return 3
	default:
		panic(fmt.Sprintf("unknown size class: %d", s))
	}
}

func (s SC) String() string {
	switch s {
	case Any:
		return "all"
	case Small:
		return fmt.Sprintf("0-%d", Small)
	case Medium:
		return fmt.Sprintf("%d-%d", Small, Medium)
	case Large:
		return fmt.Sprintf(">=%d", Medium)
	default:
		panic(fmt.Sprintf("unknown size class: %d", s))
	}
}

// TS indicates a sample and SizeClass
type TS struct {
	SizeClass SC
	Sample    string
}

type Stat struct {
	FP, FN, TP, TN int
}

func (s Stat) String() string {
	return fmt.Sprintf("precision: %.4f (%-04d / (%-04d + %-04d)) recall: %.4f (%-04d / (%-04d + %-04d))", s.Precision(), s.TP, s.TP, s.FP,
		s.Recall(), s.TP, s.TP, s.FN)
}

func (s Stat) Recall() float64 {
	return float64(s.TP) / float64(s.TP+s.FN)
}
func (s Stat) Precision() float64 {
	return float64(s.TP) / float64(s.TP+s.FP)
}

type CohortStats map[TS]*Stat

func (c CohortStats) TP(class SC) int {
	var TP int
	for cls, st := range c {
		if class != Any && cls.SizeClass != class {
			continue
		}
		TP += st.TP
	}
	return TP
}

func (c CohortStats) Tabulate() [4]Stat {
	var stats [4]Stat
	for cls, st := range c {
		stats[cls.SizeClass.Order()].TP += st.TP
		stats[cls.SizeClass.Order()].FP += st.FP
		stats[cls.SizeClass.Order()].TN += st.TN
		stats[cls.SizeClass.Order()].FN += st.FN
	}
	for _, cl := range []SC{Small, Medium, Large} {
		stats[Any.Order()].TP += stats[cl.Order()].TP
		stats[Any.Order()].FP += stats[cl.Order()].FP
		stats[Any.Order()].TN += stats[cl.Order()].TN
		stats[Any.Order()].FN += stats[cl.Order()].FN
	}
	return stats
}

func (c CohortStats) Precision(class SC) float64 {
	var TP, FP float64
	for cls, st := range c {
		if class != Any && cls.SizeClass != class {
			continue
		}
		TP += float64(st.TP)
		FP += float64(st.FP)
	}
	//log.Println("TP, FP:", TP, FP)
	return TP / (TP + FP)
}

func (c CohortStats) Recall(class SC) float64 {
	var TP, FN, FP, TN float64
	for cls, st := range c {
		if class != Any && cls.SizeClass != class {
			continue
		}
		TP += float64(st.TP)
		FN += float64(st.FN)
		FP += float64(st.FP)
		TN += float64(st.TN)
	}
	//log.Println("TP, FN, FP, TN:", TP, FN, FP, TN)
	return TP / (TP + FN)
}

func Evaluate(cnvs []CNV, truths []Truth, po float64) CohortStats {
	stat := make(map[TS]*Stat, 200)
	allSamples := make(map[string]bool)
	for _, t := range truths {
		for _, s := range t.Samples {
			allSamples[s] = true
		}
	}
	for _, c := range cnvs {
		allSamples[c.Sample] = true
	}
	sampleList := fromMap(allSamples)
	truthBySample := make(map[string][]Truth)
	truthWithoutSample := make(map[string][]Truth)
	for _, t := range truths {
		for _, s := range t.Samples {
			truthBySample[s] = append(truthBySample[s], t)
		}
		for _, s := range sampleList {
			if notin(s, t.Samples) {
				truthWithoutSample[s] = append(truthWithoutSample[s], t)
			}
		}
	}

	cnvBySample := make(map[string][]CNV)
	for _, c := range cnvs {
		cnvBySample[c.Sample] = append(cnvBySample[c.Sample], c)
	}

	for _, sample := range sampleList {
		truths, _ := truthBySample[sample]
		sort.Slice(truths, func(i, j int) bool {
			return truths[i].Chrom < truths[j].Chrom || (truths[i].Chrom == truths[j].Chrom && truths[i].Start < truths[j].Start)
		})
		cnvs, _ := cnvBySample[sample]
		sort.Slice(cnvs, func(i, j int) bool {
			return cnvs[i].Chrom < cnvs[j].Chrom || (cnvs[i].Chrom == cnvs[j].Chrom && cnvs[i].Start < cnvs[j].Start)
		})

		updatePositive(stat, truths, cnvs, po)
		others, _ := truthWithoutSample[sample]
		sort.Slice(others, func(i, j int) bool {
			return others[i].Chrom < others[j].Chrom || (others[i].Chrom == others[j].Chrom && others[i].Start < others[j].Start)
		})
		updateFP(stat, others, cnvs, truths, po)
	}

	return stat
}

func notin(a string, bs []string) bool {
	for _, b := range bs {
		if a == b {
			return false
		}
	}
	return true
}

func fromMap(samples map[string]bool) []string {
	m := make([]string, 0, len(samples))
	for k := range samples {
		m = append(m, k)
	}
	return m
}

func updateFP(stat map[TS]*Stat, others []Truth, cnvs []CNV, truths []Truth, po float64) {
	if len(cnvs) == 0 || len(others) == 0 {
		return
	}

	var i int
	for _, o := range others {
		found := false
		ts := TS{Sample: cnvs[0].Sample, SizeClass: sizeClass(o)}
		val := stat[ts]
		if val == nil {
			val = &Stat{}
			stat[ts] = val
		}
		// don't need to reset i because others is sorted.
		for ; i < len(cnvs) && (cnvs[i].Chrom < o.Chrom || (cnvs[i].Chrom == o.Chrom && cnvs[i].End < o.Start)); i++ {
		}
		if i > 0 {
			i--
		}
		for _, cnv := range cnvs[i:] {
			if cnv.Chrom > o.Chrom || (o.Chrom == cnv.Chrom && cnv.Start > o.End) {
				break
			}
			if poverlap(cnv, o) >= po && sameCN(cnv.CN, o.CN) {
				// we have a putative FP, but need to check if there's something else in the truth something
				// that also fits this criteria to avoid incorrectly calling an FP.
				tfound := false
				for _, t := range truths {
					if t.Chrom != cnv.Chrom {
						continue
					}
					if poverlap(cnv, t) >= po && sameCN(cnv.CN, t.CN) {
						tfound = true
						break
					}

				}
				if !tfound {
					val.FP++
					found = true
					break
				}
			}
			if found {
				break
			}

		}
		if !found {
			val.TN++
		}
	}
}

// given a set of truths and cnvs from a sample, update the FPs.
// truths and cnvs are sorted by chrom, then by start
func updatePositive(stat map[TS]*Stat, truths []Truth, cnvs []CNV, po float64) {
	if len(cnvs) == 0 {
		return
	}
	var i int
	for _, t := range truths {
		ts := TS{Sample: cnvs[0].Sample, SizeClass: sizeClass(t)}
		found := false
		val := stat[ts]
		if val == nil {
			val = &Stat{}
			stat[ts] = val
		}
		// since truths are sorted as well, we don't need to reset i.
		for ; i < len(cnvs) && (cnvs[i].Chrom < t.Chrom || (cnvs[i].Chrom == t.Chrom && cnvs[i].End < t.Start)); i++ {
		}
		if i > 0 {
			i--
		}
		for _, cnv := range cnvs[i:] {
			if cnv.Chrom > t.Chrom || (t.Chrom == cnv.Chrom && cnv.Start > t.End) {
				break
			}
			if poverlap(cnv, t) >= po && sameCN(cnv.CN, t.CN) {
				// used map make sure we don't double count cnvs from the sameCN
				// sample that are subsets of the full truth interval.
				if _, ok := t.used[cnv.Sample]; !ok {
					val.TP++
					found = true
					if t.used == nil {
						t.used = make(map[string]bool)
					}
					t.used[cnv.Sample] = true
				}
			}
		}
		if !found {
			val.FN++
		}
	}
}

func sizeClass(t Truth) SC {
	l := t.End - t.Start
	if l < int(Small) {
		return Small
	}
	if l < int(Medium) {
		return Medium
	}
	return Large
}

func sameCN(a, b int) bool {
	if a > 2 {
		a = 3
	}
	if b > 2 {
		b = 3
	}
	return math.Abs(float64(a)-float64(b)) < 1
}

func imin(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func imax(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func poverlap(a CNV, b Truth) float64 {
	if a.Chrom != b.Chrom {
		return 0
	}
	total := math.Min(float64(a.End-a.Start), float64(b.End-b.Start))
	ovl := imin(a.End, b.End) - imax(a.Start, b.Start)
	if ovl < 0 {
		return 0
	}
	return float64(ovl) / total
}
