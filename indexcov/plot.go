package indexcov

import (
	"encoding/json"
	"fmt"
	"html/template"
	"math/rand"
	"os"

	chartjs "github.com/brentp/go-chartjs"
	"github.com/brentp/go-chartjs/types"
	"github.com/gonum/matrix/mat64"
)

type vs struct {
	xs []float64
	ys []float64
}

func (v *vs) Xs() []float64 {
	return v.xs
}

func (v *vs) Ys() []float64 {
	return v.ys
}

func (v *vs) Rs() []float64 {
	return nil
}

// truncate depth values above this to cnMax
const cnMax = 2.5

func asValues(vals []float32, multiplier float64) chartjs.Values {

	// skip until we find non-zero.
	v := vs{xs: make([]float64, 0, len(vals)), ys: make([]float64, 0, len(vals))}
	seenNonZero := false
	for i, r := range vals {
		if r == 0 && !seenNonZero {
			continue
		}
		seenNonZero = true
		v.xs = append(v.xs, float64(i)*multiplier)
		if r > cnMax {
			r = cnMax
		}
		v.ys = append(v.ys, float64(r))
	}
	return &v
}

func randomColor(s int) *types.RGBA {
	rand.Seed(int64(s))
	return &types.RGBA{
		R: uint8(rand.Intn(256)),
		G: uint8(rand.Intn(256)),
		B: uint8(rand.Intn(256)),
		A: 240}
}

func plotDepths(depths [][]float32, samples []string, chrom string, prefix string) error {
	chart := chartjs.Chart{Label: chrom}
	xa, err := chart.AddXAxis(chartjs.Axis{Type: chartjs.Linear, Position: chartjs.Bottom, ScaleLabel: &chartjs.ScaleLabel{FontSize: 16, LabelString: "position on " + chrom, Display: chartjs.True}})
	if err != nil {
		return err
	}
	ya, err := chart.AddYAxis(chartjs.Axis{Type: chartjs.Linear, Position: chartjs.Left,
		Tick:       &chartjs.Tick{Min: 0, Max: 2.5},
		ScaleLabel: &chartjs.ScaleLabel{FontSize: 16, LabelString: "scaled coverage", Display: chartjs.True}})
	if err != nil {
		return err
	}

	for i, depth := range depths {
		xys := asValues(depth, 16384)
		c := randomColor(i)
		dataset := chartjs.Dataset{Data: xys, Label: samples[i], Fill: chartjs.False, PointRadius: 0, BorderWidth: 0.5,
			BorderColor: c, BackgroundColor: c, SteppedLine: chartjs.True, PointHitRadius: 6}
		dataset.XAxisID = xa
		dataset.YAxisID = ya
		chart.AddDataset(dataset)
	}
	chart.Options.Responsive = chartjs.False
	wtr, err := os.Create(fmt.Sprintf("%s-indexcov-depth-%s.html", prefix, chrom))
	if err != nil {
		return err
	}
	if err := chart.SaveHTML(wtr, map[string]interface{}{"width": 800, "height": 800}); err != nil {
		return err
	}
	return wtr.Close()

}

func plotCounters(counts []*counter, samples []string, prefix string) {
	c := &types.RGBA{R: 110, G: 250, B: 59, A: 240}
	chart := chartjs.Chart{}
	xa, err := chart.AddXAxis(chartjs.Axis{Type: chartjs.Linear, Position: chartjs.Bottom, ScaleLabel: &chartjs.ScaleLabel{FontSize: 16,
		LabelString: "total bins with depth < 0.15",
		Display:     chartjs.True}})
	if err != nil {
		panic(err)
	}

	ya, err := chart.AddYAxis(chartjs.Axis{Type: chartjs.Linear, Position: chartjs.Left, ScaleLabel: &chartjs.ScaleLabel{FontSize: 16,
		LabelString: "total bins with depth outside of (0.85, 1.15)",
		Display:     chartjs.True}})

	if err != nil {
		panic(err)
	}
	xys := &vs{xs: make([]float64, len(counts)), ys: make([]float64, len(counts))}
	for i, c := range counts {
		xys.xs[i] = float64(c.low)
		xys.ys[i] = float64(c.out)
	}
	dataset := chartjs.Dataset{Data: xys, Label: "samples", Fill: chartjs.False, PointHoverRadius: 6,
		PointRadius: 4, BorderWidth: 0, BorderColor: &types.RGBA{R: 150, G: 150, B: 150, A: 150},
		PointBackgroundColor: c, BackgroundColor: c, ShowLine: chartjs.False, PointHitRadius: 6}
	dataset.XAxisID = xa
	dataset.YAxisID = ya
	chart.AddDataset(dataset)
	chart.Options.Responsive = chartjs.False
	wtr, err := os.Create(fmt.Sprintf("%s-indexcov-bins.html", prefix))
	if err != nil {
		panic(err)
	}
	sjson, err := json.Marshal(samples)
	if err != nil {
		panic(err)
	}
	jsfunc := fmt.Sprintf(`
for (var i =0; i < charts.length; i++) {
    charts[i].options.tooltips.callbacks.footer = function(tts, data) {
        var names = %s
        var out = []
        tts.forEach(function(ti) {
            out.push(names[ti.index])
        })
        return out.join(",")
    }
}`, sjson)

	if err := chartjs.SaveCharts(wtr, map[string]interface{}{"custom": template.JS(jsfunc), "width": 800, "height": 800}, chart); err != nil {
		panic(err)
	}
	wtr.Close()

}

func plotPCA(mat *mat64.Dense, prefix string, samples []string, vars []float64) {

	var charts []chartjs.Chart
	c := &types.RGBA{R: 110, G: 250, B: 59, A: 240}

	for _, pc := range []int{2, 3} {

		c1 := chartjs.Chart{}
		xa, err := c1.AddXAxis(chartjs.Axis{Type: chartjs.Linear, Position: chartjs.Bottom, ScaleLabel: &chartjs.ScaleLabel{FontSize: 16,
			LabelString: fmt.Sprintf("PC1 (variance explained: %.2f%%)", 100*vars[0]),
			Display:     chartjs.True}})
		if err != nil {
			panic(err)
		}

		ya, err := c1.AddYAxis(chartjs.Axis{Type: chartjs.Linear, Position: chartjs.Left, ScaleLabel: &chartjs.ScaleLabel{FontSize: 16,
			LabelString: fmt.Sprintf("PC%d (variance explained: %.2f%%)", pc, 100*vars[pc-1]),
			Display:     chartjs.True}})

		if err != nil {
			panic(err)
		}
		xys := &vs{xs: mat64.Col(nil, 0, mat), ys: mat64.Col(nil, pc-1, mat)}
		dataset := chartjs.Dataset{Data: xys, Label: "samples", Fill: chartjs.False, PointHoverRadius: 6,
			PointRadius: 4,
			BorderWidth: 0, BorderColor: &types.RGBA{R: 150, G: 150, B: 150, A: 150}, PointBackgroundColor: c, BackgroundColor: c, ShowLine: chartjs.False, PointHitRadius: 6}
		dataset.XAxisID = xa
		dataset.YAxisID = ya
		c1.AddDataset(dataset)
		c1.Options.Responsive = chartjs.False
		charts = append(charts, c1)
	}
	wtr, err := os.Create(fmt.Sprintf("%s-indexcov-pca.html", prefix))
	if err != nil {
		panic(err)
	}
	sjson, err := json.Marshal(samples)
	if err != nil {
		panic(err)
	}
	jsfunc := fmt.Sprintf(`
for (var i =0; i < charts.length; i++) {
    charts[i].options.tooltips.callbacks.footer = function(tts, data) {
        var names = %s
        var out = []
        tts.forEach(function(ti) {
            out.push(names[ti.index])
        })
        return out.join(",")
    }
}`, sjson)

	if err := chartjs.SaveCharts(wtr, map[string]interface{}{"custom": template.JS(jsfunc), "width": 800, "height": 800}, charts...); err != nil {
		panic(err)
	}
	wtr.Close()

}

func plotROCs(rocs [][]float32, samples []string, chrom string) (chartjs.Chart, error) {

	chart := chartjs.Chart{}
	xa, err := chart.AddXAxis(chartjs.Axis{Type: chartjs.Linear, Position: chartjs.Bottom, ScaleLabel: &chartjs.ScaleLabel{FontSize: 16, LabelString: "scaled coverage for " + chrom, Display: chartjs.True}, Tick: &chartjs.Tick{Max: 1 / slotsMid}})
	if err != nil {
		return chart, err
	}
	yax := chartjs.Axis{Type: chartjs.Linear, Position: chartjs.Left,
		ScaleLabel: &chartjs.ScaleLabel{FontSize: 16, LabelString: "proportion of regions covered", Display: chartjs.True}}

	ya, err := chart.AddYAxis(yax)
	if err != nil {
		return chart, err
	}

	for i, roc := range rocs {
		xys := asValues(roc, 1/float64(slots)*1/slotsMid)
		c := randomColor(i)
		dataset := chartjs.Dataset{Data: xys, Label: samples[i], Fill: chartjs.False, PointRadius: 0.0, BorderWidth: 2, BorderColor: c, PointBackgroundColor: c, BackgroundColor: c, PointHitRadius: 8, PointHoverRadius: 3}
		dataset.XAxisID = xa
		dataset.YAxisID = ya
		chart.AddDataset(dataset)
	}
	chart.Options.Responsive = chartjs.False
	return chart, nil
}

func plotSex(sexes map[string][]float64, chroms []string, samples []string) (chartjs.Chart, string, error) {
	chart := chartjs.Chart{Label: "sex"}
	tmp := sexes["_inferred"]
	inferred := make([]int, len(tmp))
	cns := make(map[int]bool)
	for i, inf := range tmp {
		inferred[i] = int(inf)
		cns[int(inf)] = true
	}
	xa, err := chart.AddXAxis(chartjs.Axis{Type: chartjs.Linear, Position: chartjs.Bottom,
		ScaleLabel: &chartjs.ScaleLabel{FontSize: 16, LabelString: chroms[0] + " Copy Number",
			Display: chartjs.True}, Tick: &chartjs.Tick{Min: 0}})
	if err != nil {
		return chart, "", err
	}
	ya, err := chart.AddYAxis(chartjs.Axis{Type: chartjs.Linear, Position: chartjs.Left,
		ScaleLabel: &chartjs.ScaleLabel{FontSize: 16, LabelString: chroms[1] + " Copy Number",
			Display: chartjs.True}, Tick: &chartjs.Tick{Min: 0}})
	if err != nil {
		return chart, "", err
	}
	// chartjs separates into datasets so we need to track which samples in which datasets.
	jssamples := make([][]string, 0, 2)
	for cn := range cns {
		jssamples = append(jssamples, make([]string, 0))
		vals := &vs{xs: make([]float64, 0, len(inferred)),
			ys: make([]float64, 0, len(inferred))}
		for i, inf := range inferred {
			if inf != cn {
				continue
			}
			jssamples[len(jssamples)-1] = append(jssamples[len(jssamples)-1], samples[i])
			vals.xs = append(vals.xs, sexes[chroms[0]][i])
			vals.ys = append(vals.ys, sexes[chroms[1]][i])
		}
		c := randomColor(cn)
		dataset := chartjs.Dataset{Data: vals, Label: fmt.Sprintf("Inferred CN for %s: %d", chroms[0], cn), Fill: chartjs.False, PointRadius: 6, BorderWidth: 0,
			BorderColor: &types.RGBA{R: 150, G: 150, B: 150, A: 150}, PointBackgroundColor: c, BackgroundColor: c, ShowLine: chartjs.False, PointHitRadius: 6}
		dataset.XAxisID = xa
		dataset.YAxisID = ya
		chart.AddDataset(dataset)
	}
	sjson, err := json.Marshal(jssamples)
	if err != nil {
		panic(err)
	}
	jsfunc := fmt.Sprintf(`
	chart.options.tooltips.callbacks.footer = function(tts, data) {
		var names = %s
		var out = []
		tts.forEach(function(ti) {
			out.push(names[ti.datasetIndex][ti.index])
		})
		return out.join(",")
	}`, sjson)
	return chart, jsfunc, nil
}
