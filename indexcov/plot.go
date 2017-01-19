package indexcov

import (
	"encoding/json"
	"fmt"
	"html/template"
	"image/color"
	"math"
	"math/rand"
	"os"

	chartjs "github.com/brentp/go-chartjs"
	"github.com/brentp/go-chartjs/types"
	"github.com/gonum/matrix/mat64"
	"github.com/gonum/plot"
	"github.com/gonum/plot/plotter"
	"github.com/gonum/plot/vg"
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

func (v *vs) Len() int {
	return len(v.xs)
}

// make it meet gonum/plot plotter.XYer

func (v *vs) XY(i int) (x, y float64) {
	return v.xs[i], v.ys[i]
}

func (v *vs) Sample(nth int) *vs {
	o := &vs{xs: make([]float64, 0, 10+len(v.xs)/nth),
		ys: make([]float64, 0, 10+len(v.xs)/nth)}
	for i, x := range v.xs {
		if i%nth == 0 {
			o.xs = append(o.xs, x)
			o.ys = append(o.ys, v.ys[i])
		}
	}
	return o
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

func plotDepths(depths [][]float32, samples []string, chrom string, base string, writeHTML bool) error {
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
	chart.Options.Tooltip = &chartjs.Tooltip{Mode: "nearest"}
	if writeHTML {
		wtr, err := os.Create(fmt.Sprintf("%s-depth-%s.html", base, chrom))
		if err != nil {
			return err
		}
		link := template.HTML(`<a href="index.html">back to index</a>`)
		if err := chart.SaveHTML(wtr, map[string]interface{}{"width": 850, "height": 550, "customHTML": link}); err != nil {
			return err
		}
		if err := wtr.Close(); err != nil {
			return err
		}
	}
	asPng(fmt.Sprintf("%s-depth-%s.png", base, chrom), chart, 4, 3)

	return nil
}

func plotBins(counts []*counter, samples []string) (chartjs.Chart, string) {
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
	chart.Options.Tooltip = &chartjs.Tooltip{Mode: "nearest"}
	chart.Options.Legend = &chartjs.Legend{Display: chartjs.False}
	sjson, err := json.Marshal(samples)
	if err != nil {
		panic(err)
	}
	jsfunc := fmt.Sprintf(`
    bin_chart.options.tooltips.callbacks.title = function(tts, data) {
        var names = %s
        var out = []
        tts.forEach(function(ti) {
            out.push(names[ti.index])
        })
        return out.join(",")
    }`, sjson)
	return chart, jsfunc

}

func plotPCA(mat *mat64.Dense, samples []string, vars []float64) ([]chartjs.Chart, string) {

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
		c1.Options.Legend = &chartjs.Legend{Display: chartjs.False}
		c1.Options.Tooltip = &chartjs.Tooltip{Mode: "nearest"}
		charts = append(charts, c1)
	}
	sjson, err := json.Marshal(samples)
	if err != nil {
		panic(err)
	}
	jsfunc := fmt.Sprintf(`
	chart.options.hover.mode = 'index';
    chart.options.tooltips.callbacks.title = function(tts, data) {
        var names = %s
        var out = []
        tts.forEach(function(ti) {
            out.push(names[ti.index])
        })
        return out.join(",")
    }`, sjson)

	return charts, jsfunc
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
	chart.Options.Tooltip = &chartjs.Tooltip{Mode: "nearest"}
	return chart, nil
}

func plotSex(sexes map[string][]float64, chroms []string, samples []string) (*chartjs.Chart, string, error) {
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
		return nil, "", err
	}
	ya, err := chart.AddYAxis(chartjs.Axis{Type: chartjs.Linear, Position: chartjs.Left,
		ScaleLabel: &chartjs.ScaleLabel{FontSize: 16, LabelString: chroms[1] + " Copy Number",
			Display: chartjs.True}, Tick: &chartjs.Tick{Min: 0}})
	if err != nil {
		return nil, "", err
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
			BorderColor: &types.RGBA{R: 90, G: 90, B: 90, A: 150}, PointBackgroundColor: c, BackgroundColor: c, ShowLine: chartjs.False, PointHitRadius: 6}
		dataset.XAxisID = xa
		dataset.YAxisID = ya
		chart.AddDataset(dataset)
	}
	sjson, err := json.Marshal(jssamples)
	if err != nil {
		panic(err)
	}
	jsfunc := fmt.Sprintf(`
	chart.options.hover.mode = 'index'
	chart.options.tooltips.callbacks.title = function(tts, data) {
		var names = %s
		var out = []
		tts.forEach(function(ti) {
			out.push(names[ti.datasetIndex][ti.index])
		})
		return out.join(",")
	}`, sjson)
	chart.Options.Responsive = chartjs.False
	chart.Options.Tooltip = &chartjs.Tooltip{Mode: "nearest"}
	return &chart, jsfunc, nil
}

func asPng(path string, chart chartjs.Chart, wInches float64, hInches float64) {
	p, err := plot.New()
	if err != nil {
		panic(err)
	}
	p.X.Label.Text = chart.Options.Scales.XAxes[0].ScaleLabel.LabelString
	p.Y.Label.Text = chart.Options.Scales.YAxes[0].ScaleLabel.LabelString
	for i := range chart.Data.Datasets {
		ds := chart.Data.Datasets[len(chart.Data.Datasets)-i-1]
		data := ds.Data
		// gonum plotting is a significant portion of the runtime so we sample datasets.
		if data.(*vs).Len() > 2000 {
			data = data.(*vs).Sample(10)
		} else if data.(*vs).Len() > 1000 {
			data = data.(*vs).Sample(5)
		} else {
			// no data for some chromosome.
			bad := false
			for _, d := range data.Ys() {
				if math.IsNaN(d) {
					bad = true
					break
				}
			}
			if bad {
				continue
			}
		}

		l, err := plotter.NewLine(data.(*vs))
		if err != nil {
			panic(err)
		}
		c := color.RGBA(*ds.BorderColor)
		c.A = 255
		l.LineStyle.Width = vg.Points(0.8)

		l.Color = c
		p.Add(l)
	}
	if err := p.Save(vg.Length(wInches)*vg.Inch, vg.Length(hInches)*vg.Inch, path); err != nil {
		panic(err)
	}
}
