package indexcov

import (
	"fmt"
	"os"

	chartjs "github.com/brentp/go-chartjs"
	"github.com/brentp/go-chartjs/types"
)

type vs struct {
	xs   []float64
	rocs []float64
}

func (v *vs) Xs() []float64 {
	return v.xs
}

func (v *vs) Ys() []float64 {
	return v.rocs
}

func (v *vs) Rs() []float64 {
	return nil
}

func asValues(roc []float32) chartjs.Values {

	v := vs{xs: make([]float64, len(roc)), rocs: make([]float64, len(roc))}
	for i, r := range roc {
		v.xs[i] = float64(i)
		v.rocs[i] = float64(r)
	}
	return &v

}

func plotROCs(rocs [][]float32, samples []string, prefix string) error {

	set1 := []*types.RGBA{
		&types.RGBA{166, 206, 227, 240},
		&types.RGBA{31, 120, 180, 240},
		&types.RGBA{178, 223, 138, 240},
		&types.RGBA{51, 160, 44, 240},
		&types.RGBA{251, 154, 153, 240},
		&types.RGBA{227, 26, 28, 240},
		&types.RGBA{253, 191, 111, 240},
		&types.RGBA{255, 127, 0, 240},
		&types.RGBA{202, 178, 214, 240},
		&types.RGBA{106, 61, 154, 240},
		&types.RGBA{255, 255, 153, 240},
		&types.RGBA{177, 89, 40, 240},
	}

	chart := chartjs.Chart{Label: "ROC"}
	xa, err := chart.AddXAxis(chartjs.Axis{Type: chartjs.Linear, Position: chartjs.Bottom, ScaleLabel: &chartjs.ScaleLabel{FontSize: 16, LabelString: "scaled coverage", Display: chartjs.True}})
	if err != nil {
		return err
	}
	ya, err := chart.AddYAxis(chartjs.Axis{Type: chartjs.Linear, Position: chartjs.Left,
		ScaleLabel: &chartjs.ScaleLabel{FontSize: 16, LabelString: "proportion of regions covered", Display: chartjs.True}})

	for i, roc := range rocs {
		xys := asValues(roc)
		c := set1[i%len(set1)]
		dataset := chartjs.Dataset{Data: xys, Label: samples[i], Fill: chartjs.False, PointRadius: 0.01, BorderWidth: 2, BorderColor: c, PointBackgroundColor: c, BackgroundColor: c}
		dataset.XAxisID = xa
		dataset.YAxisID = ya
		chart.AddDataset(dataset)
	}
	chart.Options.Responsive = chartjs.False
	wtr, err := os.Create(fmt.Sprintf("%s-indexcov-roc.html", prefix))
	if err != nil {
		return err
	}
	if err := chart.SaveHTML(wtr, map[string]interface{}{"height": 800, "width": 800}); err != nil {
		return err
	}
	return wtr.Close()

}
