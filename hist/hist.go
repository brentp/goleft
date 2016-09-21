// Package hist accepts data on Stdin and makes a histogram.
package hist

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"regexp"
	"sort"
	"strconv"

	arg "github.com/alexflint/go-arg"
	"github.com/gonum/plot"
	"github.com/gonum/plot/plotter"
	"github.com/gonum/plot/plotutil"
	"github.com/gonum/plot/vg"
)

type dargs struct {
	Col   int    `arg:"-c,required,help:1-based column number to use for histogram"`
	Group int    `arg:"-g,help:optional 1-based column number to use grouping"`
	Sep   string `arg:"-s,help:optional sep for columns default is '\t'"`
	Bins  int    `arg:"-b,help:optional number of bins."`
	Path  string `arg:"-p,help:optional path to save plot."`
}

func pcheck(e error) {
	if e != nil {
		panic(e)
	}
}

func MustAtof(s string) float64 {
	v, err := strconv.ParseFloat(s, 64)
	pcheck(err)
	return v
}

// Main is run from the dispatcher
func Main() {

	args := dargs{Sep: "\t", Bins: 20, Path: "hist.png"}
	arg.MustParse(&args)
	args.Group--
	args.Col--
	run(args)
}

func read(args dargs) map[string][]float64 {
	sep := regexp.MustCompile(args.Sep)
	scanner := bufio.NewScanner(os.Stdin)
	scanner.Buffer(make([]byte, 0, 16384), 5e9)

	grouped := make(map[string][]float64)

	nErr := 0
	for scanner.Scan() {
		line := scanner.Text()
		toks := sep.Split(line, -1)
		v, err := strconv.ParseFloat(toks[args.Col], 64)
		if err != nil {
			if nErr < 5 {
				log.Println(err)
			}
			nErr++
			continue
		}
		var g string
		if args.Group > -1 {
			g = toks[args.Group]
		} else {
			g = "default"
		}
		if _, ok := grouped[g]; !ok {
			grouped[g] = make([]float64, 0, 256)
		}
		grouped[g] = append(grouped[g], v)

	}
	if err := scanner.Err(); err != nil {
		fmt.Fprintln(os.Stderr, "reading standard input:", err)
	}
	return grouped
}

func mapkeys(m map[string][]float64) []string {
	var ks []string
	for k := range m {
		ks = append(ks, k)
	}
	sort.Strings(ks)
	return ks
}

func run(args dargs) {
	grouped := read(args)
	keys := mapkeys(grouped)

	p, err := plot.New()
	if err != nil {
		panic(err)
	}
	p.Y.Label.Text = "Count"

	w := float64(30/len(grouped)) * float64(20) / float64(args.Bins)
	var bars []plot.Plotter

	for i, k := range keys {
		tmp, err := plotter.NewHist(plotter.Values(grouped[k]), args.Bins)
		pcheck(err)
		vals := make([]float64, len(tmp.Bins))
		for i, b := range tmp.Bins {
			vals[i] = b.Weight
		}

		bar, err := plotter.NewBarChart(plotter.Values(vals), vg.Points(w+0.01))
		pcheck(err)

		bar.LineStyle.Width = vg.Length(0.1)
		bar.Color = plotutil.Color(i)

		bar.Offset = vg.Points(float64(i) * w)
		p.Legend.Add(k, bar)
		//p.Add(bar)
		bars = append(bars, bar)
	}
	p.Add(bars...)
	log.Println(p.X.Min, p.X.Max)

	p.Legend.Top = true
	if err := p.Save(10*vg.Inch, 3*vg.Inch, args.Path); err != nil {
		panic(err)
	}
}
