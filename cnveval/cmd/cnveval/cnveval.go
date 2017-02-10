package main

import (
	"fmt"
	"io"
	"log"
	"os"
	"runtime/pprof"
	"strconv"
	"strings"

	arg "github.com/alexflint/go-arg"
	"github.com/brentp/goleft/cnveval"
	"github.com/brentp/xopen"
)

type cliargs struct {
	MinOverlap   float64 `arg:"-m,help:minimum overlap proportion for intervals to be considered as overlapping"`
	LimitSamples bool    `arg:"-s,help:only included sites in the truth-set that have matching samples in the test-set"`
	Truth        string  `arg:"positional,required,help:bed file of truth-set"`
	Test         string  `arg:"positional,required,help:bed file of test-set"`
}

func main() {

	cli := cliargs{MinOverlap: 0.4}
	if p := arg.MustParse(&cli); cli.MinOverlap > 1 || cli.MinOverlap <= 0 {
		p.Fail("minoverlap must be between 0 and 1")
	}
	tcnvs := parseTruth(cli.Test, false, nil)
	var samples map[string]bool
	if cli.LimitSamples {
		samples = make(map[string]bool)
		for _, t := range tcnvs {
			samples[t.Samples[0]] = true
		}
	}

	truths := parseTruth(cli.Truth, true, samples)
	cnvs := asCNV(tcnvs)
	f, err := os.Create("x.cpu.pprof")
	if err != nil {
		panic(err)
	}
	pprof.StartCPUProfile(f)
	defer pprof.StopCPUProfile()

	ch := cnveval.Evaluate(cnvs, truths, cli.MinOverlap)
	tabs := ch.Tabulate()
	for _, cl := range []cnveval.SC{cnveval.Small, cnveval.Medium, cnveval.Large, cnveval.Any} {
		t := tabs[cl.Order()]
		fmt.Printf("size-class: %-12s | %s\n", cl, t)
	}
}

func asCNV(ts []cnveval.Truth) []cnveval.CNV {
	cnvs := make([]cnveval.CNV, len(ts))
	for i, t := range ts {
		cnvs[i] = cnveval.CNV{Chrom: t.Chrom, Start: t.Start, End: t.End,
			CN: t.CN, Sample: t.Samples[0]}
	}
	return cnvs
}

func parseTruth(p string, cnt bool, samples map[string]bool) []cnveval.Truth {
	truths := make([]cnveval.Truth, 0, 1000)
	fh, err := xopen.Ropen(p)
	if err != nil {
		panic(err)
	}
	n := 0
	for {
		line, err := fh.ReadString('\n')
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		if line[0] == '#' {
			continue
		}
		t := parseLine(line)
		if samples != nil {
			found := false
			for _, s := range t.Samples {
				if _, ok := samples[s]; ok {
					found = true
					break
				}
			}
			if !found {
				continue
			}
		}
		truths = append(truths, parseLine(line))
		if cnt {
			n += len(truths[len(truths)-1].Samples)
		}
	}
	if cnt {
		log.Printf("%d total cnvs in truth-set (%d)", n, len(truths))
	} else {
		log.Printf("cnvs in test-set (%d)", len(truths))
	}
	return truths
}

func mustAtoi(s string) int {
	i, err := strconv.Atoi(s)
	if err != nil {
		panic(err)
	}
	return i
}

func parseLine(l string) cnveval.Truth {
	toks := strings.Split(strings.TrimRight(l, "\r\n"), "\t")
	if len(toks) < 5 {
		panic("expected five fields for CNVs")
	}
	return cnveval.Truth{Chrom: toks[0], Start: mustAtoi(toks[1]), End: mustAtoi(toks[2]),
		Samples: strings.Split(toks[4], ","), CN: mustAtoi(toks[3])}
}
