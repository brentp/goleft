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
	MinOverlap float64 `arg:"-m,help:minimum overlap proportion for intervals to be considered as overlapping"`
	Truth      string  `arg:"positional,required,help:bed file of truth-set"`
	Test       string  `arg:"positional,required,help:bed file of test-set"`
}

func main() {

	cli := cliargs{MinOverlap: 0.4}
	if p := arg.MustParse(&cli); cli.MinOverlap > 1 || cli.MinOverlap <= 0 {
		p.Fail("minoverlap must be between 0 and 1")
	}

	truths := parseTruth(cli.Truth, true)
	tcnvs := parseTruth(cli.Test, false)
	cnvs := asCNV(tcnvs)
	f, err := os.Create("x.cpu.pprof")
	if err != nil {
		panic(err)
	}
	pprof.StartCPUProfile(f)
	defer pprof.StopCPUProfile()

	ch := cnveval.Evaluate(cnvs, truths, cli.MinOverlap)
	fmt.Printf("precision: %.4f, recall: %.4f\n", ch.Precision(cnveval.Any), ch.Recall(cnveval.Any))

}

func asCNV(ts []cnveval.Truth) []cnveval.CNV {
	cnvs := make([]cnveval.CNV, len(ts))
	for i, t := range ts {
		cnvs[i] = cnveval.CNV{Chrom: t.Chrom, Start: t.Start, End: t.End,
			CN: t.CN, Sample: t.Samples[0]}
	}
	return cnvs
}

func parseTruth(p string, cnt bool) []cnveval.Truth {
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
		truths = append(truths, parseLine(line))
		if cnt {
			n += len(truths[len(truths)-1].Samples)
		}
	}
	if cnt {
		log.Printf("%d total cnvs (%d)", n, len(truths))
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
