package main

import (
	"bufio"
	"fmt"
	"os"
	"os/exec"
	"strconv"
	"strings"

	arg "github.com/alexflint/go-arg"
)

type dargs struct {
	Q          int      `arg:"-Q,help:mapping quality cutoff"`
	Chrom      string   `arg:"required,-c,help:optional chromosome to limit analysis"`
	MinCov     int      `arg:"help:minimum depth considered callable"`
	MaxCov     int      `arg:"help:maximum depth considered callable"`
	MinSamples float64  `arg:"help:proportion of samples with mincov coverage for a region to be reported"`
	Bams       []string `arg:"positional,required,help:bams for which to calculate depth"`
	minSamples int      `arg:"-"`
}

const command = "samtools depth -Q %d -d %d -r '%s' %s"

func main() {
	args := dargs{Q: 10, MinCov: 7, MaxCov: 1000, MinSamples: 0.5}
	if p := arg.MustParse(&args); len(args.Bams) == 0 {
		p.Fail("specify > 1 bam")
	}

	args.minSamples = int(0.5 + args.MinSamples*float64(len(args.Bams)))
	bams := make([]string, len(args.Bams))
	for i, b := range args.Bams {
		bams[i] = fmt.Sprintf("'%s'", b)
	}
	cmd := exec.Command("bash", "-c", fmt.Sprintf(command, args.Q, args.MaxCov, args.Chrom, strings.Join(bams, "\t")))
	pipeout, err := cmd.StdoutPipe()
	if err != nil {
		panic(err)
	}

	if err := cmd.Start(); err != nil {
		panic(err)
	}

	aggregate(bufio.NewReader(pipeout), &args, args.Chrom)
}

type site struct {
	pos0   int
	depths []int
}

func parseLine(line string) site {
	toks := strings.Split(line, "\t")
	toks[len(toks)-1] = strings.TrimSpace(toks[len(toks)-1])

	pos, err := strconv.Atoi(toks[1])
	if err != nil {
		panic(err)
	}
	s := site{pos0: pos - 1, depths: make([]int, len(toks)-2)}
	for i, str := range toks[2:] {
		s.depths[i], err = strconv.Atoi(str)
		if err != nil {
			panic(err)
		}
	}

	return s
}

func sufficientDepth(s site, a *dargs) bool {
	var n int
	for _, d := range s.depths {
		if d >= a.MinCov {
			n++
		}
	}
	return n > a.minSamples
}

func aggregate(rdr *bufio.Reader, a *dargs, chrom string) {
	stdout := bufio.NewWriter(os.Stdout)
	defer stdout.Flush()
	cache := make([]site, 0, 2000)

	for line, err := rdr.ReadString('\n'); err == nil; line, err = rdr.ReadString('\n') {
		s := parseLine(line)
		if (len(cache) == 0 || cache[len(cache)-1].pos0+1 == s.pos0) && sufficientDepth(s, a) {
			cache = append(cache, s)
		} else if len(cache) > 0 {
			dps := means(cache)
			fmt.Fprintf(stdout, "%s\t%d\t%d\t%s\n", chrom, cache[0].pos0, cache[len(cache)-1].pos0+1, strings.Join(dps, "\t"))
			cache = cache[:0]
		}
	}
	if len(cache) > 0 {
		dps := means(cache)
		fmt.Fprintf(stdout, "%s\t%d\t%d\t%s\n", chrom, cache[0].pos0, cache[len(cache)-1].pos0+1, strings.Join(dps, "\t"))
		cache = cache[:0]
	}
}

func means(sites []site) []string {
	dps := make([]float64, len(sites[0].depths))
	for _, s := range sites {
		for i, d := range s.depths {
			dps[i] += float64(d) / 1000.
		}
	}
	sdps := make([]string, len(dps))
	for i, d := range dps {
		sdps[i] = fmt.Sprintf("%.2f", d/float64(len(sites))*1000)
	}
	return sdps
}
