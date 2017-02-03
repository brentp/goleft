package main

import (
	"bufio"
	"fmt"
	"os"
	"os/exec"
	"strconv"
	"strings"

	arg "github.com/alexflint/go-arg"
	"github.com/brentp/goleft/indexcov"
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

func main() {
	args := dargs{Q: 10, MinCov: 7, MaxCov: 1000, MinSamples: 0.5}
	if p := arg.MustParse(&args); len(args.Bams) == 0 {
		p.Fail("specify > 1 bam")
	}
	hdr := make([]string, 1, len(args.Bams)+1)
	hdr[0] = "#chrom\tstart\tend"

	args.minSamples = int(0.5 + args.MinSamples*float64(len(args.Bams)))
	bams := make([]string, len(args.Bams))
	for i, b := range args.Bams {
		bams[i] = b
		nm, err := indexcov.GetShortName(b)
		if err != nil {
			panic(err)
		}
		hdr = append(hdr, nm)
	}
	cargs := append([]string{"depth", "-Q", strconv.Itoa(args.Q), "-d", strconv.Itoa(args.MaxCov), "-r", args.Chrom}, bams...)
	cmd := exec.Command("samtools", cargs...)
	cmd.Stderr = os.Stderr
	pipeout, err := cmd.StdoutPipe()
	if err != nil {
		panic(err)
	}

	if err := cmd.Start(); err != nil {
		panic(err)
	}

	fmt.Println(strings.Join(hdr, "\t"))
	aggregate(bufio.NewReader(pipeout), &args, args.Chrom)
	if err := cmd.Wait(); err != nil {
		panic(err)
	}
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
