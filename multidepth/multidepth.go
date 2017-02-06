package main

import (
	"bufio"
	"fmt"
	"os"
	"os/exec"
	"runtime"
	"strconv"
	"strings"
	"sync"

	arg "github.com/alexflint/go-arg"
	"github.com/biogo/hts/bam"
	"github.com/brentp/goleft/indexcov"
)

type dargs struct {
	Q          int      `arg:"-Q,help:mapping quality cutoff"`
	Chrom      string   `arg:"required,-c,help:optional chromosome to limit analysis"`
	MinCov     int      `arg:"help:minimum depth considered callable"`
	MaxCov     int      `arg:"help:maximum depth considered callable"`
	MaxSkip    int      `arg:"-k,help:skip this many uncovered bases before forcing a new block"`
	MinSize    int      `arg:"-m,help:only report blocks of at least this length"`
	Window     int      `arg:"-w,help:discretize into windows of this size"`
	Processes  int      `arg:"-p,help:number of processors to use"`
	MinSamples float64  `arg:"help:proportion of samples with mincov coverage for a region to be reported"`
	Bams       []string `arg:"positional,required,help:bams for which to calculate depth"`
	minSamples int      `arg:"-"`
}

func chromSize(path string, chrom string) int {
	f, err := os.Open(path)
	if err != nil {
		panic(err)
	}
	b, err := bam.NewReader(f, 1)
	if err != nil {
		panic(err)
	}
	for _, ref := range b.Header().Refs() {
		if ref.Name() == chrom {
			return ref.Len()
		}
	}
	panic(fmt.Sprintf("multidepth: chromosome %s not found in %s", chrom, path))
}

func main() {
	args := dargs{Q: 10, MinCov: 7, MaxCov: 1000, MinSamples: 0.5, MinSize: 15, MaxSkip: 10, Window: 1e7}
	if p := arg.MustParse(&args); len(args.Bams) == 0 {
		p.Fail("specify > 1 bam")
	}
	if args.Processes <= 0 {
		args.Processes = runtime.NumCPU()
	}
	runtime.GOMAXPROCS(args.Processes)
	hdr := make([]string, 1, len(args.Bams)+1)
	hdr[0] = "#chrom\tstart\tend"

	args.minSamples = int(0.5 + args.MinSamples*float64(len(args.Bams)))
	for _, b := range args.Bams {
		nm, err := indexcov.GetShortName(b)
		if err != nil {
			panic(err)
		}
		hdr = append(hdr, nm)
	}
	ch := genRegions(chromSize(args.Bams[0], args.Chrom), args.Chrom)
	out := make(chan []block, 1)

	wg := &sync.WaitGroup{}
	wg.Add(args.Processes)
	for i := 0; i < args.Processes; i++ {
		go func() {
			for reg := range ch {
				aggregate(&args, reg, out)
			}
			wg.Done()
		}()
	}

	fmt.Println(strings.Join(hdr, "\t"))
	pwg := writeOut(out)

	wg.Wait()
	close(out)
	pwg.Wait()
}

func writeOut(ch chan []block) *sync.WaitGroup {
	var wg sync.WaitGroup
	wg.Add(1)
	go func() {
		stdout := bufio.NewWriter(os.Stdout)
		for blocks := range ch {
			for _, bl := range blocks {
				stdout.WriteString(bl.String())
				stdout.WriteByte('\n')
			}
		}
		stdout.Flush()
		wg.Done()
	}()
	return &wg

}

const chunkSize = 5000000

type region struct {
	chrom string
	start int // 1-based start
	i     int
}

func (r region) String() string {
	return r.chrom + ":" + strconv.Itoa(r.start)
}

func genRegions(l int, chrom string) chan region {
	ch := make(chan region, 2)
	go func() {
		var k int
		for i := 0; i < l; i += chunkSize {
			ch <- region{chrom: chrom, start: i + 1, i: k}
			k++
		}
		close(ch)
	}()
	return ch
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

type block struct {
	chrom  string
	start  int // 0-based
	end    int // 1-based
	depths []string
}

func (b block) String() string {
	return fmt.Sprintf("%s\t%d\t%d\t%s", b.chrom, b.start, b.end, strings.Join(b.depths, "\t"))
}

func splitBlocks(chrom string, cache []site, a *dargs) []block {
	blocks := make([]block, 0, 2)

	var i, lasti int
	for i < len(cache) {
		blk := block{chrom: chrom, start: cache[i].pos0}
		i++
		for ; i < len(cache) && cache[i].pos0-blk.start < a.Window; i++ {
		}
		blk.end = cache[i-1].pos0 + 1
		blk.depths = means(cache[lasti:i])
		blocks = append(blocks, blk)
		lasti = i
	}
	return blocks
}

func aggregate(a *dargs, r region, out chan []block) {
	cache := make([]site, 0, 2000)
	blocks := make([]block, 0, 256)

	cargs := append([]string{"depth", "-q", "0", "-Q", strconv.Itoa(a.Q), "-d", strconv.Itoa(a.MaxCov), "-r", r.String()}, a.Bams...)
	cmd := exec.Command("samtools", cargs...)
	cmd.Stderr = os.Stderr
	pipeout, err := cmd.StdoutPipe()
	if err != nil {
		panic(err)
	}

	if err := cmd.Start(); err != nil {
		panic(err)
	}
	rdr := bufio.NewReader(pipeout)
	// since we are running multiple adjacent blocks in parallel:
	// only start after we've seen empty
	// only end after we've seen empty
	seen0 := false
	for line, err := rdr.ReadString('\n'); err == nil; line, err = rdr.ReadString('\n') {
		s := parseLine(line)
		suf := sufficientDepth(s, a)
		// check that we've encountered ...
		if !suf {
			seen0 = true
			// once we complete a block, we bail if we're outside of the chunk.
			if s.pos0 > r.start+chunkSize {
				if len(cache) == 0 || s.pos0-cache[len(cache)-1].pos0 >= a.MaxSkip {
					if err := cmd.Process.Kill(); err != nil {
						panic(err)
					}
					break
				}
			}

		}
		// and are past a 0 region.
		if !seen0 {
			continue
		}
		if (len(cache) == 0 || s.pos0-(cache[len(cache)-1].pos0+1) <= a.MaxSkip) && suf {
			cache = append(cache, s)
		} else if len(cache) > 0 && s.pos0-(cache[len(cache)-1].pos0+1) > a.MaxSkip {
			if len(cache) >= a.MinSize {
				blocks = append(blocks, splitBlocks(a.Chrom, cache, a)...)
				//dps := means(cache)
				//blocks = append(blocks, block{chrom: a.Chrom, start: cache[0].pos0, end: cache[len(cache)-1].pos0 + 1, depths: dps})
			}
			cache = cache[:0]
			if suf {
				cache = append(cache, s)
			}
		}
	}
	if len(cache) > 0 {
		//dps := means(cache)
		//blocks = append(blocks, block{chrom: a.Chrom, start: cache[0].pos0, end: cache[len(cache)-1].pos0 + 1, depths: dps})
		blocks = append(blocks, splitBlocks(a.Chrom, cache, a)...)
		cache = cache[:0]
	}
	if err := cmd.Wait(); err != nil {
		if _, ok := err.(*exec.ExitError); !ok {
			panic(err)
		}
	}
	out <- blocks
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
