// Package depth parallelizes calls to samtools depths and outputs:
// 1) $prefix.callable.bed that contains collapsed per-base regions of NO/LOW/or CALLABLE coverage.
// where low is < MinCov.
// 2) $prefix.depth.bed that contains the average depth for each window interval specified by WindowSize.
// TODO: output gc-content in depth windows.
package depth

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"log"
	"os"
	"regexp"
	"runtime"
	"strconv"
	"strings"

	arg "github.com/alexflint/go-arg"
	"github.com/brentp/faidx"
	"github.com/brentp/gargs/process"
	"github.com/brentp/xopen"
	"github.com/fatih/color"
)

type dargs struct {
	WindowSize   int    `arg:"-w,help:window size in which to calculate high-depth regions"`
	MaxMeanDepth int    `arg:"-m,help:windows with depth > than this are high-depth. The default reports the depth of all regions."`
	Ordered      bool   `args:"-o,help:force output to be in same order as input even with -p."`
	Q            int    `arg:"-Q,help:mapping quality cutoff"`
	Chrom        string `arg:"-c,help:optional chromosome to limit analysis"`
	MinCov       int    `arg:"help:minimum depth considered callable"`
	Stats        bool   `arg:"-s,help:report sequence stats [GC CpG masked] for each window"`
	Reference    string `arg:"-r,required,help:path to reference fasta"`
	Processes    int    `arg:"-p,help:number of processors to parallelize."`
	Bed          string `arg:"-b,help:optional file of positions or regions to restrict depth calculations."`
	Prefix       string `arg:"required,help:prefix for output files depth.bed and callable.bed"`
	Bam          string `arg:"positional,required,help:bam for which to calculate depth"`
}

// we echo the region first so the callback knows the full extents even if there is NOTE
// coverage for part of it.
const command = "echo %s; samtools depth -Q %d -d %d -r %s %s"

// this is the size in basepairs of the genomic chunks for parallelization.
const step = 10000000

var exitCode = 0

func pcheck(e error) {
	if e != nil {
		log.Fatal(e)
	}
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

// match chrom:start-end and chrom\tstart\tend
var re = regexp.MustCompile("(.+?)[:\t](\\d+)([\\-\t])(\\d+).*?")

func chromStartEndFromLine(line []byte) (string, int, int) {
	ret := re.FindSubmatch(line)
	if len(ret) != 5 {
		log.Fatal("couldn't get region from line", string(line))
	}
	chrom, start, isep, end := ret[1], ret[2], ret[3], ret[4]
	// convert from bed to chrom:start-end region so add 1 to start
	istart, err := strconv.Atoi(string(start))
	if err != nil {
		log.Fatal(err)
	}
	if bytes.Equal(isep, []byte{'-'}) {
		istart--
	}
	iend, err := strconv.Atoi(string(end))
	if err != nil {
		log.Fatal(err)
	}
	return string(chrom), max(istart, 0), iend
}

func regionFromLine(line []byte) string {
	chrom, start, end := chromStartEndFromLine(line)
	// convert from bed to chrom:start-end region so add 1 to start
	return fmt.Sprintf("%s:%d-%d", chrom, start+1, end)
}

// when the user specified a Bed file of regions for coverage, this is used.
func genFromBed(ch chan string, args dargs) {
	rdr, err := xopen.Ropen(args.Bed)
	pcheck(err)
	for {
		line, err := rdr.ReadBytes('\n')
		if err == io.EOF {
			break
		}
		pcheck(err)
		if len(line) == 0 {
			continue
		}
		region := regionFromLine(line)
		ch <- fmt.Sprintf(command, region, args.Q, args.MaxMeanDepth+1000,
			region, args.Bam)
	}
	close(ch)
}

func genCommands(args dargs) chan string {
	ch := make(chan string)
	if args.Bed != "" {
		go genFromBed(ch, args)
		return ch
	}
	go func() {

		rdr, err := xopen.Ropen(args.Reference + ".fai")
		pcheck(err)
		for {
			line, err := rdr.ReadString('\n')
			if err == io.EOF {
				break
			}
			pcheck(err)
			toks := strings.Split(line, "\t")

			chrom := toks[0]
			if args.Chrom != "" && chrom != args.Chrom {
				continue
			}
			length, err := strconv.Atoi(toks[1])
			pcheck(err)
			for i := 0; i < length; i += step {
				region := fmt.Sprintf("%s:%d-%d", chrom, i, min(i+step, length))
				ch <- fmt.Sprintf(command, region, args.Q, args.MaxMeanDepth+1000,
					region, args.Bam)
			}
		}
		close(ch)
	}()
	return ch
}

// Main is run from the dispatcher
func Main() {

	args := dargs{WindowSize: 250,
		MaxMeanDepth: 0,
		MinCov:       4,
		Q:            1}
	p := arg.MustParse(&args)
	if args.Prefix == "" {
		p.Fail("you must specify an output prefix")
	}
	runtime.GOMAXPROCS(args.Processes)
	run(args)
	os.Exit(exitCode)
}

type ipos struct {
	start int
}

func mean(sl []int) float64 {
	if len(sl) == 0 {
		return 0
	}
	l := float64(len(sl))
	avg := float64(0)
	for _, v := range sl {
		avg += float64(v) / l
	}
	return avg
}

func getStats(fa *faidx.Faidx, chrom string, start, end int) string {
	if fa == nil {
		return ""
	}
	st, err := fa.Stats(chrom, start, end)
	if err != nil {
		log.Println(err)
	}
	return fmt.Sprintf("\t%.3g\t%.3g\t%.3g", st.GC, st.CpG, st.Masked)
}

func run(args dargs) {

	callback := func(r io.Reader, w io.WriteCloser) error {
		rdr := bufio.NewReader(r)
		wtr := bufio.NewWriter(w)

		var fa *faidx.Faidx
		var err error
		if args.Stats {
			fa, err = faidx.New(args.Reference)
			defer fa.Close()
		}
		if err != nil {
			return err
		}

		defer w.Close()
		defer wtr.Flush()
		depthCache := make([]int, 0, args.WindowSize)
		var depth, pos int
		var lastCovClass string

		region, err := rdr.ReadBytes('\n')
		if err != nil {
			return err
		}
		// we echo the bounds of the region before the samtools depth call.
		chrom, regionStart, regionEnd := chromStartEndFromLine(region)
		lastWindow := max(regionStart/args.WindowSize, 0)

		cache := make([]ipos, 2)
		cache[0].start = regionStart - 1
		cache[1].start = regionStart - 1

		hdPath := fmt.Sprintf("%s.%s-%d-%d.tmp.depth.bed", args.Prefix, chrom, regionStart, regionEnd)
		fhHD, ferr := xopen.Wopen(hdPath)
		if ferr != nil {
			return ferr
		}
		caPath := fmt.Sprintf("%s.%s-%d-%d.tmp.callable.bed", args.Prefix, chrom, regionStart, regionEnd)
		fhCA, ferr := xopen.Wopen(caPath)
		if ferr != nil {
			return ferr
		}
		defer fhCA.Close()
		defer fhHD.Close()

		line, err := rdr.ReadString('\n')
		for err == nil {

			chrom = line[:strings.Index(line, "\t")]

			// toks starts after chrom. so [0] is pos and [1] is depth.
			toks := strings.SplitN(line[len(chrom)+1:], "\t", 2)
			pos, err = strconv.Atoi(toks[0])
			if err != nil {
				break
			}
			pos--

			if len(toks[1]) > 1 && toks[1][len(toks[1])-1] == '\n' {
				toks[1] = toks[1][:len(toks[1])-1]
			}

			depth, err = strconv.Atoi(toks[1])
			if err != nil {
				break
			}
			// if we have a full window...
			if pos/args.WindowSize != lastWindow {
				thisWindow := pos / args.WindowSize
				// print lastWindow along with any windows without any coverage.
				for iwindow := lastWindow; iwindow < thisWindow; iwindow++ {
					s := iwindow * args.WindowSize
					if s+args.WindowSize > regionStart {
						st := max(s, regionStart)
						e := min(regionEnd, s+args.WindowSize)
						stats := getStats(fa, chrom, s, e)
						// only the 1st loop of this will have a value. others will have 0.
						fhHD.WriteString(fmt.Sprintf("%s\t%d\t%d\t%.4g%s\n", chrom, st, e, mean(depthCache), stats))
						depthCache = depthCache[:0]
					}
				}
				lastWindow = thisWindow
			}
			depthCache = append(depthCache, depth)

			covClass := "CALLABLE"
			if depth == 0 {
				covClass = "NO_COVERAGE"
			} else if depth < args.MinCov {
				covClass = "LOW_COVERAGE"
			} else if args.MaxMeanDepth > 0 && depth >= args.MaxMeanDepth {
				covClass = "EXCESSIVE_COVERAGE"
			}
			// check for a gap or a change in the coverage class.
			if covClass != lastCovClass || pos != cache[1].start+1 {
				if lastCovClass != "" {
					fhCA.WriteString(fmt.Sprintf("%s\t%d\t%d\t%s\n", chrom, cache[0].start, cache[1].start+1, lastCovClass))
				}
				// also fill in block without any coverage.
				if pos != cache[1].start+1 {
					fhCA.WriteString(fmt.Sprintf("%s\t%d\t%d\t%s\n", chrom, cache[1].start+1, pos, "NO_COVERAGE"))
				}
				lastCovClass = covClass
				cache = cache[:0]
				// append twice so we can safely use cache[1] above.
				cache = append(cache, ipos{pos})
				cache = append(cache, ipos{pos})
			} else {
				cache[1].start = pos
			}
			line, err = rdr.ReadString('\n')
		}
		if cache[0].start != -1 && lastCovClass != "" {
			fhCA.WriteString(fmt.Sprintf("%s\t%d\t%d\t%s\n", chrom, cache[0].start, cache[1].start+1, lastCovClass))
		}
		if len(depthCache) > 0 {
			s := pos / args.WindowSize * args.WindowSize
			if s < regionEnd {
				st := max(s, regionStart)
				e := min(regionEnd, s+args.WindowSize)
				stats := getStats(fa, chrom, st, e)
				fhHD.WriteString(fmt.Sprintf("%s\t%d\t%d\t%.4g%s\n", chrom, st, e, mean(depthCache), stats))
				//depthCache = depthCache[:0]
			}

		}
		// we didn't get data for the full region, so it must end in no-coverage.
		if cache[1].start+1 < regionEnd {
			// If we had regions within section
			if cache[1].start != -1 {
				fhCA.WriteString(fmt.Sprintf("%s\t%d\t%d\tNO_COVERAGE\n", chrom, cache[1].start+1, regionEnd))
				// otherwise the whole region is NO_COVERAGE
			} else {
				fhCA.WriteString(fmt.Sprintf("%s\t%d\t%d\tNO_COVERAGE\n", chrom, regionStart, regionEnd))
			}
			for ds := max(regionStart, cache[1].start) / args.WindowSize * args.WindowSize; ds < regionEnd; ds += args.WindowSize {
				// keep de calc first.
				de := min(regionEnd, ds+args.WindowSize)
				s := max(ds, regionStart)
				stats := getStats(fa, chrom, s, de)
				fhHD.WriteString(fmt.Sprintf("%s\t%d\t%d\t%.4g%s\n", chrom, s, de, mean(depthCache), stats))
				depthCache = depthCache[:0]
			}
		}
		wtr.WriteString(caPath + "\n")
		wtr.WriteString(hdPath + "\n")
		wtr.Flush()
		return w.Close()
	}

	cancel := make(chan bool)
	defer close(cancel)
	stdout := bufio.NewWriter(os.Stdout)
	defer stdout.Flush()

	chrom := ""
	if args.Chrom != "" {
		chrom = "." + args.Chrom
	}
	fhca, err := xopen.Wopen(fmt.Sprintf("%s%s.callable.bed", args.Prefix, chrom))
	pcheck(err)
	fhhd, err := xopen.Wopen(fmt.Sprintf("%s%s.depth.bed", args.Prefix, chrom))
	pcheck(err)
	defer fhca.Flush()
	defer fhhd.Flush()
	opts := process.Options{Retries: 1, CallBack: callback, Ordered: args.Ordered}

	for cmd := range process.Runner(genCommands(args), cancel, &opts) {
		if ex := cmd.ExitCode(); ex != 0 && cmd.Err != io.EOF {
			c := color.New(color.BgRed).Add(color.Bold)
			fmt.Fprintf(os.Stderr, "%s\n", c.SprintFunc()(fmt.Sprintf("ERROR with command: %s", cmd)))
			exitCode = max(exitCode, ex)
		}
		if cmd.Err == io.EOF {
			continue
		}
		caPath, err := cmd.ReadString('\n')
		if err != nil {
			log.Println(cmd.CmdStr, err, cmd.Err)
		}
		caSrc, err := xopen.Ropen(strings.TrimSpace(caPath))
		pcheck(err)
		io.Copy(fhca, caSrc)
		os.Remove(strings.TrimSpace(caPath))

		hdPath, err := cmd.ReadString('\n')
		if err != nil {
			log.Println(err)
		}
		hdSrc, err := xopen.Ropen(strings.TrimSpace(hdPath))
		pcheck(err)
		io.Copy(fhhd, hdSrc)
		os.Remove(strings.TrimSpace(hdPath))
		cmd.Cleanup()
	}
	fhca.Flush()
	fhca.Close()
	fhhd.Flush()
	fhhd.Close()
}
