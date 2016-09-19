// Package depth parallelizes calls to samtools depths and outputs:
// 1) $prefix.callable.bed that contains collapsed per-base regions of NO/LOW/or CALLABLE coverage.
// where low is < MinCov.
// 2) $prefix.depth.bed that contains the average depth for each window interval specified by WindowSize.
// TODO: output gc-content in depth windows.
package depth

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"
	"path/filepath"
	"runtime"
	"runtime/pprof"
	"strconv"
	"strings"

	arg "github.com/alexflint/go-arg"
	"github.com/brentp/faidx"
	"github.com/brentp/gargs/process"
	"github.com/brentp/xopen"
	"github.com/fatih/color"
)

type dargs struct {
	WindowSize   int     `arg:"-w,help:window size in which to calculate high-depth regions"`
	MaxMeanDepth float64 `arg:"-m,help:windows with depth > than this are high-depth. The default reports the depth of all regions."`
	Q            int     `arg:"-Q,help:mapping quality cutoff"`
	Chrom        string  `arg:"-c,help:optional chromosome to limit analysis"`
	MinCov       int     `arg:"help:minimum depth considered callable"`
	Stats        bool    `arg:"-s,help:report sequence stats [GC CpG masked] for each window"`
	Reference    string  `arg:"-r,required,help:path to reference fasta"`
	Processes    int     `arg:"-p,help:number of processors to parallelize."`
	Prefix       string  `arg:"positional,required,help:prefix for output files [\-depth.bed, \-callable.bed]"`
	Bam          string  `arg:"positional,required,help:bam for which to calculate depth"`
}

const command = "samtools depth -Q %d -r %s --reference %s %s"

var exitCode = 0

func pcheck(e error) {
	if e != nil {
		panic(e)
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

func genCommands(args dargs) chan string {
	ch := make(chan string, 20)
	go func() {
		step := 5000000

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
				ch <- fmt.Sprintf(command, args.Q, region, args.Reference, args.Bam)
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
	arg.MustParse(&args)
	runtime.GOMAXPROCS(args.Processes)
	run(args)
	os.Exit(exitCode)
}

type ipos struct {
	chrom string
	pos   int
	depth int
}

func mean(sl []int) float64 {
	sum := 0
	for _, v := range sl {
		sum += v
	}
	return float64(sum) / float64(len(sl))
}

func getStats(fa *faidx.Faidx, chrom string, start, end int) string {
	if fa == nil {
		return ""
	}
	st, err := fa.Stats(chrom, start, end)
	if err != nil {
		log.Println(err)
	}
	return fmt.Sprintf("\t%.2f\t%.2f\t%.2f", st.GC, st.CpG, st.Masked)
}

func run(args dargs) {

	f, err := os.Create("depth.pprof")
	if err != nil {
		log.Fatal(err)
	}
	pprof.StartCPUProfile(f)
	defer pprof.StopCPUProfile()

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
		cache := make([]ipos, 2)
		depthCache := make([]int, 0, args.WindowSize)
		var depth, pos int
		var lastCov, lastChrom string
		lastWindow := -1

		line, err := rdr.ReadString('\n')
		if err != nil {
			return err
		}
		toks := strings.SplitN(line, "\t", 3)
		hdPath := filepath.Join(os.TempDir(), fmt.Sprintf("%s-%s.depth.bed", toks[0], toks[1]))
		fhHD, ferr := xopen.Wopen(hdPath)
		if ferr != nil {
			return ferr
		}
		caPath := filepath.Join(os.TempDir(), fmt.Sprintf("%s-%s.callable.bed", toks[0], toks[1]))
		fhCA, ferr := xopen.Wopen(caPath)
		if ferr != nil {
			return ferr
		}
		defer fhCA.Close()
		defer fhHD.Close()
		var chrom string

		for err == nil {

			chrom = line[:strings.Index(line, "\t")]

			toks = strings.SplitN(line[len(chrom)+1:], "\t", 2)
			pos, err = strconv.Atoi(toks[0])
			if err != nil {
				break
			}

			if len(toks[1]) > 1 && toks[1][len(toks[1])-1] == '\n' {
				toks[1] = toks[1][:len(toks[1])-1]
			}

			depth, err = strconv.Atoi(toks[1])
			if err != nil {
				break
			}
			// if we have a full window...
			if pos/args.WindowSize != lastWindow {
				s := lastWindow * args.WindowSize
				if d := mean(depthCache); d >= args.MaxMeanDepth && len(depthCache) > 0 {
					stats := getStats(fa, chrom, s, s+args.WindowSize)
					fhHD.WriteString(fmt.Sprintf("%s\t%d\t%d\t%.2f%s\n", chrom, s, s+args.WindowSize, mean(depthCache), stats))
				}
				depthCache = depthCache[:0]
				lastWindow = pos / args.WindowSize
			}
			depthCache = append(depthCache, depth)

			cov := "CALLABLE"
			if depth == 0 {
				cov = "NO_COVERAGE"
			} else if depth < args.MinCov {
				cov = "LOW_COVERAGE"
			}
			// check for a gap or a change in the coverage class.
			if cov != lastCov || pos != cache[1].pos+1 {
				if lastChrom != "" {
					stats := getStats(fa, lastChrom, cache[0].pos-1, cache[1].pos)
					fhCA.WriteString(fmt.Sprintf("%s\t%d\t%d\t%s%s\n", lastChrom, cache[0].pos-1, cache[1].pos, lastCov, stats))
				}
				lastCov = cov
				lastChrom = chrom
				cache = cache[:0]
				// append twice so we can safely use cache[1] above.
				cache = append(cache, ipos{chrom, pos, depth})
				cache = append(cache, ipos{chrom, pos, depth})
			} else {
				cache[1].pos = pos
			}

			line, err = rdr.ReadString('\n')
		}
		if len(cache) > 1 {
			fhCA.WriteString(fmt.Sprintf("%s\t%d\t%d\t%s\n", lastChrom, cache[0].pos-1, cache[1].pos, lastCov))
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

	for cmd := range process.Runner(genCommands(args), 1, cancel, callback) {
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
	}
	fhca.Flush()
	fhca.Close()
	fhhd.Flush()
	fhhd.Close()

}
