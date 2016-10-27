// Package depthwed combines files from goleft depth into matrices.
package depthwed

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"log"
	"os"
	"strconv"
	"strings"

	arg "github.com/alexflint/go-arg"
	"github.com/brentp/xopen"
)

type cliargs struct {
	Size int      `arg:"-s,required,help:sizes of windows to aggregate to must be >= window in input files."`
	Beds []string `arg:"positional,required,help:depth.bed files from goleft depth"`
}

func pcheck(e error) {
	if e != nil {
		log.Fatal(e)
	}
}

// Main is run from the dispatcher
func Main() {

	cli := cliargs{}
	arg.MustParse(&cli)
	run(cli)
}

func getNameFromFile(f string) string {
	tmp := strings.Split(f, "/")
	tmpn := tmp[len(tmp)-1]
	for _, suff := range []string{".gz", ".bed", ".depth"} {
		if strings.HasSuffix(tmpn, suff) {
			tmpn = tmpn[:len(tmpn)-len(suff)]
		}
	}
	return strings.TrimSuffix(tmpn, "\n")
}

func run(args cliargs) {

	stdout := bufio.NewWriter(os.Stdout)
	defer stdout.Flush()

	beds := make([]*xopen.Reader, len(args.Beds))
	names := make([]string, 3+len(args.Beds))
	names[0], names[1], names[2] = "#chrom", "start", "end"
	var err error
	for i, f := range args.Beds {
		beds[i], err = xopen.Ropen(f)
		pcheck(err)
		names[i+3] = getNameFromFile(f)
	}
	stdout.WriteString(strings.Join(names, "\t") + "\n")

	depths, eof := next(beds, args.Size)
	for ; !eof; depths, eof = next(beds, args.Size) {
		fmt.Fprintf(stdout, "%s\t%d\t%d", depths[0].chrom, depths[0].start, depths[0].end)
		for _, d := range depths {
			fmt.Fprintf(stdout, "\t%d", d.depth)
		}
		stdout.Write([]byte{'\n'})
	}

}

type depth struct {
	chrom string
	start int
	end   int
	depth int
}

func mustAtoi(s string) int {
	v, err := strconv.Atoi(s)
	pcheck(err)
	return v
}
func mustAtof(s string) float64 {
	v, err := strconv.ParseFloat(s, 64)
	pcheck(err)
	return v
}

func sFromLine(l string) depth {

	toks := strings.Split(strings.TrimSuffix(l, "\n"), "\t")
	dep := mustAtof(toks[3])
	d := depth{
		chrom: toks[0],
		start: mustAtoi(toks[1]),
		end:   mustAtoi(toks[2]),
	}

	d.depth = int(0.5 + dep*float64(d.end-d.start))
	return d

}

func getNextChrom(r *bufio.Reader) string {
	chromLine, _ := r.Peek(100)
	var chrom string
	if len(chromLine) > 0 {
		chrom = string(chromLine[:bytes.Index(chromLine, []byte("\t"))])
	}
	return chrom
}

func next(beds []*xopen.Reader, size int) (depths []depth, eof bool) {
	depths = make([]depth, len(beds))
	eof = false
	k := 0
	// endSeen makes sure we only print the indivisible message once per chrom
	endSeen := false

	chrom := getNextChrom(beds[0].Reader)

	for !eof && depths[0].end-depths[0].start < size && chrom == getNextChrom(beds[0].Reader) {

		for i, bed := range beds {
			line, err := bed.ReadString('\n')
			if err == io.EOF {
				if i > 0 && !eof {
					panic("not all files have same number of records")
				}
				eof = true
				continue
			} else if err != nil {
				panic(err)
			}
			if k == 0 {
				depths[i] = sFromLine(line)
				if depths[i].chrom != chrom {
					log.Fatalf("got unexpected chromosome from %s: %s", bed, depths[i].chrom)
				}
				if size%(depths[i].end-depths[i].start) != 0 && !endSeen {
					endSeen = true
					log.Printf("size %d indivisible by interval in line: %s likely chromosome change.", size, line)
				}
			} else {
				tmp := sFromLine(line)
				depths[i].end = tmp.end
				depths[i].depth += tmp.depth
			}

		}
		k++
	}
	return depths, eof
}
