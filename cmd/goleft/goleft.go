package main

import (
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"

	"github.com/brentp/goleft"
	"github.com/brentp/goleft/covstats"
	"github.com/brentp/goleft/depth"
	"github.com/brentp/goleft/depthwed"
	"github.com/brentp/goleft/indexcov"
	"github.com/brentp/goleft/indexsplit"
	"github.com/brentp/goleft/samplename"
)

type progPair struct {
	help string
	main func()
}

var progs = map[string]progPair{
	"depth":      progPair{"parallelize calls to samtools in user-defined windows", depth.Main},
	"depthwed":   progPair{"matricize output from depth to n-sites * n-samples", depthwed.Main},
	"covstats":   progPair{"coverage stats across bams by sampling", covstats.Main},
	"indexcov":   progPair{"quick coverage estimate using only the bam index", indexcov.Main},
	"indexsplit": progPair{"create regions of even coverage across bams/crams", indexsplit.Main},
	"samplename": progPair{"report samplename(s) from a bam's SM tag", samplename.Main},
}

func printProgs(ret int) {

	var wtr io.Writer = os.Stdout

	fmt.Fprintf(wtr, "goleft Version: %s\n\n", goleft.Version)
	var keys []string
	l := 5
	for k := range progs {
		keys = append(keys, k)
		if len(k) > l {
			l = len(k)
		}
	}
	fmtr := "%-" + strconv.Itoa(l) + "s : %s\n"
	sort.Strings(keys)
	for _, k := range keys {
		fmt.Fprintf(wtr, fmtr, k, progs[k].help)

	}
	os.Exit(ret)

}

func main() {

	if len(os.Args) < 2 {
		printProgs(0)
	}
	var p progPair
	var ok bool
	if p, ok = progs[os.Args[1]]; !ok {
		printProgs(1)
	}
	// remove the prog name from the call
	os.Args = append(os.Args[:1], os.Args[2:]...)
	p.main()
}
