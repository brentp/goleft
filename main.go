package main

import (
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"

	"github.com/brentp/goleft/covmed"
	"github.com/brentp/goleft/depth"
	"github.com/brentp/goleft/depthwed"
	"github.com/brentp/goleft/indexcov"
)

const Version = "0.1.7"

type progPair struct {
	help string
	main func()
}

var progs = map[string]progPair{
	"depth":    progPair{"parallelize calls to samtools in user-defined windows", depth.Main},
	"depthwed": progPair{"matricize output from depth to n-sites * n-samples", depthwed.Main},
	"covmed":   progPair{"calculate median coverage on a bam by sampling", covmed.Main},
	"indexcov": progPair{"quick coverage estimate using only the bam index", indexcov.Main},
}

func printProgs() {

	var wtr io.Writer = os.Stdout

	fmt.Fprintf(wtr, "goleft Version: %s\n\n", Version)
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
	os.Exit(1)

}

func main() {

	if len(os.Args) < 2 {
		printProgs()
	}
	var p progPair
	var ok bool
	if p, ok = progs[os.Args[1]]; !ok {
		printProgs()
	}
	// remove the prog name from the call
	os.Args = append(os.Args[:1], os.Args[2:]...)
	p.main()
}
