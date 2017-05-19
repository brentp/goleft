package main

import (
	"fmt"
	"log"
	"os"
	"os/exec"
	"strconv"
	"time"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/brentp/xopen"
)

func main() {
	if len(os.Args) < 3 {
		fmt.Println("\nUsage: anonymize-for-indexcov name *.bam.\n\nThis will create files like sample_0001.bam through sample_$name_$n.bam with read-groups and file names changed as well.\n")
		os.Exit(1)
	}
	if xopen.Exists(os.Args[1]) {
		log.Fatal("send (arbitrary) name as 1st argument. must not be an exisiting file.")
	}
	for i := 2; i < len(os.Args); i++ {
		f, err := os.Open(os.Args[i])
		if err != nil {
			panic(err)
		}
		name := fmt.Sprintf("sample_%s_%04d", os.Args[1], i-1)

		br, err := bam.NewReader(f, 1)
		if err != nil {
			panic(err)
		}
		defer f.Close()
		defer br.Close()

		rgs := br.Header().RGs()
		hdr := br.Header().Clone()
		if len(rgs) == 0 {
			fmt.Fprintf(os.Stderr, "no readgroups in %s\n", os.Args[i])
		}
		rg, err := sam.NewReadGroup(name, "", strconv.Itoa(i-1), "XX", "indexcov-anon", "illumina", "", name, "", "", time.Now(), 1000)
		if err != nil {
			panic(err)
		}
		if err := hdr.AddReadGroup(rg); err != nil {
			panic(err)
		}

		fo, err := os.Create(fmt.Sprintf("%s.bam", name))
		if err != nil {
			panic(err)
		}
		defer fo.Close()

		o, err := bam.NewWriter(fo, hdr, 1)
		if err != nil {
			panic(err)
		}

		if err := o.Close(); err != nil {
			panic(err)
		}
		bai := ""
		if _, err := os.Stat(f.Name() + ".bai"); err == nil {
			bai = f.Name() + ".bai"
		} else if _, err := os.Stat(f.Name()[:len(f.Name())-4] + ".bai"); err == nil {
			bai = f.Name()[:len(f.Name())-4] + ".bai"
		}
		if bai == "" {
			panic(fmt.Sprintf("unable to find bam index for %s", f.Name()))
		}
		if err := exec.Command("cp", "-f", bai, fmt.Sprintf("%s.bam.bai", name)).Run(); err != nil {
			panic(err)
		}
		fmt.Printf("wrote: %s.bam\n", name)
	}
}
