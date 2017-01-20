package main

import (
	"fmt"
	"os"
	"os/exec"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

func main() {
	if len(os.Args) < 2 {
		fmt.Println("\nUsage: anonymize-for-indexcov *.bam.\n\nThis will create files like sample_0001.bam through sample_$n.bam with read-groups and file names anonymized\n")
		os.Exit(1)
	}
	for i := 1; i < len(os.Args); i++ {
		f, err := os.Open(os.Args[i])
		if err != nil {
			panic(err)
		}

		br, err := bam.NewReader(f, 1)
		if err != nil {
			panic(err)
		}
		defer f.Close()
		defer br.Close()

		rgs := br.Header().RGs()
		hdr, err := sam.NewHeader(nil, br.Header().Refs())
		if err != nil {
			panic(err)
		}
		if len(rgs) == 0 {
			fmt.Fprintf(os.Stderr, "no readgroups in %s\n", os.Args[i])
		}
		for _, rg := range rgs {
			rg.Set(sam.Tag([2]byte{'S', 'M'}), "ACDSDF")
		}
		hdr.AddReadGroup(rgs[0])

		fo, err := os.Create(fmt.Sprintf("sample_%04d.bam", i))
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
		if err := exec.Command("cp", "-f", bai, fmt.Sprintf("sample_%04d.bam.bai", i)).Run(); err != nil {
			panic(err)
		}
		fmt.Printf("wrote: sample_%04d.bam\n", i)
	}
}
