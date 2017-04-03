package samplename

import (
	"fmt"
	"os"
	"strings"

	arg "github.com/alexflint/go-arg"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/brentp/goleft"
)

func Names(h *sam.Header) []string {
	rgs := h.RGs()
	if len(rgs) == 1 {
		v := rgs[0].Get(sam.Tag([2]byte{'S', 'M'}))
		if v == "" {
			return nil
		}
		return []string{v}
	}
	m := make(map[string]bool)
	tag := sam.Tag([2]byte{'S', 'M'})
	for _, rg := range rgs {
		v := rg.Get(tag)
		if v == "" {
			continue
		}
		m[v] = true
	}
	names := make([]string, 0, len(m))
	for sm := range m {
		names = append(names, sm)
	}
	return names
}

type cliargs struct {
	Bam        string `arg:"positional,required,help:bam for to get sample name(s)"`
	ErrorMulti bool   `arg:"-e,help:return an error if there is not exactly 1 sample in the bam."`
}

func (c cliargs) Version() string {
	return fmt.Sprintf("samplename %s", goleft.Version)
}

func Main() {
	cli := &cliargs{}
	arg.MustParse(cli)

	f, err := os.Open(cli.Bam)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	b, err := bam.NewReader(f, 1)
	if err != nil {
		panic(err)
	}
	defer b.Close()

	names := Names(b.Header())
	if cli.ErrorMulti && len(names) != 1 {
		panic(fmt.Sprintf("goleft/samplename: found multiple samples in %s", cli.Bam))
	}
	fmt.Println(strings.Join(names, "\n"))
}
