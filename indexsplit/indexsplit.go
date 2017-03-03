// Package indexsplit is used to quickly generate evenly sized (by amount of data) regions across
// a cohort. It does this by reading the bam (or cram) index and using the file offsets as proxies
// for the amount of data. It sums the values in these bins across all samples. This gives a good
// estimate for actual reads in the region but without having to parse the bam file.
//
// A common use of this will be to generate regions to be use to parallelize variant calling fairly
// by splitting in to `N` regions with approximately equal amounts of data **across the cohort**.
package indexsplit

import (
	"fmt"
	"strings"

	arg "github.com/alexflint/go-arg"
	"github.com/biogo/hts/sam"
	"github.com/brentp/goleft/indexcov"
	"github.com/gonum/floats"
)

type cliargs struct {
	N       int      `arg:"-n,required,help:number of regions to split to."`
	Fai     string   `arg:"--fai,required,help:fasta index file."`
	Indexes []string `arg:"positional,required,help:bai's/crais to use for splitting genome."`
}

func imin(a, b int) int {
	if a < b {
		return a
	}
	return b
}

// return the proportion of data in each chromosome.
func getPercents(sizes [][]float64) ([]float64, []float64) {
	var tot float64
	pcts := make([]float64, len(sizes))
	sums := make([]float64, len(sizes))
	for i, s := range sizes {
		sums[i] = floats.Sum(s)
		tot += sums[i]
	}
	for i := range pcts {
		pcts[i] = sums[i] / tot
	}

	return pcts, sums
}

// Chunk is a region of the genome create by `Split`.
type Chunk struct {
	Chrom string
	Start int
	End   int
	Sum   float64 // amount of data in this Chunk
}

func (c Chunk) String() string {
	return fmt.Sprintf("%s\t%d\t%d\t%.2f", c.Chrom, c.Start, c.End, c.Sum)
}

// Split takes paths of bams or crais and generates `N` `Chunks`
func Split(paths []string, refs []*sam.Reference, N int) chan Chunk {

	ch := make(chan Chunk)
	go func() {

		// sizes will be the sum of values for all samples
		// use a float since we divide by a large number to avoid overflow
		var sizes [][]float64
		scalar := float64(1000000000)

		for _, path := range paths {
			osz := indexcov.ReadIndex(path).Sizes()
			for i := range refs {
				for i >= len(sizes) {
					sizes = append(sizes, make([]float64, 0))
				}
				if i > len(osz) {
					break
				}
				s, o := sizes[i], osz[i]
				m := imin(len(s), len(o))
				var j int
				// we add for as long as we have data from both...
				for j = 0; j < m; j++ {
					s[j] += float64(o[j]) / scalar
				}
				for ; j < len(o); j++ {
					s = append(s, float64(o[j])/scalar)
				}
				sizes[i] = s
			}
		}

		percents, sums := getPercents(sizes)

		for ri, ref := range refs {
			if ri >= len(sizes) || len(sizes[ri]) == 0 {
				continue
			}
			n := int(percents[ri] * float64(N))
			if n == 0 && percents[ri] > 0 {
				n = 1
			}
			// we get `chunk` as a sum and then we know we have enough data.
			chunk := sums[ri] / float64(n)
			size := sizes[ri]
			ni := 1

			for i := 0; i < len(size); {
				var sum float64
				var j int
				for j = i; j < len(size); j++ {
					sum += size[j]
					if sum > chunk && (ni < n || j == len(size)-1) {
						if ni == n {
							ch <- Chunk{Chrom: ref.Name(), Start: i * indexcov.TileWidth, End: ref.Len(), Sum: sum}
						} else {
							ni++
							ch <- Chunk{Chrom: ref.Name(), Start: i * indexcov.TileWidth, End: j * indexcov.TileWidth, Sum: sum}
						}
						break
					}
				}
				i = j
			}

		}
		close(ch)

	}()
	return ch
}

// Main is called from the goleft dispatcher.
func Main() {

	cli := &cliargs{}
	arg.MustParse(cli)

	var refs []*sam.Reference
	if strings.HasSuffix(cli.Indexes[0], ".bam") {
		refs = indexcov.RefsFromBam(cli.Indexes[0], "")
	} else {
		refs = indexcov.ReadFai(cli.Fai, "")
	}
	for chunk := range Split(cli.Indexes, refs, cli.N) {
		fmt.Println(chunk)
	}
}
