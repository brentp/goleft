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

	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/stat"

	arg "github.com/alexflint/go-arg"
	"github.com/biogo/hts/sam"
	"github.com/biogo/store/interval"
	"github.com/brentp/goleft/depth"
	"github.com/brentp/goleft/indexcov"
)

type cliargs struct {
	N           int      `arg:"-n,required,help:number of regions to split to."`
	Fai         string   `arg:"--fai,help:fasta index file."`
	Problematic string   `arg:"-p,help:pipe-delimited list of regions to split small."`
	Indexes     []string `arg:"positional,required,help:bai's/crais to use for splitting genome."`
}

func imin(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func chop(sizes [][]float64) {
	for k, size := range sizes {
		m, std := stat.MeanStdDev(size, nil)
		max := m + 3*std
		for i, s := range size {
			if s > max {
				size[i] = 8 * m
			}
		}
		sizes[k] = size
	}
}

// return the proportion of data in each chromosome.
func getPercents(sizes [][]float64) ([]float64, []float64) {
	chop(sizes)
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
	Chrom  string
	Start  int
	End    int
	Sum    float64 // amount of data in this Chunk
	Splits int     // number of splits
}

func (c Chunk) String() string {
	return fmt.Sprintf("%s\t%d\t%d\t%.2f\t%d", c.Chrom, c.Start, c.End, c.Sum, c.Splits)
}

// Split takes paths of bams or crais and generates `N` `Chunks`
func Split(paths []string, refs []*sam.Reference, N int, probs map[string]*interval.IntTree) chan Chunk {

	ch := make(chan Chunk)
	go func() {

		// sizes will be the sum of values for all samples
		// use a float since we divide by a large number to avoid overflow
		var sizes [][]float64
		scalar := float64(1000000000)

		for _, path := range paths {
			osz := indexcov.ReadIndex(path).Sizes()
			for _, ref := range refs {
				i := ref.ID()
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

		for _, ref := range refs {
			ri := ref.ID()
			if ri >= len(sizes) || len(sizes[ri]) == 0 {
				// output the empty chrom with a sum of 0 the user isn't.
				ch <- Chunk{Chrom: ref.Name(), Start: 0, End: ref.Len(), Sum: 0, Splits: 0}
				continue
			}
			n := int(percents[ri] * float64(N))
			if n == 0 && percents[ri] > 0 {
				n = 1
			} else if n == 0 {
				ch <- Chunk{Chrom: ref.Name(), Start: 0, End: ref.Len(), Sum: 0, Splits: 0}
				continue
			}
			// we get `chunk` as a sum and then we know we have enough data.
			chunk := sums[ri] / float64(n)
			size := sizes[ri]

			var sum float64
			var lasti int
			var tree *interval.IntTree
			if probs != nil {
				if t, ok := probs[ref.Name()]; ok {
					tree = t
				}
			}
			// loop over the tiles and yield regions as soon as each is > chunk.
			for i := 0; i < len(size); i++ {
				// for single Tiles > chunk, we split into smaller regions.
				// 120 heuristic is arbitrary, may need to be tuned.
				ovl := depth.Overlaps(tree, i*indexcov.TileWidth, (i+1)*indexcov.TileWidth)
				if size[i] > chunk || (size[i] >= 0.05*chunk && ovl) {
					//if sum >= 0 && i > lasti {
					if i > lasti {
						ch <- Chunk{Chrom: ref.Name(), Start: lasti * indexcov.TileWidth, End: i * indexcov.TileWidth, Sum: sum, Splits: 1}
					}
					sum = size[i]
					nsplits := int(0.5 + (sum / (chunk / 2)))
					if nsplits > 8 {
						nsplits = 8
					} else if nsplits < 1 {
						nsplits = 1
						if ovl {
							nsplits = 3
						}
					}
					start := i * indexcov.TileWidth
					l := int(float64(indexcov.TileWidth)/float64(nsplits) + 1)
					for k := 0; k < nsplits; k++ {
						if i+k == len(size)+1 {
							ch <- Chunk{Chrom: ref.Name(), Start: start, End: ref.Len(), Sum: sum / float64(nsplits), Splits: nsplits}
						} else {
							ch <- Chunk{Chrom: ref.Name(), Start: start, End: imin(start+l, (i+1)*indexcov.TileWidth), Sum: sum / float64(nsplits), Splits: nsplits}
						}
						start += l
					}

					lasti, sum = i+1, 0
					continue
				}
				sum += size[i]
				if sum >= chunk || i == len(size)-1 || (sum >= 0.2*chunk && ovl) {
					if i == len(size)-1 {
						ch <- Chunk{Chrom: ref.Name(), Start: lasti * indexcov.TileWidth, End: ref.Len(), Sum: sum, Splits: 1}
					} else {
						ch <- Chunk{Chrom: ref.Name(), Start: lasti * indexcov.TileWidth, End: (i + 1) * indexcov.TileWidth, Sum: sum, Splits: 1}
					}
					lasti = i + 1
					sum = 0
				}
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

	var probs map[string]*interval.IntTree
	if cli.Problematic != "" {
		probs = depth.ReadTree(cli.Problematic)
	}

	var refs []*sam.Reference
	if strings.HasSuffix(cli.Indexes[0], ".bam") {
		refs = indexcov.RefsFromBam(cli.Indexes[0], "")
	} else {
		refs = indexcov.ReadFai(cli.Fai, "")
	}
	for chunk := range Split(cli.Indexes, refs, cli.N, probs) {
		fmt.Println(chunk)
	}
}
