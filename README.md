<!--
#CGO_ENABLED=0 GOARCH=amd64 go build -o goleft_linux64 --ldflags '-extldflags "-static"' main.go
#GOOS=darwin GOARCH=amd64 CGO_ENABLED=0 go build -o goleft_osx --ldflags '-extldflags "-static"' main.go
-->
# goleft

goleft is a collection of bioinformatics tools written in
[go](https://gitub.com/golang/go) distributed together
as a single binary under a liberal (MIT) license.

Running the binary `goleft` will give a list of subcommands
with a short description. Running any subcommand without
arguments will give a full help for that command.

# Installation

The easiest way to install goleft is to download the latest binary from
the [releases](https://github.com/brentp/goleft/releases) and make sure to chmod +x the resulting binary.

If you are using [go](https://github.com/golang/go), you can build from source with:
```
go get -u github.com/brentp/goleft/...
```

If you want to install the binary to your machine:

```
go install github.com/brentp/goleft/cmd/goleft
```


`goleft` is also available in [bioconda](https://bioconda.github.io)

# Commands

### depth

depth parallelizes calls to [samtools](https://samtools.github.io) in user-defined windows.

##### Usage 

```
goleft depth -Q 1 --reference $fasta --prefix t $bam -p 32 -w 50 --stats
```
will use 32 cpus to parallelize the depth coverage counting only reads
with a mapping quality (-Q) of 1 or greater. The output bed files
will have 3 additional columns for the GC content, CpG content, and fraction
of masked (lower-case) bases in the reference.


### depthwed

`depthwed` takes output from `depth` and makes a matrix -file of n-sites * n-samples

### covmed

covmed calculates median coverage by reading the bam index and getting mean read length.
It outputs median coverage, mean insert-size, sd of insert-size, mean of template length, sd of template length
to stdout.

##### Usage 

```
goleft covmed $bam
```
This will output an estimate of median coverage to stdout.

### indexcov

quickly estimate coverage from the bam index.

##### Usage 

```
goleft indexcov -c $chrom *.bam > depth.bed
```

This will create a bed file where each additional column is the normalized, estimated depth for each
sample.

