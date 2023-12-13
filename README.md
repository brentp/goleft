<!--
#CGO_ENABLED=0 GOARCH=amd64 go build -o goleft_linux64 --ldflags '-extldflags "-static"' cmd/goleft/goleft.go
#GOOS=darwin GOARCH=amd64 CGO_ENABLED=0 go build -o goleft_osx --ldflags '-extldflags "-static"' cmd/goleft/goleft.go
-->
# goleft

[![Build Status](https://travis-ci.org/brentp/goleft.svg)](https://travis-ci.org/brentp/goleft)


goleft is a collection of bioinformatics tools written in
[go](https://github.com/golang/go) distributed together
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
go install github.com/brentp/goleft/cmd/goleft
```

`goleft` is also available in [bioconda](https://bioconda.github.io)

# Commands

+ [covstats](https://github.com/brentp/goleft/tree/master/covstats#covstats)   : estimate coverage and insert-size statistics on bams by sampling
+ [depth](https://github.com/brentp/goleft/tree/master/depth#depth)    : parallelize calls to samtools in user-defined windows
+ depthwed : matricize output from depth to n-sites * n-samples
+ [indexcov](https://github.com/brentp/goleft/tree/master/indexcov#indexcov) : quick coverage estimate using only the bam index
+ [indexsplit](https://github.com/brentp/goleft/tree/master/indexsplit#indexsplit) : generate regions of even data across a cohort (for parallelization)
+ [samplename](https://github.com/brentp/goleft/tree/master/samplename#samplename): report samplename(s) from a bam's SM tag
