<!--
#CGO_ENABLED=0 GOARCH=amd64 go build -o goleft_linux64 --ldflags '-extldflags "-static"' main.go
#GOOS=darwin GOARCH=amd64 CGO_ENABLED=0 go build -o goleft_osx --ldflags '-extldflags "-static"' main.go
-->
# goleft

goleft is a collection of bioinformatics tools written in
[go](https://gitub.com/golang.org) distributed together
as a single binary under a liberal (MIT) license.

Running the binary `goleft` will give a list of subcommands
with a short description. Running any subcommand without
arguments will give a full help for that command.

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

##### Usage 

```
goleft covmed $bam
```
This will output an estimate of median coverage to stdout.
