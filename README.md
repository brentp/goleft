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
