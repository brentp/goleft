this is a small tool to create a new stub bam file with **all real readgroups and data removed** 
so users can share the full bam index but a (small) bam file containing only a spartan header.

The index is left unchanged and just copied to a new file with an unformative name.

usage is like:

```
anonymize-for-indexcov *.bam
```

this will create sample\_0001.bam, sample\_0001.bam.bai ... sample\_$n.bam, sample\_$n.bam.bai 
which can be shared.
