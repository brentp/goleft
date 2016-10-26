[![GoDoc] (https://godoc.org/github.com/brentp/goleft/emdepth?status.png)](https://godoc.org/github.com/brentp/goleft/emdepth)

# emdepth
--
    import "github.com/brentp/goleft/emdepth"

emdepth uses a simplified EM algorithm to assign copy-number given a set of depths.
it is based off the implementation in cn.mops.
Like cn.mops, it works best with more samples and it assumes that most samples will
have copy-number 2.
emdepth consists of a single function EMDepth that iteratively assigns depths to copy-numbers, adjusts
the center of each copy-number bin. and re-assigns...
This package does no normalization and therefore expects incoming data to be normalized.
