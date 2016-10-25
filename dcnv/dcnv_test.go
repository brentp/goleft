package main

import (
	"log"
	"testing"
)

func TestMedian(t *testing.T) {
	m := []*interval{
		&interval{log2: 22},
		&interval{log2: 42},
		&interval{log2: 32},
		&interval{log2: 92},
		&interval{log2: 12},
		&interval{log2: 12},
		&interval{log2: 12},
		&interval{log2: 12},
		&interval{log2: 12},
		&interval{log2: 92},
		&interval{log2: 12},
		&interval{log2: 13},
		&interval{log2: 14},
		&interval{log2: 100},
		&interval{log2: 12},
		&interval{log2: 13},
		&interval{log2: 14},
	}

	log.Println(m)
	correctByMovingMedian(m, 3)
	log.Println(m)

}
