package emdepth

import (
	"fmt"
	"reflect"
	"testing"
)

var p = Position{123, 456}

func TestDepth(t *testing.T) {
	v := []float32{1, 8, 33, 34, 35, 37, 31, 22, 66}
	cns := EMDepth(v, p).CN()
	exp := []int{0, 1, 2, 2, 2, 2, 2, 2, 4}
	if !reflect.DeepEqual(cns, exp) {
		t.Errorf("expected: %v, got: %v", exp, cns)
	}

	exp = []int{2, 2, 2, 2, 2, 2, 2, 2, 2}
	v = []float32{30, 28, 33, 34, 35, 37, 31, 22, 38}
	cns = EMDepth(v, p).CN()
	if !reflect.DeepEqual(cns, exp) {
		t.Errorf("expected: %v, got: %v", exp, cns)
	}

}

func TestBig(t *testing.T) {
	v := []float32{296.6, 16.7, 17.0, 3019.2, 14.4, 16.5, 14.2, 26, 7}
	cns := EMDepth(v, p)
	exp := []int{8, 2, 2, 8, 2, 2, 2, 3, 1}
	_ = exp

	if !reflect.DeepEqual(cns.CN(), exp) {
		t.Errorf("expected: %v, got: %v", exp, cns.CN())
	}

	v2 := []float32{96.6, 16.7, 17.0, 319.2, 14.4, 16.5, 14.2, 7, 16}
	em2 := EMDepth(v2, p)
	idx, changed, pct := em2.Same(cns)
	if pct != 7./9. {
		t.Errorf("percent of sample with same state... expected: %v, got: %v", 7./9., pct)
	}
	if !reflect.DeepEqual(idx, []int{0, 3}) {
		t.Errorf("sample with non CN2 state... expected: 0, 3, got: %v", idx)
	}
	if !reflect.DeepEqual(changed, []int{7, 8}) {
		t.Errorf("expected change in 7,8, got: %v", changed)
	}

	fmt.Println(changed)
	fmt.Println(cns.Log2FC())
	fmt.Println(em2.Log2FC())
}

func BenchmarkEMDepth(b *testing.B) {
	v := []float32{296.6, 16.7, 17.0, 319.2, 14.4, 16.5, 14.2, 22, 33, 44, 66, 22, 33, 11, 15, 18, 22, 22, 44, 31, 22, 66, 22, 21, 23, 16, 17, 19}
	v = append(v, v...)

	s := 0
	for i := 0; i < b.N; i++ {
		cns := EMDepth(v, p).CN()
		s += len(cns)
	}

}
