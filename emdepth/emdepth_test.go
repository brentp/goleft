package emdepth

import (
	"reflect"
	"testing"
)

func TestDepth(t *testing.T) {
	v := []float32{1, 8, 33, 34, 35, 37, 31, 22, 66}
	cns := EMDepth(v)
	exp := []int{0, 1, 2, 2, 2, 2, 2, 2, 4}
	if !reflect.DeepEqual(cns, exp) {
		t.Errorf("expected: %v, got: %v", exp, cns)
	}

	exp = []int{2, 2, 2, 2, 2, 2, 2, 2, 2}
	v = []float32{30, 28, 33, 34, 35, 37, 31, 22, 38}
	cns = EMDepth(v)
	if !reflect.DeepEqual(cns, exp) {
		t.Errorf("expected: %v, got: %v", exp, cns)
	}

}

func TestBig(t *testing.T) {
	v := []float32{296.6, 16.7, 17.0, 319.2, 14.4, 16.5, 14.2}
	cns := EMDepth(v)
	exp := []int{8, 2, 2, 8, 2, 2, 2}

	if !reflect.DeepEqual(cns, exp) {
		t.Errorf("expected: %v, got: %v", exp, cns)
	}
}

func BenchmarkEMDepth(b *testing.B) {
	v := []float32{296.6, 16.7, 17.0, 319.2, 14.4, 16.5, 14.2, 22, 33, 44, 66, 22, 33, 11, 15, 18, 22, 22, 44, 31, 22, 66, 22, 21, 23, 16, 17, 19}

	s := 0
	for i := 0; i < b.N; i++ {
		cns := EMDepth(v)
		s += len(cns)
	}

}
