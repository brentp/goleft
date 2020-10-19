package indexcov

import (
	"testing"
)

func TestShortName(t *testing.T) {

	r, err := GetShortName("asdf.cram", true)
	if r != "asdf" || err != nil {
		t.Errorf("expected: 'asdf', got: %s", r)
	}

	r, err = GetShortName("/path/to/v1/asdf.123.cram", true)
	if r != "asdf-123" || err != nil {
		t.Errorf("expected: 'asdf-123', got: %s", r)
	}

}
