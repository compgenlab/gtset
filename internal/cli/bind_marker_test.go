package cli

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// bind names the archive marker gtset.json.gz. A gtset bound before that
// carries gtstore.json.gz: it must still answer queries, and the next bind
// moves it to the current name, so nobody reconverts for a file name.
func TestBindWritesGtsetMarkerAndAdoptsTheOldName(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	runVcf(t, "convert", "--out", filepath.Join(base, "chr1"), "testdata/multi_chr1.vcf")
	runVcf(t, "convert", "--out", filepath.Join(base, "chr2"), "testdata/multi_chr2.vcf")
	runVcf(t, "bind", base)

	current := filepath.Join(base, "gtset.json.gz")
	previous := filepath.Join(base, "gtstore.json.gz")
	if _, err := os.Stat(current); err != nil {
		t.Fatalf("bind did not write gtset.json.gz: %v", err)
	}
	if _, err := os.Stat(previous); !os.IsNotExist(err) {
		t.Fatalf("bind wrote gtstore.json.gz as well: %v", err)
	}
	want := dataRowsOnly(runVcf(t, "query", "--variant", "chr2", base))

	// A gtset from before the rename.
	if err := os.Rename(current, previous); err != nil {
		t.Fatal(err)
	}
	if got := dataRowsOnly(runVcf(t, "query", "--variant", "chr2", base)); got != want {
		t.Errorf("query against the old marker name differs:\n%s\nwant\n%s", got, want)
	}
	if out := runVcf(t, "summary", base); !strings.Contains(out, "volumes") {
		t.Errorf("summary on the old marker name:\n%s", out)
	}

	runVcf(t, "bind", base)
	if _, err := os.Stat(current); err != nil {
		t.Errorf("rebind did not write gtset.json.gz: %v", err)
	}
	if _, err := os.Stat(previous); !os.IsNotExist(err) {
		t.Errorf("rebind left gtstore.json.gz behind: %v", err)
	}
}
