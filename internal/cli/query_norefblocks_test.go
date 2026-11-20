package cli

import (
	"path/filepath"
	"strings"
	"testing"
)

// --no-ref-blocks answers from one store as though its gVCF blocks had never
// been attached, so the two can be compared without keeping two stores.
func TestQueryNoRefBlocks(t *testing.T) {
	roster, s1 := refBlocksFixture(t, "S1")
	frag := filepath.Join(t.TempDir(), "S1")
	runVcf(t, "gvcf-extract", "--roster", roster, "--out", frag, s1)
	runVcf(t, "gvcf-attach", roster, frag)

	const off = "chr1:150:C:T" // between catalog sites, inside S1's block
	if got := homRefSamplesAt(t, roster, off); len(got) == 0 {
		t.Fatalf("fixture: %s is not answered from blocks", off)
	}
	out := runVcf(t, "query", "--hom-ref", "--variant", off, "--no-ref-blocks", roster)
	if rows := tsvDataRows(out); len(rows) != 0 {
		t.Errorf("--no-ref-blocks still answered %s:\n%s", off, rows)
	}
	// And it says so rather than claiming blocks will answer.
	if strings.Contains(out, "samples are reported from reference blocks") ||
		!strings.Contains(out, "warning: "+off+" is not in the source") {
		t.Errorf("off-catalog note under --no-ref-blocks:\n%s", out)
	}
}
