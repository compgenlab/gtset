package cli

import (
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"
)

// fragmentsFor writes one fragment per sample for roster and returns their
// directories, in sample order.
func fragmentsFor(t *testing.T, roster string, samples ...string) []string {
	t.Helper()
	dir := t.TempDir()
	var out []string
	for _, s := range samples {
		_, g := refBlocksFixture(t, s)
		f := filepath.Join(dir, s)
		runVcf(t, "gvcf-extract", "--roster", roster, "--out", f, g)
		out = append(out, f)
	}
	return out
}

// --append on a store with no blocks is the first attach; after that it adds
// only samples not yet attached, so it can be rerun with every fragment there
// is and a growing fragment directory needs no bookkeeping.
func TestAttachAppendAddsOnlyNewSamples(t *testing.T) {
	roster, _ := refBlocksFixture(t, "S1")
	frags := fragmentsFor(t, roster, "S1", "S2")

	runVcf(t, "gvcf-attach", "--append", roster, frags[0])
	if got := homRefSamplesAt(t, roster, "chr1:150:C:T"); !reflect.DeepEqual(got, []string{"S1"}) {
		t.Fatalf("after the first --append: %v", got)
	}

	out := runVcf(t, "gvcf-attach", "--append", roster, frags[0], frags[1])
	if !strings.Contains(out, "skipping 1 sample(s) already attached: S1") ||
		!strings.Contains(out, "appending 1 sample(s)") {
		t.Errorf("second --append did not say what it did:\n%s", out)
	}
	if got := homRefSamplesAt(t, roster, "chr1:150:C:T"); !reflect.DeepEqual(got, []string{"S1", "S2"}) {
		t.Errorf("after appending S2: %v", got)
	}
	if s := runVcf(t, "summary", roster); !strings.Contains(s, "from all 2 samples") {
		t.Errorf("summary after --append:\n%s", s)
	}

	// Everything already there: a no-op, not an error.
	if out := runVcf(t, "gvcf-attach", "--append", roster, frags[0], frags[1]); !strings.Contains(out, "nothing to append") {
		t.Errorf("rerun:\n%s", out)
	}
	if got := homRefSamplesAt(t, roster, "chr1:150:C:T"); !reflect.DeepEqual(got, []string{"S1", "S2"}) {
		t.Errorf("after a no-op --append: %v", got)
	}

	if err := runVcfErr(t, "gvcf-attach", "--append", "--rebuild", roster, frags[0]); err == nil ||
		!strings.Contains(err.Error(), "alternatives") {
		t.Errorf("--append with --rebuild: %v", err)
	}
}

// On an archive "already attached" means attached in every volume. A run that
// died partway leaves a sample in some volumes only; --append attaches it again
// everywhere, so rerunning it finishes the job instead of skipping the sample.
func TestAttachAppendOnAnArchiveFinishesAPartialAttach(t *testing.T) {
	base := filepath.Join(t.TempDir(), "cohort")
	runVcf(t, "convert", "--out", filepath.Join(base, "chr1"), "testdata/multi_chr1.vcf")
	runVcf(t, "convert", "--out", filepath.Join(base, "chr2"), "testdata/multi_chr2.vcf")
	runVcf(t, "create", base)
	frags := fragmentsFor(t, base, "S1", "S2")

	runVcf(t, "gvcf-attach", "--append", base, frags[0])
	// S2 reaches chr1's volume only, as if the run had died after it. (A
	// chr1-only fragment, since attaching straight to one volume would make
	// it take chr2 as well.)
	var chr1Only []string
	for _, l := range strings.SplitAfter(refBlocksGVCF("S2"), "\n") {
		if !strings.HasPrefix(l, "chr2\t") {
			chr1Only = append(chr1Only, l)
		}
	}
	g := filepath.Join(t.TempDir(), "S2.chr1.g.vcf")
	if err := os.WriteFile(g, []byte(strings.Join(chr1Only, "")), 0o644); err != nil {
		t.Fatal(err)
	}
	partial := filepath.Join(t.TempDir(), "S2")
	runVcf(t, "gvcf-extract", "--roster", base, "--out", partial, g)
	runVcf(t, "gvcf-attach", "--replace", "S2", filepath.Join(base, "chr1"), partial)
	if got := homRefSamplesAt(t, base, "chr2:150:C:T"); !reflect.DeepEqual(got, []string{"S1"}) {
		t.Fatalf("setup: chr2:150 should have only S1, got %v", got)
	}

	out := runVcf(t, "gvcf-attach", "--append", base, frags[0], frags[1])
	if !strings.Contains(out, "skipping 1 sample(s) already attached: S1") ||
		!strings.Contains(out, "appending 1 sample(s)") {
		t.Errorf("--append over a partial attach:\n%s", out)
	}
	if s := runVcf(t, "summary", base); !strings.Contains(s, "from all 2 samples") || strings.Contains(s, "differ by volume") {
		t.Errorf("summary after finishing the attach:\n%s", s)
	}
	for _, locus := range []string{"chr1:150:C:T", "chr2:150:C:T"} {
		if got := homRefSamplesAt(t, base, locus); !reflect.DeepEqual(got, []string{"S1", "S2"}) {
			t.Errorf("%s: %v", locus, got)
		}
	}
}
