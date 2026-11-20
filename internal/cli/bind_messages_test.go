package cli

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// create names what it adds, as it adds it, and names every subdirectory it
// passed over with the reason -- an unfinished conversion must not silently
// leave a chromosome out of the store.
func TestCreateReportsVolumes(t *testing.T) {
	dir := t.TempDir()
	store := filepath.Join(dir, "cohort")
	for _, c := range []string{"chr1", "chr2"} {
		src := filepath.Join(dir, c+".vcf")
		body := strings.ReplaceAll(extrasVCF, "chr1", c)
		if err := os.WriteFile(src, []byte(body), 0o644); err != nil {
			t.Fatal(err)
		}
		runVcf(t, "convert", "--out", filepath.Join(store, c), src)
	}
	// A conversion that has not finished: tables, no manifest.
	src := filepath.Join(dir, "chr3.vcf")
	os.WriteFile(src, []byte(strings.ReplaceAll(extrasVCF, "chr1", "chr3")), 0o644)
	runVcf(t, "convert", "--out", filepath.Join(store, "chr3"), src)
	if err := os.Remove(filepath.Join(store, "chr3", "manifest.json.gz")); err != nil {
		t.Fatal(err)
	}
	// And something that was never a volume.
	os.MkdirAll(filepath.Join(store, "logs"), 0o755)

	out := runVcf(t, "bind", "-v", store)
	for _, want := range []string{
		"skipping chr3: no manifest.json.gz -- a conversion that is still running or did not finish",
		"skipping logs: not a volume",
		"found chr1: chr1",
		"adding 2 volume(s): chr1, chr2",
		"[1/2] chr1 ... chr1, 4 sites, ",
		"[2/2] chr2 ... chr2, 4 sites, ",
		"checked 2 volume(s) in ",
	} {
		if !strings.Contains(out, want) {
			t.Errorf("missing %q in:\n%s", want, out)
		}
	}

	// Naming an unfinished volume is an error that names it.
	if err := runVcfErr(t, "bind", store, "chr1", "chr3"); err == nil || !strings.Contains(err.Error(), "volume chr3") {
		t.Errorf("unfinished named volume: %v", err)
	}
}
