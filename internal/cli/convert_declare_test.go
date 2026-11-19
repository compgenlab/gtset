package cli

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// A chrX file that was not imputed: no IMPUTED or DS in its header.
const chrXNotImputed = `##fileformat=VCFv4.2
##contig=<ID=chrX,length=1000>
##INFO=<ID=R2,Number=1,Type=Float,Description="r2">
##INFO=<ID=AF,Number=A,Type=Float,Description="af">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chrX	100	.	A	G	.	PASS	R2=0.9;AF=0.2	GT:DP	0/1:30	0/0:30	0/0:30
`

// The same flags convert every chromosome: the declarations add what chrX's
// header lacks and agree with what chr1's has, so both volumes capture the
// same fields and a filter on IMPUTED works across the whole store.
func TestDeclaredFieldsMakeVolumesAgree(t *testing.T) {
	dir := t.TempDir()
	store := filepath.Join(dir, "cohort")
	flags := []string{"--info", "IMPUTED,R2", "--format", "DS",
		"--declare-info", "IMPUTED:Flag", "--declare-format", "DS:Float:A"}

	chr1 := filepath.Join(dir, "chr1.vcf")
	os.WriteFile(chr1, []byte(extrasVCF), 0o644)
	chrX := filepath.Join(dir, "chrX.vcf")
	os.WriteFile(chrX, []byte(chrXNotImputed), 0o644)

	// Without the declaration chrX is refused, as before.
	if err := runVcfErr(t, "convert", "--out", filepath.Join(store, "chrX"), "--info", "IMPUTED", chrX); err == nil ||
		!strings.Contains(err.Error(), "declares no ##INFO field") {
		t.Fatalf("undeclared --info: %v", err)
	}

	runVcf(t, append([]string{"convert", "--out", filepath.Join(store, "chr1")}, append(flags, chr1)...)...)
	runVcf(t, append([]string{"convert", "--force", "--out", filepath.Join(store, "chrX")}, append(flags, chrX)...)...)
	runVcf(t, "create", store)

	// chrX's site carries no IMPUTED, so excluding imputed sites keeps it, and
	// chr1's imputed site at 200 is gone -- one query across both volumes.
	out := runVcf(t, "query", "--variant", "chr1", "--variant", "chrX", "--exclude-info", "IMPUTED",
		"--show-info", "IMPUTED,R2", "--show-format", "DS", store)
	rows := strings.Join(dataLines(out), "\n")
	if !strings.Contains(rows, "chrX\t100\tA\tG\tS1\t0/1") || strings.Contains(rows, "chr1\t200\t") {
		t.Errorf("--exclude-info across volumes:\n%s", rows)
	}
	// Shown as absent on chrX, present on chr1.
	if !strings.Contains(rows, "chrX\t100\tA\tG\tS1\t0/1\t30\t30\t.\t.\t.\t.\t.\t0.9\t.") {
		t.Errorf("chrX values:\n%s", rows)
	}
}

// A declaration that disagrees with the header is refused: the same key typed
// two ways would read differently in different volumes.
func TestDeclarationMustMatchHeader(t *testing.T) {
	dir := t.TempDir()
	chr1 := filepath.Join(dir, "chr1.vcf")
	os.WriteFile(chr1, []byte(extrasVCF), 0o644)
	for _, tc := range []struct{ flag, spec, want string }{
		{"--declare-info", "IMPUTED:Integer", "must match"},
		{"--declare-format", "DS:Float:1", "must match"},
		{"--declare-info", "IMPUTED", "want KEY:Type"},
		{"--declare-info", "X:Bool", "is not one of"},
		{"--declare-format", "Y:Flag", "cannot be Flags"},
	} {
		err := runVcfErr(t, "convert", "--out", filepath.Join(dir, "s"+tc.flag+tc.spec), tc.flag, tc.spec, chr1)
		if err == nil || !strings.Contains(err.Error(), tc.want) {
			t.Errorf("%s %s: %v, want %q", tc.flag, tc.spec, err, tc.want)
		}
	}
}
