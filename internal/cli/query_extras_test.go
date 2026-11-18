package cli

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// A mixed cohort: a typed site, an imputed one, and a multiallelic with a
// Number=A field missing for its second allele.
const extrasVCF = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=400>
##INFO=<ID=IMPUTED,Number=0,Type=Flag,Description="imputed">
##INFO=<ID=R2,Number=1,Type=Float,Description="r2">
##INFO=<ID=AF,Number=A,Type=Float,Description="af">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=DS,Number=A,Type=Float,Description="DS">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	100	.	A	G	.	PASS	R2=0.870;AF=0.25	GT:DP:DS	0/1:30:1.000	0/0:30:0	0/0:5:0.01
chr1	200	.	C	T	.	PASS	IMPUTED;R2=0.5;AF=0.4	GT:DP:DS	1/1:30:2	0/0:30:0	0/1:30:1
chr1	300	.	G	A,C	.	PASS	AF=0.1,.;R2=0.95	GT:DP:DS	1/2:30:0.9,0.6	0/0:30:0,0	0/1:30:0.8,0
`

func extrasFixture(t *testing.T) (vcfPath, store string) {
	t.Helper()
	plain := filepath.Join(t.TempDir(), "cohort.vcf")
	if err := os.WriteFile(plain, []byte(extrasVCF), 0o644); err != nil {
		t.Fatal(err)
	}
	// Indexed, so the VCF backend can answer a region query too.
	vcfPath = plain + ".gz"
	bgzipIndex(t, plain, vcfPath)
	return vcfPath, convert(t, vcfPath, "--info", "IMPUTED,R2,AF", "--format", "DS")
}

// dataLines drops the ## provenance and returns header + rows.
func dataLines(out string) []string {
	var keep []string
	for _, l := range strings.Split(strings.TrimSpace(out), "\n") {
		if strings.HasPrefix(l, "#") || strings.HasPrefix(l, "warning:") || strings.HasPrefix(l, "note:") {
			continue
		}
		keep = append(keep, l)
	}
	return keep
}

func TestQueryExcludeAndRequireInfo(t *testing.T) {
	_, store := extrasFixture(t)
	excl := strings.Join(dataLines(runVcf(t, "query", "--variant", "chr1", "--exclude-info", "IMPUTED", store)), "\n")
	if strings.Contains(excl, "\t200\t") || !strings.Contains(excl, "\t100\t") || !strings.Contains(excl, "\t300\t") {
		t.Errorf("--exclude-info IMPUTED:\n%s", excl)
	}
	req := dataLines(runVcf(t, "query", "--variant", "chr1", "--require-info", "IMPUTED", store))
	for _, l := range req[1:] {
		if !strings.Contains(l, "\t200\t") {
			t.Errorf("--require-info IMPUTED kept %s", l)
		}
	}
	if len(req) < 2 {
		t.Errorf("--require-info IMPUTED returned nothing")
	}
	// A field conversion never captured is refused, not treated as absent.
	if err := runVcfErr(t, "query", "--variant", "chr1", "--exclude-info", "NOPE", store); err == nil ||
		!strings.Contains(err.Error(), "not captured") {
		t.Errorf("uncaptured field: %v", err)
	}
}

func TestQueryShowInfoAndFormatTSV(t *testing.T) {
	vcfPath, store := extrasFixture(t)
	for name, src := range map[string]string{"store": store, "vcf": vcfPath} {
		lines := dataLines(runVcf(t, "query", "--variant", "chr1", "--show-info", "R2,IMPUTED,AF",
			"--show-format", "DS", src))
		if got := lines[0]; !strings.HasSuffix(got, "ref_source\tinfo_R2\tinfo_IMPUTED\tinfo_AF\tfmt_DS") {
			t.Errorf("%s: header %q", name, got)
		}
		joined := strings.Join(lines, "\n")
		for _, want := range []string{
			"chr1\t100\tA\tG\tS1\t0/1\t30\t30\t.\t.\t.\t.\t0.87\t.\t0.25\t1",
			"chr1\t200\tC\tT\tS1\t1/1\t30\t30\t.\t.\t.\t.\t0.5\t1\t0.4\t2",
			"chr1\t300\tG\tC\tS1\t./1\t30\t30\t.\t.\t.\t.\t0.95\t.\t.\t0.6",
		} {
			if !strings.Contains(joined, want) {
				t.Errorf("%s: missing row\n  %s\nin\n%s", name, want, joined)
			}
		}
	}
}

func TestQueryShowInfoJSON(t *testing.T) {
	_, store := extrasFixture(t)
	out := runVcf(t, "query", "--variant", "chr1:200:C:T", "--format", "json", "--show-info", "R2",
		"--show-format", "DS", store)
	if !strings.Contains(out, `"info":{"R2":"0.5"}`) || !strings.Contains(out, `"format":{"DS":"2"}`) {
		t.Errorf("json:\n%s", out)
	}
	// A query that asked for nothing gains no new keys.
	plain := runVcf(t, "query", "--variant", "chr1:200:C:T", "--format", "json", store)
	if strings.Contains(plain, "Extra") || strings.Contains(plain, `"info"`) {
		t.Errorf("plain json changed:\n%s", plain)
	}
}

func TestQueryShowFieldsInVCF(t *testing.T) {
	_, store := extrasFixture(t)
	out := runVcf(t, "query", "--variant", "chr1", "--format", "vcf", "--show-info", "R2,IMPUTED",
		"--show-format", "DS", "--min-dp", "10", store)
	for _, want := range []string{
		`##INFO=<ID=R2,Number=1,Type=Float`,
		`##INFO=<ID=IMPUTED,Number=0,Type=Flag`,
		`##FORMAT=<ID=DS,Number=A,Type=Float`,
	} {
		if !strings.Contains(out, want) {
			t.Errorf("header lacks %s", want)
		}
	}
	var rec200 string
	for _, l := range strings.Split(out, "\n") {
		if strings.HasPrefix(l, "chr1\t200\t") {
			rec200 = l
		}
	}
	f := strings.Split(rec200, "\t")
	if len(f) < 12 {
		t.Fatalf("no chr1:200 record in:\n%s", out)
	}
	if !strings.Contains(f[7], "R2=0.5") || !strings.Contains(";"+f[7]+";", ";IMPUTED;") {
		t.Errorf("INFO %q", f[7])
	}
	if f[8] != "GT:DP:GQ:DS" || !strings.HasPrefix(f[9], "1/1:30:.:2") {
		t.Errorf("FORMAT %q, S1 %q", f[8], f[9])
	}

	// AF is recomputed by --format vcf, so the source's cannot share its key.
	if err := runVcfErr(t, "query", "--variant", "chr1", "--format", "vcf", "--show-info", "AF", store); err == nil ||
		!strings.Contains(err.Error(), "writes its own AF") {
		t.Errorf("--show-info AF with --format vcf: %v", err)
	}
	// "*" places what it can and says what it left out.
	star := runVcf(t, "query", "--variant", "chr1", "--format", "vcf", "--show-info", "*", store)
	if !strings.Contains(star, "not written to the VCF") || !strings.Contains(star, "INFO/AF") {
		t.Errorf("--show-info '*' did not note the skipped AF")
	}
}

// A reference genome whose base at 300 is not the store's REF.
func writeFasta(t *testing.T, mismatchAt int) string {
	t.Helper()
	seq := []byte(strings.Repeat("A", 400))
	seq[199] = 'C'
	seq[299] = 'g' // lowercase (soft-masked) still matches
	if mismatchAt > 0 {
		seq[mismatchAt-1] = 'T'
	}
	dir := t.TempDir()
	fa := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(fa, []byte(">chr1\n"+string(seq)+"\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	fai := fmt.Sprintf("chr1\t400\t6\t400\t401\n")
	if err := os.WriteFile(fa+".fai", []byte(fai), 0o644); err != nil {
		t.Fatal(err)
	}
	return fa
}

func TestQueryReferenceCheck(t *testing.T) {
	_, store := extrasFixture(t)
	ok := runVcf(t, "query", "--variant", "chr1", "--reference", writeFasta(t, 0), store)
	if strings.Contains(ok, "warning:") || !strings.Contains(ok, "\t300\t") {
		t.Errorf("matching reference dropped or warned:\n%s", ok)
	}
	bad := runVcf(t, "query", "--variant", "chr1", "--hom-ref", "--min-dp", "10", "--reference", writeFasta(t, 300), store)
	if strings.Contains(strings.Join(dataLines(bad), "\n"), "\t300\t") {
		t.Errorf("mismatched site still reported:\n%s", bad)
	}
	for _, want := range []string{"chr1:300 REF G, reference has T; site skipped", "1 of 3 site(s) skipped"} {
		if !strings.Contains(bad, want) {
			t.Errorf("no %q in:\n%s", want, bad)
		}
	}
	// A contig the genome lacks is the wrong genome, not a bad site.
	other := filepath.Join(t.TempDir(), "other.fa")
	os.WriteFile(other, []byte(">chrZ\nACGT\n"), 0o644)
	os.WriteFile(other+".fai", []byte("chrZ\t4\t6\t4\t5\n"), 0o644)
	if err := runVcfErr(t, "query", "--variant", "chr1", "--reference", other, store); err == nil ||
		!strings.Contains(err.Error(), "has no contig chr1") {
		t.Errorf("missing contig: %v", err)
	}
}

func TestQueryPassingTargets(t *testing.T) {
	_, store := extrasFixture(t)
	panel := filepath.Join(t.TempDir(), "panel.vcf")
	os.WriteFile(panel, []byte("##fileformat=VCFv4.2\n##FILTER=<ID=Low,Description=\"x\">\n"+
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"+
		"chr1\t100\t.\tA\tG\t.\tPASS\t.\n"+
		"chr1\t200\t.\tC\tT\t.\tLow\t.\n"), 0o644)
	all := strings.Join(dataLines(runVcf(t, "query", "--variant", panel, store)), "\n")
	passing := strings.Join(dataLines(runVcf(t, "query", "--variant", panel, "--passing", "-v", store)), "\n")
	if !strings.Contains(all, "\t200\t") {
		t.Errorf("without --passing the filtered record should be asked about:\n%s", all)
	}
	if strings.Contains(passing, "\t200\t") || !strings.Contains(passing, "\t100\t") {
		t.Errorf("--passing:\n%s", passing)
	}
}
