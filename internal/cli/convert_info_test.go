package cli

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/compgenlab/cghts/varstore"
	"github.com/compgenlab/cghts/vcf"
)

// The imputed-VCF shape --info exists for: a quality score, typed/imputed flags,
// per-ALT panel frequencies, and one multi-allelic record to prove Number=A is
// split rather than repeated.
const infoVCF = `##fileformat=VCFv4.2
##contig=<ID=chr22,length=50818468>
##INFO=<ID=R2,Number=1,Type=Float,Description="Estimated imputation accuracy">
##INFO=<ID=IMP,Number=0,Type=Flag,Description="Imputed marker">
##INFO=<ID=AF,Number=A,Type=Float,Description="Panel allele frequency">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Panel allele count">
##INFO=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">
##INFO=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##INFO=<ID=IMPUTED,Number=.,Type=String,Description="Imputation sources">
##INFO=<ID=XS,Number=.,Type=Float,Description="Variable-length scores">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=XF,Number=.,Type=Integer,Description="Variable-length per-sample values">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr22	1000	.	G	A	50	PASS	R2=0.98;AF=0.31;AC=2	GT:DP	0/1:30	0/0:30	0/1:30
chr22	2000	.	C	T	50	PASS	R2=0.42;IMP;AF=0.02;AC=1	GT:DP	0/1:30	0/0:30	0/0:30
chr22	3000	.	A	G	50	PASS	IMP;AF=0.005;AC=0	GT:DP	0/0:30	0/0:30	0/0:30
chr22	4000	.	T	C,G	50	PASS	R2=0.0;IMP;AF=0.11,0.07;AC=1,2;IMPUTED=beagle,panel5;XS=0.91,0.88	GT:DP:XF	0/1:30:3,1,4	0/2:30:5	0/2:30:.
`

func infoHeader(t *testing.T) *vcf.VcfHeader {
	t.Helper()
	p := filepath.Join(t.TempDir(), "in.vcf")
	if err := os.WriteFile(p, []byte(infoVCF), 0o644); err != nil {
		t.Fatal(err)
	}
	r, err := vcf.NewVcfFile(p)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	h, err := r.Header()
	if err != nil {
		t.Fatal(err)
	}
	return h
}

// Types and cardinality come from the header, never from the caller: the VCF
// must declare them, so asking the user to restate them would only be a way to
// get them wrong, and a mismatch lands as a column of zeros rather than an error.
func TestResolveInfoFieldsReadsTheHeader(t *testing.T) {
	h := infoHeader(t)
	got, _, err := resolveInfoFields(h, []string{"R2,IMP", "AF"})
	if err != nil {
		t.Fatal(err)
	}
	want := map[string]varstore.InfoField{
		"R2":  {Name: "R2", Column: "info_r2", Type: varstore.InfoFloat, Number: "1"},
		"IMP": {Name: "IMP", Column: "info_imp", Type: varstore.InfoFlag, Number: "0"},
		"AF":  {Name: "AF", Column: "info_af", Type: varstore.InfoFloat, Number: "A"},
	}
	if len(got) != len(want) {
		t.Fatalf("resolved %d fields, want %d: %+v", len(got), len(want), got)
	}
	for _, f := range got {
		if w, ok := want[f.Name]; !ok || f != w {
			t.Errorf("%s resolved to %+v, want %+v", f.Name, f, w)
		}
	}
}

// Number=G and Number=R carry more values than a site row can hold, and a site
// row is one ALT. Keeping "the first value" would put a number on a row that
// means something different there, which is worse than not capturing it.
func TestResolveInfoFieldsRefusesUnstorableCardinality(t *testing.T) {
	h := infoHeader(t)
	for _, id := range []string{"PL", "AD"} {
		if _, _, err := resolveInfoFields(h, []string{id}); err == nil {
			t.Errorf("%s was accepted despite Number=%s", id, mustNumber(t, h, id))
		} else if !strings.Contains(err.Error(), id) {
			t.Errorf("error for %s does not name it: %v", id, err)
		}
	}
}

func mustNumber(t *testing.T, h *vcf.VcfHeader, id string) string {
	t.Helper()
	d, ok := h.InfoDef(id)
	if !ok {
		t.Fatalf("no such INFO field %s", id)
	}
	return d.Number
}

// A misremembered field must fail loudly. Silently converting without it
// produces a store that looks complete and has quietly lost the column the run
// was for -- and the header already knows the right spelling.
func TestResolveInfoFieldsNamesTheNearMiss(t *testing.T) {
	h := infoHeader(t)
	_, _, err := resolveInfoFields(h, []string{"r2"})
	if err == nil {
		t.Fatal("a lowercase r2 was accepted")
	}
	if !strings.Contains(err.Error(), "R2") || !strings.Contains(err.Error(), "case sensitive") {
		t.Errorf("error does not offer the header's own spelling: %v", err)
	}

	if _, _, err := resolveInfoFields(h, []string{"GNOMAD_*"}); err == nil {
		t.Error("a glob matching nothing was accepted")
	}
}

// A glob is a convenience, so a member it happens to match but that cannot be
// stored is SKIPPED and reported -- otherwise `--info '*'` could never succeed
// against a real VCF, since almost every one declares a Number=R or Number=G
// field. Naming that same field explicitly still fails: see the test above.
func TestResolveInfoFieldsSkipsUnstorableGlobMatches(t *testing.T) {
	h := infoHeader(t)
	got, skipped, err := resolveInfoFields(h, []string{"A*"})
	if err != nil {
		t.Fatalf("a glob matching an unstorable field refused the run: %v", err)
	}
	names := map[string]bool{}
	for _, f := range got {
		names[f.Name] = true
	}
	if !names["AF"] || !names["AC"] {
		t.Errorf("glob A* resolved to %v, want AF and AC", names)
	}
	if names["AD"] {
		t.Error("AD is Number=R and cannot be stored, but was captured")
	}
	if len(skipped) != 1 || !strings.Contains(skipped[0], "AD") {
		t.Errorf("skipped = %v, want it to name AD so the caller is told", skipped)
	}
}

// One record with several ALTs becomes several site rows, and a Number=A field
// has one value per ALT. Attaching the record's first value to every row would
// give two alleles the same frequency -- plausible-looking and wrong.
func TestCaptureInfoSplitsPerAlt(t *testing.T) {
	h := infoHeader(t)
	fields, _, err := resolveInfoFields(h, []string{"AF", "AC", "R2", "IMP"})
	if err != nil {
		t.Fatal(err)
	}

	rec := readRecord(t, 4) // chr22:4000 T>C,G
	buf := map[string]any{}

	captureInfo(buf, rec, fields, 0)
	if buf["AF"] != 0.11 || buf["AC"] != int32(1) {
		t.Errorf("ALT 0 got AF=%v AC=%v, want 0.11 and 1", buf["AF"], buf["AC"])
	}
	captureInfo(buf, rec, fields, 1)
	if buf["AF"] != 0.07 || buf["AC"] != int32(2) {
		t.Errorf("ALT 1 got AF=%v AC=%v, want 0.07 and 2", buf["AF"], buf["AC"])
	}
	// Number=1 is the same for both rows, which is correct rather than a bug.
	if buf["R2"] != 0.0 {
		t.Errorf("R2 = %v, want 0", buf["R2"])
	}
	if buf["IMP"] != true {
		t.Errorf("IMP = %v, want true", buf["IMP"])
	}
}

// A key the record does not carry is LEFT OUT rather than zeroed, so the store
// writes a null. "The program emitted no R2 here" and "R2 is 0 here" are
// different claims and a filter for R2 >= 0.3 only treats them alike by accident.
func TestCaptureInfoLeavesAbsentKeysOut(t *testing.T) {
	h := infoHeader(t)
	fields, _, err := resolveInfoFields(h, []string{"R2", "IMP"})
	if err != nil {
		t.Fatal(err)
	}
	buf := map[string]any{}

	// chr22:3000 carries IMP but no R2.
	captureInfo(buf, readRecord(t, 3), fields, 0)
	if _, ok := buf["R2"]; ok {
		t.Errorf("R2 is absent from the record but was captured as %v", buf["R2"])
	}
	if buf["IMP"] != true {
		t.Error("IMP is present in the record but was not captured")
	}

	// And the buffer is reused, so the previous record's R2 must not survive.
	captureInfo(buf, readRecord(t, 1), fields, 0) // has R2=0.98, no IMP
	if buf["R2"] != 0.98 {
		t.Errorf("R2 = %v, want 0.98", buf["R2"])
	}
	if _, ok := buf["IMP"]; ok {
		t.Error("IMP carried over from the previous record through the reused buffer")
	}
}

// readRecord returns the nth data record (1-based) of the fixture.
func readRecord(t *testing.T, n int) *vcf.VcfRecord {
	t.Helper()
	p := filepath.Join(t.TempDir(), "in.vcf")
	if err := os.WriteFile(p, []byte(infoVCF), 0o644); err != nil {
		t.Fatal(err)
	}
	r, err := vcf.NewVcfFile(p)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	for i := 0; i < n; i++ {
		rec, err := r.NextRecord()
		if err != nil {
			t.Fatalf("reading record %d: %v", i+1, err)
		}
		if rec == nil {
			t.Fatalf("fixture has fewer than %d records", n)
		}
		if i == n-1 {
			return rec
		}
	}
	t.Fatal("unreachable")
	return nil
}

// Number=. is captured whole, as text, whatever its declared Type: XS is a Float
// in the header, but "0.91,0.88" fits no numeric column. Both rows of the
// multiallelic record get the full value -- nothing says which value belongs to
// which allele, so splitting it would invent that.
func TestCaptureInfoKeepsNumberDotWhole(t *testing.T) {
	h := infoHeader(t)
	fields, _, err := resolveInfoFields(h, []string{"IMPUTED,XS"})
	if err != nil {
		t.Fatalf("a Number=. field was refused: %v", err)
	}
	for _, f := range fields {
		if f.Type != varstore.InfoString || f.Number != "." {
			t.Errorf("%s resolved to %s/%s, want String/.", f.Name, f.Type, f.Number)
		}
	}
	rec := readRecord(t, 4) // chr22:4000 T>C,G
	buf := map[string]any{}
	for alt := 0; alt < 2; alt++ {
		captureInfo(buf, rec, fields, alt)
		if buf["IMPUTED"] != "beagle,panel5" || buf["XS"] != "0.91,0.88" {
			t.Errorf("ALT %d got IMPUTED=%v XS=%v, want the full values", alt, buf["IMPUTED"], buf["XS"])
		}
	}
	captureInfo(buf, readRecord(t, 1), fields, 0)
	if len(buf) != 0 {
		t.Errorf("a record without the fields captured %v", buf)
	}
}

// End to end: the store holds the verbatim text on every split row, and the
// manifest says the column is text from a Number=. field.
func TestConvertCapturesNumberDotInfo(t *testing.T) {
	p := filepath.Join(t.TempDir(), "in.vcf")
	if err := os.WriteFile(p, []byte(infoVCF), 0o644); err != nil {
		t.Fatal(err)
	}
	base := convert(t, p, "--info", "IMPUTED,XS")
	s, err := varstore.OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()
	for _, f := range s.InfoFields() {
		if f.Type != varstore.InfoString || f.Number != "." {
			t.Errorf("manifest records %+v, want String/.", f)
		}
	}
	rows := 0
	if err := s.SitesInfo(func(site varstore.Site, info varstore.InfoRow) bool {
		xs, ok := info.String("info_xs")
		if site.Pos == 4000 {
			rows++
			if xs != "0.91,0.88" {
				t.Errorf("chr22:4000 %s info_xs = %q", site.Alt, xs)
			}
		} else if ok {
			t.Errorf("chr22:%d carried no XS but reads %q", site.Pos, xs)
		}
		return true
	}); err != nil {
		t.Fatal(err)
	}
	if rows != 2 {
		t.Errorf("read %d rows at chr22:4000, want one per ALT", rows)
	}
}

// FORMAT follows the same rule: XF is an Integer in the header, but a call row
// keeps the sample's whole value as text, and "." stays absent.
func TestConvertCapturesNumberDotFormat(t *testing.T) {
	p := filepath.Join(t.TempDir(), "in.vcf")
	if err := os.WriteFile(p, []byte(infoVCF), 0o644); err != nil {
		t.Fatal(err)
	}
	base := convert(t, p, "--format", "XF")
	s, err := varstore.OpenParquet(base)
	if err != nil {
		t.Fatal(err)
	}
	defer s.Close()
	for _, f := range s.FormatFields() {
		if f.Type != varstore.InfoString || f.Number != "." {
			t.Errorf("manifest records %+v, want String/.", f)
		}
	}
	got := map[string]string{}
	if err := s.CallsFormat(func(c varstore.Call, f varstore.FormatRow) bool {
		if v, ok := f.String("fmt_xf"); ok {
			got[c.SampleID+":"+c.Alt] = v
		}
		return true
	}); err != nil {
		t.Fatal(err)
	}
	want := map[string]string{"S1:C": "3,1,4", "S2:G": "5"}
	if len(got) != len(want) {
		t.Errorf("captured %v, want %v", got, want)
	}
	for k, v := range want {
		if got[k] != v {
			t.Errorf("%s fmt_xf = %q, want %q", k, got[k], v)
		}
	}
}
