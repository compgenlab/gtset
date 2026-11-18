package cli

import (
	"fmt"
	"io"
	"sort"
	"strings"

	"github.com/compgenlab/cghts/gtstore"
	"github.com/compgenlab/cghts/seqio"
	"github.com/compgenlab/cghts/vcf"
)

// Query options over captured INFO/FORMAT fields, and the reference check.
var (
	vcfVarQueryExcludeInfo []string
	vcfVarQueryRequireInfo []string
	vcfVarQueryShowInfo    []string
	vcfVarQueryShowFormat  []string
	vcfVarQueryPassing     bool
	vcfVarQueryReference   string
)

// splitList flattens repeatable, comma-separated flag values.
func splitList(vals []string) []string {
	var out []string
	for _, v := range vals {
		for _, p := range strings.Split(v, ",") {
			if p = strings.TrimSpace(p); p != "" {
				out = append(out, p)
			}
		}
	}
	return out
}

// fieldLister is what every backend offers for captured fields: a store's
// captured INFO/FORMAT, or a VCF header's declarations.
type fieldLister interface {
	InfoFields() []gtstore.InfoField
	FormatFields() []gtstore.FormatField
}

// shown is the resolved --show-info/--show-format request, in output order.
type shown struct {
	info   []gtstore.InfoField
	format []gtstore.FormatField
}

func (s shown) any() bool { return len(s.info)+len(s.format) > 0 }

// resolveShown expands "*" and orders the requested fields, so columns and VCF
// header lines come out the same way every time. A name the source lacks is
// left for the query itself to refuse (ErrNotCaptured), so the message is the
// library's one.
func resolveShown(store gtstore.Store, q gtstore.Query) shown {
	fl, ok := store.(fieldLister)
	if !ok {
		return shown{}
	}
	var out shown
	infoBy := map[string]gtstore.InfoField{}
	for _, f := range fl.InfoFields() {
		infoBy[f.Name] = f
	}
	seen := map[string]bool{}
	for _, k := range q.WithInfo {
		if k == "*" {
			for _, f := range fl.InfoFields() {
				if !seen[f.Name] {
					seen[f.Name] = true
					out.info = append(out.info, f)
				}
			}
			continue
		}
		if f, ok := infoBy[k]; ok && !seen[k] {
			seen[k] = true
			out.info = append(out.info, f)
		}
	}
	fmtBy := map[string]gtstore.FormatField{}
	for _, f := range fl.FormatFields() {
		fmtBy[f.Name] = f
	}
	seen = map[string]bool{}
	for _, k := range q.WithFormat {
		if k == "*" {
			for _, f := range fl.FormatFields() {
				if !seen[f.Name] {
					seen[f.Name] = true
					out.format = append(out.format, f)
				}
			}
			continue
		}
		if f, ok := fmtBy[k]; ok && !seen[k] {
			seen[k] = true
			out.format = append(out.format, f)
		}
	}
	return out
}

// columns are the extra tabular columns, after the standard layout.
func (s shown) columns() []string {
	var out []string
	for _, f := range s.info {
		out = append(out, "info_"+f.Name)
	}
	for _, f := range s.format {
		out = append(out, "fmt_"+f.Name)
	}
	return out
}

// cells renders one row's extra columns, "." where the value is absent.
func (s shown) cells(c gtstore.Call) []string {
	var out []string
	get := func(m map[string]string, k string) string {
		if v, ok := m[k]; ok && v != "" {
			return v
		}
		return "."
	}
	var info, format map[string]string
	if c.Extra != nil {
		info, format = c.Extra.Info, c.Extra.Format
	}
	for _, f := range s.info {
		out = append(out, get(info, f.Name))
	}
	for _, f := range s.format {
		out = append(out, get(format, f.Name))
	}
	return out
}

// vcfInfoReserved are the INFO keys --format vcf recomputes over the samples in
// the file; a shown field by the same name would be a second, different claim
// under one key.
var vcfInfoReserved = map[string]bool{"AC": true, "AN": true, "AF": true, "NS": true, "nhomalt": true}

// vcfPlacement decides which shown fields --format vcf can place, refusing an
// explicitly named one that collides with a field the output writes itself and
// skipping a "*"-expanded one with a note.
func vcfPlacement(s shown, q gtstore.Query, dosage bool, warn io.Writer) (shown, error) {
	named := func(keys []string, k string) bool {
		for _, x := range keys {
			if x == k {
				return true
			}
		}
		return false
	}
	var out shown
	var skipped []string
	for _, f := range s.info {
		if !vcfInfoReserved[f.Name] {
			out.info = append(out.info, f)
			continue
		}
		if named(q.WithInfo, f.Name) {
			return out, fmt.Errorf("--show-info %s: --format vcf writes its own %s, recomputed over the "+
				"samples in the file, and one key cannot carry both; use --format tsv or json to see the source's",
				f.Name, f.Name)
		}
		skipped = append(skipped, "INFO/"+f.Name)
	}
	fmtReserved := map[string]bool{"GT": true, "DP": true, "GQ": true}
	if dosage {
		fmtReserved["DS"] = true
	}
	for _, f := range s.format {
		if !fmtReserved[f.Name] {
			out.format = append(out.format, f)
			continue
		}
		if named(q.WithFormat, f.Name) {
			return out, fmt.Errorf("--show-format %s: --format vcf already writes FORMAT/%s%s",
				f.Name, f.Name, map[bool]string{true: " (from --dosage)", false: ""}[f.Name == "DS"])
		}
		skipped = append(skipped, "FORMAT/"+f.Name)
	}
	if len(skipped) > 0 {
		fmt.Fprintf(warn, "note: not written to the VCF, which carries its own: %s\n", strings.Join(skipped, ", "))
	}
	return out, nil
}

// headerDefs declares the placed fields with the Number and Type the source
// gave them, so a reader parses the values as they were captured.
func (s shown) headerDefs(h *vcf.VcfHeader) {
	for _, f := range s.info {
		h.AddInfo(&vcf.AnnotationDef{IsInfo: true, ID: f.Name, Number: vcfNumber(f.Number, f.Type),
			Type: string(f.Type), Description: "From the source (captured " + f.Column + ")"})
	}
	for _, f := range s.format {
		h.AddFormat(&vcf.AnnotationDef{ID: f.Name, Number: vcfNumber(f.Number, f.Type),
			Type: string(f.Type), Description: "From the source (captured " + f.Column + ")"})
	}
}

// vcfNumber is the Number a placed field is declared with. Every record this
// output writes has one ALT, so a Number=A field holds that ALT's one value.
func vcfNumber(number string, t gtstore.InfoType) string {
	if t == gtstore.InfoFlag {
		return "0"
	}
	if number == "" {
		return "1"
	}
	return number
}

// refChecker drops sites whose REF does not match the reference genome.
//
// A store answers in the coordinates it was converted with, and nothing in a
// query can tell that the reference a caller has in mind is the same one: a
// GRCh37 store read as GRCh38, or a liftover gone wrong, answers confidently and
// wrongly. Checking REF against the genome catches that at the first site.
type refChecker struct {
	path  string
	ref   seqio.ReferenceReader
	names map[string]string // CanonKey -> the FASTA's own name
	seen  map[string]bool   // chrom:pos:ref -> matches
	warn  io.Writer

	nBad  int
	nSeen int
}

func openRefChecker(path string, warn io.Writer) (*refChecker, error) {
	ref, err := seqio.OpenReference(path)
	if err != nil {
		return nil, fmt.Errorf("--reference %s: %w", path, err)
	}
	rc := &refChecker{path: path, ref: ref, names: map[string]string{}, seen: map[string]bool{}, warn: warn}
	for _, n := range ref.Names() {
		rc.names[gtstore.CanonKey(n)] = n
	}
	return rc, nil
}

func (rc *refChecker) Close() error { return rc.ref.Close() }

// keep reports whether a row's REF matches the genome. Each site is checked
// once, however many samples it has. A contig the FASTA does not have is an
// error rather than a mismatch: that is the wrong genome, not a bad site.
func (rc *refChecker) keep(c gtstore.Call) (bool, error) {
	key := fmt.Sprintf("%s:%d:%s", gtstore.CanonKey(c.Chrom), c.Pos, c.Ref)
	if ok, done := rc.seen[key]; done {
		return ok, nil
	}
	name, ok := rc.names[gtstore.CanonKey(c.Chrom)]
	if !ok {
		return false, fmt.Errorf("--reference %s has no contig %s (it has %d, e.g. %s); "+
			"is it the genome this store was called against?", rc.path, c.Chrom, len(rc.names), sampleNames(rc.ref.Names()))
	}
	start := int(c.Pos) - 1
	got, err := rc.ref.GetSequenceRange(name, start, start+len(c.Ref))
	if err != nil {
		return false, fmt.Errorf("--reference %s at %s:%d: %w", rc.path, c.Chrom, c.Pos, err)
	}
	match := refMatches(c.Ref, string(got))
	rc.seen[key] = match
	rc.nSeen++
	if !match {
		rc.nBad++
		const maxListed = 10
		if rc.nBad <= maxListed {
			fmt.Fprintf(rc.warn, "warning: %s:%d REF %s, reference has %s; site skipped\n",
				c.Chrom, c.Pos, c.Ref, orDot(string(got)))
		}
	}
	return match, nil
}

// refMatches compares REF to the genome case-insensitively. An N in the
// genome matches any base -- it says nothing -- but a REF running past the
// contig's end does not match.
func refMatches(ref, genome string) bool {
	if len(genome) != len(ref) {
		return false
	}
	for i := 0; i < len(ref); i++ {
		a, b := upper(ref[i]), upper(genome[i])
		if a != b && b != 'N' {
			return false
		}
	}
	return true
}

func upper(b byte) byte {
	if b >= 'a' && b <= 'z' {
		return b - 'a' + 'A'
	}
	return b
}

func orDot(s string) string {
	if s == "" {
		return "nothing (past the contig's end)"
	}
	return s
}

func sampleNames(names []string) string {
	n := append([]string(nil), names...)
	sort.Strings(n)
	if len(n) > 3 {
		n = n[:3]
	}
	return strings.Join(n, ", ")
}

// report states how many sites failed the check -- always, not only under -v:
// a skipped site is a missing answer, and a reader must not take it for a
// negative.
func (rc *refChecker) report() {
	if rc.nBad == 0 {
		return
	}
	fmt.Fprintf(rc.warn, "warning: %d of %d site(s) skipped: REF does not match %s\n", rc.nBad, rc.nSeen, rc.path)
}
