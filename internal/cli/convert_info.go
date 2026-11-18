package cli

import (
	"fmt"
	"sort"
	"strings"

	"github.com/compgenlab/cghts/gtstore"
	"github.com/compgenlab/cghts/vcf"
)

// --info: capturing source INFO fields into the sites catalog.
//
// WHAT THIS IS FOR. Some things worth knowing about a variant are properties of
// THIS FILE rather than of the variant: imputation quality depends on which
// reference panel was used, VQSLOD on which model was trained. They cannot come
// from an annotation service, because a different VCF of the same cohort gives
// different numbers at the same locus. Without capture they are simply lost at
// conversion.
//
// THE TYPES COME FROM THE HEADER, not from the caller. A VCF must declare
// Number and Type for every INFO field it uses, so asking the user to restate
// them would only be a way to get them wrong -- and a mismatch would land as a
// column of zeros rather than as an error.

// resolveInfoFields turns the --info selectors into declared fields, reading
// each one's Type and Number from the header.
//
// Selectors may be globs, so `--info 'R2,AF_*'` works; a glob matching nothing
// is an error rather than a silent no-op, because the usual cause is a
// misremembered field name and the store would otherwise be quietly built
// without it.
func resolveInfoFields(h *vcf.VcfHeader, selectors []string) ([]gtstore.InfoField, []string, error) {
	if len(selectors) == 0 {
		return nil, nil, nil
	}
	seen := map[string]bool{}
	var out []gtstore.InfoField
	var skipped []string

	for _, raw := range selectors {
		for _, sel := range strings.Split(raw, ",") {
			sel = strings.TrimSpace(sel)
			if sel == "" {
				continue
			}
			ids := []string{sel}
			glob := strings.ContainsAny(sel, "*?")
			if glob {
				ids = h.MatchInfoIDs(sel)
				if len(ids) == 0 {
					return nil, nil, fmt.Errorf("--info %s: no INFO field in the header matches", sel)
				}
			}
			for _, id := range ids {
				if seen[id] {
					continue
				}
				def, ok := h.InfoDef(id)
				if !ok {
					return nil, nil, fmt.Errorf(
						"--info %s: the header declares no ##INFO field by that name%s",
						id, nearestInfo(h, id))
				}
				// InfoFieldFor stores a Number=. field as its verbatim text.
				f := gtstore.InfoFieldFor(id, def.Type, def.Number)
				// NAMING A FIELD IS A REQUEST; MATCHING ONE IS A SIDE EFFECT.
				// An explicit --info AD that cannot be stored must fail, or the
				// run quietly produces a store missing the column it was for.
				// A field that merely fell inside a glob is skipped instead --
				// otherwise `--info '*'` could never succeed against a real VCF,
				// since almost every one declares something with Number=R or G.
				//
				// Refused one at a time and by name either way, because the
				// remedy is to drop that field rather than abandon the run.
				if err := gtstore.ValidateInfo([]gtstore.InfoField{f}); err != nil {
					if glob {
						skipped = append(skipped, fmt.Sprintf("%s (Number=%s)", id, def.Number))
						continue
					}
					return nil, nil, err
				}
				seen[id] = true
				out = append(out, f)
			}
		}
	}
	sort.Slice(out, func(i, j int) bool { return out[i].Name < out[j].Name })
	return out, skipped, nil
}

// nearestInfo offers the header's own field names when one is not found. A
// misremembered case ("r2" for "R2") is the common mistake and the header
// already holds the answer.
func nearestInfo(h *vcf.VcfHeader, want string) string {
	var near []string
	for _, id := range h.InfoIDs() {
		if strings.EqualFold(id, want) {
			return fmt.Sprintf(" (it does declare %s -- INFO keys are case sensitive)", id)
		}
		if strings.Contains(strings.ToUpper(id), strings.ToUpper(want)) {
			near = append(near, id)
		}
	}
	if len(near) > 0 {
		if len(near) > 6 {
			near = near[:6]
		}
		return fmt.Sprintf(" (it declares %s)", strings.Join(near, ", "))
	}
	return ""
}

// captureInfo reads one record's values for the captured fields into dst, which
// is reused between records.
//
// A KEY ABSENT FROM THE RECORD IS LEFT OUT rather than zeroed: the store's
// columns are optional precisely so "this program emitted no R2 here" stays
// distinct from "R2 is 0 here". A value the header declared but that will not
// parse is also left out rather than guessed at.
//
// altIdx selects which value a Number=A field contributes, since one VCF record
// with several ALTs becomes several site rows.
func captureInfo(dst map[string]any, rec *vcf.VcfRecord, fields []gtstore.InfoField, altIdx int) {
	clear(dst)
	// cghts's InfoValue is the one reading of a captured value, shared with the
	// query side: a Flag is present or not, a Number=A field takes this ALT's
	// element (none when the record has fewer values than ALTs -- another
	// allele's number would mean something else), and "." or an unparseable
	// value is absent.
	for _, f := range fields {
		if v, ok := gtstore.InfoValue(rec, f, altIdx); ok {
			dst[f.Name] = v
		}
	}
}
