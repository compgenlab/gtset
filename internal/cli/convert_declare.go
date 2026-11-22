package cli

import (
	"fmt"
	"io"
	"strings"

	"github.com/compgenlab/cghts/gtstore"
	"github.com/compgenlab/cghts/vcf"
)

// --declare-info / --declare-format: capture a field a header does not declare.
//
// WHY. A callset split by chromosome is not always uniform: chrX and chrY are
// often not imputed, so their headers do not declare IMPUTED, TYPED or DS while
// the autosomes' do. Converting them without the field leaves their volumes
// without the column, and every query that filters on it then fails with "not
// captured" the moment it touches chrX -- where what is true is "nothing here
// was imputed". Rewriting the VCFs' headers is not always an option.
//
// A declaration states the field's Type and Number as the header would have,
// so it is captured like any other: absent in every record that does not carry
// it (a Flag reads as unset), which is exactly that answer.
//
// It is also safe to pass for every chromosome. Where the header DOES declare
// the field, the declaration must agree with it; a mismatch is an error, since
// the same key typed two ways in two volumes would read differently per
// chromosome.

var (
	vcfToGtstoreDeclareInfo   []string
	vcfToGtstoreDeclareFormat []string
)

// parseDeclaration reads KEY:Type[:Number]. Number defaults to 0 for a Flag and
// 1 otherwise.
func parseDeclaration(flag, spec string, allowFlag bool) (*vcf.AnnotationDef, error) {
	parts := strings.Split(spec, ":")
	if len(parts) < 2 || len(parts) > 3 || parts[0] == "" {
		return nil, fmt.Errorf("%s %q: want KEY:Type or KEY:Type:Number, e.g. IMPUTED:Flag or DS:Float:A", flag, spec)
	}
	key, typ := parts[0], parts[1]
	switch typ {
	case "Integer", "Float", "String":
	case "Flag":
		if !allowFlag {
			return nil, fmt.Errorf("%s %s: FORMAT fields cannot be Flags", flag, key)
		}
	default:
		return nil, fmt.Errorf("%s %s: Type %q is not one of Integer, Float, Flag, String", flag, key, typ)
	}
	number := "1"
	if typ == "Flag" {
		number = "0"
	}
	if len(parts) == 3 {
		number = parts[2]
	}
	if typ == "Flag" && number != "0" {
		return nil, fmt.Errorf("%s %s: a Flag has Number=0", flag, key)
	}
	return &vcf.AnnotationDef{ID: key, Type: typ, Number: number,
		Description: "Declared by gtset convert " + flag + "; not in the source header"}, nil
}

// applyDeclarations adds declared fields to a header that lacks them, checks
// the ones it has, and reports which were added.
func applyDeclarations(h *vcf.VcfHeader, info, format []string, warn io.Writer) error {
	type kind struct {
		flag      string
		specs     []string
		allowFlag bool
		get       func(string) (*vcf.AnnotationDef, bool)
		add       func(*vcf.AnnotationDef)
	}
	for _, k := range []kind{
		{"--declare-info", info, true, h.InfoDef, func(d *vcf.AnnotationDef) { d.IsInfo = true; h.AddInfo(d) }},
		{"--declare-format", format, false, h.FormatDef, h.AddFormat},
	} {
		for _, spec := range splitList(k.specs) {
			d, err := parseDeclaration(k.flag, spec, k.allowFlag)
			if err != nil {
				return err
			}
			if have, ok := k.get(d.ID); ok {
				// Compared as the store will record them: a Number=. field is
				// captured as String whatever its declared Type.
				if !sameCapture(have, d, k.allowFlag) {
					return fmt.Errorf("%s %s: the header declares %s as Type=%s,Number=%s; the declaration "+
						"must match it, or the field would mean different things in different volumes",
						k.flag, spec, d.ID, have.Type, have.Number)
				}
				continue
			}
			k.add(d)
			if vcfVerbose {
				fmt.Fprintf(warn, "declared %s (%s, Number=%s): not in the header, so absent in every record that lacks it\n",
					d.ID, d.Type, d.Number)
			}
		}
	}
	return nil
}

func sameCapture(have, want *vcf.AnnotationDef, info bool) bool {
	if info {
		a := gtstore.InfoFieldFor(have.ID, have.Type, have.Number)
		b := gtstore.InfoFieldFor(want.ID, want.Type, want.Number)
		return a.Type == b.Type && a.Number == b.Number
	}
	a := gtstore.FormatFieldFor(have.ID, have.Type, have.Number)
	b := gtstore.FormatFieldFor(want.ID, want.Type, want.Number)
	return a.Type == b.Type && a.Number == b.Number
}
