package cli

import (
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"time"

	"github.com/compgenlab/cghts/gtstore"
	"github.com/compgenlab/gtstore/internal/buildinfo"
	"github.com/spf13/cobra"
)

// A gtstore: several volumes, disjoint by chromosome, read as one.
//
// WHY THIS IS A SEPARATE STEP rather than something a conversion does. A
// whole-genome callset arrives as one VCF per chromosome and is converted one
// volume at a time -- often in parallel, often on different machines, often
// days apart. Combining is the statement that those volumes belong together,
// and it can only be made once they all exist, because making it means checking
// that they agree.
//
// It also means a gtstore can be rebuilt, or one volume replaced, without
// reconverting anything.
//
// Note that gtstore convert already produces a gtstore -- a complete one of a
// single volume. This does not promote anything to a different kind of object;
// it composes volumes into a larger archive of the same kind.

var bindCmd = &cobra.Command{
	Annotations: map[string]string{"since": "v0.1.0"},
	Use:         "bind <store-dir> [volume ...]",
	Short:       "Bind existing volumes into one gtstore",
	Long: `Bind existing volumes into one gtstore: write its manifest over them.
No genotype is read or rewritten; this checks that the volumes agree and
records them as one store.

  gtstore bind cohort/ chr1 chr2 chr3
  gtstore bind s3://bucket/cohort chr1 chr2
  gtstore bind cohort/

A volume named without a scheme or a leading slash is relative to the store
directory, which keeps a gtstore one thing to copy, move or delete. Named with
a locator, it is used as given, so volumes may live apart when they must.

Given no volumes, every subdirectory of the store that is a readable volume is
taken, in name order. A subdirectory that is not -- a conversion still running
or one that failed -- is named on stderr with the reason, so a chromosome cannot
go missing from the store unannounced. Each volume is then reported as it is
checked, with its chromosomes, counts and the time it took.

THE VOLUMES MUST AGREE, and this is where that is established rather than
rediscovered on every query:

  the same samples, in the same order    genotype columns are positional, and a
                                         mismatch attributes calls to the wrong
                                         person
  the same --min-dp                      otherwise they do not mean the same
                                         thing by "callable"
  the same --depth-bands                 otherwise a recorded depth describes a
                                         different span in each
  no shared chromosome                   otherwise a locus has two answers and
                                         no rule for choosing

A gtstore that cannot be built is a better outcome than one that answers with a
different population depending on which chromosome you ask about.

The manifest is written LAST and named gtstore.json.gz. Its presence is what
makes the directory an archive of several volumes, and gtstore and anything else
reading through cghts will then open it as one.

A gtstore is several volumes, disjoint by chromosome, read as one.

Whole-genome data ships one VCF per chromosome, so converting it gives one
volume per chromosome -- and twenty-five volumes is the wrong unit to manage,
register or reason about. A person has one genome; a gtstore is one object that
means it.

It is also what makes small shards affordable. A shard index lives in a manifest
read on every open, so its size is a budget: chr22 split every 500 sites carries
about 100 KB of index, while a whole genome split the same way carries 2.6 MB. A
single undivided archive forces a choice between coarse shards and an unreadable
manifest, and coarse shards are what make a locus query slow. Split by
chromosome, each index covers one chromosome and a query opens only the
chromosomes it touches.

The words, since two of them used to be one word:

  gtstore   the whole archive -- what you hold, name and query
  volume     one chromosome's worth of it, a complete gtstore on its own
  shard      a coordinate range within a volume`,
	Args: cobra.MinimumNArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		base := strings.TrimSuffix(args[0], "/")
		volumes := args[1:]
		errOut := cmd.ErrOrStderr()

		if len(volumes) == 0 {
			found, skipped, err := discoverVolumes(cmd, base)
			if err != nil {
				return err
			}
			// Said whether or not -v: a subdirectory that is not a finished
			// volume -- a conversion still running, or one that failed -- would
			// otherwise be left out of the store without a word, and the store
			// would answer "not assayed" for that whole chromosome.
			for _, sk := range skipped {
				fmt.Fprintf(errOut, "skipping %s: %s\n", sk.name, sk.why)
			}
			if len(found) == 0 {
				return fmt.Errorf(
					"%s holds no volumes; name them, or convert them into it first", base)
			}
			volumes = found
		}
		fmt.Fprintf(errOut, "adding %d volume(s): %s\n", len(volumes), strings.Join(volumes, ", "))

		// Each volume is opened and reported as it is reached, so a slow or
		// failing one is visible while it happens rather than after the whole
		// store is checked. BuildStore opens them again to establish that they
		// agree; an open reads only footers, so the second costs little.
		started := time.Now()
		for i, v := range volumes {
			t0 := time.Now()
			fmt.Fprintf(errOut, "  [%d/%d] %s ... ", i+1, len(volumes), v)
			st, err := gtstore.OpenParquetContext(cmd.Context(), volumeLocator(base, v))
			if err != nil {
				fmt.Fprintln(errOut, "failed")
				return fmt.Errorf("volume %s: %w", v, err)
			}
			var chroms []string
			var sites, calls int64
			if m := st.VolumeManifest(); m != nil {
				for _, c := range m.Chromosomes {
					chroms = append(chroms, c.Name)
				}
				sites, calls = m.Counts.Sites, m.Counts.Calls
			}
			st.Close()
			fmt.Fprintf(errOut, "%s, %s sites, %s calls (%s)\n", strings.Join(chroms, " "),
				comma(sites), comma(calls), time.Since(t0).Round(time.Millisecond))
		}
		man, err := gtstore.BuildStore(cmd.Context(), base, volumes,
			buildinfo.Program(), buildinfo.CommandLine())
		if err != nil {
			return err
		}
		if vcfVerbose {
			fmt.Fprintf(errOut, "checked %d volume(s) in %s\n", len(volumes), time.Since(started).Round(time.Millisecond))
		}
		if err := gtstore.WriteStoreManifest(base, *man); err != nil {
			return err
		}

		out := cmd.OutOrStdout()
		fmt.Fprintf(out, "\nwrote %s: %d volumes, %d samples\n\n",
			gtstore.StoreManifestPath(base), len(man.Volumes), len(man.Samples))
		for _, m := range man.Volumes {
			fmt.Fprintf(out, "  %-24s %-28s %s sites\n",
				m.Name, strings.Join(m.Chroms, " "), comma(m.Sites))
		}
		fmt.Fprintln(out)
		return nil
	},
}

// discoverVolumes lists the subdirectories of a store that are readable volumes.
//
// Only local paths: listing is not something every locator scheme offers, and a
// caller on object storage knows its own volume names. Naming them is also the
// safer habit -- a discovered archive silently gains whatever was left in the
// directory.
//
// It also reports what it passed over and why, so a chromosome whose conversion
// has not finished is visible rather than silently missing from the store.
func discoverVolumes(cmd *cobra.Command, base string) ([]string, []skippedDir, error) {
	entries, err := osReadDirNames(base)
	if err != nil {
		return nil, nil, fmt.Errorf("listing %s: %w (name the volumes instead)", base, err)
	}
	var out []string
	var skipped []skippedDir
	for _, name := range entries {
		dir := joinLocal(base, name)
		m, err := gtstore.ReadVolumeManifest(dir)
		if err != nil {
			skipped = append(skipped, skippedDir{name, whyNotVolume(dir, err)})
			continue
		}
		if vcfVerbose {
			var chroms []string
			for _, c := range m.Chromosomes {
				chroms = append(chroms, c.Name)
			}
			fmt.Fprintf(cmd.ErrOrStderr(), "found %s: %s\n", name, strings.Join(chroms, " "))
		}
		out = append(out, name)
	}
	sort.Strings(out)
	return out, skipped, nil
}

type skippedDir struct{ name, why string }

// whyNotVolume says why a subdirectory is not a usable volume, in the terms a
// person running conversions would act on.
func whyNotVolume(dir string, err error) string {
	if _, statErr := os.Stat(gtstore.VolumeManifestPath(dir)); statErr != nil {
		if fileExists(filepath.Join(dir, "calls.parquet")) || fileExists(filepath.Join(dir, "calls")) {
			return "no manifest.json.gz -- a conversion that is still running or did not finish"
		}
		return "not a volume (no manifest.json.gz)"
	}
	return "unreadable manifest: " + err.Error()
}

// volumeLocator resolves a volume name the way BuildStore does: relative to the
// store unless it carries a scheme or is absolute.
func volumeLocator(base, v string) string {
	if strings.Contains(v, "://") || strings.HasPrefix(v, "/") {
		return v
	}
	if strings.Contains(base, "://") {
		return strings.TrimSuffix(base, "/") + "/" + v
	}
	return filepath.Join(base, v)
}

func fileExists(p string) bool {
	_, err := os.Stat(p)
	return err == nil
}

func comma(n int64) string {
	s := fmt.Sprint(n)
	if len(s) <= 3 {
		return s
	}
	var b strings.Builder
	for i, r := range s {
		if i > 0 && (len(s)-i)%3 == 0 {
			b.WriteByte(',')
		}
		b.WriteRune(r)
	}
	return b.String()
}

func osReadDirNames(p string) ([]string, error) {
	es, err := os.ReadDir(p)
	if err != nil {
		return nil, err
	}
	var out []string
	for _, e := range es {
		if e.IsDir() {
			out = append(out, e.Name())
		}
	}
	return out, nil
}

func joinLocal(a, b string) string { return filepath.Join(a, b) }

func init() {
	addVerboseFlag(bindCmd, "Also list each volume found while searching the store directory, and the total time")
}
