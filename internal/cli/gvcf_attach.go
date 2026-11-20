package cli

import (
	"bufio"
	"context"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/compgenlab/cghts/gtstore"
	"github.com/compgenlab/cghts/iosource"
	"github.com/compgenlab/pendium/internal/buildinfo"
	"github.com/spf13/cobra"
)

var (
	attachRebuild bool
	attachAppend  bool
	attachReplace []string
	attachFrom    string
	attachBatch   int
	attachTempDir string
)

var gvcfAttachCmd = &cobra.Command{
	Annotations: map[string]string{"since": "v0.1.0"},
	Use:         "gvcf-attach <store> [fragment ...]",
	Short:       "Merge reference-block fragments into a store's ref_blocks",
	Long: `Merge the per-sample fragments pendium-gt gvcf-extract wrote into a store's ref_blocks
table, so the store answers at positions off its catalog.

  pendium-gt gvcf-attach cohort/ frags/
  pendium-gt gvcf-attach cohort/ --from fragments.txt
  pendium-gt gvcf-attach cohort/ frags/ --append
  pendium-gt gvcf-attach cohort/ frags/S0042 --replace S0042

A fragment argument is a fragment directory, or a local directory whose
subdirectories are fragments; --from reads one per line, for when there are
too many for a command line. Every fragment must be finished, one per sample,
for a sample on the store's roster, and built with the store's --min-dp and
depth bands. A fragment from an older gvcf-extract that merged across gaps is
refused, since its blocks may claim bases that were never assayed.

Every chromosome a fragment has blocks on becomes answerable. A store routes a
query by chromosome to the volume that holds it, so a chromosome with no
catalog site -- chrM, unplaced contigs, chrY in an all-female cohort -- gets a
BLOCKS-ONLY volume of its own on an archive (same samples and settings, empty
catalog), or joins the volume itself when the store is a single volume.

The first attach builds the table from the fragments. After that, say which:

  --append            add the samples not attached yet, and keep everyone
                      already attached as they are; a fragment for a sample
                      already in the table is skipped. Safe to rerun with the
                      same fragments, and to point at a directory that keeps
                      growing. On a store with no blocks yet it is the first
                      attach.
  --rebuild           build it again from these fragments alone; the old
                      table is ignored. The fragments are the source of truth.
  --replace S1,S2     keep the existing table for everyone else, and swap in
                      these samples' fragments -- for a patient whose gVCF was
                      redone. A replaced sample with no fragment is removed.
                      Unchanged fragments need not have been kept. Use
                      this, not --append, to redo a sample already attached.

Either way every file a changed sample appears in is rewritten, because the
table is sorted by position across all samples. What differs is what is read.

Fragments are merged in batches of --batch (default 100), through local
intermediate files, so a ten-thousand-sample cohort never holds ten thousand
files open. The output is a new generation of files and each volume's manifest
is rewritten last: until then the volume reads exactly as before, and a run
that dies leaves it as it was. On an archive every volume is attached and the
archive manifest rewritten after them; a run that dies partway leaves volumes
that disagree, which the archive refuses to open until it is rerun.`,
	Args: cobra.MinimumNArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		ctx := cmd.Context()
		base := strings.TrimSuffix(args[0], "/")

		dirs, err := expandFragments(args[1:], attachFrom)
		if err != nil {
			return err
		}
		if len(dirs) == 0 && len(attachReplace) == 0 {
			return fmt.Errorf("no fragments given")
		}
		if attachAppend && (attachRebuild || len(attachReplace) > 0) {
			return fmt.Errorf("--append, --rebuild and --replace are alternatives; give one")
		}
		frags, err := gtstore.LoadFragments(ctx, dirs)
		if err != nil {
			return err
		}
		opts := gtstore.AttachOpts{
			Rebuild: attachRebuild,
			Replace: attachReplace,
			Batch:   attachBatch,
			TempDir: attachTempDir,
		}
		if attachAppend {
			var skipped []string
			frags, opts.Replace, skipped, err = planAppend(ctx, base, frags)
			if err != nil {
				return err
			}
			warn := cmd.ErrOrStderr()
			if len(skipped) > 0 {
				fmt.Fprintf(warn, "skipping %d sample(s) already attached: %s\n",
					len(skipped), listSome(skipped, 10))
			}
			if len(frags) == 0 {
				fmt.Fprintln(cmd.OutOrStdout(), "nothing to append: every fragment's sample is already attached")
				return nil
			}
			fmt.Fprintf(warn, "appending %d sample(s)\n", len(frags))
		}

		res, err := gtstore.AttachRefBlocksToStore(ctx, base, frags, opts,
			buildinfo.Program(), buildinfo.CommandLine())
		if err != nil {
			return err
		}
		out := cmd.OutOrStdout()
		names := make([]string, 0, len(res.Volumes))
		for n := range res.Volumes {
			names = append(names, n)
		}
		sort.Strings(names)
		fmt.Fprintf(out, "attached %d fragment(s) to %d volume(s)\n", len(frags), len(names))
		for _, n := range names {
			r := res.Volumes[n]
			fmt.Fprintf(out, "  %-32s %12s blocks in %d file(s), %d sample(s)\n",
				n, comma(r.Rows), r.Shards, r.Samples)
		}
		for _, n := range res.Created {
			fmt.Fprintf(out, "  created blocks-only volume %s: a chromosome no volume held\n", n)
		}
		return nil
	},
}

// expandFragments turns arguments and a --from list into fragment directories.
// A local directory without a fragment.json is taken as a parent of fragments,
// which is how ten thousand of them arrive.
func expandFragments(args []string, from string) ([]string, error) {
	var out []string
	add := func(a string) error {
		if iosource.IsRemote(a) || hasFragmentMeta(a) {
			out = append(out, a)
			return nil
		}
		entries, err := os.ReadDir(a)
		if err != nil {
			return fmt.Errorf("%s is neither a fragment nor a directory of them: %w", a, err)
		}
		var found []string
		for _, e := range entries {
			if e.IsDir() && hasFragmentMeta(filepath.Join(a, e.Name())) {
				found = append(found, filepath.Join(a, e.Name()))
			}
		}
		if len(found) == 0 {
			return fmt.Errorf("%s holds no finished fragments", a)
		}
		sort.Strings(found)
		out = append(out, found...)
		return nil
	}
	for _, a := range args {
		if err := add(a); err != nil {
			return nil, err
		}
	}
	if from != "" {
		f, err := os.Open(from)
		if err != nil {
			return nil, err
		}
		defer f.Close()
		sc := bufio.NewScanner(f)
		for sc.Scan() {
			line := strings.TrimSpace(sc.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}
			out = append(out, line)
		}
		if err := sc.Err(); err != nil {
			return nil, err
		}
	}
	return out, nil
}

// planAppend splits fragments into those to attach and those whose sample is
// already in the table, and says which samples the attach replaces.
//
// "Already attached" is read from each volume's own manifest, not the archive's,
// and means attached in EVERY volume. An attach that died partway leaves volumes
// that disagree; a sample present in some of them is attached again everywhere,
// as a replace, which for the volumes that already had it swaps in the same
// blocks. So rerunning an interrupted --append finishes it.
//
// A table attached before the sample list was recorded cannot say who is in it,
// and is refused rather than guessed at: appending a sample already there would
// count its coverage twice.
func planAppend(ctx context.Context, base string, frags []gtstore.Fragment) (
	attach []gtstore.Fragment, replace, skipped []string, err error) {
	vols, err := volumeManifests(ctx, base)
	if err != nil {
		return nil, nil, nil, err
	}
	withBlocks := 0
	for _, v := range vols {
		if v.man.Params.RefBlocks {
			withBlocks++
		}
	}
	if withBlocks == 0 {
		return frags, nil, nil, nil // the first attach
	}
	if withBlocks < len(vols) {
		return nil, nil, nil, fmt.Errorf("%s: some volumes have ref_blocks and some do not; "+
			"--append cannot tell who is attached, so --rebuild from every fragment", base)
	}
	var attached map[string]bool
	for _, v := range vols {
		if v.man.RefBlockSamples == nil {
			return nil, nil, nil, fmt.Errorf("%s: ref_blocks were attached before the attached samples "+
				"were recorded, so --append cannot tell who is already in them; --rebuild once from "+
				"every fragment, and --append works from then on", v.name)
		}
		in := map[string]bool{}
		for _, s := range *v.man.RefBlockSamples {
			if attached == nil || attached[s] {
				in[s] = true
			}
		}
		attached = in
	}
	for _, f := range frags {
		if attached[f.Meta.SampleID] {
			skipped = append(skipped, f.Meta.SampleID)
			continue
		}
		attach = append(attach, f)
		replace = append(replace, f.Meta.SampleID)
	}
	return attach, replace, skipped, nil
}

type namedManifest struct {
	name string
	man  *gtstore.VolumeManifest
}

// volumeManifests reads the manifest of a single volume, or of every volume an
// archive lists.
func volumeManifests(ctx context.Context, base string) ([]namedManifest, error) {
	base = gtstore.TrimStoreSuffix(base)
	if !gtstore.IsStore(ctx, base) {
		m, err := gtstore.ReadVolumeManifestContext(ctx, base)
		if err != nil {
			return nil, err
		}
		return []namedManifest{{base, m}}, nil
	}
	st, err := gtstore.OpenStore(ctx, base)
	if err != nil {
		return nil, err
	}
	var out []namedManifest
	for _, v := range st.StoreManifest().Volumes {
		loc := v.Name
		if !strings.Contains(loc, "://") && !strings.HasPrefix(loc, "/") {
			loc = strings.TrimSuffix(base, "/") + "/" + loc
		}
		m, err := gtstore.ReadVolumeManifestContext(ctx, loc)
		if err != nil {
			return nil, fmt.Errorf("volume %s: %w", v.Name, err)
		}
		out = append(out, namedManifest{v.Name, m})
	}
	return out, nil
}

func listSome(xs []string, n int) string {
	if len(xs) <= n {
		return strings.Join(xs, ", ")
	}
	return strings.Join(xs[:n], ", ") + fmt.Sprintf(", ... (%d more)", len(xs)-n)
}

func hasFragmentMeta(dir string) bool {
	_, err := os.Stat(filepath.Join(dir, gtstore.FragmentMetaFile))
	return err == nil
}

func init() {
	f := gvcfAttachCmd.Flags()
	f.BoolVar(&attachRebuild, "rebuild", false, "Rebuild ref_blocks from these fragments alone")
	f.BoolVar(&attachAppend, "append", false, "Add samples not attached yet; skip fragments for samples already attached")
	f.StringSliceVar(&attachReplace, "replace", nil, "Samples to replace in the existing ref_blocks")
	f.StringVar(&attachFrom, "from", "", "Read fragment directories from FILE, one per line")
	f.IntVar(&attachBatch, "batch", gtstore.DefaultAttachBatch, "Fragments merged at once")
	f.StringVar(&attachTempDir, "tmp", "", "Directory for intermediate merge files")
}
