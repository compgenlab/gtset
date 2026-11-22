# Quick start: from a cohort VCF (and, optionally, gVCFs) to a gtset

This walks through building one gtset from a joint-called cohort VCF, then
optionally adding each sample's gVCF reference blocks. Every command here is
real; `gtset <command> --help` has the full flags, and
[convert.md](convert.md) is the reference for what a store is and why.

- **Part 1** builds the store from the cohort VCF. That is a complete, usable
  store on its own.
- **Part 2 is optional.** It adds per-sample gVCF reference blocks so the store
  can also answer at positions the cohort VCF never mentions. You can skip it,
  or come back to it at any time later without reconverting anything.

## What you need

- The joint-called cohort VCF (pVCF), bgzipped. One file per chromosome is
  normal and is what the commands below assume; tabix indexes are recommended.
- *(Part 2 only)* one gVCF per sample, for samples in the cohort VCF.
- A `gtset` binary: `make build` writes `bin/gtset.<os>_<arch>`.

---

## Part 1: build the store

### 1. Convert each chromosome to a volume

```bash
for c in chr{1..22} chrX chrY; do
  gtset convert -v --out cohort/$c \
      --meta-dataset "MyCohort release 3" --meta-reference GRCh38 \
      calls/$c.vcf.gz 2> cohort/$c.log
done
```

With `-v`, a conversion reports on stderr about every five minutes (when it
has moved), so a cluster job's log shows where it is and roughly how long is
left:

```
  chr2: starting at 10,177 of 242,193,529
  chr2:121,450,332 / 242,193,529 (50.1%)  records 61,204 (34/s)  genotypes 1.84G (1.0M/s)  calls 2,310,442  elapsed 30m0s  eta ~30m
```

When the input has a tabix index that records how many records each
chromosome holds (`tabix` and `bcftools index` write this), the percentage and
estimate count records instead, which tracks the actual work more closely:

```
  chr2:121,450,332  record 61,204 of 122,310 (50.0%)  34/s  genotypes 1.84G (1.0M/s)  calls 2,310,442  elapsed 30m0s  eta ~30m
```

Otherwise they come from the position against the chromosome's length, which
assumes the rest of the chromosome is as dense as what came before, so treat it
as rough. The summary at the end gives the overall rate,
which is the number to plan the next batch with.

Each chromosome becomes a **volume**: a directory holding the `calls`,
`catalog` and `catalog_runs` tables and a `manifest.json.gz` written last, so a
volume either finished or is unreadable. The chromosomes are independent jobs,
so run them in parallel on a cluster if you have one.

Worth deciding now, because every volume (and every gVCF fragment later) must
agree on them:

- `--min-dp` (default 10): the depth at which a sample counts as called. A
  `0/0` the store reports later is one observed at or above this depth.
- Depth bands (default: geometric from `--min-dp`, ratio 1.5). Leave the
  default unless you have a reason.

Also useful:

- `--meta-dataset`, `--meta-reference`: what the store *is*. Nothing else in a
  store records the assembly, and a store read against the wrong one does not
  fail, it answers with the wrong coordinates.
- `--shard-sites N`: split each volume into coordinate shards, so a locus
  query reads one small file instead of a large footer. Worth it for big
  cohorts (e.g. `--shard-sites 100000`).
- `--info AF,DR2`: keep INFO fields as columns on the site catalog.

Memory grows with the number of samples (the width of a VCF line), not with
the size of the file.

### 2. Combine the volumes into one store

```bash
gtset bind cohort/
gtset summary cohort
```

`bind` writes `cohort/gtset.json.gz` over every volume it finds there; no
genotype is read or rewritten, so it takes seconds.
It refuses volumes that disagree on samples (or their order), `--min-dp` or
depth bands, or that share a chromosome. `summary` reads only manifests, so it
is instant; check that every chromosome is there.

### 3. Query it

```bash
# Who carries a variant
gtset query cohort --variant chr17:43094464:G:A

# ...and who is confidently reference (0/0) there
gtset query cohort --variant chr17:43094464:G:A --hom-ref

# What one subject carries in a region
gtset query cohort --sample S0001 --variant chr17:43000000-43200000

# Many sites at once (VCF, BED or "chrom pos [ref alt]" list): one pass, not one per site
gtset query cohort --variant panel.bed --sample S0001 --sample S0002

# A genotype matrix as VCF
gtset query cohort --variant panel.bed --format vcf -o panel.vcf.gz
```

Every row has a `ref_source` column. For a `0/0` it says what the reference
call rests on: `catalog_run` means the cohort VCF called the sample there.

**Without Part 2, a position that is not in the cohort VCF is not assayed for
everyone.** A plain VCF only reports variant sites, so the store cannot claim
a sample was reference anywhere else. That is the gap Part 2 fills.

---

## Part 2 (optional): add gVCF reference blocks

Do this if you want answers at positions the cohort VCF does not contain: a
variant private to another batch, one found after the cohort was called, a
clinical panel position nobody in this cohort happened to carry. Skipping it
changes nothing at sites the cohort VCF does contain, which always answer from
the cohort VCF. You can add gVCFs in batches, as they arrive.

### 4. Extract each sample's reference blocks into a fragment

```bash
gtset gvcf-extract gvcf/S0001.g.vcf.gz --roster cohort --out frags/S0001
gtset gvcf-extract gvcf/S0002.g.vcf.gz --roster cohort --out frags/S0002
# ... one job per sample, in parallel
```

A **fragment** holds one sample's reference blocks: one Parquet file per
chromosome plus `fragment.json`, written last as the completion marker.
`--roster` is checked first: the sample must be in the store, and `--min-dp`
and the depth bands default to the store's and are refused if they differ, so
an unusable fragment fails in seconds rather than at attach time.

Blocks below `--min-dp` are dropped, and abutting blocks in the same depth band
are merged. Nothing is merged **across** a gap or a variant record: a gap is
bases the gVCF did not vouch for, and bridging it would report `0/0` where
nobody looked. A finished fragment is kept on a rerun unless `--force`.

### 5. Attach the fragments to the store

```bash
gtset gvcf-attach cohort/ frags/
```

`frags/` is a directory of fragments; you can also name fragments one by one,
or list them in a file with `--from fragments.txt` when there are thousands.
The merge runs in batches of `--batch` (default 100) through intermediate files
in `--tmp`, so point `--tmp` at a disk with room. Each volume's manifest is
rewritten last, so an interrupted attach leaves the store as it was; rerun it.

Chromosomes with blocks but no cohort VCF site (chrM, unplaced contigs) get a
blocks-only volume, so they are answerable too.

### 6. Add the next batch of gVCFs later

Once a store has blocks, an attach must say what it is doing. To add samples
as their gVCFs arrive, use `--append` with the same fragment directory:

```bash
gtset gvcf-extract gvcf/S0011.g.vcf.gz --roster cohort --out frags/S0011
# ... S0012 to S0020
gtset gvcf-attach --append cohort/ frags/
```

`--append` attaches only the samples not in the store yet and skips fragments
for samples already attached, saying which, so you never have to list who is
new. Rerunning it is safe, and it finishes an attach that died partway. On a
store with no blocks yet it is the first attach, so a job script can use
`--append` every time.

To redo a sample already attached (its gVCF changed), name it with
`--replace`, which swaps in the new fragment and keeps everyone else:

```bash
gtset gvcf-attach cohort/ frags/S0003 --replace S0003
```

To rebuild the whole table from fragments instead, pass every fragment with
`--rebuild`.

Each attach rewrites the block files, since they are sorted by position across
all samples, so attaching in larger batches is cheaper than one sample at a
time.

### 7. Check and query

```bash
gtset summary cohort
gtset query cohort --variant chr1:123456:C:T --hom-ref
```

`summary` now reports `ref blocks yes` and **whose gVCFs are attached**
(`from 9,800 of 10,000 samples; not attached: ...`). Read that line: a sample
whose gVCF was never attached looks exactly like one whose gVCF did not cover a
position. Both are not assayed.

Off-catalog answers come only for an **exact locus** (`chrom:pos:ref:alt`)
with `--hom-ref`, and each such `0/0` has `ref_source` `ref_block` and
`min_dp` set to the block's depth floor:

```
note: chr1:123456:C:T is not in the catalog; samples are reported from reference blocks
      where a block covers it, and are otherwise not assayed.
chrom  pos     ref alt sample gt  dp min_dp ad_ref ad_alt gq ref_source
chr1   123456  C   T   S0001  0/0 .  25     .      .      .  ref_block
chr1   123456  C   T   S0002  0/0 .  25     .      .      .  ref_block
```

A region query never lists off-catalog positions base by base.

---

## What you end up with

```
cohort/
  gtset.json.gz            the gtset: roster, shared settings, volume list
  chr1/
    manifest.json.gz       the volume's completion marker and row counts
    calls.parquet          one row per ALT-carrying genotype
    catalog.parquet        one row per site in the cohort VCF
    catalog_runs.parquet   runs of sites each sample was called at
    ref_blocks/...         Part 2 only: the attached gVCF blocks
  chr2/
  ...
frags/                     Part 2 only; keep it if you may --rebuild
  S0001/
    fragment.json
    00000.parquet ...
```

## See also

- [convert.md](convert.md): the format, its guarantees, and every flag.
- `gtset <command> --help` for each step.
