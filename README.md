# gtset

A sparse, queryable genotype store for cohort VCFs, and the command-line tool
that builds and queries it. A **gtset** is a set of
**volumes** (normally one chromosome each), each optionally split into
coordinate **shards**.

A store keeps only ALT-carrying genotypes plus enough context to tell a
confidently-called reference genotype (`0/0`) from a position that was never
assayed. Queries against a store and against the VCF it came from return the same
answer; the store is just faster to ask. gVCF reference blocks can be attached
per sample, so a store can also answer at positions no variant record mentions.

The format and its query API live in
[`cghts/gtstore`](https://github.com/compgenlab/cghts); this repo is the CLI.
The project was called gtstore (and briefly Pendium), and the format kept the
gtstore name on disk: a
gtset's top-level marker is `gtstore.json.gz` and its metadata keys are
`gtstore.*`, so stores written before the rename still open.
General VCF tooling is in [`cgkit`](https://github.com/compgenlab/cgkit).

## Building

```bash
make build     # bin/gtset.{darwin,linux}_{arm64,amd64}, windows_amd64.exe
make test
```

## Pipeline

Step by step, with the gVCF part optional and adding gVCFs in batches:
[docs/quickstart.md](docs/quickstart.md).

```bash
# 1. Convert the joint-called pVCF, one volume per chromosome
gtset convert --out cohort/chr1 chr1.vcf.gz
gtset convert --out cohort/chr2 chr2.vcf.gz

# 2. Combine the volumes into one store
gtset bind cohort/

# 3. (optional) Reference blocks from each sample's gVCF -- one job per sample
gtset gvcf-extract --roster cohort --out frags/S0001 S0001.g.vcf.gz

# 4. (optional) Attach them; rerun with --append as more gVCFs arrive,
#    and --replace S0001 swaps in a redone sample
gtset gvcf-attach --append cohort/ frags/

# 5. Query: sites, samples, or both
gtset query cohort --variant chr1:1000:A:T
gtset query cohort --variant panel.bed --sample S0001 --hom-ref
gtset query cohort.vcf.gz --variant chr1:1000:A:T   # same answer, from the VCF

# Describe a store or a VCF without reading genotypes
gtset summary cohort
```

Run `gtset <command> --help` for flags. `--profile=FILE` writes a CPU profile
from anywhere on the command line.

## Commands

| Command | Description |
|---------|-------------|
| `convert` | Convert one or more VCFs into a store volume ([format docs](docs/convert.md)) |
| `bind` | Bind volumes, disjoint by chromosome, into one store |
| `gvcf-extract` | Extract one sample's gVCF reference blocks into a fragment |
| `gvcf-attach` | Merge fragments into a store's `ref_blocks` table |
| `query` | Query genotypes by site, by sample, or both, from a store or a VCF |
| `summary` | Describe a store or VCF: samples, contigs, provenance, census |

## Remote inputs

Inputs may be local paths, `http(s)://` URLs or `s3://` objects, and index
sidecars resolve over the same transport, so an indexed query still seeks.
Outputs are local, except that `convert --out s3://bucket/store` writes a store
straight to S3. S3 credentials come from the standard AWS chain;
`AWS_ENDPOINT_URL` points at an S3-compatible gateway.
