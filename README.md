# FinnGen_LDpred-2

An end-to-end **LDpred2** pipeline (WDL + [**bigsnpr**](https://privefl.github.io/bigsnpr/)): from GWAS
summary statistics to **per-individual polygenic scores (PGS)** on imputed `bgen` genotypes.

Two stages, each a standalone WDL, plus a top-level workflow that runs them together:

| Stage | Workflow | Does |
|---|---|---|
| **Infer** | `wdl/ldpred2_infer.wdl` | GWAS sumstats + LD reference → LDpred2-auto posterior weights (per-SNP β) |
| **Score** | `wdl/ldpred2_score.wdl` | weights + `bgen` → per-individual PGS (`Σ_j dosage_ij · β_j`) |
| **End-to-end** | `wdl/ldpred2_pipeline.wdl` | sumstats → weights → PGS in one submission |

```mermaid
flowchart LR
    SS["GWAS summary statistics<br/>(per trait)"] --> INF
    subgraph LDref ["LD reference"]
        direction TB
        LDa["in-sample: build_ld_reference.R<br/>(snp_cor over cohort bgen)"]
        LDb["or precomputed parts"]
    end
    LDref --> INF["ldpred2_infer<br/>snp_ldsc h2 → snp_ldpred2_auto"]
    INF --> W["harmonised weights<br/>bgen_id·chr·pos·a1_effect·a0·weight"]
    W --> SC["ldpred2_score<br/>scatter-gather over bgen chunks"]
    G["target cohort bgen"] --> SC
    SC --> PGS["per-individual PGS<br/>.sscore"]
```

The scoring stage is also useful on its own for any externally-inferred weights (LDpred2, LDAK
MegaPRS, PRS-CS, …). Both stages run from a single Docker image and are generic across chunked `bgen`
releases (built for FinnGen R13).

---

## Stage 1 — Inference (LDpred2-auto)

`ldpred2_infer.wdl` turns GWAS summary statistics into posterior per-SNP weights with
`bigsnpr::snp_ldpred2_auto` — the parameter-free "auto" model (no validation set needed). Per trait:
match sumstats to the LD reference, optional SD-based QC, estimate h² by LD-score regression
(`snp_ldsc`), run a multi-chain `snp_ldpred2_auto`, drop diverged chains, and average the kept chains'
effects into the final weights.

```mermaid
flowchart TD
    SS["sumstats (chr,pos,a1,a0,beta,beta_se,n_eff[,freq,info])"] --> M
    LD["LD reference parts<br/>(corr + map + LD scores, per chr)"] --> M["snp_match + subset<br/>assemble SFBM aligned to betas"]
    M --> QC["sd-based QC (if freq given)"]
    QC --> H["snp_ldsc → h² init"]
    H --> A["snp_ldpred2_auto<br/>(30 chains, multi p_init)"]
    A --> F["filter diverged chains → average β"]
    F --> W["weights.tsv<br/>bgen_id·chr·pos·a1_effect·a0·weight"]
```

### LD reference — two options

The LD reference is one or more **per-chromosome parts**, each an `.rds` holding
`list(corr = <sparse correlation>, map = data.frame(chr,pos,a0,a1), ld = <LD scores>)`. Parts store
plain sparse matrices (portable); the on-disk SFBM that LDpred2 needs is rebuilt inside the inference
task, so nothing with an absolute backing path crosses machines.

- **In-sample (default, recommended for a specific cohort):** the workflow scatters over the cohort's
  per-chromosome `bgen`s and runs `build_ld_reference.R` (`snp_cor` on a subsample of ~`n_ld`
  individuals, restricted to a supplied variant list, e.g. HapMap3+). Most accurate for the target
  ancestry (e.g. Finnish) and matches the scoring genotypes.
- **Precomputed:** pass `precomputed_ld_parts` (per-chr parts, or a one-file bundle from
  `combine_ld_reference.R` as a single-element array) and the build is skipped. Quicker/portable; use
  a reference matched to your ancestry.

### Inference inputs / outputs

Configure [`wdl/ldpred2_infer.inputs.json`](wdl/ldpred2_infer.inputs.json): `sumstats_files` (one per
trait), the LD-reference inputs, and the sumstats **column mapping** (`*_col` on the `infer` task;
defaults `chr,pos,a1,a0,beta,beta_se,n_eff`, with `a1` = effect allele; optional `freq_col` enables
SD-QC, `info_col` enables an INFO filter). Provide `n_eff` per variant (`n_col`) or a constant. Output:
one **`<trait>.weights.tsv`** (harmonised, feeds scoring) and a **`<trait>.report.tsv`** (h², p, chains
kept) per trait.

---

## Stage 2 — Scoring: how it works

For each `bgen` chunk (a scatter shard), `score_chunk` reads only the variants that overlap the
weight files, aligns alleles, computes a **partial** PGS per weight file, and writes one gzipped
partial. `gather_scores` then **sums the partials across all chunks** into one `.sscore` per weight
file. Because scoring is additive over disjoint variant sets, the sum of per-chunk partials equals the
genome-wide score.

```mermaid
flowchart TD
    subgraph inputs [Inputs]
        W["Harmonised weight files<br/>(bgen_id · chr · pos · a1_effect · a0 · weight)"]
        B["bgen chunk list<br/>(one gs:// path per line)"]
    end

    B --> SCAT{{"scatter over N bgen chunks"}}
    W --> SCAT

    subgraph shard ["score_chunk (per chunk, runs in parallel)"]
        R1["read .bgi index → variants in this chunk"]
        R2["canonical-key match weights ↔ chunk<br/>(chr:pos:sorted-alleles)"]
        R3["snp_readBGEN → dosage matrix<br/>(individuals × matched variants)"]
        R4["snp_match → align alleles, flip β"]
        R5["big_prodVec → partial PGS per weight file"]
        R1 --> R2 --> R3 --> R4 --> R5
    end

    SCAT --> shard
    R5 --> P["partial scores<br/>scorebase@@chunk.partial.tsv.gz"]
    P --> GATH["gather_scores<br/>sum partials across all chunks"]
    GATH --> OUT["Final per-individual PGS<br/>{lab}.{model}_{component}.sscore<br/>+ pgs_qc.tsv"]
```

### WDL structure

```mermaid
flowchart LR
    CL["bgen_chunk_list"] -->|read_lines| S
    SF["score_files[]"] --> S
    subgraph wf [workflow ldpred2_score]
        S["scatter: task score_chunk<br/>(cpu, mem, disk, preemptible)"] -->|"flatten(partial_scores)"| G["task gather_scores<br/>(preemptible = 0)"]
    end
    G --> O1["Array[File] sscores"]
    G --> O2["File qc (pgs_qc.tsv)"]
```

---

## Inputs

Edit [`wdl/ldpred2_score.inputs.json`](wdl/ldpred2_score.inputs.json):

| Input | Description |
|---|---|
| `ldpred2_score.docker` | The scoring image (built from `docker/`). |
| `ldpred2_score.bgen_chunk_list` | A text file with **one `gs://…bgen` path per line**. Each chunk must have sibling `.bgen.bgi` and `.bgen.sample` files (localized automatically). |
| `ldpred2_score.score_files` | Array of **harmonised weight files** (see format below), one per PGS to compute. |

### Harmonised weight-file format

Tab-separated, one row per SNP, with a header. One file per PGS.

```
bgen_id          chr   pos      a1_effect   a0   weight
1_785910_G_C     1     785910   C           G    1.8585e-04
```

| Column | Meaning |
|---|---|
| `bgen_id` | variant id (informational; matching is done on chr/pos/alleles) |
| `chr`, `pos` | chromosome (int or `chr`-prefixed) and base-pair position, on the **same build as the bgen** |
| `a1_effect` | effect allele (the `weight`/β is expressed w.r.t. this allele) |
| `a0` | other allele |
| `weight` | per-SNP weight / posterior effect (β) |

Producing these harmonised files (rsID → bgen coordinate/allele mapping) is an **upstream step**,
outside this repo.

## Outputs

- `sscores` — one `.sscore` per weight file, columns `FINNGENID` and `PGS` (the genome-wide score).
- `qc` — `pgs_qc.tsv`: per weight file, `n_individuals, mean, sd, n_chunks_summed`.

---

## Repository layout

```
FinnGen_LDpred-2/
├── wdl/
│   ├── ldpred2_pipeline.wdl             # end-to-end: infer → score (imports the two below)
│   ├── ldpred2_infer.wdl                # inference: (build LD) → LDpred2-auto weights
│   ├── ldpred2_score.wdl               # scoring: scatter score_chunk → gather_scores
│   ├── ldpred2_pipeline.inputs.json     # inputs templates …
│   ├── ldpred2_infer.inputs.json
│   ├── ldpred2_score.inputs.json
│   └── cromwell_workflow_options.json   # labels + call-caching
├── scripts/
│   ├── build_ld_reference.R             # in-sample LD (snp_cor) for one chromosome
│   ├── combine_ld_reference.R           # (optional) bundle per-chr LD parts into one file
│   ├── ldpred2_infer.R                  # snp_ldsc h2 + snp_ldpred2_auto → weights
│   ├── score_chunk.R                    # per-chunk bigsnpr scoring
│   └── gather_scores.R                  # sum partials → final .sscore
├── docker/
│   ├── Dockerfile                       # R + bigsnpr + deps (one image for all stages)
│   └── cloudbuild.yaml                  # Cloud Build config
├── submit_ldpred2_score.sh              # example Cromwell submission helper
├── requirements.txt
└── LICENSE
```

---

## Usage

### 1. Build the image

No local Docker is required — build with Cloud Build (see [`docker/cloudbuild.yaml`](docker/cloudbuild.yaml)):

```bash
gcloud builds submit . \
  --config docker/cloudbuild.yaml \
  --substitutions=_TAG=0.1,_IMAGE=REGION-docker.pkg.dev/PROJECT/REGISTRY/ldpred2 \
  --project YOUR_GCP_PROJECT --region YOUR_REGION \
  --gcs-source-staging-dir gs://YOUR_STAGING_BUCKET/cloud_build_staging
```

> If your GCP org enforces `constraints/gcp.resourceLocations`, the `--gcs-source-staging-dir` flag
> (pointing at a bucket in an allowed region) is **required** — otherwise Cloud Build's default US
> staging bucket triggers `HTTPError 412: 'us' violates constraint`.

### 2. Configure inputs

Fill the relevant inputs template (image tag, buckets, columns) and your labels in
`wdl/cromwell_workflow_options.json`.

### 3. Submit

The image is shared by all stages. Run whichever entry point you need:

```bash
# End-to-end: sumstats → weights → PGS (imports the two sub-workflows, so zip them as deps)
( cd wdl && zip ldpred2_deps.zip ldpred2_infer.wdl ldpred2_score.wdl )
cromwell --port 5000 submit --wdl wdl/ldpred2_pipeline.wdl \
  --inputs wdl/ldpred2_pipeline.inputs.json --deps wdl/ldpred2_deps.zip \
  --options wdl/cromwell_workflow_options.json

# Or a single stage on its own (no --deps needed):
cromwell --port 5000 submit --wdl wdl/ldpred2_infer.wdl --inputs wdl/ldpred2_infer.inputs.json
./submit_ldpred2_score.sh          # thin wrapper for the scoring stage
```

Adapt to your Cromwell submitter (Cromshell, REST API, Terra, or
[CromwellInteract](https://github.com/FINNGEN/CromwellInteract)).

### 4. Collect outputs

Retrieve the per-trait `weights.tsv` / `report.tsv` and the per-individual `.sscore` files from your
Cromwell execution bucket once the workflow succeeds.

---

## ⚠️ Validation before trusting scores

This pipeline is numerically silent on a few points — validate them once for your genotype release
before using the scores downstream:

1. **Allele/dosage orientation (highest priority).** The score sign depends on which allele
   `snp_readBGEN` counts. The code assumes it returns the dosage of the **second** allele of each
   `<chr>_<pos>_<a1>_<a2>` id (the documented `bgen` behaviour) and orients `snp_match` accordingly. If
   your `bgen`/bigsnpr version counts the first allele, every PGS would be **sign-inverted with no
   error**. Confirm on ~50 individuals and one weight file against an independent scorer, e.g.
   `plink2 --bgen … ref-first --score <weights> …`, and check both **sign and magnitude** match.
2. **Non-autosomal chunks.** `chr_to_int` maps X/Y/XY/MT → 23–26, but bigsnpr's `format_snp_id` may
   reject non-1..22 codes. If your chunk list includes sex/MT chromosomes, test one such chunk (it
   would crash, not mis-score).
3. **Sample-id column.** The `FINNGENID` output is read from column 2 (`ID_2`) of the Oxford
   `.sample` file. Confirm that is your intended sample identifier.

The built-in guards that already prevent *silent* wrongness: `snp_match(strand_flip = FALSE)` keeps
ambiguous SNPs (harmonised data is same-strand); `gather_scores` **refuses to write** if any weight
file is missing chunks (would under-sum); and `score_chunk` aborts if two weight files share a
basename (would collide).

## Implementation notes & gotchas

These are the non-obvious pitfalls this pipeline explicitly handles — useful if you adapt it:

- **`bgen` variant-id / chromosome convention.** bigsnpr `snp_readBGEN` matches `list_snp_id` in the
  form `<chr>_<pos>_<a1>_<a2>` and its `format_snp_id()` requires a **numeric** chromosome. Some `bgen`
  indexes store `chr1` (string). `score_chunk.R` rewrites the localized `.bgi` `chromosome` column to a
  bare integer in place (only the text column; byte offsets untouched) and builds matching numeric
  ids. Variant matching itself uses an order-independent canonical key (`chr:pos:{sorted alleles}`) so
  it is robust to ref/alt ordering and `chr`-prefix differences.
- **Allele orientation.** `snp_match` reconciles effect-allele orientation and flips β for
  reversed/complementary alleles, so the score is not silently inverted; ambiguous (strand) SNPs are
  dropped.
- **Single-threaded BLAS.** bigstatsr aborts with *"Two levels of parallelism are used"* when
  `ncores > 1` runs on a multi-threaded BLAS. The workflow exports `OPENBLAS_NUM_THREADS=1` (etc.) and
  runs `--ncores 1`; per-chunk scoring is cheap, so this is not a bottleneck.
- **Compression.** The image's `data.table` is built without zlib, so partials are written plain and
  compressed with the `gzip` binary; `gather_scores` reads them via `gzip -dc`.
- **Dependencies.** `bigsnpr::snp_readBGEN` needs `dbplyr` at runtime — it is installed in the image.
- **`gather` is not preemptible.** It is a single reduce over all shards; `preemptible = 0` avoids
  losing it late.

---

## Requirements

See [`requirements.txt`](requirements.txt). Everything runs in the Docker image; the submitting host
only needs the Google Cloud SDK and a WDL runner.

## License

[MIT](LICENSE).
