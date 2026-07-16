## ldpred2_score.wdl
##
## Pure SCORING application of precomputed per-SNP PGS weight files to imputed
## bgen genotypes (e.g. FinnGen R13), using the bigsnpr / LDpred2 toolchain as a
## dot-product scoring engine only. NO weight re-inference, NO LD matrix.
##
## WDL syntax: draft-2 (matches regenie-pipelines/wdl/gwas/regenie_sub.wdl -
## no `version` line, `${var}` interpolation, `command <<< >>>`, task-level
## input declarations). Do NOT add a `version 1.0` header or this file will be
## parsed under the wrong WDL version.
##
## Scatter model: one shard per R13 bgen chunk (535 chunks). Each shard scores
## every weight file against the variants present in that chunk and emits one
## partial-PGS file per weight file. gather_scores sums the partials across all
## chunks per weight file to produce the final per-individual .sscore files.

task score_chunk {

    # One bgen chunk plus its sibling index and sample files. Localized the same
    # way regenie_sub.wdl does (File bgi = bgen + ".bgi", sample = bgen + ".sample").
    File bgen
    File bgi = bgen + ".bgi"
    File sample = bgen + ".sample"

    # All harmonised score files (columns: bgen_id chr pos a1_effect a0 weight).
    Array[File] score_files

    # Column names in the harmonised score files (override for the raw-fallback path).
    String id_col = "bgen_id"
    String a1_col = "a1_effect"
    String a0_col = "a0"
    String weight_col = "weight"

    # Chunk tag used to make partial-score filenames unique across shards.
    String chunk_tag = basename(bgen, ".bgen")

    String docker
    Int cpu = 2
    String mem = "13 GB"
    # NOTE: this Cromwell's WDL `size()` accepts only a single File, not Array[File],
    # so we size off the bgen chunk (~2 GB) + a fixed buffer covering the bigsnpr FBM
    # backing file, the (small) score files, and the gzipped partial outputs.
    Int disk_gb = ceil(size(bgen, "G")) + 25
    Int preemptible = 2

    command <<<
        set -euxo pipefail

        # Force single-threaded BLAS at the process level. bigstatsr's big_prodVec
        # aborts with "Two levels of parallelism are used" when ncores>1 runs on top
        # of a multi-threaded BLAS; the in-script set_blas_ncores(1) proved unreliable
        # across VMs (run 6: 18/535 shards passed, the rest died here). Pinning BLAS
        # threads via env vars is deterministic, and --ncores 1 makes assert_cores a
        # no-op regardless. Scoring one ~520k x ~2k chunk is trivial single-threaded.
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
        export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OPENBLAS_MAIN_FREE=1

        Rscript /scripts/score_chunk.R \
            --bgen ${bgen} \
            --sample ${sample} \
            --score_files ${sep="," score_files} \
            --id_col ${id_col} \
            --a1_col ${a1_col} \
            --a0_col ${a0_col} \
            --weight_col ${weight_col} \
            --chunk_tag ${chunk_tag} \
            --ncores 1

        echo "partial score files produced:"
        ls -1 *.partial.tsv.gz
    >>>

    output {
        # One gzipped partial-score file per weight file, named
        # <scorebase>@@<chunk_tag>.partial.tsv.gz (FINNGENID<TAB>partial_PGS).
        # Gzipped to keep the gather-stage localization tractable (~520k rows/file
        # x 32 files x 535 chunks).
        Array[File] partial_scores = glob("*.partial.tsv.gz")
    }

    runtime {
        docker: "${docker}"
        cpu: cpu
        memory: mem
        disks: "local-disk " + disk_gb + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: preemptible
        noAddress: true
    }
}

task gather_scores {

    # All partial-score files from all chunks, flattened.
    Array[File] partials
    # Manifest of the partial paths (one per line); passed to the R script so we
    # avoid a giant command line for 535 chunks x 32 weight files.
    File partial_manifest = write_lines(partials)

    String docker
    Int cpu = 1
    String mem = "8 GB"
    # `size(Array[File])` is unsupported here; the gzipped partials total ~60-85 GB
    # (weight_files x chunks x ~5 MB). Fixed generous default (overridable via inputs).
    Int disk_gb = 150
    # Single reduce over all shards - not worth losing to preemption.
    Int preemptible = 0

    command <<<
        set -euxo pipefail

        mkdir -p sscore

        Rscript /scripts/gather_scores.R \
            --partials ${partial_manifest} \
            --out_dir sscore

        echo "final sscore files:"
        ls -1 sscore/*.sscore
        echo "QC:"
        cat sscore/pgs_qc.tsv
    >>>

    output {
        # 32 final per-individual PGS files (FINNGENID<TAB>PGS).
        Array[File] sscores = glob("sscore/*.sscore")
        File qc = "sscore/pgs_qc.tsv"
    }

    runtime {
        docker: "${docker}"
        cpu: cpu
        memory: mem
        disks: "local-disk " + disk_gb + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: preemptible
        noAddress: true
    }
}

workflow ldpred2_score {

    String docker

    # File with one gs:// bgen chunk path per line (e.g. the R13 535-chunk list).
    File bgen_chunk_list
    Array[String] chunks = read_lines(bgen_chunk_list)

    # Harmonised weight files (bgen_id chr pos a1_effect a0 weight), one per PGS.
    Array[File] score_files

    scatter (chunk in chunks) {
        call score_chunk {
            input:
                docker = docker,
                bgen = chunk,
                score_files = score_files
        }
    }

    # score_chunk.partial_scores is Array[Array[File]] (chunks x weight files);
    # flatten to a single Array[File] and sum per weight file in one gather task.
    Array[File] all_partials = flatten(score_chunk.partial_scores)

    call gather_scores {
        input:
            docker = docker,
            partials = all_partials
    }

    output {
        Array[File] sscores = gather_scores.sscores
        File qc = gather_scores.qc
    }
}
