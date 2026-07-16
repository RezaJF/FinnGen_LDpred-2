#!/usr/bin/env bash
# Submit the ldpred2_score PGS-scoring workflow to a Cromwell server.
#
# This is an example helper for a Cromwell deployment driven by CromwellInteract
# (https://github.com/FINNGEN/CromwellInteract) over an SSH SOCKS tunnel. Adapt
# it to your own Cromwell submission mechanism (Cromshell, the REST API, Terra, …).
#
# Prerequisites:
#   1. Build + push the scoring image (see docker/cloudbuild.yaml) and set the
#      image tag in wdl/ldpred2_score.inputs.json ("ldpred2_score.docker").
#   2. Fill wdl/ldpred2_score.inputs.json with your bgen chunk list + harmonised
#      score-file paths, and wdl/cromwell_workflow_options.json with your labels.
#   3. Open a tunnel / authenticate to your Cromwell server.
#
# Configurable via environment:
#   CROMWELL_INTERACT   path to cromwell_interact.py   (default: on PATH)
#   CROMWELL_PORT       local SOCKS/API port           (default: 5000)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WDL="${SCRIPT_DIR}/wdl/ldpred2_score.wdl"
INPUTS="${SCRIPT_DIR}/wdl/ldpred2_score.inputs.json"
OPTIONS="${SCRIPT_DIR}/wdl/cromwell_workflow_options.json"
PORT="${CROMWELL_PORT:-5000}"
CROMWELL_INTERACT="${CROMWELL_INTERACT:-cromwell_interact.py}"

echo "Submitting ldpred2_score (PGS scoring)"
echo "  WDL:     ${WDL}"
echo "  Inputs:  ${INPUTS}"
echo "  Options: ${OPTIONS}"
echo "  Port:    ${PORT}"
echo ""

# --options carries the (optional) google_labels; this workflow has no
# subworkflow imports, so no --deps zip is required.
python3 "${CROMWELL_INTERACT}" \
    --port "${PORT}" submit \
    --wdl "${WDL}" \
    --inputs "${INPUTS}" \
    --options "${OPTIONS}"

echo ""
echo "Track:    ${CROMWELL_INTERACT} --port ${PORT} meta <workflow-id> -s"
echo "Outputs:  ${CROMWELL_INTERACT} --port ${PORT} outfiles <workflow-id>"
