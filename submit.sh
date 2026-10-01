#!/bin/bash
#
# Submit a pipeline step to SLURM from the project root: logs go to logs/<step>.<jobid>.{o,e},
# the root is passed to the job as $BUNDLE (so it can source config.sh), and the submission is
# appended to logs/submissions.tsv. Extra args pass through to sbatch.
#
#   ./submit.sh scripts/04_relatedness.sh --dependency=afterok:12345678
#
# Human-run on biowulf (login + 2FA are yours).

set -euo pipefail

BUNDLE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <scripts/NN_step.sh> [sbatch args...]" >&2
    exit 2
fi

SCRIPT="$1"; shift
[[ -f "$SCRIPT" ]] || SCRIPT="${BUNDLE}/${SCRIPT}"
[[ -f "$SCRIPT" ]] || { echo "No such script: $1" >&2; exit 1; }

STEP="$(basename "$SCRIPT")"; STEP="${STEP%.sh}"
LOG_DIR="${BUNDLE}/logs"
mkdir -p "$LOG_DIR"

# --job-name  names the logs too (steps 1-2 run once per callset and would otherwise collide).
# --export    sbatch keeps only the LAST --export, so a caller's would drop ALL,BUNDLE=; merge.
EXPORTS="ALL,BUNDLE=${BUNDLE}"
ARGS=()
for a in "$@"; do
    case "$a" in
        --job-name=*)
            STEP="${a#--job-name=}"
            ARGS+=("$a")
            ;;
        --export=*)
            v="${a#--export=}"
            v="${v#ALL,}"; [[ "$v" == "ALL" ]] && v=""      # ALL is already in EXPORTS
            [[ -n "$v" ]] && EXPORTS="${EXPORTS},${v}"
            ;;
        *)
            ARGS+=("$a")
            ;;
    esac
done

# %j keeps every attempt's logs; without it a retry overwrites the evidence of the last failure.
JOBID=$(sbatch --parsable \
    --job-name="$STEP" \
    --output="${LOG_DIR}/${STEP}.%j.o" \
    --error="${LOG_DIR}/${STEP}.%j.e" \
    --export="$EXPORTS" \
    ${ARGS+"${ARGS[@]}"} "$SCRIPT")
JOBID="${JOBID%%;*}"

# ── the run log ──────────────────────────────────────────────────────────────
# Append-only: sacct has each job's outcome but not its --export values, and it ages out.
# scripts/runlog.sh joins the two.
SUBLOG="${LOG_DIR}/submissions.tsv"
if [[ ! -f "$SUBLOG" ]]; then
    printf 'submitted_utc\tjobid\tstep\tscript\tsbatch_args\texports\tsubmitted_by\n' > "$SUBLOG"
fi
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$JOBID" "$STEP" "${SCRIPT#$BUNDLE/}" \
    "${ARGS[*]:-}" "${EXPORTS#ALL,BUNDLE=$BUNDLE,}" "${USER:-unknown}" >> "$SUBLOG"

echo "Submitted ${STEP} as job ${JOBID}"
echo
echo "  watch:  squeue -j ${JOBID}"
echo "  out:    tail -f ${LOG_DIR}/${STEP}.${JOBID}.o"
echo "  err:    tail -f ${LOG_DIR}/${STEP}.${JOBID}.e"
echo "  chain:  ./submit.sh scripts/<next>.sh --dependency=afterok:${JOBID}"
