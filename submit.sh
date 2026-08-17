#!/bin/bash
#
# Submit a pipeline step to SLURM from the bundle root.
#
# This exists for one reason: so no script has to hardcode a log path. It resolves the
# bundle location, creates logs/, points --output/--error at logs/<step>.<jobid>.{o,e}, and
# hands the bundle root to the job as $BUNDLE so the script can source config.sh.
#
#   ./submit.sh scripts/03_merge.sh
#   ./submit.sh scripts/02_normalize.sh --job-name=norm_wgs --export=TAG=wgs,...
#   ./submit.sh scripts/04_relatedness.sh --dependency=afterok:12345678
#
# Anything after the script name is passed through to sbatch verbatim, so an explicit
# --output=/--error= of your own still wins.
#
# You run this, not the AI (login + 2FA are yours). Plain `sbatch` works too — see README §3.

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

# Walk the passthrough args once:
#   --job-name  steps 1 and 2 run once per callset, so a per-step log name would collide across
#               the three jobs. When you pass --job-name, the logs take that name instead.
#   --export    sbatch takes the LAST --export and discards earlier ones, so a caller's
#               --export=PFILE=...  would silently drop ALL,BUNDLE= and the job would lose both
#               its environment and the path to config.sh. Merge into one directive instead.
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

# %j is the job id, and it is the whole point of this naming. Without it every retry of a
# step overwrites the previous attempt's logs, so a failure that is only reproduced on some
# nodes destroys its own evidence each time you chase it. That is exactly what happened to
# genotools_br_dsnwgs: three failures (26843275, 27057269, 27103380) left nothing to read,
# MaxRSS was blank, and the cause was invisible until a fourth run was finally caught. Logs
# are the cheapest thing in the project; keep all of them.
JOBID=$(sbatch --parsable \
    --job-name="$STEP" \
    --output="${LOG_DIR}/${STEP}.%j.o" \
    --error="${LOG_DIR}/${STEP}.%j.e" \
    --export="$EXPORTS" \
    ${ARGS+"${ARGS[@]}"} "$SCRIPT")
JOBID="${JOBID%%;*}"

# ── the run log ──────────────────────────────────────────────────────────────
# Append-only record of what was submitted, with the arguments. sacct knows the OUTCOME of
# every job but not the --export values that produced it, and it ages out; this keeps the
# intent. Written here rather than maintained by hand so it cannot drift from what ran —
# a log that is only sometimes right is worse than none. scripts/runlog.sh joins it to sacct.
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
