#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --job-name=mkcomplex

# Assemble one design's multi-chain complex sequence per array task, and write a
# per-design TSV (name, status, sequence, n_chains, chain_lens,
# design_chain_index, total_len) for collect_mkcomplex.py to merge back.
#
# Pure string work -- fixed chains from literals or a FASTA, the design's own
# sequence in the middle, joined by the chainbreak separator -- so there is no
# worker module and no dependency beyond bash + awk. CPU-only.

set -euo pipefail

# Shared scaffolding: sets MANIFEST/OUT_DIR/SAPIA_TASK_ID/SAPIA_LINE.
source "${SAPIA_PRELUDE:?}"

# Site-specific activation — see docs/configuration.md. This tool needs nothing
# but bash and awk, so the hook is OPTIONAL here (unlike every other tool's .sh):
# it is sourced when SAPIA_ACTIVATE_MKCOMPLEX is set, and skipped when it is not,
# rather than failing a task that has no environment to activate. Under modal
# sapia_activate is a no-op anyway.
if [[ -n "${SAPIA_ACTIVATE_MKCOMPLEX:-}" ]]; then
    sapia_activate SAPIA_ACTIVATE_MKCOMPLEX
fi

NAME=$(echo "$SAPIA_LINE" | cut -f1)
SEQ=$(echo "$SAPIA_LINE" | cut -f2)
PREPEND_SEQS=$(echo "$SAPIA_LINE" | cut -f3)
PREPEND_FASTA=$(echo "$SAPIA_LINE" | cut -f4)
APPEND_SEQS=$(echo "$SAPIA_LINE" | cut -f5)
APPEND_FASTA=$(echo "$SAPIA_LINE" | cut -f6)
SEP=$(echo "$SAPIA_LINE" | cut -f7)
REPEAT=$(echo "$SAPIA_LINE" | cut -f8)

RESULT_TSV="$OUT_DIR/${NAME}.tsv"

# The 20 standard one-letter amino-acid codes, uppercase. A stray FASTA header or
# a lowercase/masked residue slipping into a target chain would silently corrupt
# every downstream prediction, so it is refused here instead.
AA_RE='^[ACDEFGHIKLMNPQRSTVWY]+$'

write_result() { # status sequence n_chains chain_lens design_chain_index total_len
    printf 'name\tstatus\tsequence\tn_chains\tchain_lens\tdesign_chain_index\ttotal_len\n' \
        >"$RESULT_TSV"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$NAME" "$1" "$2" "$3" "$4" "$5" "$6" >>"$RESULT_TSV"
}

# Errors are recorded as DATA (status 'error: ...') and the task still exits 0, so
# one bad row never fails the array and partial runs still collect.
fail() {
    echo "$NAME: $1" >&2
    write_result "error: $1" "" "" "" "" ""
    exit 0
}

# One sequence per record, whitespace stripped, in file order.
fasta_seqs() {
    awk '
        /^>/ { if (seen) print seq; seq = ""; seen = 1; next }
        { gsub(/[ \t\r]/, "", $0); seq = seq $0 }
        END { if (seen) print seq }
    ' "$1"
}

CHAINS=()
N_FIXED=0

add_fixed() { # one fixed chain -> REPEAT consecutive copies
    local s=$1
    s="${s//[[:space:]]/}"
    [[ -z "$s" ]] && return 0
    if [[ ! "$s" =~ $AA_RE ]]; then
        fail "invalid residue in fixed chain '${s:0:40}' (only the 20 standard uppercase one-letter codes)"
    fi
    local i
    for ((i = 0; i < REPEAT; i++)); do
        CHAINS+=("$s")
    done
    N_FIXED=$((N_FIXED + 1))
}

add_side() { # csv fasta -> the side's fixed chains, in order
    local csv=$1 fasta=$2 chain
    local raw=()
    if [[ -n "$csv" ]]; then
        mapfile -t raw < <(printf '%s' "$csv" | tr ',' '\n')
    elif [[ -n "$fasta" ]]; then
        [[ -f "$fasta" ]] || fail "missing fasta ($fasta)"
        mapfile -t raw < <(fasta_seqs "$fasta")
        if ((${#raw[@]} == 0)); then
            fail "no records in fasta ($fasta)"
        fi
    else
        return 0
    fi
    for chain in ${raw[@]+"${raw[@]}"}; do
        add_fixed "$chain"
    done
}

echo "[$(date +%T)] task $SAPIA_TASK_ID: mkcomplex for $NAME"

add_side "$PREPEND_SEQS" "$PREPEND_FASTA"
DESIGN_INDEX=${#CHAINS[@]}

SEQ="${SEQ//[[:space:]]/}"
case "$SEQ" in
    "" | NA | NaN | nan | None | none | null)
        fail "empty input sequence" ;;
esac
if [[ "$SEQ" == *"$SEP"* ]]; then
    fail "input sequence already contains the separator '$SEP' (already a multi-chain string)"
fi
CHAINS+=("$SEQ")

add_side "$APPEND_SEQS" "$APPEND_FASTA"

if ((N_FIXED == 0)); then
    fail "no fixed chains to add (empty --prepend-*/--append-*)"
fi

ASSEMBLED=""
CHAIN_LENS=""
TOTAL_LEN=0
for chain in "${CHAINS[@]}"; do
    ASSEMBLED="${ASSEMBLED:+${ASSEMBLED}${SEP}}${chain}"
    CHAIN_LENS="${CHAIN_LENS:+${CHAIN_LENS},}${#chain}"
    TOTAL_LEN=$((TOTAL_LEN + ${#chain}))
done

write_result "OK" "$ASSEMBLED" "${#CHAINS[@]}" "$CHAIN_LENS" "$DESIGN_INDEX" "$TOTAL_LEN"
echo "[$(date +%T)] $NAME: ${#CHAINS[@]} chains, $TOTAL_LEN residues, design at index $DESIGN_INDEX"
