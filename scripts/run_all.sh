#!/bin/bash

run_all_main() (
    set -eo pipefail

    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
    cd "$SCRIPT_DIR"

    if [[ ! -f "$SCRIPT_DIR/set_ifx" ]]; then
        echo "ERROR: Intel environment setup script not found: $SCRIPT_DIR/set_ifx" >&2
        return 1
    fi

    source "$SCRIPT_DIR/set_ifx"
    set -u

    # Usage:
    #   ./run_all.sh [file_type] [first] [last]
    #
    # file_type:
    #   1 = Sim, runs Sims/Sim_XX.dat
    #   2 = DB,  runs Sims/DB_XX.dat
    FILE_TYPE="${1:-1}"
    FIRST="${2:-1}"
    LAST="${3:-2}"
    EXECUTABLE="$PROJECT_ROOT/build/T3ST"

    if [[ "$FILE_TYPE" != "1" && "$FILE_TYPE" != "2" ]]; then
        echo "ERROR: file_type must be 1 (Sim) or 2 (DB), got: $FILE_TYPE" >&2
        return 1
    fi

    if ! [[ "$FIRST" =~ ^[0-9]+$ && "$LAST" =~ ^[0-9]+$ ]]; then
        echo "ERROR: first and last must be positive integers." >&2
        return 1
    fi

    if (( FIRST < 1 || LAST < FIRST || LAST > 99 )); then
        echo "ERROR: invalid simulation range: $FIRST..$LAST" >&2
        return 1
    fi

    if [[ ! -x "$EXECUTABLE" ]]; then
        echo "ERROR: executable not found or not executable: $EXECUTABLE" >&2
        echo "Run from the project root: source scripts/set_ifx && source scripts/compile" >&2
        return 1
    fi

    for ((val = FIRST; val <= LAST; val++)); do
        printf '=====================================\n'
        printf ' RUNNING %s NO. : %02d\n' "$([[ "$FILE_TYPE" == "1" ]] && echo Sim || echo DB)" "$val"
        printf ' Using script: %s\n' "$0"
        printf ' Started at: %s\n' "$(date)"
        printf '=====================================\n'

        "$EXECUTABLE" "$FILE_TYPE" "$val"

        printf ' Finished at: %s\n' "$(date)"
        printf '%s\n' '----------------------'
    done
)

run_all_main "$@"
status=$?
unset -f run_all_main
return "$status" 2>/dev/null || exit "$status"
