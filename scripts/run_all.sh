#!/bin/bash

run_all_main() (
    set -o pipefail
    child_pid=""

    handle_interrupt() {
        printf '\nInterrupted. Stopping run_all.sh and returning to the shell.\n' >&2
        if [[ -n "${child_pid:-}" ]] && kill -0 "$child_pid" 2>/dev/null; then
            kill -INT "$child_pid" 2>/dev/null || true
            wait "$child_pid" 2>/dev/null || true
        fi
        exit 130
    }

    trap handle_interrupt INT TERM

    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
    cd "$SCRIPT_DIR"

    if [[ ! -f "$SCRIPT_DIR/set_ifx" ]]; then
        echo "ERROR: Intel environment setup script not found: $SCRIPT_DIR/set_ifx" >&2
        return 1
    fi

    original_args=("$@")
    set +u
    source "$SCRIPT_DIR/set_ifx"
    set -- "${original_args[@]}"
    set -u

    # Usage:
    #   ./run_all.sh [file_type] [first] [last]
    #
    # file_type:
    #   1 = Sim, runs Sims/Sim_XXX.dat
    #   2 = DB,  runs Sims/DB_XXX.dat
    if (( $# > 3 )); then
        echo "ERROR: too many arguments." >&2
        echo "Usage: ./run_all.sh [file_type] [first] [last]" >&2
        return 1
    fi

    FILE_TYPE="${1:-1}"
    FIRST="${2:-1}"
    LAST="${3:-$FIRST}"
    EXECUTABLE="$PROJECT_ROOT/build/T3ST"

    if [[ "$FILE_TYPE" != "1" && "$FILE_TYPE" != "2" ]]; then
        echo "ERROR: file_type must be 1 (Sim) or 2 (DB), got: $FILE_TYPE" >&2
        return 1
    fi

    if ! [[ "$FIRST" =~ ^[0-9]+$ && "$LAST" =~ ^[0-9]+$ ]]; then
        echo "ERROR: first and last must be positive integers." >&2
        return 1
    fi

    if (( FIRST < 1 || LAST < FIRST || LAST > 999 )); then
        echo "ERROR: invalid simulation range: $FIRST..$LAST" >&2
        return 1
    fi

    if [[ ! -x "$EXECUTABLE" ]]; then
        echo "ERROR: executable not found or not executable: $EXECUTABLE" >&2
        echo "Run from the project root: source scripts/set_ifx && source scripts/compile" >&2
        return 1
    fi

    printf 'Requested file_type=%s, range=%03d..%03d\n' "$FILE_TYPE" "$FIRST" "$LAST"

    for ((val = FIRST; val <= LAST; val++)); do
        printf '=====================================\n'
        printf ' RUNNING %s NO. : %03d\n' "$([[ "$FILE_TYPE" == "1" ]] && echo Sim || echo DB)" "$val"
        printf ' Using script: %s\n' "$0"
        printf ' Started at: %s\n' "$(date)"
        printf '=====================================\n'

        "$EXECUTABLE" "$FILE_TYPE" "$val" &
        child_pid=$!
        wait "$child_pid"
        status=$?
        child_pid=""

        if (( status == 130 || status == 143 )); then
            handle_interrupt
        elif (( status != 0 )); then
            printf 'ERROR: %s exited with status %d\n' "$EXECUTABLE" "$status" >&2
            return "$status"
        fi

        printf ' Finished at: %s\n' "$(date)"
        printf '%s\n' '----------------------'
    done
)

run_all_main "$@"
status=$?
unset -f run_all_main
return "$status" 2>/dev/null || exit "$status"
