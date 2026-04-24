#!/bin/bash
set -euo pipefail

tmp=$(mktemp)
cleanup() {
  rm -f "$tmp"
}
trap cleanup EXIT

cargo run --release -- --note "autoresearch-checks" >"$tmp" 2>&1 || {
  tail -80 "$tmp"
  exit 1
}

# Qubit budget guard: user-requested cap is 2800 for this Toffoli session.
qubits=$(awk -F: '/qubits/{gsub(/ /, "", $2); print $2; exit}' "$tmp")
if [[ -n "${qubits:-}" ]] && (( qubits > 2800 )); then
  echo "CHECKS FAIL: peak qubits $qubits exceeds session cap 2800"
  exit 1
fi
