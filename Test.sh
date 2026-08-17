#!/usr/bin/env bash

cd /home/pierbacq/Softwares/FreeFem-sources || exit 1
LOGDIR=$(mktemp -d /tmp/ff-dtest-XXXXXX); NRUN=0; NBAD=0

ffrun() {  # ffrun <np> <edp> [args...]
  local np=$1 edp=$2; shift 2
  local tag; tag=$(printf '%s np=%s %s' "$(basename "$edp" .edp)" "$np" "$*")
  local log="$LOGDIR/$(echo "$tag" | tr -c 'A-Za-z0-9' '_').log"
  NRUN=$((NRUN+1)); printf '%-62s ' "$tag"
  if timeout -k 10 600 src/mpi/ff-mpirun -np "$np" "$edp" -ns "$@" >"$log" 2>&1; then
    if grep -q '>>> ALL PASS' "$log"; then
      echo "OK    $(grep -c '=> PASS' "$log") criteres"
    else
      echo "FAIL  -> $log"; NBAD=$((NBAD+1))
    fi
  else
    local rc=$?
    [ $rc -eq 124 ] && echo "TIMEOUT (blocage : collectif non uniforme ?)  -> $log" \
                    || echo "ERREUR rc=$rc  -> $log"
    NBAD=$((NBAD+1))
  fi
}

echo "=== scatter (S7 = seul controle de PoU possible en mode scatter) ==="
ffrun 2 examples/mpi/test_distribute_scatter.edp
ffrun 3 examples/mpi/test_distribute_scatter.edp
ffrun 4 examples/mpi/test_distribute_scatter.edp
ffrun 3 examples/mpi/test_distribute_scatter.edp -ov 2
ffrun 3 examples/mpi/test_distribute_scatter.edp -case 1 -strip 1 -nc 8
ffrun 3 examples/mpi/test_distribute_scatter.edp -case 1 -strip 2 -nc 8
ffrun 3 examples/mpi/test_distribute_scatter.edp -case 1 -strip 3 -nc 8
ffrun 3 examples/mpi/test_distribute_scatter.edp -v 5

echo "=== multi (V6, W6, X1) -- exige Scotch ET ParMETIS ==="
ffrun 3 examples/mpi/test_distribute_multi.edp -wg
ffrun 4 examples/mpi/test_distribute_multi.edp -wg
ffrun 3 examples/mpi/test_distribute_multi.edp -wg -ov 2
ffrun 3 examples/mpi/test_distribute_multi.edp -wg -case 1 -strip 1 -nc 8
ffrun 3 examples/mpi/test_distribute_multi.edp -wg -case 1 -strip 2 -nc 8
ffrun 3 examples/mpi/test_distribute_multi.edp -wg -case 1 -strip 3 -nc 8
ffrun 3 examples/mpi/test_distribute_multi.edp -wg -v 5

echo
echo "=== BILAN : $((NRUN-NBAD)) / $NRUN runs OK   (logs: $LOGDIR) ==="
echo "--- valeurs des nouveaux criteres ---"
grep -h -o 'V6=[^ ]*\|W6=[^ ]*\|S7=[^ ]*\|X1=[^ ]*' "$LOGDIR"/*.log | sort | uniq -c | sort -rn
echo "--- diagnostic gratuit (runs -v 5) ---"
grep -h 'max|sum PoU\|intersections symetriques' "$LOGDIR"/*v_5*.log | sort -u

