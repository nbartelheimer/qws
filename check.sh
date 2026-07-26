#!/bin/bash
#
# check.sh - Validate a QWS run against a reference and FAIL on any anomaly.
#
# Usage: ./check.sh CASE data/CASE
#          CASE       : output produced by ./main
#          data/CASE  : reference output shipped with the repo
#
# Exit status:
#   0  all checked values are finite, aligned with the reference, and within
#      tolerance.
#   1  at least one anomaly was detected (nan/inf, missing/extra/reordered row,
#      or an out-of-tolerance value).
#   2  usage / unreadable-file error.
#
# Rationale: the previous implementation dropped any result line whose value was
# not a leading digit (e.g. "-nan"). That silently removed failed rows, which
# misaligned the paste against the reference and made awk abort on a divide-by-
# zero *after* it had already accumulated a clean error count - so a numerically
# broken run (rnorm^2 = -nan) was reported as passing. This version retains such
# rows, gates nan/inf explicitly, and verifies row-by-row alignment so that a
# "weird" run always errors out.
#
# On failure it also reports WHY, so a CI log is enough to diagnose without
# re-running: which rows failed and by how much relative to the tolerance that
# was applied, where a truncated run stopped, and any crash evidence found in
# the run output or in a companion "<CASE>.err" stderr file.
#----------------------------------------------------------------
readonly program=$(basename "$0")

print_usage_and_exit() {
  echo >&2 "Usage: ./${program} CASE data/CASE"
  exit 2
}

if [ $# -ne 2 ]; then
  print_usage_and_exit
fi

run=$1
ref=$2

for f in "$run" "$ref"; do
  if [ ! -r "$f" ]; then
    echo >&2 "${program}: cannot read file '$f'"
    exit 2
  fi
done

#----------------------------------------------------------------
# Post-mortem context, printed whenever we fail. Keeps the "why" in the log.
#----------------------------------------------------------------
report_context() {
  echo "----------------------------------------------------------------"
  echo "context: run='$run' ($(wc -l < "$run") lines), ref='$ref' ($(wc -l < "$ref") lines)"

  # Crash evidence in the run output and in its companion stderr file, if any.
  local errfile="${run}.err"
  local pattern='Assertion|Segmentation fault|signal|Aborted|terminate called|corrupt|malloc|free\(\)|stack smashing|MPI_ABORT|Killed|out of memory'
  local hits
  hits=$(grep -Ein "$pattern" "$run" | head -5)
  if [ -n "$hits" ]; then
    echo "crash evidence in '$run':"
    printf '  %s\n' "$hits"
  fi
  if [ -r "$errfile" ] && [ -s "$errfile" ]; then
    echo "stderr file '$errfile' (last 10 lines):"
    tail -10 "$errfile" | sed 's/^/  /'
  else
    echo "hint: no non-empty '$errfile'; capture stderr as well"
    echo "      (./main ... > $run 2> $errfile) so aborts/assertions are visible here"
  fi

  echo "last 5 lines of '$run':"
  tail -5 "$run" | sed 's/^/  /'
  echo "----------------------------------------------------------------"
}

#----------------------------------------------------------------
# 1) Hard gate: any nan/inf immediately after a ':' or '=' on a result line is
#    an unconditional failure, independent of the reference comparison below.
#----------------------------------------------------------------
if grep -Ev 'git commit hash' "$run" \
     | grep -Eiq '[:=][[:space:]]*[-+]?(nan|inf)([^a-zA-Z0-9]|$)'; then
  echo "check: FAIL - non-finite value(s) (nan/inf) in '$run':"
  grep -Ev 'git commit hash' "$run" \
    | grep -Ein '[:=][[:space:]]*[-+]?(nan|inf)([^a-zA-Z0-9]|$)' >&2
  # Continue so the row table below is still printed for context, but remember
  # that we have already failed.
  hard_fail=1
else
  hard_fail=0
fi

#----------------------------------------------------------------
# 2) Extract the checked rows, in file order, as "<label>\t<value>".
#    Two kinds of rows are checked, exactly as before:
#      - kernel self-checks:   "<label> :  <value>"
#      - solver norms:         "<bnorm|rnorm|xnorm>^2 = <value>"
#    Non-numeric values (nan/inf) are RETAINED here on purpose.
#----------------------------------------------------------------
select_rows() {
  grep -Ev 'git commit hash' "$1" \
  | awk '
      / : +[^[:space:]]/                        { print $1 "\t" $NF; next }
      /^[[:space:]]*(bnorm|rnorm|xnorm)\^2[[:space:]]*=/ { print $1 "\t" $NF; next }
    '
}

runrows=$(select_rows "$run")
refrows=$(select_rows "$ref")

nrun=$(printf '%s\n' "$runrows" | grep -c .)
nref=$(printf '%s\n' "$refrows" | grep -c .)

if [ "$nref" -eq 0 ]; then
  echo >&2 "check: FAIL - reference '$ref' contains no comparable rows"
  exit 1
fi

if [ "$nrun" -ne "$nref" ]; then
  echo "check: FAIL - row count mismatch: run has $nrun checked rows, reference has $nref"
  echo "       (a missing row usually means the run crashed or diverged early)"

  if [ "$nrun" -lt "$nref" ]; then
    last_ok=$(printf '%s\n' "$runrows"  | grep -c . >/dev/null; printf '%s\n' "$runrows" | awk -F'\t' 'NF{l=$1} END{print l}')
    next_exp=$(printf '%s\n' "$refrows" | awk -F'\t' -v n="$nrun" 'NF{i++; if (i==n+1) {print $1; exit}}')
    echo "       run stopped after row $nrun ('${last_ok}'); the next expected row is '${next_exp}'"
    echo "       -> the run died or was cut short there; the $((nref-nrun)) remaining rows were never produced"
  else
    echo "       run produced $((nrun-nref)) MORE rows than the reference"
    echo "       -> output format changed, or two runs were concatenated into '$run'"
  fi

  # Which labels differ, accounting for duplicates (bnorm^2 etc. appear twice).
  label_counts() { printf '%s\n' "$1" | awk -F'\t' 'NF{print $1}' | sort | uniq -c \
                   | awk '{printf "%s (x%s)\n", $2, $1}'; }
  missing=$(comm -13 <(label_counts "$runrows") <(label_counts "$refrows"))
  extra=$(comm -23 <(label_counts "$runrows") <(label_counts "$refrows"))
  [ -n "$missing" ] && { echo "       rows expected from the reference but missing (or fewer) in the run:"
                         printf '%s\n' "$missing" | sed 's/^/         /'; }
  [ -n "$extra" ]   && { echo "       rows present in the run but not in the reference:"
                         printf '%s\n' "$extra" | sed 's/^/         /'; }

  report_context
  exit 1
fi

#----------------------------------------------------------------
# 3) Row-by-row comparison. Alignment (label match) is verified explicitly so a
#    reordered / substituted row cannot slip through. Divisions are guarded.
#      rnorm^2 : converged residual must be < 1e-15 (absolute)
#      *_s / *_s_ (single precision) : relative diff < 3e-6
#      everything else (double)      : relative diff < 1e-14
#    Failing rows are repeated in a summary at the end, with the tolerance that
#    was applied and how far past it the value is.
#----------------------------------------------------------------
table=$(paste <(printf '%s\n' "$runrows") <(printf '%s\n' "$refrows") \
| awk -F'\t' -v hard_fail="$hard_fail" '
    BEGIN { err = 0; nbad = 0 }
    {
      rl = $1; rv = $2; fl = $3; fv = $4; row++

      if (rl != fl) {
        printf("%-27s  run=%-22s  ref=%-22s  MISALIGNED(ref label=%s)\n", rl, rv, "", fl)
        bad[++nbad] = sprintf("  row %d: label mismatch - run has \"%s\", reference has \"%s\"\n" \
                              "         -> the run printed a different set/order of checks", row, rl, fl)
        err++; next
      }

      low = tolower(rv)
      if ((rv + 0) != rv || index(low, "nan") || index(low, "inf")) {
        printf("%-27s  run=%-22s  ref=%-22s  NON-FINITE\n", rl, rv, fv)
        bad[++nbad] = sprintf("  row %d: %s = %s is not a finite number (ref %s)\n" \
                              "         -> the solver diverged or read uninitialised data", row, rl, rv, fv)
        err++; next
      }

      d = rv - fv; if (d < 0) d = -d
      denom = (rv < 0) ? -rv : rv

      if (rl == "rnorm^2") {
        tol = 1e-15; kind = "absolute (converged residual)"
        metric = rv + 0; ok = (metric < tol)
      } else if (rl ~ /_s_$/ || rl ~ /_s$/) {
        tol = 3e-6; kind = "relative (single precision)"
        metric = (denom > 0) ? d / denom : d; ok = (denom > 0) ? (metric < tol) : (d == 0)
      } else {
        tol = 1e-14; kind = "relative (double precision)"
        metric = (denom > 0) ? d / denom : d; ok = (denom > 0) ? (metric < tol) : (d == 0)
      }

      e = ok ? 0 : 1
      err += e
      if (e) {
        over = (tol > 0) ? metric / tol : 0
        bad[++nbad] = sprintf("  row %d: %s = %s vs ref %s\n" \
                              "         |diff| = %.3e, %s = %.3e, tolerance = %.1e -> exceeded by %.1fx",
                              row, rl, rv, fv, d, kind, metric, tol, over)
      }
      printf("%-27s  run=%-22s  ref=%-22s  d=%.3e  %s\n", rl, rv, fv, d, e ? "ERR" : "ok")
    }
    END {
      if (err > 0 || hard_fail == "1") {
        printf("check: FAIL - %d row error(s)%s\n", err, (hard_fail == "1") ? " + nan/inf gate" : "")
        if (nbad > 0) {
          print  "why it failed:"
          for (i = 1; i <= nbad; i++) print bad[i]
        }
        if (hard_fail == "1")
          print "  plus: the nan/inf gate above matched - see the lines it printed"
        exit 1
      }
      print "check: OK - all rows finite and within tolerance"
      exit 0
    }
  ')
status=$?
printf '%s\n' "$table"
[ $status -ne 0 ] && report_context
exit $status
