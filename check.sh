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
  exit 1
fi

#----------------------------------------------------------------
# 3) Row-by-row comparison. Alignment (label match) is verified explicitly so a
#    reordered / substituted row cannot slip through. Divisions are guarded.
#      rnorm^2 : converged residual must be < 1e-15 (absolute)
#      *_s / *_s_ (single precision) : relative diff < 3e-6
#      everything else (double)      : relative diff < 1e-14
#----------------------------------------------------------------
paste <(printf '%s\n' "$runrows") <(printf '%s\n' "$refrows") \
| awk -F'\t' '
    BEGIN { err = 0 }
    {
      rl = $1; rv = $2; fl = $3; fv = $4

      if (rl != fl) {
        printf("%-27s  run=%-22s  ref=%-22s  MISALIGNED(ref label=%s)\n", rl, rv, "", fl)
        err++; next
      }

      low = tolower(rv)
      if ((rv + 0) != rv || index(low, "nan") || index(low, "inf")) {
        printf("%-27s  run=%-22s  ref=%-22s  NON-FINITE\n", rl, rv, fv)
        err++; next
      }

      d = rv - fv; if (d < 0) d = -d
      denom = (rv < 0) ? -rv : rv

      if (rl == "rnorm^2") {
        ok = ((rv + 0) < 1e-15)
      } else if (rl ~ /_s_$/ || rl ~ /_s$/) {
        ok = (denom > 0) ? (d / denom < 3e-6) : (d == 0)
      } else {
        ok = (denom > 0) ? (d / denom < 1e-14) : (d == 0)
      }

      e = ok ? 0 : 1
      err += e
      printf("%-27s  run=%-22s  ref=%-22s  d=%.3e  %s\n", rl, rv, fv, d, e ? "ERR" : "ok")
    }
    END {
      if (err > 0 || "'"$hard_fail"'" == "1") {
        printf("check: FAIL - %d row error(s)%s\n", err, ("'"$hard_fail"'" == "1") ? " + nan/inf gate" : "")
        exit 1
      }
      print "check: OK - all rows finite and within tolerance"
      exit 0
    }
  '
