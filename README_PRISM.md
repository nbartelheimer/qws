# No-sync PRISM backend on current master

This branch starts from master commit 478cbdb and retains its numerical kernels,
solver algorithms, profiling targets, and timer collector. Halo persistent
requests use OMPI_Psend_init / OMPI_Precv_init, OMPI_Start_no_sync for the
double-buffered single-precision SAP halos, OMPI_Wait, and
OMPI_Prequest_free from Open MPI's persistent extension. Collectives and
communicators continue to use ordinary MPI. The MPI fallbacks in the RDMA and
half-precision adapters use the extension as well.

qws_prism.h selects the start API. This no-sync branch changes the synchronized
port's two start aliases to OMPI_Start_no_sync / OMPI_Startall_no_sync; allocation,
buffer parity, completion waits, cleanup, and timing code remain identical.
Single-buffered double-precision and half-precision exchanges explicitly retain
OMPI_Start. Removing their receiver-ready synchronization allowed a fast sender
to overwrite unread data: the initial all-no-sync CASE8 run produced incorrect
residuals despite normal process exit. No-sync is therefore restricted to the
single-precision SAP path with alternating receive buffers. It does not remove
QWS's receive waits or send-completion waits.

## Build

Use an Open MPI installation that exports the persistent extension, including
OMPI_Prequest_free. With GNU compilers:

```sh
make clean
make -j4 fugaku_benchmark= compiler=openmpi-gnu arch=skylake \
  mpi=1 omp=1 rdma= powerapi= profiler=timing target=send \
  CFLAGS='-O2 -g -fopenmp'
```

The Makefile automatically defines _NO_OMP_SINGLE for MPI builds, ensuring
master-thread communication and consistent thread-0 send timing. Profiling
target defines are in MYFLAGS, so overriding CFLAGS does not silently discard
target=send. This is important when selecting architecture optimization flags.
The driver requests MPI_THREAD_FUNNELED (MPI_THREAD_MULTIPLE for threaded RDMA).
Library callers must initialize MPI with the corresponding thread support.

## Timers and overlap probes

Ordinary targets such as send, recv, send_post, all_calc, all, and reduc1/2/3
retain master's instrumentation. timing.c includes timing.h so the schema guard
uses real target constants, rather than comparing two undefined names as zero.
The collector retains master's thread-0 check before list initialization.

Master's polling probe uses MPI_Test. The extension has no OMPI_Test, and
extension handles cannot be passed to standard MPI request routines. Thus
polling probe modes fail compilation with a diagnostic. The explicit
target=overlap_probe_full_wait_only mode is supported: parent/vector timers,
overlap-window timings, directional completion-at-wait observations, and send
wait timings remain available, but hidden completion is not observed. Its
schema reports probe_mode:wait_only and hidden_recv_completion:not_observed.

All extension calls are checked. Before cleanup, neighboring ranks exchange the
active receive masks and drain the still-preposted receives using the existing
send buffers. This finalization-only traffic runs after timing and cannot alter
the numerical result. It makes buffered requests inactive without MPI_Cancel.
Extension requests are then released through OMPI_Prequest_free before their
buffers are freed, including receive slots that were never used numerically.
A barrier before boundary receive waits prevents overlap with send initiation
when a caller has no intervening bulk calculation.

## Validation

Run local build/instrumentation checks with:

```sh
python3 tests/check_prism_build.py
```

Numerical validation uses README CASE0, two-rank CASE1, and eight-rank CASE7 with
check.sh against the existing reference data. Do not update references to hide
numerical discrepancies. Tests on the InfiniBand hosts and branch-specific
results are recorded in the accompanying validation/prism-master directory of
the original checkout. Hardware-specific uTofu/RDMA and half-precision runtime
paths require their respective toolchains and are not covered by those tests.

The synchronized parent port passed CASE0, CASE1, CASE7, CASE8, and the short
wait-only probe on 2026-11-09. This branch is validated separately with no-sync
starts; see the results in the original checkout's validation/prism-master.

No-sync validation on 2026-11-09: CASE0 (one rank/one thread), CASE1 (two
ranks/eight threads), CASE7 (eight ranks/eight threads), and CASE8 (16 ranks/four
threads) all pass the supplied references with GNU 11.5.0, Open MPI 5.0.9,
arch=skylake, -O2 -g -fopenmp, profiler=timing target=send. Timer sample counts
are 10,000 for CASE7 and 11,000 for the other grids. The short two-rank wait-only
probe also completes with the correct schema and finite numerical results.
