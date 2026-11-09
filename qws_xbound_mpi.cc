//****************************************************************************************
//
//  Copyright (c) 2015-2020, Yoshifumi Nakamura <nakamura@riken.jp>
//  Copyright (c) 2015-2020, Yuta Mukai         <mukai.yuta@fujitsu.com>
//  Copyright (c) 2018-2020, Ken-Ichi Ishikawa  <ishikawa@theo.phys.sci.hirosima-u.ac.jp>
//  Copyright (c) 2019-2020, Issaku Kanamori    <kanamori-i@riken.jp>
//
//
//  All rights reserved.
//
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are
//  met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer. 
//
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer listed
//    in this license in the documentation and/or other materials
//    provided with the distribution.
//
//  * Neither the name of the copyright holders nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT  
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
//  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  
//  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
//
//----------------------------------------------------------------------------------------
//  ACKNOWLEDGMENT
//
//  This software has been developed in a co-design working group for the lattice QCD
//  supported by MEXT's programs for the Development and Improvement for the Next
//  Generation Ultra High-Speed Computer System, under its Subsidies for Operating the
//  Specific Advanced Large Research Facilities, and Priority Issue 9 
//  (Elucidation of the Fundamental Laws and Evolution of the Universe) to be tackled by
//  using the Supercomputer Fugaku.
//
//****************************************************************************************
#include "qws.h"
#include "timing.h"
#if defined(PROF_TARGET) && PROF_TARGET == TARGET_OVERLAP_PROBE && !defined(OVERLAP_PROBE_WAIT_ONLY)
#error "PRISM has no OMPI_Test; use target=overlap_probe_full_wait_only for completion-at-wait observations"
#endif
#include <stdio.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
//#define min(a,b) (a)>(b)?(b):(a) 

//#include "timing.h"
//#include "eml_lib.h"

#ifdef _MPI_


#ifdef __cplusplus
extern "C"{
#endif

#ifdef _USE_RANKMAP
#include "rankmap_lib.h"
#endif

  extern   projscd1_t *xfd_send;
  extern   projscd1_t *xfd_recv;
  extern   projscd1_t *xbd_send;
  extern   projscd1_t *xbd_recv;
  extern   projscd_t *yfd_send;
  extern   projscd_t *yfd_recv;
  extern   projscd_t *ybd_send;
  extern   projscd_t *ybd_recv;
  extern   projscd_t *zfd_send;
  extern   projscd_t *zfd_recv;
  extern   projscd_t *zbd_send;
  extern   projscd_t *zbd_recv;
  extern   projscd_t *tfd_send;
  extern   projscd_t *tfd_recv;
  extern   projscd_t *tbd_send;
  extern   projscd_t *tbd_recv;

  extern   projscs1_t *xfs_send;
  extern   projscs1_t *xfs_recv;
  extern   projscs1_t *xbs_send;
  extern   projscs1_t *xbs_recv;
  extern   projscs_t *yfs_send;
  extern   projscs_t *yfs_recv;
  extern   projscs_t *ybs_send;
  extern   projscs_t *ybs_recv;
  extern   projscs_t *zfs_send;
  extern   projscs_t *zfs_recv;
  extern   projscs_t *zbs_send;
  extern   projscs_t *zbs_recv;
  extern   projscs_t *tfs_send;
  extern   projscs_t *tfs_recv;
  extern   projscs_t *tbs_send;
  extern   projscs_t *tbs_recv;

  extern   int npe[4];
  extern   int nx;
  extern   int ny;
  extern   int nz;
  extern   int nt;
  extern   int nxh;
  extern   int nxd;
  extern   int nxs;

  projscs1_t *xfs_recv_array[2];
  projscs1_t *xbs_recv_array[2];
  projscs_t *yfs_recv_array[2];
  projscs_t *ybs_recv_array[2];
  projscs_t *zfs_recv_array[2];
  projscs_t *zbs_recv_array[2];
  projscs_t *tfs_recv_array[2];
  projscs_t *tbs_recv_array[2];
  int bufs_parity;


  MPI_Request sd_req[8];
  MPI_Request rd_req[8];
  MPI_Request *ss_req;
  MPI_Request *rs_req;
  MPI_Request *rs_req_opposite;
  MPI_Request ss_req_array[2][8];
  MPI_Request rs_req_array[2][8];

  int recv_started_array[2][8];
  int *recv_started;
  int *recv_started_opposite;
  static int send_peer[8], recv_peer[8];

  projscd1_t *xfd_recv0;
  projscd1_t *xbd_recv0;

#if defined(PROF_TARGET) && PROF_TARGET == TARGET_OVERLAP_PROBE && defined(_CHECK_TIMING)
  static int overlap_probe_active = 0;
  static int overlap_probe_recv_done = 0;
  static int overlap_probe_exposed_wait_active = 0;
  static int overlap_probe_post_recv_window_active = 0;
  static int overlap_probe_req_active[8];
  static int overlap_probe_req_done[8];
  static int overlap_probe_req_post_recorded[8];
  static int overlap_probe_send_label_valid = 0;
  static int overlap_probe_recv_active_count = 0;
  static const char *overlap_probe_next_label = "unlabeled";
  static const char *overlap_probe_active_label = "unlabeled";
  static const char *overlap_probe_send_label = "unlabeled";
  static const char *overlap_probe_req_name[8] = {
    "xf", "xb", "yf", "yb", "zf", "zb", "tf", "tb"
  };
  static double overlap_probe_begin_clock = 0.0;
  static double overlap_probe_done_clock = 0.0;
  static double overlap_probe_calc_end_clock = 0.0;
  static double overlap_probe_exposed_wait_start_clock = 0.0;
  static double overlap_probe_last_phase_clock = 0.0;
  static double overlap_probe_req_done_clock[8];
  static double overlap_probe_req_pending_elapsed[8];
  static double overlap_probe_req_post_elapsed[8];
  static double overlap_probe_req_exposed_elapsed[8];

  static void overlap_probe_add_timing(const char *id, double elapsed) {
    if (elapsed < 0.0) elapsed = 0.0;
    add_timing_delta_(id, &elapsed);
  }

  static void overlap_probe_add_labeled_timing(const char *metric,
                                               const char *label,
                                               const char *phase,
                                               double elapsed) {
    char id[160];
    if (phase == NULL) {
      snprintf(id, sizeof(id), "%s.%s", metric, label);
    } else {
      snprintf(id, sizeof(id), "%s.%s.%s", metric, label, phase);
    }
    overlap_probe_add_timing(id, elapsed);
  }

  static void overlap_probe_add_metric(const char *metric, double elapsed) {
    overlap_probe_add_timing(metric, elapsed);
    overlap_probe_add_labeled_timing(metric, overlap_probe_active_label, NULL, elapsed);
  }

  static int overlap_probe_is_unlabeled(const char *label) {
    return label == NULL || strcmp(label, "unlabeled") == 0;
  }

  static void overlap_probe_add_phase_metric(const char *metric,
                                             const char *phase,
                                             double elapsed) {
    overlap_probe_add_timing(metric, elapsed);
    overlap_probe_add_labeled_timing(metric, overlap_probe_active_label, NULL, elapsed);
    overlap_probe_add_labeled_timing(metric, overlap_probe_active_label, phase, elapsed);
  }

  static void overlap_probe_add_dir_metric(const char *metric,
                                           int req,
                                           double elapsed) {
    char id[160];
    if (req < 0 || req >= 8) return;
    snprintf(id, sizeof(id), "%s.%s.%s",
             metric, overlap_probe_active_label, overlap_probe_req_name[req]);
    overlap_probe_add_timing(id, elapsed);
  }

  static void overlap_probe_add_dir_phase_metric(const char *metric,
                                                 int req,
                                                 const char *phase,
                                                 double elapsed) {
    char id[160];
    if (req < 0 || req >= 8 || phase == NULL) return;
    snprintf(id, sizeof(id), "%s.%s.%s.%s",
             metric, overlap_probe_active_label, phase, overlap_probe_req_name[req]);
    overlap_probe_add_timing(id, elapsed);
  }

  static void overlap_probe_add_dir_check_metric(const char *metric,
                                                 int req,
                                                 double value) {
    char id[160];
    double abs_value = value < 0.0 ? -value : value;
    if (req < 0 || req >= 8) return;
    snprintf(id, sizeof(id), "%s.%s.%s",
             metric, overlap_probe_active_label, overlap_probe_req_name[req]);
    overlap_probe_add_timing(id, abs_value);
  }

  static void overlap_probe_add_send_metric(const char *metric,
                                            double elapsed) {
    overlap_probe_add_timing(metric, elapsed);
    if (overlap_probe_send_label_valid) {
      overlap_probe_add_labeled_timing(metric, overlap_probe_send_label, NULL, elapsed);
    }
  }

  static void overlap_probe_add_send_dir_metric(const char *metric,
                                                int req,
                                                double elapsed) {
    char id[160];
    if (req < 0 || req >= 8) return;
    if (!overlap_probe_send_label_valid) return;
    snprintf(id, sizeof(id), "%s.%s.%s",
             metric, overlap_probe_send_label, overlap_probe_req_name[req]);
    overlap_probe_add_timing(id, elapsed);
  }

  static void overlap_probe_add_send_check_metric(const char *metric,
                                                  double value) {
    double abs_value = value < 0.0 ? -value : value;
    overlap_probe_add_send_metric(metric, abs_value);
  }

  static double overlap_probe_wait_send_request(MPI_Request *request,
                                                MPI_Status *status) {
    double start_clock = get_clock();
    CHECK(OMPI_Wait(request, status));
    return get_clock() - start_clock;
  }

  static void overlap_probe_record_send_wait(int req,
                                             double elapsed,
                                             double *sum_elapsed,
                                             double *max_elapsed,
                                             int *max_req,
                                             int *wait_count) {
    overlap_probe_add_send_dir_metric("overlap_probe_send_dir_wait", req, elapsed);
    *sum_elapsed += elapsed;
    if (*wait_count == 0 || elapsed > *max_elapsed) {
      *max_elapsed = elapsed;
      *max_req = req;
    }
    (*wait_count)++;
  }

  static void overlap_probe_record_phase(const char *phase) {
    double now_clock;
    if (phase == NULL || !overlap_probe_active) return;
    now_clock = get_clock();
    overlap_probe_add_phase_metric("overlap_probe_phase",
                                   phase,
                                   now_clock - overlap_probe_last_phase_clock);
    overlap_probe_last_phase_clock = now_clock;
  }

  static void overlap_probe_record_window_begin(void) {
    overlap_probe_add_metric("overlap_probe_window_begin", 0.0);
    overlap_probe_add_metric("overlap_probe_recv_active_count",
                             (double)overlap_probe_recv_active_count);
    if (overlap_probe_is_unlabeled(overlap_probe_active_label)) {
      overlap_probe_add_metric("overlap_probe_check_unlabeled_window", 1.0);
    }
    if (overlap_probe_send_label_valid) {
      overlap_probe_add_timing("overlap_probe_check_send_label_not_consumed", 1.0);
      overlap_probe_add_labeled_timing("overlap_probe_check_send_label_not_consumed",
                                       overlap_probe_send_label,
                                       NULL,
                                       1.0);
    }
    for (int req=0; req<8; req++) {
      if (!overlap_probe_req_active[req]) continue;
      overlap_probe_add_dir_metric("overlap_probe_recv_active_dir", req, 0.0);
    }
  }

  static int overlap_probe_all_receives_done(void) {
    for (int req=0; req<8; req++) {
      if (overlap_probe_req_active[req] && !overlap_probe_req_done[req]) return 0;
    }
    return 1;
  }

  static void overlap_probe_record_all_done(double done_clock,
                                            const char *phase) {
    double max_done_clock = overlap_probe_begin_clock;
    if (overlap_probe_recv_done) return;
    overlap_probe_done_clock = done_clock;
    if (overlap_probe_exposed_wait_active) {
      overlap_probe_add_metric("overlap_probe_recv_exposed_wait",
                               overlap_probe_done_clock - overlap_probe_exposed_wait_start_clock);
      overlap_probe_exposed_wait_active = 0;
    } else {
      overlap_probe_add_metric("overlap_probe_recv_pending_window",
                               overlap_probe_done_clock - overlap_probe_begin_clock);
      overlap_probe_post_recv_window_active = 1;
    }
    overlap_probe_add_metric("overlap_probe_recv_observed",
                             overlap_probe_done_clock - overlap_probe_begin_clock);
    if (phase != NULL) {
      overlap_probe_add_phase_metric("overlap_probe_recv_observed_at",
                                     phase,
                                     overlap_probe_done_clock - overlap_probe_begin_clock);
    }
    overlap_probe_recv_done = 1;
    for (int req=0; req<8; req++) {
      if (!overlap_probe_req_active[req]) continue;
      if (overlap_probe_req_done_clock[req] > max_done_clock) {
        max_done_clock = overlap_probe_req_done_clock[req];
      }
    }
    overlap_probe_add_metric("overlap_probe_check_recv_observed_vs_dir_max",
                             overlap_probe_done_clock - max_done_clock);
  }

  static void overlap_probe_record_dir_done(int req,
                                            const char *phase,
                                            double done_clock) {
    if (!overlap_probe_req_active[req] || overlap_probe_req_done[req]) return;
    overlap_probe_req_done[req] = 1;
    overlap_probe_req_done_clock[req] = done_clock;
    if (overlap_probe_exposed_wait_active) {
      overlap_probe_req_exposed_elapsed[req] =
        done_clock - overlap_probe_exposed_wait_start_clock;
      overlap_probe_add_dir_metric("overlap_probe_recv_dir_exposed_wait",
                                   req,
                                   overlap_probe_req_exposed_elapsed[req]);
    } else {
      overlap_probe_req_pending_elapsed[req] =
        done_clock - overlap_probe_begin_clock;
      overlap_probe_add_dir_metric("overlap_probe_recv_dir_pending_window",
                                   req,
                                   overlap_probe_req_pending_elapsed[req]);
    }
    overlap_probe_add_dir_metric("overlap_probe_recv_dir_observed",
                                 req,
                                 done_clock - overlap_probe_begin_clock);
    overlap_probe_add_dir_phase_metric("overlap_probe_recv_dir_observed_at",
                                       req,
                                       phase,
                                       done_clock - overlap_probe_begin_clock);
    if (overlap_probe_all_receives_done()) {
      overlap_probe_record_all_done(done_clock, phase);
    }
  }

  static void overlap_probe_record_dir_post_windows(double calc_end_clock) {
    for (int req=0; req<8; req++) {
      if (!overlap_probe_req_active[req]) continue;
      if (!overlap_probe_req_done[req]) continue;
      if (overlap_probe_req_post_recorded[req]) continue;
      if (overlap_probe_req_done_clock[req] > calc_end_clock) continue;
      overlap_probe_req_post_elapsed[req] =
        calc_end_clock - overlap_probe_req_done_clock[req];
      overlap_probe_add_dir_metric("overlap_probe_recv_dir_post_window",
                                   req,
                                   overlap_probe_req_post_elapsed[req]);
      overlap_probe_add_dir_check_metric("overlap_probe_check_dir_window_balance",
                                         req,
                                         overlap_probe_req_pending_elapsed[req]
                                         + overlap_probe_req_post_elapsed[req]
                                         - (calc_end_clock - overlap_probe_begin_clock));
      overlap_probe_req_post_recorded[req] = 1;
    }
  }

  static void overlap_probe_record_dir_observed_checks(void) {
    for (int req=0; req<8; req++) {
      if (!overlap_probe_req_active[req]) continue;
      if (!overlap_probe_req_done[req]) continue;
      overlap_probe_add_dir_check_metric("overlap_probe_check_dir_observed_balance",
                                         req,
                                         overlap_probe_req_pending_elapsed[req]
                                         + overlap_probe_req_exposed_elapsed[req]
                                         - (overlap_probe_req_done_clock[req]
                                            - overlap_probe_begin_clock));
    }
  }

  static int overlap_probe_on_master_thread(void) {
#ifdef _OPENMP
    return omp_get_thread_num() == 0;
#else
    return 1;
#endif
  }

  static int overlap_probe_calc_end_only_mode(void) {
#if defined(OVERLAP_PROBE_CALC_END_ONLY) || defined(OVERLAP_PROBE_WAIT_ONLY)
    return 1;
#else
    return 0;
#endif
  }

  static int overlap_probe_wait_only_mode(void) {
#ifdef OVERLAP_PROBE_WAIT_ONLY
    return 1;
#else
    return 0;
#endif
  }

  static int overlap_probe_coarse_phase_mode(void) {
#ifdef OVERLAP_PROBE_COARSE_PHASE
    return 1;
#else
    return 0;
#endif
  }

  static int overlap_probe_phase_is_calc_end(const char *phase) {
    return phase != NULL && strcmp(phase, "calc_end") == 0;
  }

  static int overlap_probe_active_label_has_accum(void) {
    return overlap_probe_active_label != NULL
      && strstr(overlap_probe_active_label, ".ddd_in_accum") != NULL;
  }

  static int overlap_probe_should_coalesce_poll(const char *phase) {
    return overlap_probe_coarse_phase_mode()
      && phase != NULL
      && strcmp(phase, "ddd_in") == 0;
  }

  static const char *overlap_probe_completion_phase(const char *phase) {
    if (overlap_probe_coarse_phase_mode()
        && overlap_probe_active_label_has_accum()
        && phase != NULL
        && (strcmp(phase, "accum_addsub") == 0
            || strcmp(phase, "accum_add") == 0)) {
      return "ddd_in_accum";
    }
    return phase;
  }

  static int overlap_probe_test_receives(const char *phase) {
    // PRISM exposes Wait but no Test. In the explicit wait-only probe,
    // completion observations are recorded by xbound_wait instead of polling.
    (void)phase;
    return overlap_probe_all_receives_done();
  }
#endif

#ifdef _USE_RANKMAP
  int rankmap_id;
#endif

  void xbound_set_parity(int parity, int prec){
    if(prec == 8){
      return;
    }
    bufs_parity=parity;
    ss_req = ss_req_array[parity];
    rs_req = rs_req_array[parity];
    rs_req_opposite = rs_req_array[1-parity];
    recv_started = &recv_started_array[parity][0];
    recv_started_opposite = &recv_started_array[1-parity][0];
  }

  void xbound_flip_parity(int prec){
    if(prec == 8){
      return;
    }
    xbound_set_parity(1-bufs_parity, prec);
  }

  void xbound_recv_updateall(int prec){
    if(prec==4){
      xfs_recv=xfs_recv_array[bufs_parity];
      xbs_recv=xbs_recv_array[bufs_parity];
      yfs_recv=yfs_recv_array[bufs_parity];
      ybs_recv=ybs_recv_array[bufs_parity];
      zfs_recv=zfs_recv_array[bufs_parity];
      zbs_recv=zbs_recv_array[bufs_parity];
      tfs_recv=tfs_recv_array[bufs_parity];
      tbs_recv=tbs_recv_array[bufs_parity];
    } else {
      xfd_recv=xfd_recv0;
      xbd_recv=xbd_recv0;
    }
  }

 //---------------------------------------------------------------------------------------- get rankmap
  int xbound_get_rankmap(int myrank, int *rank_coords, int *rank_size,
			 int *rankf, int *rankb) {
#ifdef _USE_RANKMAP
    static int dim=4;
    int mapid=rankmap_lib_set_rankmap4d();
    rankmap_id=mapid;
    if(mapid<0){ // no rank map is available
      return mapid;
    }
    int neighbor_ranks[2*dim];
    rankmap_lib_get_rankmap(rank_coords, neighbor_ranks, rank_size);
    for(int dir=0; dir<dim; dir++){
      rankf[dir]=neighbor_ranks[2*dir];
      rankb[dir]=neighbor_ranks[2*dir+1];
    }
    return mapid;
#else
    return -1;  // empty
#endif

  }

  //---------------------------------------------------------------------------------------- init communication
  void xbound_init(int rank,
		   int pxf, int pyf, int pzf, int ptf,
		   int pxb, int pyb, int pzb, int ptb)
  {

    const int send_ranks[8] = {pxb, pxf, pyb, pyf, pzb, pzf, ptb, ptf};
    const int recv_ranks[8] = {pxf, pxb, pyf, pyb, pzf, pzb, ptf, ptb};
    memcpy(send_peer, send_ranks, sizeof(send_peer));
    memcpy(recv_peer, recv_ranks, sizeof(recv_peer));

    int nxh=nx/2;
    int nxd=nx/2/VLEND;
    int nxs=nx/2/VLENS;

    // allocate communication buffers (double prec.)
    xfd_send = (projscd1_t*)malloc( sizeof(projscd1_t) * ny *nz*nt);
    xbd_send = (projscd1_t*)malloc( sizeof(projscd1_t) * ny *nz*nt);
    yfd_send = (projscd_t*)malloc( sizeof(projscd_t) * nxd*nz*nt);
    ybd_send = (projscd_t*)malloc( sizeof(projscd_t) * nxd*nz*nt);
    zfd_send = (projscd_t*)malloc( sizeof(projscd_t) * nxd*ny*nt);
    zbd_send = (projscd_t*)malloc( sizeof(projscd_t) * nxd*ny*nt);
    tfd_send = (projscd_t*)malloc( sizeof(projscd_t) * nxd*ny*nz);
    tbd_send = (projscd_t*)malloc( sizeof(projscd_t) * nxd*ny*nz);
    xfd_recv = (projscd1_t*)malloc( sizeof(projscd1_t) * ny *nz*nt);
    xbd_recv = (projscd1_t*)malloc( sizeof(projscd1_t) * ny *nz*nt);
    yfd_recv = (projscd_t*)malloc( sizeof(projscd_t) * nxd*nz*nt);
    ybd_recv = (projscd_t*)malloc( sizeof(projscd_t) * nxd*nz*nt);
    zfd_recv = (projscd_t*)malloc( sizeof(projscd_t) * nxd*ny*nt);
    zbd_recv = (projscd_t*)malloc( sizeof(projscd_t) * nxd*ny*nt);
    tfd_recv = (projscd_t*)malloc( sizeof(projscd_t) * nxd*ny*nz);
    tbd_recv = (projscd_t*)malloc( sizeof(projscd_t) * nxd*ny*nz);
    xfd_recv0 = xfd_recv; // keep the original value
    xbd_recv0 = xbd_recv;

    // allocate communication buffers (single prec.)
    void *tmp;
    posix_memalign(&tmp, CLS, sizeof(projscs1_t)*ny*nz*nt);
    xfs_send = (projscs1_t*)tmp;
    posix_memalign(&tmp, CLS, sizeof(projscs1_t)*ny*nz*nt);
    xbs_send = (projscs1_t*)tmp;
    posix_memalign(&tmp, CLS, sizeof(projscs_t)*nxs*nz*nt);
    yfs_send = (projscs_t*)tmp;
    posix_memalign(&tmp, CLS, sizeof(projscs_t)*nxs*nz*nt);
    ybs_send = (projscs_t*)tmp;
    posix_memalign(&tmp, CLS, sizeof(projscs_t)*nxs*ny*nt);
    zfs_send = (projscs_t*)tmp;
    posix_memalign(&tmp, CLS, sizeof(projscs_t)*nxs*ny*nt);
    zbs_send = (projscs_t*)tmp;
    posix_memalign(&tmp, CLS, sizeof(projscs_t)*nxs*ny*nz);
    tfs_send = (projscs_t*)tmp;
    posix_memalign(&tmp, CLS, sizeof(projscs_t)*nxs*ny*nz);
    tbs_send = (projscs_t*)tmp;
    for(int parity=0; parity<2; parity++){
      posix_memalign(&tmp, CLS, sizeof(projscs1_t)*ny*nz*nt);
      xfs_recv_array[parity] = (projscs1_t*)tmp;
      posix_memalign(&tmp, CLS, sizeof(projscs1_t)*ny*nz*nt);
      xbs_recv_array[parity] = (projscs1_t*)tmp;
      posix_memalign(&tmp, CLS, sizeof(projscs_t)*nxs*nz*nt);
      yfs_recv_array[parity] = (projscs_t*)tmp;
      posix_memalign(&tmp, CLS, sizeof(projscs_t)*nxs*nz*nt);
      ybs_recv_array[parity] = (projscs_t*)tmp;
      posix_memalign(&tmp, CLS, sizeof(projscs_t)*nxs*ny*nt);
      zfs_recv_array[parity] = (projscs_t*)tmp;
      posix_memalign(&tmp, CLS, sizeof(projscs_t)*nxs*ny*nt);
      zbs_recv_array[parity] = (projscs_t*)tmp;
      posix_memalign(&tmp, CLS, sizeof(projscs_t)*nxs*ny*nz);
      tfs_recv_array[parity] = (projscs_t*)tmp;
      posix_memalign(&tmp, CLS, sizeof(projscs_t)*nxs*ny*nz);
      tbs_recv_array[parity] = (projscs_t*)tmp;
    }

    // initializing communications: double prec.
    CHECK(OMPI_Psend_init(xfd_send, 12*ny *nz*nt, MPI_DOUBLE_PRECISION, pxb, 0, MPI_COMM_WORLD, &sd_req[0]));
    CHECK(OMPI_Psend_init(xbd_send, 12*ny *nz*nt, MPI_DOUBLE_PRECISION, pxf, 1, MPI_COMM_WORLD, &sd_req[1]));
    CHECK(OMPI_Psend_init(yfd_send, 12*nxh*nz*nt, MPI_DOUBLE_PRECISION, pyb, 2, MPI_COMM_WORLD, &sd_req[2]));
    CHECK(OMPI_Psend_init(ybd_send, 12*nxh*nz*nt, MPI_DOUBLE_PRECISION, pyf, 3, MPI_COMM_WORLD, &sd_req[3]));
    CHECK(OMPI_Psend_init(zfd_send, 12*nxh*ny*nt, MPI_DOUBLE_PRECISION, pzb, 4, MPI_COMM_WORLD, &sd_req[4]));
    CHECK(OMPI_Psend_init(zbd_send, 12*nxh*ny*nt, MPI_DOUBLE_PRECISION, pzf, 5, MPI_COMM_WORLD, &sd_req[5]));
    CHECK(OMPI_Psend_init(tfd_send, 12*nxh*ny*nz, MPI_DOUBLE_PRECISION, ptb, 6, MPI_COMM_WORLD, &sd_req[6]));
    CHECK(OMPI_Psend_init(tbd_send, 12*nxh*ny*nz, MPI_DOUBLE_PRECISION, ptf, 7, MPI_COMM_WORLD, &sd_req[7]));
    CHECK(OMPI_Precv_init(xfd_recv, 12*ny *nz*nt, MPI_DOUBLE_PRECISION, pxf, 0, MPI_COMM_WORLD, &rd_req[0]));
    CHECK(OMPI_Precv_init(xbd_recv, 12*ny *nz*nt, MPI_DOUBLE_PRECISION, pxb, 1, MPI_COMM_WORLD, &rd_req[1]));
    CHECK(OMPI_Precv_init(yfd_recv, 12*nxh*nz*nt, MPI_DOUBLE_PRECISION, pyf, 2, MPI_COMM_WORLD, &rd_req[2]));
    CHECK(OMPI_Precv_init(ybd_recv, 12*nxh*nz*nt, MPI_DOUBLE_PRECISION, pyb, 3, MPI_COMM_WORLD, &rd_req[3]));
    CHECK(OMPI_Precv_init(zfd_recv, 12*nxh*ny*nt, MPI_DOUBLE_PRECISION, pzf, 4, MPI_COMM_WORLD, &rd_req[4]));
    CHECK(OMPI_Precv_init(zbd_recv, 12*nxh*ny*nt, MPI_DOUBLE_PRECISION, pzb, 5, MPI_COMM_WORLD, &rd_req[5]));
    CHECK(OMPI_Precv_init(tfd_recv, 12*nxh*ny*nz, MPI_DOUBLE_PRECISION, ptf, 6, MPI_COMM_WORLD, &rd_req[6]));
    CHECK(OMPI_Precv_init(tbd_recv, 12*nxh*ny*nz, MPI_DOUBLE_PRECISION, ptb, 7, MPI_COMM_WORLD, &rd_req[7]));

    // initializing communications: single prec.
    for(int parity=0; parity<2; parity++){
      // send buffers are common for both parity
      int tag0=parity*8+8;
      CHECK(OMPI_Psend_init(xfs_send, 12*ny *nz*nt, MPI_REAL, pxb, tag0+0, MPI_COMM_WORLD, &ss_req_array[parity][0]));
      CHECK(OMPI_Psend_init(xbs_send, 12*ny *nz*nt, MPI_REAL, pxf, tag0+1, MPI_COMM_WORLD, &ss_req_array[parity][1]));
      CHECK(OMPI_Psend_init(yfs_send, 12*nxh*nz*nt, MPI_REAL, pyb, tag0+2, MPI_COMM_WORLD, &ss_req_array[parity][2]));
      CHECK(OMPI_Psend_init(ybs_send, 12*nxh*nz*nt, MPI_REAL, pyf, tag0+3, MPI_COMM_WORLD, &ss_req_array[parity][3]));
      CHECK(OMPI_Psend_init(zfs_send, 12*nxh*ny*nt, MPI_REAL, pzb, tag0+4, MPI_COMM_WORLD, &ss_req_array[parity][4]));
      CHECK(OMPI_Psend_init(zbs_send, 12*nxh*ny*nt, MPI_REAL, pzf, tag0+5, MPI_COMM_WORLD, &ss_req_array[parity][5]));
      CHECK(OMPI_Psend_init(tfs_send, 12*nxh*ny*nz, MPI_REAL, ptb, tag0+6, MPI_COMM_WORLD, &ss_req_array[parity][6]));
      CHECK(OMPI_Psend_init(tbs_send, 12*nxh*ny*nz, MPI_REAL, ptf, tag0+7, MPI_COMM_WORLD, &ss_req_array[parity][7]));
      CHECK(OMPI_Precv_init(xfs_recv_array[parity], 12*ny *nz*nt, MPI_REAL, pxf, tag0+0, MPI_COMM_WORLD, &rs_req_array[parity][0]));
      CHECK(OMPI_Precv_init(xbs_recv_array[parity], 12*ny *nz*nt, MPI_REAL, pxb, tag0+1, MPI_COMM_WORLD, &rs_req_array[parity][1]));
      CHECK(OMPI_Precv_init(yfs_recv_array[parity], 12*nxh*nz*nt, MPI_REAL, pyf, tag0+2, MPI_COMM_WORLD, &rs_req_array[parity][2]));
      CHECK(OMPI_Precv_init(ybs_recv_array[parity], 12*nxh*nz*nt, MPI_REAL, pyb, tag0+3, MPI_COMM_WORLD, &rs_req_array[parity][3]));
      CHECK(OMPI_Precv_init(zfs_recv_array[parity], 12*nxh*ny*nt, MPI_REAL, pzf, tag0+4, MPI_COMM_WORLD, &rs_req_array[parity][4]));
      CHECK(OMPI_Precv_init(zbs_recv_array[parity], 12*nxh*ny*nt, MPI_REAL, pzb, tag0+5, MPI_COMM_WORLD, &rs_req_array[parity][5]));
      CHECK(OMPI_Precv_init(tfs_recv_array[parity], 12*nxh*ny*nz, MPI_REAL, ptf, tag0+6, MPI_COMM_WORLD, &rs_req_array[parity][6]));
      CHECK(OMPI_Precv_init(tbs_recv_array[parity], 12*nxh*ny*nz, MPI_REAL, ptb, tag0+7, MPI_COMM_WORLD, &rs_req_array[parity][7]));
      for(int i=0; i<8; i++){
	recv_started_array[parity][i]=0;
      }
    }

    xbound_set_parity(0, 4);
    xbound_recv_updateall(4);
    for(int req=0; req<8; req++){ // start recieving
      if (npe[req/2] != 1) {
	CHECK(QWS_PRISM_START(&rs_req[req]));
	recv_started[req]=1;
      }
    }
  }

  void xbound_finalize() { 

    // make sure all the communications are finished
    MPI_Barrier(MPI_COMM_WORLD);

    // Future receives remain preposted after the final numerical exchange.
    // PRISM cannot cancel an active buffered receive. Tell each sender which
    // receives remain active, then complete them using the existing buffers.
    // This drain runs after solver timing and never contributes to the result.
    for (int parity=0; parity<2; ++parity) {
      int drain_send[8] = {0};
      for (int req=0; req<8; ++req) {
        if (npe[req/2] == 1) continue;
        CHECK(MPI_Sendrecv(&recv_started_array[parity][req], 1, MPI_INT,
                           recv_peer[req], 32+8*parity+req,
                           &drain_send[req], 1, MPI_INT,
                           send_peer[req], 32+8*parity+req,
                           MPI_COMM_WORLD, MPI_STATUS_IGNORE));
      }
      for (int req=0; req<8; ++req) {
        if (drain_send[req]) CHECK(QWS_PRISM_START(&ss_req_array[parity][req]));
      }
      for (int req=0; req<8; ++req) {
        if (recv_started_array[parity][req]) {
          CHECK(OMPI_Wait(&rs_req_array[parity][req], MPI_STATUS_IGNORE));
          recv_started_array[parity][req] = 0;
        }
      }
      for (int req=0; req<8; ++req) {
        if (drain_send[req]) CHECK(OMPI_Wait(&ss_req_array[parity][req], MPI_STATUS_IGNORE));
      }
    }

    // Extension handles must never reach standard MPI request routines.
    for (int req=0; req<8; ++req) {
      CHECK(OMPI_Prequest_free(&sd_req[req]));
      CHECK(OMPI_Prequest_free(&rd_req[req]));
      for (int parity=0; parity<2; ++parity) {
        CHECK(OMPI_Prequest_free(&ss_req_array[parity][req]));
        CHECK(OMPI_Prequest_free(&rs_req_array[parity][req]));
      }
    }

    // send/recv buffers: double precision
    xfd_recv=xfd_recv0;  // recover from the original value
    xbd_recv=xbd_recv0;

    if(xfd_send) free(xfd_send);
    if(xfd_recv) free(xfd_recv);
    if(xbd_send) free(xbd_send);
    if(xbd_recv) free(xbd_recv);

    if(yfd_send) free(yfd_send);
    if(yfd_recv) free(yfd_recv);
    if(ybd_send) free(ybd_send);
    if(ybd_recv) free(ybd_recv);

    if(zfd_send) free(zfd_send);
    if(zfd_recv) free(zfd_recv);
    if(zbd_send) free(zbd_send);
    if(zbd_recv) free(zbd_recv);

    if(tfd_send) free(tfd_send);
    if(tfd_recv) free(tfd_recv);
    if(tbd_send) free(tbd_send);
    if(tbd_recv) free(tbd_recv);

  // for double buffering (single precision)
    if(xfs_send) free(xfs_send);
    if(xbs_send) free(xbs_send);
    if(yfs_send) free(yfs_send);
    if(ybs_send) free(ybs_send);
    if(zfs_send) free(zfs_send);
    if(zbs_send) free(zbs_send);
    if(tfs_send) free(tfs_send);
    if(tbs_send) free(tbs_send);

    for(int parity=0; parity<2; parity++){
      if(xfs_recv_array[parity]) free(xfs_recv_array[parity]);
      if(xbs_recv_array[parity]) free(xbs_recv_array[parity]);
      if(yfs_recv_array[parity]) free(yfs_recv_array[parity]);
      if(ybs_recv_array[parity]) free(ybs_recv_array[parity]);
      if(zfs_recv_array[parity]) free(zfs_recv_array[parity]);
      if(zbs_recv_array[parity]) free(zbs_recv_array[parity]);
      if(tfs_recv_array[parity]) free(tfs_recv_array[parity]);
      if(tbs_recv_array[parity]) free(tbs_recv_array[parity]);
    }

    // wait for all cleaning is done
    MPI_Barrier(MPI_COMM_WORLD);
  }


  //---------------------------------------------------------------------------------------- COMM
  void xbound_start(int req, int prec) {
    if (npe[req/2] != 1) {
      if (prec == 8) {
	CHECK(QWS_PRISM_START(&rd_req[req]));
      }
    }
  }

  //---------------------------------------------------------------------------------------- COMM
  void xbound(int req, int prec) {
    if (npe[req/2] != 1) {
      if (prec == 8) {
	CHECK(QWS_PRISM_START(&sd_req[req]));
      } else {
	CHECK(QWS_PRISM_START(&ss_req[req]));
      }
    } else {
      switch (req) {
      case 0 :
	if (prec == 8) {
	  memcpy(xfd_recv, xfd_send, sizeof(double)*12*ny*nz*nt);
	} else {
	  memcpy(xfs_recv, xfs_send, sizeof(float)*12*ny*nz*nt);
	}
	break;
      case 1 :
	if (prec == 8) {
	  memcpy(xbd_recv, xbd_send, sizeof(double)*12*ny*nz*nt);
	} else {
	  memcpy(xbs_recv, xbs_send, sizeof(float)*12*ny*nz*nt);
	}
	break;
      case 2 :
	if (prec == 8) {
	  memcpy(yfd_recv, yfd_send, sizeof(projscd_t)*nxd*nz*nt);
	} else {
	  memcpy(yfs_recv, yfs_send, sizeof(projscs_t)*nxs*nz*nt);
	}
	break;
      case 3 :
	if (prec == 8) {
	  memcpy(ybd_recv, ybd_send, sizeof(projscd_t)*nxd*nz*nt);
	} else {
	  memcpy(ybs_recv, ybs_send, sizeof(projscs_t)*nxs*nz*nt);
	}
	break;
      case 4 :
	if (prec == 8) {
	  memcpy(zfd_recv, zfd_send, sizeof(projscd_t)*nxd*ny*nt);
	} else {
	  memcpy(zfs_recv, zfs_send, sizeof(projscs_t)*nxs*ny*nt);
	}
	break;
      case 5 :
	if (prec == 8) {
	  memcpy(zbd_recv, zbd_send, sizeof(projscd_t)*nxd*ny*nt);
	} else {
	  memcpy(zbs_recv, zbs_send, sizeof(projscs_t)*nxs*ny*nt);
	}
	break;
      case 6 :
	if (prec == 8) {
	  memcpy(tfd_recv, tfd_send, sizeof(projscd_t)*nxd*ny*nz);
	} else {
	  memcpy(tfs_recv, tfs_send, sizeof(projscs_t)*nxs*ny*nz);
	}
	break;
      case 7 :
	if (prec == 8) {
	  memcpy(tbd_recv, tbd_send, sizeof(projscd_t)*nxd*ny*nz);
	} else {
	  memcpy(tbs_recv, tbs_send, sizeof(projscs_t)*nxs*ny*nz);
	}
	break;
      }
    }
  }
  //---------------------------------------------------------------------------------------- wait
  void xbound_wait(int req, int prec) {
    MPI_Status status;
    if (npe[req/2] != 1) {
      if (prec == 8) {
	CHECK(OMPI_Wait(&rd_req[req], &status));
      } else {
	CHECK(OMPI_Wait(&rs_req[req], &status));
#if defined(PROF_TARGET) && PROF_TARGET == TARGET_OVERLAP_PROBE && defined(_CHECK_TIMING)
        xbound_overlap_probe_after_req_wait(req);
#endif
	recv_started[req]=0;
	if(!recv_started_opposite[req]){
	  CHECK(QWS_PRISM_START(&rs_req_opposite[req]));
	  recv_started_opposite[req]=1;
	}
      }
    }
  }
  //---------------------------------------------------------------------------------------- reset
  void xbound_reset_comm(int req, int prec) {
    if (prec == 4) {
      if (npe[req/2] != 1){
	if(!recv_started_opposite[req]){
	  CHECK(QWS_PRISM_START(&rs_req_opposite[req]));
	  recv_started_opposite[req]=1;
	}
      }
    }
  }
  //---------------------------------------------------------------------------------------- send
  void xbound_send_waitall(int prec) {
    int i;
    MPI_Status status;
#if defined(PROF_TARGET) && PROF_TARGET == TARGET_OVERLAP_PROBE && defined(_CHECK_TIMING)
    double send_waitall_start_clock = get_clock();
    double send_dir_wait_sum = 0.0;
    double send_dir_wait_max = 0.0;
    int send_label_valid_at_start = overlap_probe_send_label_valid;
    int send_wait_max_req = -1;
    int send_wait_count = 0;
#endif
    if (prec == 8) {
      for (i=0;i<4;i++) {
		if (npe[i] != 1) {
#if defined(PROF_TARGET) && PROF_TARGET == TARGET_OVERLAP_PROBE && defined(_CHECK_TIMING)
		  double elapsed = overlap_probe_wait_send_request(&sd_req[0+2*i], &status);
		  overlap_probe_record_send_wait(0+2*i, elapsed,
                                 &send_dir_wait_sum, &send_dir_wait_max,
                                 &send_wait_max_req, &send_wait_count);
		  elapsed = overlap_probe_wait_send_request(&sd_req[1+2*i], &status);
		  overlap_probe_record_send_wait(1+2*i, elapsed,
                                 &send_dir_wait_sum, &send_dir_wait_max,
                                 &send_wait_max_req, &send_wait_count);
#else
		  CHECK(OMPI_Wait(&sd_req[0+2*i], &status));
		  CHECK(OMPI_Wait(&sd_req[1+2*i], &status));
#endif
  }
      }
    } else {
      for (i=0;i<4;i++) {
		if (npe[i] != 1) {
#if defined(PROF_TARGET) && PROF_TARGET == TARGET_OVERLAP_PROBE && defined(_CHECK_TIMING)
		  double elapsed = overlap_probe_wait_send_request(&ss_req[0+2*i], &status);
		  overlap_probe_record_send_wait(0+2*i, elapsed,
                                 &send_dir_wait_sum, &send_dir_wait_max,
                                 &send_wait_max_req, &send_wait_count);
		  elapsed = overlap_probe_wait_send_request(&ss_req[1+2*i], &status);
		  overlap_probe_record_send_wait(1+2*i, elapsed,
                                 &send_dir_wait_sum, &send_dir_wait_max,
                                 &send_wait_max_req, &send_wait_count);
#else
		  CHECK(OMPI_Wait(&ss_req[0+2*i], &status));
		  CHECK(OMPI_Wait(&ss_req[1+2*i], &status));
#endif
		}
      }
    }
#if defined(PROF_TARGET) && PROF_TARGET == TARGET_OVERLAP_PROBE && defined(_CHECK_TIMING)
    if (send_wait_count > 0) {
      double send_waitall_elapsed = get_clock() - send_waitall_start_clock;
      if (!send_label_valid_at_start) {
        overlap_probe_add_timing("overlap_probe_send_waitall_unlabeled",
                                 send_waitall_elapsed);
      }
      overlap_probe_add_send_metric("overlap_probe_send_waitall",
                                    send_waitall_elapsed);
      overlap_probe_add_send_metric("overlap_probe_send_waitall_dir_max",
                                    send_dir_wait_max);
      overlap_probe_add_send_dir_metric("overlap_probe_send_dir_wait_max",
                                        send_wait_max_req,
                                        send_dir_wait_max);
      overlap_probe_add_send_check_metric("overlap_probe_check_send_waitall_vs_dir_sum",
                                          send_waitall_elapsed - send_dir_wait_sum);
      overlap_probe_send_label_valid = 0;
      overlap_probe_send_label = "unlabeled";
    }
#endif
  }

  void xbound_recv_okall(int prec) {
    return;
  }

  void xbound_overlap_probe_set_label(const char *label) {
#if defined(PROF_TARGET) && PROF_TARGET == TARGET_OVERLAP_PROBE && defined(_CHECK_TIMING)
    if (!overlap_probe_on_master_thread()) return;
    overlap_probe_next_label = label == NULL ? "unlabeled" : label;
#endif
  }

  void xbound_overlap_probe_begin(void) {
#if defined(PROF_TARGET) && PROF_TARGET == TARGET_OVERLAP_PROBE && defined(_CHECK_TIMING)
    if (!overlap_probe_on_master_thread()) return;
    overlap_probe_active = 1;
    overlap_probe_recv_done = 0;
    overlap_probe_exposed_wait_active = 0;
    overlap_probe_post_recv_window_active = 0;
    overlap_probe_active_label = overlap_probe_next_label;
    overlap_probe_begin_clock = get_clock();
    overlap_probe_recv_active_count = 0;
    for (int req=0; req<8; req++) {
      overlap_probe_req_active[req] = (npe[req/2] != 1 && recv_started[req]);
      overlap_probe_req_done[req] = overlap_probe_req_active[req] ? 0 : 1;
      if (overlap_probe_req_active[req]) overlap_probe_recv_active_count++;
      overlap_probe_req_post_recorded[req] = 0;
      overlap_probe_req_done_clock[req] = 0.0;
      overlap_probe_req_pending_elapsed[req] = 0.0;
      overlap_probe_req_post_elapsed[req] = 0.0;
      overlap_probe_req_exposed_elapsed[req] = 0.0;
    }
    overlap_probe_done_clock = 0.0;
    overlap_probe_calc_end_clock = 0.0;
    overlap_probe_exposed_wait_start_clock = 0.0;
    overlap_probe_last_phase_clock = overlap_probe_begin_clock;
    if (overlap_probe_recv_active_count == 0) {
      overlap_probe_recv_done = 1;
    }
    overlap_probe_record_window_begin();
    if (overlap_probe_wait_only_mode()) {
      overlap_probe_add_metric("overlap_probe_mode_wait_only", 0.0);
    } else if (overlap_probe_coarse_phase_mode()) {
      overlap_probe_add_metric("overlap_probe_mode_coarse_phase", 0.0);
    } else if (overlap_probe_calc_end_only_mode()) {
      overlap_probe_add_metric("overlap_probe_mode_calc_end_only", 0.0);
    }
#endif
  }

  void xbound_overlap_probe_poll_labeled(const char *phase) {
#if defined(PROF_TARGET) && PROF_TARGET == TARGET_OVERLAP_PROBE && defined(_CHECK_TIMING)
    if (!overlap_probe_on_master_thread()) return;
    if (!overlap_probe_active) return;
    overlap_probe_record_phase(phase);
    if (overlap_probe_should_coalesce_poll(phase)) {
      overlap_probe_add_phase_metric("overlap_probe_poll_coalesced", phase, 0.0);
      return;
    }
    if (overlap_probe_calc_end_only_mode()) {
      overlap_probe_add_phase_metric("overlap_probe_poll_skipped", phase, 0.0);
      return;
    }
    if (overlap_probe_recv_done) return;
    overlap_probe_test_receives(overlap_probe_completion_phase(phase));
#endif
  }

  void xbound_overlap_probe_calc_end(void) {
#if defined(PROF_TARGET) && PROF_TARGET == TARGET_OVERLAP_PROBE && defined(_CHECK_TIMING)
    if (!overlap_probe_on_master_thread()) return;
    if (!overlap_probe_active) return;
    overlap_probe_calc_end_clock = get_clock();
    if (!overlap_probe_recv_done && !overlap_probe_wait_only_mode()) {
      overlap_probe_test_receives("calc_end");
    }
    if (!overlap_probe_recv_done) {
      for (int req=0; req<8; req++) {
        if (!overlap_probe_req_active[req] || overlap_probe_req_done[req]) continue;
        overlap_probe_req_pending_elapsed[req] =
          overlap_probe_calc_end_clock - overlap_probe_begin_clock;
        overlap_probe_add_dir_metric("overlap_probe_recv_dir_pending_window",
                                     req,
                                     overlap_probe_req_pending_elapsed[req]);
      }
      overlap_probe_add_metric("overlap_probe_recv_pending_window",
                               overlap_probe_calc_end_clock - overlap_probe_begin_clock);
      overlap_probe_exposed_wait_start_clock = get_clock();
      overlap_probe_exposed_wait_active = 1;
    }
    if (overlap_probe_post_recv_window_active) {
      overlap_probe_add_metric("overlap_probe_post_recv_window",
                               overlap_probe_calc_end_clock - overlap_probe_done_clock);
      overlap_probe_post_recv_window_active = 0;
    }
    overlap_probe_record_dir_post_windows(overlap_probe_calc_end_clock);
    overlap_probe_record_dir_observed_checks();
    overlap_probe_add_metric("overlap_probe_window",
                             overlap_probe_calc_end_clock - overlap_probe_begin_clock);
#endif
  }

  void xbound_overlap_probe_after_req_wait(int req) {
#if defined(PROF_TARGET) && PROF_TARGET == TARGET_OVERLAP_PROBE && defined(_CHECK_TIMING)
    if (!overlap_probe_on_master_thread()) return;
    if (!overlap_probe_active || overlap_probe_recv_done) return;
    overlap_probe_record_dir_done(req, "wait", get_clock());
#endif
  }

  void xbound_overlap_probe_after_wait(void) {
#if defined(PROF_TARGET) && PROF_TARGET == TARGET_OVERLAP_PROBE && defined(_CHECK_TIMING)
    if (!overlap_probe_on_master_thread()) return;
    if (!overlap_probe_active) return;
    if (!overlap_probe_recv_done) {
      overlap_probe_done_clock = get_clock();
      for (int req=0; req<8; req++) {
        overlap_probe_record_dir_done(req, "after_wait", overlap_probe_done_clock);
      }
      if (!overlap_probe_recv_done) {
        overlap_probe_record_all_done(overlap_probe_done_clock, "after_wait");
      }
    }
    overlap_probe_record_dir_observed_checks();
    overlap_probe_send_label = overlap_probe_active_label;
    overlap_probe_send_label_valid = 1;
    overlap_probe_active = 0;
    overlap_probe_recv_done = 0;
    overlap_probe_post_recv_window_active = 0;
    for (int req=0; req<8; req++) {
      overlap_probe_req_active[req] = 0;
      overlap_probe_req_done[req] = 0;
      overlap_probe_req_post_recorded[req] = 0;
      overlap_probe_req_done_clock[req] = 0.0;
      overlap_probe_req_pending_elapsed[req] = 0.0;
      overlap_probe_req_post_elapsed[req] = 0.0;
      overlap_probe_req_exposed_elapsed[req] = 0.0;
    }
    overlap_probe_active_label = "unlabeled";
#endif
  }

#ifdef __cplusplus
}
#endif

#else // _MPI_
#error _MPI_ is not defeind

#endif // _MPI_
