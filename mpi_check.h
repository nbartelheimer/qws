#ifndef QWS_MPI_CHECK_H
#define QWS_MPI_CHECK_H

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

static inline void qws_mpi_check(int error, const char *call,
                                 const char *file, int line)
{
  int initialized = 0;
  int finalized = 0;
  int rank = -1;
  int length = 0;
  char message[MPI_MAX_ERROR_STRING];

  if (error == MPI_SUCCESS) return;

  fprintf(stderr, "%s:%d: %s failed with MPI error %d\n",
          file, line, call, error);
  fflush(stderr);

  /* Check diagnostic calls directly to avoid recursive failures. Only use
     communicators and error strings while MPI is active (Init/Finalize can fail). */
  if (MPI_Initialized(&initialized) == MPI_SUCCESS && initialized &&
      MPI_Finalized(&finalized) == MPI_SUCCESS && !finalized) {
    if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS) rank = -1;
    if (MPI_Error_string(error, message, &length) == MPI_SUCCESS) {
      fprintf(stderr, "MPI rank %d: %.*s\n", rank, length, message);
    }
    fflush(stderr);
    /* Abort the whole job so peers do not hang in communication. */
    if (MPI_Abort(MPI_COMM_WORLD, error) != MPI_SUCCESS) {
      fprintf(stderr, "MPI_Abort failed; terminating locally\n");
      fflush(stderr);
    }
  }
  /* Also covers errors before initialization, after finalization, or if abort returns. */
  exit(EXIT_FAILURE);
}

#define CHECK(call) do { qws_mpi_check((call), #call, __FILE__, __LINE__); } while (0)

#endif
