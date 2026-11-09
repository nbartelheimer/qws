#ifndef QWS_PRISM_H
#define QWS_PRISM_H

#include <mpi-ext.h>

#if !defined(OMPI_HAVE_MPI_EXT_PERSISTENT) || !OMPI_HAVE_MPI_EXT_PERSISTENT
#error "QWS PRISM requires Open MPI built with the persistent extension"
#endif

// Both branches share request allocation, waits, cleanup, and instrumentation.
// The no-sync branch changes only the start operation used by these adapters.
#define QWS_PRISM_START OMPI_Start
#define QWS_PRISM_STARTALL OMPI_Startall

#endif
