#ifndef MPIROUNTINES
#define MPIROUNTINES

#include <mpi.h>

namespace mpimod{
  extern int myid_w, nprocs_w;
  
  constexpr int mreq = 300;
  extern MPI_Status stat[mreq];
  extern MPI_Request req[mreq];
  
  extern MPI_Comm mpi_comm_hyd;
  extern int myid_hyd, nprocs_hyd;
  
  extern MPI_Comm comm3d;
  extern int myid, nprocs;
  
  constexpr int dimx=0, dimy=1, dimz=2;
  extern int periodic[3];
  extern int ntiles[3];
  extern int coords[3];
  
  extern bool reorder;
  extern int n1m, n1p, n2m, n2p, n3m, n3p;
  extern int nreq, nsub;
}

int InitializeMPI();


#endif //MPIROUNTINES


