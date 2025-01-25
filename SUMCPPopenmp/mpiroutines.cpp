#include <mpi.h>
#include <iostream>

#include "mpiroutines.hpp"

namespace mpimod{
  int myid_w, nprocs_w;
  
  MPI_Status stat[mreq];
  MPI_Request req[mreq];
  
  MPI_Comm mpi_comm_hyd;
  int myid_hyd, nprocs_hyd;
  
  MPI_Comm comm3d;
  int myid, nprocs;
  
  int periodic[3];
  int ntiles[3];
  int coords[3];
  
  bool reorder;
  int n1m, n1p, n2m, n2p, n3m, n3p;
  int nreq, nsub;
}


int InitializeMPI(){
  using namespace std;
  using namespace mpimod;
  int key,color;
  int np_hyd;
  int argc;
  char** argv;

  // Initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs_w);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid_w  );
  // call setMPI
  periodic[dimx]= 1;//true
  periodic[dimy]= 1;//true
  periodic[dimz]= 1;//true
  if(myid_w == 0) {
    cout <<  "MPI process=" << nprocs_w << endl;
    cout << "decomposition=" << ntiles[dimx] <<" "<< ntiles[dimy] <<" "<< ntiles[dimz] << endl;
  }

  MPI_Bcast(ntiles,3,MPI_INTEGER,0,MPI_COMM_WORLD);
  MPI_Bcast(periodic,3,MPI_LOGICAL,0,MPI_COMM_WORLD);

  // Making 3D strucure
  np_hyd = ntiles[dimx]*ntiles[dimy]*ntiles[dimz];
  color = int(myid_w/np_hyd);
  key   = myid_w ;
  MPI_Comm_split(MPI_COMM_WORLD,color,key,&mpi_comm_hyd);
  MPI_Comm_size( mpi_comm_hyd, &nprocs_hyd);
  MPI_Comm_rank( mpi_comm_hyd, &myid_hyd)  ;   
  
  // Create a virtual Cartesian topology for the domain decomposition.
  MPI_Cart_create( mpi_comm_hyd, 3, ntiles, periodic, reorder, &comm3d);
  MPI_Comm_rank( comm3d, &myid);
  MPI_Comm_size( comm3d, &nprocs);

  // Find the ranks of my neighbors; find my virtual Cartesian coords.
  MPI_Cart_shift( comm3d, dimx, 1, &n1m, &n1p);
  MPI_Cart_shift( comm3d, dimy, 1, &n2m, &n2p);
  MPI_Cart_shift( comm3d, dimz, 1, &n3m, &n3p);
  
  MPI_Cart_coords( comm3d, myid, 3, coords);
      
  return 1;
  }
