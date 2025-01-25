
#include <omp.h>
#include <mpi.h>
#include <iostream>
#include "mpiroutines.hpp"

namespace basicmod{
  const int nhymax=2000;
  const int ngridx=90;
  const int ngridy=180;
  const int ngridz=1;

  const int in = ngridx;
  const int jn = ngridy;
  const int kn = ngridz;

  const int is = 0;
  const int js = 0;
  const int ks = 0;
    
  const int ie = ngridx-1;
  const int je = ngridy-1;
  const int ke = ngridz-1;

  double S[kn][jn][in];
  
  double sumall;

}

int setmpi(){
  using namespace mpimod;
  ntiles[dimx] = 4;
  ntiles[dimy] = 2;
  ntiles[dimz] = 1;
  return 1;
}

int InitializeVariables(){
  using namespace basicmod;
 
    for (int k=ks; k<= ke;k++){
    for (int j=js; j<= je;j++){
    for (int i=is; i<= ie;i++){
      S[k][j][i] = 1.0;
    }
    }
    }

 
  return 1;
}

int SumVariables(){
  using namespace basicmod;
  using namespace mpimod;
  double sum;
  double sumsend,sumrecv;
  
  sum = 0.0;
#pragma omp parallel for collapse(3) reduction(+:sum)
  for (int k=ks; k<= ke;k++){
  for (int j=js; j<= je;j++){
  for (int i=is; i<= ie;i++){
    sum = sum + S[k][j][i];
  }
  }
  }

  sumsend = sum;
  MPI_Allreduce(&sumsend, &sumrecv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm3d);
  sumall = sumrecv;
  return 1;
}

int main(){
  using namespace std;
  using namespace basicmod;
  double time_begin, time_end;
  int threadsnum;
  InitializeMPI();
  cout<< "grid size for x y "<< ngridx << " " << ngridy << endl;

  threadsnum = omp_get_max_threads();
  cout<< "threads="<<threadsnum<< endl;
  
  InitializeVariables();

  cout<< "entering main loop" << endl;
  time_begin = omp_get_wtime();
  for (int nhy=0;nhy<nhymax;nhy++){
    SumVariables();
  }
  time_end = omp_get_wtime();

  cout<< "sum=" << sumall << endl;
  cout<< "sim time [s] = " << time_end-time_begin << endl;
  cout<< "time/count/cell = " << (time_end-time_begin)/(ngridx*ngridy)/nhymax << endl;

  
  cout<< "program has been finished"<<endl;
  return 1;
}
