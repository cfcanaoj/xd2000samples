#include <omp.h>

namespace basicmod{
  const int nhymax=2000;
  const int ngridx=90;
  const int ngridy=180;
  const int ngridz=1;

  const int in = ngridx;
  const int jn = ngridy;
  const int kn = ngridz;

  const int is = 1;
  const int js = 1;
  const int ks = 1;
    
  const int ie = ngridx;
  const int je = ngridy;
  const int ke = ngridz;

  double S[kn][jn][in];


}


int InitializeVariables(){
  using namespace basicmod;
 
    for (int k=0; k< kn;k++){
    for (int j=0; j< jn;j++){
    for (int i=0; i< in;i++){
      S[k][j][i] = 1.0;
    }
    }
    }

 
  return 1;
}

int SumVariables(){
  using namespace basicmod;
  double sum;
  sum = 0.0;
#pragma omp parallel for collapse(3) reduction(+:sum)
  for (int k=0; k< kn;k++){
  for (int j=0; j< jn;j++){
  for (int i=0; i< in;i++){
    sum = sum + S[k][j][i];
  }
  }
  }
  
 
  return 1;
}

int main(){
  using namespace basicmod;
  InitializeVariables();
  for (int nhy=0;nhy<nhymax;nhy++){
    SumVariables();
  }
    
  return 1;
}
