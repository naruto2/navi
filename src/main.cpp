#include <cstdio>
#include <vector>
#include <map>
#include "op.hpp"
#include "matrix.hpp"
#include "navi.hpp"
using namespace std;


double u(double x, double y)
{
  if ( islabel("v0")||islabel("v1")||islabel("v2")||islabel("v3") )
    return 0.0;

  if ( islabel("e0")||islabel("e1")||islabel("e3") )
    return 0.0;

  if ( islabel("e2") )
    return 1.0;

  fprintf(stderr,"A wrong label is exist.\n");
  abort();
}


double v(double x, double y)
{
  if ( islabel("v0")||islabel("v1")||islabel("v2")||islabel("v3") )
    return 0.0;

  if ( islabel("e0")||islabel("e1")||islabel("e2")||islabel("e3") )
    return 0.0;
  
  fprintf(stderr,"A wrong label is exist.\n");
  abort();
}


int main(int argc, char **argv)
{
  double Re, dt, res;

  initop(argc,argv);
  
  if(initnavi("cavity32.mesh",Re=5000,dt=0.001, u, v)) return 0;

  matrix<double> A; vector<double> U, b;

  for ( int T=0; T<= 36000000; T++) {
    fprintf(stderr,"T");
    navi(A,U,b); A[0][0] = 1.0;

    fprintf(stderr,"= %05d\n",T);
    solver(A,U,b);

    if ( defop("-test") ) {
      if ( T > 5 ) return 0;
      fprintuv(U); plotuv(U);
    } else if (!(T%100)) fprintuv(U);
  }
  return 0;
}
