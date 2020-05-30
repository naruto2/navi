#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <unistd.h>
#include "foreach.hpp"
#include "matrix.hpp"
#include "femesh.hpp"


static double B1=0.0, B2=0.0, C1=0.0, C2=0.0;

void setB1B2C1C2(double b1, double b2, double c1, double c2) {
  B1 = b1; B2 = b2; C1 = c1; C2 = c2;
}


static double a(long i){
  switch(i) {
  case 1: return  0.0;
  case 2: return  0.0;
  case 3: return  1.0;
  case 4: return  0.0;
  case 5: return  0.0;
  case 6: return  0.0;
  default: abort();
  }
  abort();
  return NAN;
}


static double b(long i){
  switch(i) {
  case 1: return -1.0;
  case 2: return  0.0;
  case 3: return -3.0;
  case 4: return  0.0;
  case 5: return  4.0;
  case 6: return  0.0;
  default: abort();
  }
  abort();
  return NAN;
}


static double c(long i){
  switch(i) {
  case 1: return  0.0;
  case 2: return -1.0;
  case 3: return -3.0;
  case 4: return  4.0;
  case 5: return  0.0;
  case 6: return  0.0;
  default: abort();
  }
  abort();
  return NAN;
}


static double d(long i){
  switch(i) {
  case 1: return  0.0;
  case 2: return  0.0;
  case 3: return  4.0;
  case 4: return -4.0;
  case 5: return -4.0;
  case 6: return  4.0;
  default: abort();
  }
  abort();
  return NAN;
}


static double e(long i){
  switch(i) {
  case 1: return  2.0;
  case 2: return  0.0;
  case 3: return  2.0;
  case 4: return  0.0;
  case 5: return -4.0;
  case 6: return  0.0;
  default: abort();
  }
  abort();
  return NAN;
}


static double f(long i){
  switch(i) {
  case 1: return  0.0;
  case 2: return  2.0;
  case 3: return  2.0;
  case 4: return -4.0;
  case 5: return  0.0;
  case 6: return  0.0;
  default: abort();
  }
  abort();
  return NAN;
}


static double alphaB(long j)
{
  switch(j){
  case 1:
  case 2:
  case 3:
  case 4:
  case 5:
  case 6: return b(j)*B1 + c(j)*B2;
  default: abort();
  }
  abort();
  return NAN;
}


static double betaB(long j)
{
  switch(j){
  case 1:
  case 2:
  case 3:
  case 4:
  case 5:
  case 6: return 2.0*e(j)*B1 + d(j)*B2;
  default: abort();
  }
  abort();
  return NAN;
}


static double gammaB(long j)
{
  switch(j){
  case 1:
  case 2:
  case 3:
  case 4:
  case 5:
  case 6: return d(j)*B1 + 2.0*f(j)*B2;
  default: abort();
  }
  abort();
  return NAN;
}


static double alphaC(long j)
{
  switch(j){
  case 1:
  case 2:
  case 3:
  case 4:
  case 5:
  case 6: return b(j)*C1 + c(j)*C2;
  default: abort();
  }
  abort();
  return NAN;
}


static double betaC(long j)
{
  switch(j){
  case 1:
  case 2:
  case 3:
  case 4:
  case 5:
  case 6: return 2.0*e(j)*C1 + d(j)*C2;
  default: abort();
  }
  abort();
  return NAN;
}


static double gammaC(long j)
{
  switch(j){
  case 1:
  case 2:
  case 3:
  case 4:
  case 5:
  case 6: return d(j)*C1 + 2.0*f(j)*C2;
  default: abort();
  }
  abort();
  return NAN;
}


static double ad(long j){
  switch(j) {
  case 1: return  0.0;
  case 2: return  0.0;
  case 3: return  1.0;
  default: abort();
  }
  abort();
  return NAN;
}


static double bd(long j){
  switch(j) {
  case 1: return  1.0;
  case 2: return  0.0;
  case 3: return -1.0;
  default: abort();
  }
  abort();
  return NAN;
}


static double cd(long j){
  switch(j) {
  case 1: return  0.0;
  case 2: return  1.0;
  case 3: return -1.0;
  default: abort();
  }
  abort();
  return NAN;
}

static double mij(long i, long j, double delta)
{
  return  (delta/180.0)*
    (
     180.0 * (  a(i)*a(j)                                                                                       )
     +60.0  * (  a(i)*b(j) + b(i)*a(j) + a(i)*c(j) + c(i)*a(j)                                                  )
     +15.0  * (  a(i)*d(j) + d(i)*a(j) + b(i)*c(j) + c(i)*b(j)                                                  )
     +30.0  * (  b(i)*b(j) + c(i)*c(j) + a(i)*e(j) + e(i)*a(j) + a(i)*f(j) + f(i)*a(j)                          )
     +18.0  * (  b(i)*e(j) + e(i)*b(j) + c(i)*f(j) + f(i)*c(j)                                                  )
     + 6.0  * (  b(i)*d(j) + d(i)*b(j) + c(i)*e(j) + e(i)*c(j) + b(i)*f(j) + f(i)*b(j) + c(i)*d(j) + d(i)*c(j)  )
     + 2.0  * (  d(i)*d(j) + e(i)*f(j) + f(i)*e(j)                                                              )
     + 3.0  * (  d(i)*e(j) + e(i)*d(j) + d(i)*f(j) + f(i)*d(j)                                                  )
     +12.0  * (  e(i)*e(j) + f(i)*f(j)                                                                          )
     );
}


static double axij(long i, long j, double *u, double B1, double B2, double B3)
{
  double Au=0.0, Bu=0.0, Cu=0.0, Du=0.0, Eu=0.0, Fu=0.0,
         alphaBj, betaBj, gammaBj, axij;

  for (long k = 1; k <= 6; k++) {
    Au += a(k)*u[k];
    Bu += b(k)*u[k];
    Cu += c(k)*u[k];
    Du += d(k)*u[k];
    Eu += e(k)*u[k];
    Fu += f(k)*u[k];
  }
  alphaBj =     b(j)*B1 +     c(j)*B2;
  betaBj  = 2.0*e(j)*B1 +     d(j)*B2;
  gammaBj =     d(j)*B1 + 2.0*f(j)*B2;

  axij =
  + 420.0 * (a(i)*Au                     ) * (3.0*alphaBj +     betaBj +     gammaBj)
  + 105.0 * (a(i)*Bu + b(i)*Au           ) * (4.0*alphaBj + 2.0*betaBj +     gammaBj)
  + 105.0 * (a(i)*Cu + c(i)*Au           ) * (4.0*alphaBj +     betaBj + 2.0*gammaBj)
  +  21.0 * (a(i)*Du + d(i)*Au           ) * (5.0*alphaBj + 2.0*betaBj + 2.0*gammaBj)
  +  21.0 * (b(i)*Cu + c(i)*Bu           ) * (5.0*alphaBj + 2.0*betaBj + 2.0*gammaBj)
  +  42.0 * (a(i)*Eu + e(i)*Au + b(i)*Bu ) * (5.0*alphaBj + 3.0*betaBj +     gammaBj)
  +  42.0 * (a(i)*Fu + f(i)*Au + c(i)*Cu ) * (5.0*alphaBj +     betaBj + 3.0*gammaBj)
  +   7.0 * (b(i)*Du + d(i)*Bu           ) * (6.0*alphaBj + 3.0*betaBj + 2.0*gammaBj)
  +   7.0 * (c(i)*Eu + e(i)*Cu           ) * (6.0*alphaBj + 3.0*betaBj + 2.0*gammaBj)
  +   7.0 * (b(i)*Fu + f(i)*Bu           ) * (6.0*alphaBj + 2.0*betaBj + 3.0*gammaBj)
  +   7.0 * (c(i)*Du + d(i)*Cu           ) * (6.0*alphaBj + 2.0*betaBj + 3.0*gammaBj)
  +  21.0 * (b(i)*Eu + e(i)*Bu           ) * (6.0*alphaBj + 4.0*betaBj +     gammaBj)
  +  21.0 * (c(i)*Fu + f(i)*Cu           ) * (6.0*alphaBj +     betaBj + 4.0*gammaBj)
  +   3.0 * (d(i)*Eu + e(i)*Du           ) * (7.0*alphaBj + 4.0*betaBj + 2.0*gammaBj)
  +   3.0 * (d(i)*Fu + f(i)*Du           ) * (7.0*alphaBj + 2.0*betaBj + 4.0*gammaBj)
  +   2.0 * (e(i)*Fu + f(i)*Eu + d(i)*Du ) * (7.0*alphaBj + 3.0*betaBj + 3.0*gammaBj)
  +  12.0 * (e(i)*Eu                     ) * (7.0*alphaBj + 5.0*betaBj +     gammaBj)
  +  12.0 * (f(i)*Fu                     ) * (7.0*alphaBj +     betaBj + 5.0*gammaBj)
  ;
  return axij/1260.0;
}



static double ayij(long i, long j, double *v, double C1, double C2, double C3)
{
  double Av=0.0, Bv=0.0, Cv=0.0, Dv=0.0, Ev=0.0, Fv=0.0,
    alphaCj, betaCj, gammaCj, ayij;

  for (long k = 1; k <= 6; k++) {
    Av += a(k)*v[k];
    Bv += b(k)*v[k];
    Cv += c(k)*v[k];
    Dv += d(k)*v[k];
    Ev += e(k)*v[k];
    Fv += f(k)*v[k];
  }
  alphaCj =     b(j)*C1 +     c(j)*C2;
  betaCj  = 2.0*e(j)*C1 +     d(j)*C2;
  gammaCj =     d(j)*C1 + 2.0*f(j)*C2;

  ayij =
    420.0 * (a(i)*Av                    ) * (3.0*alphaCj +     betaCj +     gammaCj)
  + 105.0 * (a(i)*Bv + b(i)*Av          ) * (4.0*alphaCj + 2.0*betaCj +     gammaCj)
  + 105.0 * (a(i)*Cv + c(i)*Av          ) * (4.0*alphaCj +     betaCj + 2.0*gammaCj)
  +  21.0 * (a(i)*Dv + d(i)*Av          ) * (5.0*alphaCj + 2.0*betaCj + 2.0*gammaCj)
  +  21.0 * (b(i)*Cv + c(i)*Bv          ) * (5.0*alphaCj + 2.0*betaCj + 2.0*gammaCj)
  +  42.0 * (a(i)*Ev + e(i)*Av + b(i)*Bv) * (5.0*alphaCj + 3.0*betaCj +     gammaCj)
  +  42.0 * (a(i)*Fv + f(i)*Av + c(i)*Cv) * (5.0*alphaCj +     betaCj + 3.0*gammaCj)
  +   7.0 * (b(i)*Dv + d(i)*Bv          ) * (6.0*alphaCj + 3.0*betaCj + 2.0*gammaCj)
  +   7.0 * (c(i)*Ev + e(i)*Cv          ) * (6.0*alphaCj + 3.0*betaCj + 2.0*gammaCj)
  +   7.0 * (b(i)*Fv + f(i)*Bv          ) * (6.0*alphaCj + 2.0*betaCj + 3.0*gammaCj)
  +   7.0 * (c(i)*Dv + d(i)*Cv          ) * (6.0*alphaCj + 2.0*betaCj + 3.0*gammaCj)
  +  21.0 * (b(i)*Ev + e(i)*Bv          ) * (6.0*alphaCj + 4.0*betaCj +     gammaCj)
  +  21.0 * (c(i)*Fv + f(i)*Cv          ) * (6.0*alphaCj +     betaCj + 4.0*gammaCj)
  +   3.0 * (d(i)*Ev + e(i)*Dv          ) * (7.0*alphaCj + 4.0*betaCj + 2.0*gammaCj)
  +   3.0 * (d(i)*Fv + f(i)*Dv          ) * (7.0*alphaCj + 2.0*betaCj + 4.0*gammaCj)
  +   2.0 * (e(i)*Fv + f(i)*Ev + d(i)*Dv) * (7.0*alphaCj + 3.0*betaCj + 3.0*gammaCj)
  +  12.0 * (e(i)*Ev                    ) * (7.0*alphaCj + 5.0*betaCj +     gammaCj)
  +  12.0 * (f(i)*Fv                    ) * (7.0*alphaCj +     betaCj + 5.0*gammaCj)
  ;
  return ayij/1260.0;
}


static double dij(long i, long j)
{
  return (1.0/12.0) * ( 12.0  *(alphaB(i)*alphaB(j) + alphaC(i)*alphaC(j))
                          + 4.0 *(alphaB(i)* betaB(j) +  betaB(i)*alphaB(j) + alphaB(i)*gammaB(j) + gammaB(i)*alphaB(j))
                          + 4.0 *(alphaC(i)* betaC(j) +  betaC(i)*alphaC(j) + alphaC(i)*gammaC(j) + gammaC(i)*alphaC(j))
                          + 1.0 *( betaB(i)*gammaB(j) + gammaB(i)* betaB(j) +  betaC(i)*gammaC(j) + gammaC(i)* betaC(j))
                          + 2.0 *( betaB(i)* betaB(j) + gammaB(i)*gammaB(j) +  betaC(i)* betaC(j) + gammaC(i)*gammaC(j)) );
}


static double hxij(long i, long j)
{
  return (1.0/12.0) * (12.0*(alphaB(i)*ad(j))
                         +4.0*( betaB(i)*ad(j) + alphaB(i)*bd(j) + gammaB(i)*ad(j) + alphaB(i)*cd(j))
                         +1.0*( betaB(i)*cd(j) + gammaB(i)*bd(j))
                         +2.0*( betaB(i)*bd(j) + gammaB(i)*cd(j))                                     );
}


static double hyij(long i, long j)
{
  return (1.0/12.0) * (12.0*(alphaC(i)*ad(j))
                         +4.0*( betaC(i)*ad(j) + alphaC(i)*bd(j) + gammaC(i)*ad(j) + alphaC(i)*cd(j))
                         +1.0*( betaC(i)*cd(j) + gammaC(i)*bd(j))
                         +2.0*( betaC(i)*bd(j) + gammaC(i)*cd(j))                                     );
}


extern void f2mesh(FILE*, std::vector<xyc>&, std::vector<nde>&);
extern long dimp2(std::vector<nde>&N);
extern double delta(int, std::vector<xyc>&, std::vector<nde>&);

using namespace std;

void printmx(matrix<double>&A){
  for (int i=1; i<A.size(); i++) for ( auto it: A[i]) { int j = it.first;
      printf("%d %d %e\n",i,j,A[i][j]);
    }
  exit(0);
}


void makeM(matrix<double>&M,vector<xyc>&Z,vector<nde>&N)
{
  int e, m=dimp2(N), n=N.size(), i, j, I, J;
  int a, b, c, A, B, C;
  double s;
  
  M.clear();
  M.resize(m+1);
  
  for (e=1; e<n; e++) {
    a=N[e].a; b=N[e].b; c=N[e].c; A=N[e].A; B=N[e].B; C=N[e].C;
    s=delta(e,Z,N);
    i=1;
    foreach(I,&a,&b,&c,&A,&B,&C){
      j=1;
      foreach(J,&a,&b,&c,&A,&B,&C){
	M[I][J] += mij(i,j,s);
	j++;
      }
      i++;
    }
  }
}


void makeAx(matrix<double>&Ax,vector<double>&U,vector<xyc>&Z,vector<nde>&N)
{
  double del, B1, B2, B3, u[7];
  int e, m, n, i, j, I, J, a, b, c, A, B, C;
  n = N.size();
  m = dimp2(N);

  Ax.clear();
  Ax.resize(m+1);
  
  for (e=1; e<n; e++) {
    a = N[e].a; b = N[e].b; c = N[e].c; A = N[e].A; B = N[e].B; C = N[e].C;

    del = delta(e,Z,N);
    i = 0;
    foreach(I, &a, &b, &c, &A, &B, &C){
      ++i; j=0;
      foreach(J, &a, &b, &c, &A, &B, &C) {
	B1 = Z[b].y - Z[c].y; B1 /= 2.0*del;
	B2 = Z[c].y - Z[a].y; B2 /= 2.0*del;
	B3 = Z[a].y - Z[b].y; B3 /= 2.0*del;
	u[1] = U[a];
	u[2] = U[b];
	u[3] = U[c];
	u[4] = U[A];
	u[5] = U[B];
	u[6] = U[C];
	Ax[I][J] += del*axij(i,++j, u, B1, B2, B3);
      }
    }
  }
}


void makeAy(matrix<double>&Ay,vector<double>&V,vector<xyc>&Z,vector<nde>&N)
{
  double del, C1, C2, C3, v[7];
  int e, m, n, i, j, I, J, a, b, c, A, B, C;
  n = N.size();
  m = dimp2(N);

  Ay.clear();
  Ay.resize(m+1);
  
  for (e=1; e<n; e++) {
    a = N[e].a; b = N[e].b; c = N[e].c; A = N[e].A; B = N[e].B; C = N[e].C;
    del = delta(e,Z,N);
    i = 0;
    foreach(I, &a, &b, &c, &A, &B, &C) {
      ++i; j=0;
      foreach(J, &a, &b, &c, &A, &B, &C) {
        C1 = Z[c].x - Z[b].x; C1 /= 2.0*del;
        C2 = Z[a].x - Z[c].x; C2 /= 2.0*del;
        C3 = Z[b].x - Z[a].x; C3 /= 2.0*del;
        v[1] = V[a+m];
        v[2] = V[b+m];
        v[3] = V[c+m];
        v[4] = V[A+m];
        v[5] = V[B+m];
        v[6] = V[C+m];
        Ay[I][J] += del*ayij(i,++j, v, C1, C2, C3);
      }
    }
  }
}


void makeD(matrix<double>&D,vector<xyc>&Z,vector<nde>&N)
{
  int e, m, n, i, j, I, J;
  int a, b, c, A, B, C;
  double del;
  
  m = dimp2(N);
  n = N.size();

  D.clear();
  D.resize(m+1);

  
  for (e=1; e<n; e++) {
    a=N[e].a; b=N[e].b; c=N[e].c; A=N[e].A; B=N[e].B; C=N[e].C;
    del = delta(e,Z,N);
    setB1B2C1C2((Z[b].y-Z[c].y)/(2.0*del), (Z[c].y-Z[a].y)/(2.0*del),
		(Z[c].x-Z[b].x)/(2.0*del), (Z[a].x-Z[c].x)/(2.0*del));

    i=1;
    foreach(I,&a,&b,&c,&A,&B,&C) {
      j=1;
      foreach(J,&a,&b,&c,&A,&B,&C) {
        D[I][J] += del*dij(i,j);
        j++;
      }
      i++;
    }
  }
}


void makeHx(matrix<double>&Hx,vector<xyc>&Z,vector<nde>&N)
{
  int e, m, n, i, j, I, J;
  int a, b, c, A, B, C;
  double del;
  
  m = dimp2(N);
  n = N.size();

  Hx.clear();
  Hx.resize(m+1);
  
  for (e=1; e<n; e++) {
    a=N[e].a; b=N[e].b; c=N[e].c; A=N[e].A; B=N[e].B; C=N[e].C;
    del = delta(e,Z,N);
    setB1B2C1C2((Z[b].y-Z[c].y)/(2.0*del), (Z[c].y-Z[a].y)/(2.0*del),
		(Z[c].x-Z[b].x)/(2.0*del), (Z[a].x-Z[c].x)/(2.0*del));

    i=1;
    foreach(I,&a,&b,&c,&A,&B,&C) {
      j=1;
      foreach(J,&a,&b,&c) {
        Hx[I][J] += del*hxij(i,j);
        j++;
      }
      i++;
    }
  }
}


void makeHy(matrix<double>&Hy,vector<xyc>&Z,vector<nde>&N)
{
  int e, m, n, i, j, I, J;
  int a, b, c, A, B, C;
  double del;
  
  m = dimp2(N);
  n = N.size();

  Hy.clear();
  Hy.resize(m+1);
  
  for (e=1; e<n; e++) {
    a=N[e].a; b=N[e].b; c=N[e].c; A=N[e].A; B=N[e].B; C=N[e].C;
    del = delta(e,Z,N);
    setB1B2C1C2((Z[b].y-Z[c].y)/(2.0*del), (Z[c].y-Z[a].y)/(2.0*del),
		(Z[c].x-Z[b].x)/(2.0*del), (Z[a].x-Z[c].x)/(2.0*del));

    i=1;
    foreach(I,&a,&b,&c,&A,&B,&C) {
      j=1;
      foreach(J,&a,&b,&c) {
        Hy[I][J] += del*hyij(i,j);
        j++;
      }
      i++;
    }
  }
}

static double staticdt=0.001;
static double staticRe=5000.0;

double tau(void)
{
  return staticdt;
}


double Re(void)
{
  return staticRe;
}


static char* border(char *s, char *t)
{
  return (strcmp(s,t)<0?s:t);
}


void makeMid(vector<xyc>&Mid,vector<xyc>&Z,vector<nde>&N) {
  int e, n=N.size();


  Mid.clear();
  Mid.resize(dimp2(N)+1);

  for(e=1; e<n; e++) {
    Mid[N[e].A].x = (Z[N[e].b].x + Z[N[e].c].x)/2.0;
    Mid[N[e].A].y = (Z[N[e].b].y + Z[N[e].c].y)/2.0;
    Mid[N[e].A].label = border(Z[N[e].b].label,Z[N[e].c].label);

    Mid[N[e].B].x = (Z[N[e].c].x + Z[N[e].a].x)/2.0;
    Mid[N[e].B].y = (Z[N[e].c].y + Z[N[e].a].y)/2.0;
    Mid[N[e].B].label = border(Z[N[e].c].label,Z[N[e].a].label);

    Mid[N[e].C].x = (Z[N[e].a].x + Z[N[e].b].x)/2.0;
    Mid[N[e].C].y = (Z[N[e].a].y + Z[N[e].b].y)/2.0;
    Mid[N[e].C].label = border(Z[N[e].a].label,Z[N[e].b].label);
  }

  for(int i=1;i<Z.size();i++){
    Mid[i].x = Z[i].x;
    Mid[i].y = Z[i].y;
    Mid[i].label = Z[i].label;
  }
}  


static double (*u)(double x, double y);
static double (*v)(double x, double y);

static char *label;

int islabel(const char *str)
{
  return !strcmp(label,str);
}

void makeA(matrix<double>&A,vector<double>&U,vector<double>&b,vector<xyc>&Z,vector<nde>&N,vector<xyc>&Mid)
{
  static matrix<double> Aa, M, Ax, Ay, D, Hx, Hy;
  int i, j, num = 2*dimp2(N)+Z.size(), m = dimp2(N);

  U.resize(num);
  A.clear(); A.resize(num);
  b.clear(); b.resize(num);

  static int init=0;
  if(!init){
    makeM(M,Z,N);
    makeD(D,Z,N);
    makeHx(Hx,Z,N);
    makeHy(Hy,Z,N);
    init = 1;
    Aa.resize(num);
    for (i=1; i<=m; i++) for (auto it : M[i]) { j = it.first;
	Aa[  i][  j] = M[i][j]/tau();
	Aa[m+i][m+j] = M[i][j]/tau();
      }

    for (i=1; i<=m; i++) for (auto it : D[i]) { j = it.first;
	Aa[  i][  j] += D[i][j]/Re();
	Aa[m+i][m+j] += D[i][j]/Re();
      }

    for (i=1; i<=m; i++) for (auto it : Hx[i]) { j = it.first;
	Aa[    i][2*m+j] = -Hx[i][j];
	Aa[2*m+j][    i] = -Hx[i][j];
      }
    
    for (i=1; i<=m; i++) for (auto it : Hy[i]) { j = it.first;
	Aa[  m+i][2*m+j] = -Hy[i][j];
	Aa[2*m+j][  m+i] = -Hy[i][j];
      }
  }

  for (i=1; i<num; i++) for (auto it : Aa[i]) { j = it.first;
      A[i][j] = Aa[i][j];
    }

  makeAx(Ax,U,Z,N);
  makeAy(Ay,U,Z,N);

  for (i=1; i<=m; i++) for (auto it : Ax[i]) { j = it.first;
      A[  i][  j] += Ax[i][j];
      A[m+i][m+j] += Ax[i][j];
    }

  for (i=1; i<=m; i++) for (auto it : Ay[i]) { j = it.first;
      A[  i][  j] += Ay[i][j];
      A[m+i][m+j] += Ay[i][j];
    }

  for (i=1; i<=m; i++) {
    b[  i] = 0.0;
    b[m+i] = 0.0;
    for (auto it:M[i] ) { j = it.first;
      b[  i] += M[i][j]*U[  j]/tau();
      b[m+i] += M[i][j]*U[m+j]/tau();
    }
  }
#if 0
  for(i=1;i<=m;i++) if(!strcmp(Mid[i].label,"v0")) {
      A[i].clear();
      A[i+m].clear();
      A[i][i] = 1.0;
      A[i+m][i+m] = 1.0;
      b[i] = 0.0;
      b[i+m] = 0.0;
    }

  for(i=1;i<=m;i++) if(!strcmp(Mid[i].label,"v1")) {
      A[i].clear();
      A[i+m].clear();
      A[i][i] = 1.0;
      A[i+m][i+m] = 1.0;
      b[i] = 0.0;
      b[i+m] = 0.0;
    }

  for(i=1;i<=m;i++) if(!strcmp(Mid[i].label,"v2")) {
      A[i].clear();
      A[i+m].clear();
      A[i][i] = 1.0;
      A[i+m][i+m] = 1.0;
      b[i] = 0.0;
      b[i+m] = 0.0;
    }

  for(i=1;i<=m;i++) if(!strcmp(Mid[i].label,"v3")) {
      A[i].clear();
      A[i+m].clear();
      A[i][i] = 1.0;
      A[i+m][i+m] = 1.0;
      b[i] = 0.0;
      b[i+m] = 0.0;
    }

  for(i=1;i<=m;i++) if(!strcmp(Mid[i].label,"e0")) {
      A[i].clear();
      A[i+m].clear();
      A[i][i] = 1.0;
      A[i+m][i+m] = 1.0;
      b[i] = 0.0;
      b[i+m] = 0.0;
    }

  for(i=1;i<=m;i++) if(!strcmp(Mid[i].label,"e1")) {
      A[i].clear();
      A[i+m].clear();
      A[i][i] = 1.0;
      A[i+m][i+m] = 1.0;
      b[i] = 0.0;
      b[i+m] = 0.0;
    }

  for(i=1;i<=m;i++) if(!strcmp(Mid[i].label,"e2")) {
      A[i].clear();
      A[i+m].clear();
      A[i][i] = 1.0;
      A[i+m][i+m] = 1.0;
      b[i] = 1.0;
      b[i+m] = 0.0;
    }

  for(i=1;i<=m;i++) if(!strcmp(Mid[i].label,"e3")) {
      A[i].clear();
      A[i+m].clear();
      A[i][i] = 1.0;
      A[i+m][i+m] = 1.0;
      b[i] = 0.0;
      b[i+m] = 0.0;
    }
#endif

  for(i=1;i<=m;i++) if (strcmp(Mid[i].label,"")){
      label=Mid[i].label;
      A[i].clear();
      A[i+m].clear();
      A[i][i] = 1.0;
      A[i+m][i+m] = 1.0;
      b[i] = u(Mid[i].x,Mid[i].y);
      b[i+m] = v(Mid[i].y,Mid[i].y);
    }
  
  A[2*m+1].clear();
  A[2*m+1][2*m+1] = 1.0;
  b[2*m+1] = 0.0;
}


static void internal_plotuv(vector<double>&U,vector<xyc>&Z,vector<nde>&N,vector<xyc>&Mid)
{
  static FILE *pp = NULL;

  if ( pp == NULL ) pp = popen("export QT_LOGGING_RULES='*=false'; gnuplot","w");
  
  double scale=0.4;
  long arrow =1 ;
  int i, m=Mid.size()-1;
  
  for (i=1; i<=m; i++)
    fprintf(pp, "set arrow %ld from %f,%f to %f,%f\n",
	    arrow++,Mid[i].x,Mid[i].y,Mid[i].x+U[i]*scale,Mid[i].y+U[i+m]*scale);
  fprintf(pp,"set size square\n");
  fprintf(pp,"set xrange [0:1]\n");
  fprintf(pp,"set yrange [0:1]\n");
  fprintf(pp,"plot '-' w l\n");
  for(int e=1;e<N.size();e++){
    fprintf(pp,"%f %f\n",Z[N[e].a].x,Z[N[e].a].y);
    fprintf(pp,"%f %f\n",Z[N[e].b].x,Z[N[e].b].y);
    fprintf(pp,"%f %f\n",Z[N[e].c].x,Z[N[e].c].y);
    fprintf(pp,"%f %f\n\n",Z[N[e].a].x,Z[N[e].a].y);
  }
  fprintf(pp,"e\n");
  fflush(pp);
}


void saveuv(vector<double>&U,vector<xyc>&Z,vector<nde>&N,vector<xyc>&Mid)
{
  static int T = 0;
  static FILE *pp = NULL;
  static char str[1000];

  sprintf(str,"gzip -f>%05d.gz",T++);
  if ( pp == NULL ) pp = popen(str,"w");
  
  double scale=0.4;
  long arrow =1 ;
  int i, m=Mid.size()-1;
  
  for (i=1; i<=m; i++)
    fprintf(pp, "set arrow %ld from %f,%f to %f,%f\n",
	    arrow++,Mid[i].x,Mid[i].y,Mid[i].x+U[i]*scale,Mid[i].y+U[i+m]*scale);
  fprintf(pp,"set size square\n");
  fprintf(pp,"set xrange [0:1]\n");
  fprintf(pp,"set yrange [0:1]\n");
  fprintf(pp,"plot '-' w l\n");
  for(int e=1;e<N.size();e++){
    fprintf(pp,"%f %f\n",Z[N[e].a].x,Z[N[e].a].y);
    fprintf(pp,"%f %f\n",Z[N[e].b].x,Z[N[e].b].y);
    fprintf(pp,"%f %f\n",Z[N[e].c].x,Z[N[e].c].y);
    fprintf(pp,"%f %f\n\n",Z[N[e].a].x,Z[N[e].a].y);
  }
  fprintf(pp,"e\n");
  fflush(pp);
  pclose(pp);
  pp = NULL;
}



vector<double> sparse__bicgstab(matrix<double>&, vector<double>&);
  

void swapcolumn(matrix<double>&A,int p, int q)
{
  for (int j=0; j<A.size(); j++) swap(A[j][p],A[j][q]);
}


void matrixreorder(map<int,int>&Aindex,matrix<double>&A)
{
  int m, n, i, j, k;

  /* A[k][k]以降はゼロ成分 */
  for ( k=1;A[k][k]!=0.0;k++ );

  set<int> NG;
  NG.clear();
  double eps = 0.001, max = 0.0;
  n = A.size();
  int Notfound;

  for( ; k<n; k++) 
    {
      for ( auto it : NG ) printf("%d ",it);
      printf("\n");
      
      max = -10.0;
      for ( auto it : A[k] ) {
	i = it.first;
	if ( NG.find(i) == NG.end() ) if ( max < A[k][i] ) max = A[k][i];
      }
      for ( auto it: A[k] ) {
	Notfound = 0;
	i = it.first;
	if ( NG.find(i) == NG.end() ) {
	  printf("Not found %d\n",i);
	  if ( max * 0.4 < A[k][i] ) 
	    {
	      if ( A[i][k] > eps ) 
		{
		  NG.insert(i);
		  Aindex[i] = k;
		  printf("max = %f %d %d\n",max,i,k);
		  break;
		}
	    }
	}
      }
      printf("Notfound = %d\n",Notfound);
    }
  NG.clear();

  for(auto it: Aindex ) { i = it.first; swapcolumn(A,Aindex[i],i);}
}

void matrixreversereorder(map<int,int>&Aindex,matrix<double>&A)
{
  int m, n, i, j, k, p;

  /* A[k][k]以降はゼロ成分 */
  for ( k=1;A[k][k]!=0.0;k++ );

  set<int> NG;
  NG.clear();
  double eps = 0.001, max = 0.0;
  int binx = 0, biny = 0;
  
  m = k/2;
  p = k;
  p = A.size()-4;
  for(k=A.size()-1; p <= k; k--) 
    {
      for ( auto it : NG ) 
      max = 0.0;
      for ( j = k-1; 0 <= j; j--) {
	if ( NG.find(j) == NG.end() ) if ( max < (A[k][j]) )  max = (A[k][j]);

      }
      for ( j = k-1; m <= j; j--) {

	if ( NG.find(j) == NG.end() ) {
	  if ( max * 0.5 < (A[k][j]) ) 
	    {
	      if ( (A[j][k]) > max * 0.5 ) 
		{
		  printf("swap %f %f\n",A[k][j],A[j][k]);
		  NG.insert(j);
		  Aindex[j] = k;
		  break;
		}
	    }
	}
      }
    }
  NG.clear();

  for(auto it: Aindex ) {
    i = it.first; swapcolumn(A,Aindex[i],i);
    printf("%d %d\n",Aindex[i],i);
  }
}


void printdiag(matrix<double>&A)
{
  int k;
  for (k=0; k<A.size(); k++) if(A[k][k] == 0.0)
			       printf("A[%d][%d] = 0.0\n",k,k);
}


void fprintAindex(FILE *fp, map<int,int> Bindex)
{
  for ( auto it : Bindex )
    fprintf(fp,"%d %d\n",it.first,it.second);
  fclose(fp);
}


#include <f2c.h>
extern "C" {
int forful_(integer *ia, integer *l, integer *ma, integer *m,
	    integer *ip, integer *jp, integer *ir, integer *ic, integer *kerns,
	    integer *ier);
}


matrix<double> T(matrix<double>&A)
{
  int i, j, n = A.size();
  matrix<double> AT(n);

  for ( i = 1; i < n; i++) for ( auto it : A[i] ) {
      j = it.first;
      if ( A[i][j] != 0.0 ) AT[j][i] = A[i][j];
    }
  return AT;
}


int stwart(matrix<double>&A);
int GLU1(matrix<double>&A);
int GSLV1(matrix<double>&A, vector<double>&b);
int estiva__bicgstab(matrix<double>&,vector<double>&);

static vector<xyc> staticZ;
static vector<nde> staticN;
static vector<xyc> staticMid;

int navierstokes_init(const char *filename, double Re, double dt,
		      double (*upointer)(double x, double y),
		      double (*vpointer)(double x, double y))
{
  FILE *fp;
  u = upointer; v = vpointer;
  if ( NULL == (fp=fopen(filename,"r"))) {
      fprintf(stderr,"Can't open the mesh file %s\n",filename);
      return -1;
    }

  f2mesh(fp,staticZ,staticN);
  makeMid(staticMid,staticZ,staticN);
  staticRe = Re; staticdt = dt;
  return 0;
}


void navierstokes(matrix<double>&A, vector<double>&U, vector<double>&b)
{
    makeA(A,U,b,staticZ,staticN,staticMid);
}


void plotuv(vector<double>&U)
{
  internal_plotuv(U,staticZ,staticN,staticMid);
}


void fprintuv(vector<double>&U)
{
  saveuv(U,staticZ,staticN,staticMid);
}
