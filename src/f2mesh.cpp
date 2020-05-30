#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "femesh.hpp"
using namespace std;


double delta(int i, vector<xyc>&Z, vector<nde>&N)
{
  double xi, xj, xk, yi, yj, yk;

  xi = Z[N[i].a].x;
  xj = Z[N[i].b].x;
  xk = Z[N[i].c].x;
  yi = Z[N[i].a].y;
  yj = Z[N[i].b].y;
  yk = Z[N[i].c].y;
  return (xi*yj+xj*yk+xk*yi-yi*xj-yj*xk-yk*xi)/2.0;
}



long dimp2(vector<nde>&N)
{
  long e, E;
  int n=1;
  E = N.size();
  for(e=1; e<E; e++) {
    n = max(n,N[e].A);
    n = max(n,N[e].B);
    n = max(n,N[e].C);
  }
  return n;
}


static FILE *fcopy(void){
  int c;
  FILE *fp = tmpfile();

  while ( EOF != ( c = getchar() )) putc(c,fp);
  rewind(fp);
  return fp;
}


void f2mesh(FILE *fp, vector<xyc>&Z, vector<nde>&N){
  int m=0, n=0, i;
  char *buf, *label;
  double x, y;
  int a,b,c,A,B,C;

  if ( fp == stdin ) fp = fcopy();

  buf   = (char*)calloc(sizeof(char),1000);
  label = (char*)calloc(sizeof(char),1000);

  while(strncmp(buf,"<xyc>",5)) { fgets(buf,999,fp); }

  while(1) {
    fgets(buf,999,fp);
    if(!strncmp(buf,"<nde>",5)) break;
    sscanf(buf,"%d",&i);
    n = max(n,i);
  } 
  while(1) {
    fgets(buf,999,fp);
    if (feof(fp)) break;
    sscanf(buf,"%d",&i);
    m = max(m,i);
  } 
  rewind(fp);
  Z.resize(n+1); N.resize(m+1);

  while(strncmp(buf,"<xyc>",5)) { fgets(buf,999,fp); }

  while(1) {
    fgets(buf,999,fp);
    if(!strncmp(buf,"<nde>",5)) break;
    label[0] = 0;
    sscanf(buf,"%d %lf %lf %s",&i,&x,&y,label);
    Z[i].label = strdup(label);
    Z[i].x = x;
    Z[i].y = y;
  } 
  while(1) {
    fgets(buf,999,fp);
    if (feof(fp)) break;
    sscanf(buf,"%d %ld %ld %ld %ld %ld %ld",&i,
	   &a,&b,&c,
	   &A,&B,&C);
    N[i].a = a;
    N[i].b = b;
    N[i].c = c;
    N[i].A = A;
    N[i].B = B;
    N[i].C = C;
  } 
  free(buf);
  free(label);
  fclose(fp);
}
