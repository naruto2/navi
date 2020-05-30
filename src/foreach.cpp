#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <cstring>
#include <vector>
#include "femesh.hpp"


#ifndef _ESTIVA_ARY_H_
#define _ESTIVA_ARY_H_

static void    estiva_ary1(void** v,long n_1, size_t o);
static void    estiva_ary2(void** v,long m_1, long n_1, size_t o);
static size_t  estiva_dim0(void* v);
static long    estiva_dim1(void* v);
static long    estiva_dim2(void* v);

#define ary1(b,n_1)        estiva_ary1((void **)&b,n_1,sizeof(*b))
#define ary2(A,m_1,n_1)    estiva_ary2((void *)&A,m_1,n_1,sizeof(**A))
#define dim0(v)            estiva_dim0(v)
#define dim1(v)            estiva_dim1(v)
#define dim2(v)            estiva_dim2(v)

#endif

#ifndef _ESTIVA_STD_H_
#define _ESTIVA_STD_H_

static void  estiva_swap(void *A, void *B, int size1, int size2);
static void  estiva_cp(void *A, void *B, int size1, int size2);
static void  estiva_std_R(void *x, void *y);
static void *estiva_std_f(void *x);
static void  estiva_std_Rdestroy(void *x);
static void *estiva_std_f2(long size, void *x);

#define forall(m,i,n) for(i=m;i<=n;i++)
#define  max(x,y)     (x>y?(x):(y))
#define  min(x,y)     (x<y?(x):(y))
#define swap(a,b)     estiva_swap(&(a),&(b),sizeof(a),sizeof(b))
#define   cp(a,b)     estiva_cp(&(a),&(b),sizeof(a),sizeof(b))
#define    R(x,y)     estiva_std_R(x,y)
#define Rdestroy(x)   estiva_std_Rdestroy(x)
#define static_new(type,x)   Rnew(x,type)
#define static_free(x)       Rdestroy(x)
#define static_bind(type,x) (*(type*)estiva_std_f2(sizeof(type),x))
#endif 
using namespace std;


typedef struct{
  size_t dim0;
  long dim1, dim2;
  long padding;
} dim;

static dim* r;

static void check(void* v)
{
  if(v != NULL){
    r = (dim*)v;
    r--;
    return;
  }
  fprintf(stderr,"dim(): You need alloc!\n");
  abort();
}

static size_t estiva_dim0(void* v)
{
  check(v);
  return r->dim0;
}

static long estiva_dim1(void* v)
{
  check(v);
  return r->dim1;
}

static long estiva_dim2(void* v)
{
  check(v);
  return r->dim2;
}

static void *pointer_array[1024];
static long array_index;

static void freefunc(void)
{
  long i;
  for (i=0; i<array_index; i++) free(pointer_array[i]);
}

static void atexit_free(void)
{
  static int init;
  if ( init == 0 ) {
    init = 1;    
    atexit(freefunc);
  }
  if ( array_index > 1000 ) { printf("Too many ary\n"); abort(); }
}

static void* alloc(size_t n)
{
  void *p;
  p = calloc(1,n);
  if(p != NULL) {
    pointer_array[array_index++] = p;
    atexit_free();
    return p;
  }
  else{
    fprintf(stderr,"ary(): Can't alloc memory!\n");
    abort();
  }
  return NULL;
}

static void new_ary1(void** v, long n_1, size_t o)
{
  long n;
  dim* r;
  
  n = n_1-1;
  
  r = (dim*)alloc(sizeof(dim)+n_1*o);
  
  r->dim2 = 0;   r->dim1 = n;   r->dim0 = o;
  
  r++;
  
  *v = r;
}

static void del_ary1(void** v)
{
  dim* r;

  r = (dim*)*v;
  r--;
  free(r);

  *v = NULL;
}


static void estiva_ary1(void** v, long n_1, size_t o)
{
  if(*v == NULL){ new_ary1(v,n_1,o); return;}

  if(n_1 == dim1(*v)+1) return;
  
  del_ary1(v);

  if(n_1 != 0) new_ary1(v,n_1,o);
}


static void new_ary2(void** v, long m_1, long n_1, size_t o)
{
  long i, m, n;
  dim* r;
  char** a;

  m = m_1-1;    n = n_1-1;

  r = (dim*)alloc(sizeof(dim) + m_1*sizeof(void *));
  
  r->dim2 = m;    r->dim1 = n;    r->dim0 = o;
  
  r++;
  
  *v = r;
  a = (char**)*v;
  
  a[0] = (char*)alloc(m_1*n_1*o);
  for(i=1;i<=m;i++) a[i] = &a[i-1][n_1*o];
}


static void del_ary2(void** v)
{
  dim* r;
  char** a;
  long i, m;

  a = (char**)*v;
  m = dim2(a);
  for(i=0;i<=m;i++) free(a[i]);

  r = (dim*)*v;
  r--;
  free(r);

  *v = NULL;
}


static void estiva_ary2(void** v, long m_1, long n_1, size_t o)
{
  if(*v == NULL){ new_ary2(v,m_1,n_1,o); return;}
  
  if(m_1==dim2(*v)+1 && n_1==dim1(*v)+1) return;
  
  del_ary2(v); 

  if(m_1!=0 && n_1!=0) new_ary2(v,m_1,n_1,o); 
}

#define n static_bind(long,x)

static int top=0, limit = 1000;
static void *x_array[1000], *y_array[1000];

static void set(int i, void *x, void *y)
{
  x_array[i] = x; y_array[i] = y;
}

static int where(void *x)
{ 
  int i; 
  forall(0,i,top) if( x_array[i] == x ) break; 
  return i;
}


static void estiva_std_R(void *x, void *y)
{ 
  int i;

  if ( y == NULL ) { 
    i = where(x);
    set(i,NULL,NULL); 
    if ( i == top ) top--;
  }
  else{
    R(x,NULL);
    i = where(NULL);
    set(i,x,y); 
    if ( i == top ) {
      top++;
      if ( top >= limit ) {
        static void **tx, **ty;
        
        ary1(tx,limit); ary1(ty,limit);
        
        forall(0,i,top) tx[i] = x_array[i];
        forall(0,i,top) ty[i] = y_array[i];

        limit++;
        
        ary1(x_array,limit); ary1(y_array,limit);

        forall(0,i,top) x_array[i] = tx[i];
        forall(0,i,top) y_array[i] = ty[i];
      }
    }
  }
}


static void *estiva_std_f(void *x)
{ 
  return y_array[where(x)];
}


static void estiva_std_Rnew(void *x,size_t size)
{
  if( estiva_std_f(x) == NULL )
    R(x,calloc(1,size));
  if( estiva_std_f(x) == NULL )
    abort();
}

static void estiva_std_Rdestroy(void *x)
{
  free(estiva_std_f(x));
  R(x, NULL);
}


static void *estiva_std_f2(long size, void *x)
{
  estiva_std_Rnew(x, size);
  return estiva_std_f(x);
}


void estiva_forgammap1(long *x)
{
  n = 1;
}

int estiva_forgammap1_loop(long *x, const char *name, vector<xyc> &Z)
{
  while ( n < (long )Z.size() ) {
    if (Z[n].label && !strcmp(Z[n].label, name) ) {
      *x = n;
      n++;
      return 1;
    }
    else {
      n++;
    }
  }
  static_free(x);
  return 0;
}


int estiva_foreach(int xsize, void *x, ...)
{ 
  va_list ap;
  void *xi, *xn;
  int i; 

  va_start(ap,x); 
  xi = x;
  forall(1,i,n) 
    xi = va_arg(ap, void *);
  xn = va_arg(ap,void *);
  va_end(ap); 

  if ( xi != x ) memcpy(xi,x,xsize);
  
  if ( xn != NULL ) { 
    memcpy(x,xn,xsize);
    n++; 
    return 1;
  }
  static_free(x);
  return 0;
}

void *estiva_foreachend()
{ 
  return NULL;
}
