#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cfloat>

#define Np 1024
#define rho 0.8
#define Nn 500
#define sqrt(Np/rho)
#define dtbdhs 0.1
#define dim 2
#define temp 0.2
#define cut 2.5
#define skin 1.0
#define delta 0.15
#define polydispersity 0.1
using namespace std;

double unif_rand(double left, double right)
{
  return left + (right - left)*rand()/RAND_MAX;
}

double gaussian_rand(void)
{
  static double iset = 0;
  static double gset;
  double fac, rsq, v1, v2;

  if (iset == 0) {
    do {
      v1 = unif_rand(-1, 1);
      v2 = unif_rand(-1, 1);
      rsq = v1*v1 + v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0*log(rsq)/rsq);

    gset = v1*fac;
    iset = 0.50;
    return v2*fac;
  } else {
    iset = 0;
    return gset;
  }
}

void set_diameter(double *a){
  for(int i=0;i<Np;i++)
    a[i]=1.0+polydispersity*gaussian_rand();
}

void ini_coord_rand(double (*x)[dim]){
  for(int i = 0 ; i < Np; i++){
    x[i][0] = L*unif_rand(0.,1.);
    x[i][1] = L*unif_rand(0.,1.);
    x[i][2] = L*unif_rand(0.,1.);
  }
}

void p_boundary(double (*x)[dim]){
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++)
      x[i][j]-=L*floor(x[i][j]/L);
}

void p_boundary_mc(double (*x)[dim],int i){
  for(int j=0;j<dim;j++)
    x[i][j]-=L*floor(x[i][j]/L);
}

void ini_matrix(double (*x)[dim]){
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++)
      x[i][j]=0.0;
}

void ini_array(double *x){
  for(int i=0;i<Np;i++)
    x[i]=0.0;
}

int f(int i,int M)
{
  int k;

  k=i;

  if(k<0)
    k+=M;
  if(k>=M)
    k-=M;

  return k;
}

void cell_list(int (*list)[Nn],double (*x)[dim],int M)
{
  int i,j,k;
  int nx,ny;
  int l,m;
  double dx,dy,r2;
  double thresh=cut*7./6.+skin;

  int (*map)[Np]=new int[M*M][Np];

  for(i=0;i<M;i++)
    for(j=0;j<M;j++)
      map[i+M*j][0]=0;

  for(i=0;i<Np;i++){
    nx=f((int)(x[i][0]*M/L),M);
    ny=f((int)(x[i][1]*M/L),M);

    for(m=ny-1;m<=ny+1;m++){
      for(l=nx-1;l<=nx+1;l++){
        map[f(l,M)+M*f(m,M)][map[f(l,M)+M*f(m,M)][0] +1]=i;
        map[f(l,M)+M*f(m,M)][0]++;
      }
    }
  }

  for(i=0;i<Np;i++){
    list[i][0]=0;
    nx = f((int)(x[i][0]*M/L),M);
    ny = f((int)(x[i][1]*M/L),M);

    for (k=1; k<=(map[nx+M*ny][0]); k++){
      j = map[nx+M*ny][k];
      if(j>i){
        dx =x[i][0] - x[j][0];
        dy =x[i][1] - x[j][1];
        
        dx-=L*floor((dx+0.5*L)/L);
        dy-=L*floor((dy+0.5*L)/L);
        r2 = dx*dx + dy*dy;
        
        if(r2<thresh*thresh){
          list[i][0]++;
          list[i][list[i][0]]=j;
        }
      }
    }
  }
  delete []map;
}

void calc_force_hs(double (*x)[dim],double (*f)[dim],double *a,double *U,int (*list)[Nn]){
  double dx,dy,dy_temp,dr2,dUr,aij,eij,dUrcut,Ucut,dr,t;
  ini_matrix(f);
  *U=0.0;
	
  for(int i=0;i<Np;i++){
    for(int j=1;j<=list[i][0];j++){
      dx=x[i][0]-x[list[i][j]][0];
      dy=x[i][1]-x[list[i][j]][1];
      dy-=L*floor((dy+0.5*L)/L);
      dx-=L*floor((dx+0.5*L)/L);
      dr2=dx*dx+dy*dy;
      aij = (a[i]+a[list[i][j]])/2.0;
      
      if(dr2<aij*aij){
        dr = sqrt(dr2);
        t=dr/aij;
        dUr=-(1.0-t)/aij;
        f[i][0]-=dUr*dx/dr;
        f[list[i][j]][0]+=dUr*dx/dr;
        f[i][1]-=dUr*dy/dr;
        f[list[i][j]][1]+=dUr*dy/dr;
        f[i][2]-=dUr*dz/dr;
        f[list[i][j]][2]+=dUr*dz/dr;
        *U += 0.5*(1.0-t)*(1.0-t);
      }
    }
  }
}

void calc_force_LJ(double (*x)[dim],double (*f)[dim],double *a,double *U,int (*list)[Nn]){
  double dx,dy,dr2,dUr,w2,w6,w12,w2cut,w6cut,w12cut,aij,dUrcut,Ucut,dr;
  ini_array(f);
  *U=0;
  for(int i=0;i<Np;i++){
    for(int j=1;j<=list[i][0];j++){
      dx=x[i][0]-x[list[i][j]][0];
      dy=x[i][1]-x[list[i][j]][1];
      dx-=L*floor((dx+0.5*L)/L);
      dy-=L*floor((dy+0.5*L)/L);
      dr2=dx*dx+dy*dy;
      aij=0.5*(a[i]+a[list[i][j]]);
      if(dr2<cut*cut*aij*aij){
        dr=sqrt(dr2);
        w2=aij*aij/dr2;
        w6=w2*w2*w2;
        w12=w6*w6;

        w2cut=1./cut/cut;
        w6cut=w2cut*w2cut*w2cut;
        w12cut=w6cut*w6cut;
        dUrcut=-48.*w12cut/(cut*aij)+24.*w6cut/(cut*aij);
        Ucut=4.*w12cut-4.*w6cut;

        dUr=(-48.*w12+24*w6)/dr2-dUrcut/dr;
        f[i][0]-=dUr*dx;
        f[list[i][j]][0]+=dUr*dx;
        f[i][1]-=dUr*dy;
        f[list[i][j]][1]+=dUr*dy;
        *U+=4.*w12-4.*w6-Ucut-dUrcut*(dr-cut*aij);
      }
    }
  }
}

void eom_langevin_hs(double (*v)[dim],double (*x)[dim],double (*f)[dim],double *a,double *U,double dt,int (*list)[Nn],double *kine){
  //double zeta=1.0;
  //*kine=0.0;
  calc_force_hs(x,f,a,&(*U),list);
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++){
      v[i][j]f[i][j];
      x[i][j]+=v[i][j]*dt;
      //*kine +=v[i][j]*v[i][j]/2.0/Np;
    }
  p_boundary(x);
}

void mc(double (*x)[dim],double *a,double *U,double temp0,int *count,double (*f)[dim],int (*list)[Nn],int *trial_count){
  double dr[2];
  double U0;
  double p;
  double *x,*y,*sigma;
  int  i = (int)(Np*unif_rand(0,1.0));
  int j = (int)(Np*unif_rand(0,1.0));
  while(j=i){
  int  j = (int)(Np*unif_rand(0,1.0));
  if(j!=i)
    break;
  }
  *x = x[i][0];
  *y = x[i][1];
  *sigma = a[i];
  U0=*U;
  if(*trial_count<5){
    for(int k=0;k<dim;k++){
      dr[k]=delta*unif_rand(-1.0,1.0);
      x[i][k]+=dr[k];
    }
    calc_force_LJ(x,f,a,&(*U),list);
    p=unif_rand(0,1.0);
    if(p > 1./exp((*U-U0)/temp0)){
      *count+=1;
      for(int k=0;k<dim;k++)
        x[i][k]-=dr[k];
      *U=U0;
    }
    *trial_count++;
  }
  if(*trial_count==5){
    for(int k=0;k<dim;k++){
      x[i][k] = x[j][k];
      a[i] = a[j];
      x[j][0] = *x;
      x[j][1] = *y;
      a[j] = *sigma;
    }
    calc_force_LJ(x,f,a,&(*U),list);
    p=unif_rand(0,1.0);
    if(p > 1./exp((*U-U0)/temp0)){
      *count+=1;
      for(int k=0;k<dim;k++){
	x[j][k] = x[i][k];
	a[i] = a[j];
	x[i][0] = *x;
        x[i][1] = *y;
        a[i] = *sigma;
      }
    *trial_count = 1;
  }
  p_boundary_mc(x,i);
}

void output_coord(double (*x)[dim],double *a){
  char filename[128];
  ofstream file;
  static int j=0; 
  sprintf(filename,"coord_mc_T%.3f_%d.dat",temp,j);
  file.open(filename); 
  for(int i=0;i<Np;i++)
    file << a[i] << " " << x[i][0] << " " << x[i][1] << endl;
  file.close();
  j++;
}

void update(double (*x_update)[dim],double (*x)[dim])
{
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++)
      x_update[i][j]=x[i][j];
}

void calc_disp_max(double *disp_max,double (*x)[dim],double (*x_update)[dim])
{
  int i;
  double dx,dy;
  double disp,disp_max0;
  disp_max0 = *disp_max;
  for(i=0;i<Np;i++){
    dx=x[i][0]-x_update[i][0];
    dy=x[i][1]-x_update[i][1];
    dx-=L*floor((dx+0.5*L)/L);
    dy-=L*floor((dy+0.5*L)/L);
    disp = dx*dx+dy*dy;
    if(disp > disp_max0)
      disp_max0 =disp;
  }
  *disp_max = disp_max0;
}

void auto_list_update(double *disp_max,double (*x)[dim],double (*x_update)[dim],int (*list)[Nn],int M)
{
  static int count=0;
  count++;
  calc_disp_max(&(*disp_max),x,x_update);
  if(*disp_max > skin*skin*0.25){
    cell_list(list,x,M);
    update(x_update,x);
    *disp_max=0.0;
    count=0;
  }
}

int main(){
  double x[Np][dim],x_update[Np][dim],v[Np][dim],f[Np][dim],kine,U,a[Np],disp_max=0.0,teq=1.e+3,t;
  int j=0,out=0,count=0,trial_count=1;
  int list[Np][Nn];
  int M=(int)(L/(cut+skin));
  
  set_diameter(a);
  ini_coord_rand(x);
  ini_matrix(v);
  cell_list(list,x,M);
  
  for(t=0.0;t<teq;t+=dtbdhs){
    auto_list_update(&disp_max,x,x_update,list,M);
    eom_langevin_hs(v,x,f,a,&U,dtbdhs,list,&kine);
  }
  
  while(j<mcstep_max*Np){
    j++;
    auto_list_update(&disp_max,x,x_update,list,M);
    mc(x,a,&U,temp,&count,f,list,&trial_count);
    if(j>out){
      output(x,a);
      out+=1000*Np;   
    }
  }
  
  return 0;
}
