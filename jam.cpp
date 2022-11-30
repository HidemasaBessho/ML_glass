include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <time.h>
#include <iomanip>
#include <float.h>
using namespace std;
#define Np 1024
#define Nn 200
#define dim 2
#define skin 0.1

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

void ini_coord_rand(double (*x)[dim],double L){
  for(int i = 0 ; i < Np; i++){
    x[i][0] = L*unif_rand(0.,1.);
    x[i][1] = L*unif_rand(0.,1.);
  }
}

void set_diameter(double *a,double r1,double r2){
  int i,p=0;
  for(i=0;i<Np;i++){
    if(p==0){
      a[i] = r1;
      p = p+1;
    }
    else{
      a[i] = r2;
      p = p-1;
    }
  }
}

void p_boundary(double (*x)[dim],double L){
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++)
      x[i][j]-=L*floor(x[i][j]/L);
}

void ini_array(double (*x)[dim]){
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++)
      x[i][j]=0.0;
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

void cell_list(int (*list)[Nn],double (*x)[dim],int M,double Rcell,double L){
  int i,j,k;
  int nx,ny;
  int l,m;
  double dx,dy,r2;

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

        if(r2<Rcell*Rcell){
          list[i][0]++;
          list[i][list[i][0]]=j;
        }
      }
    }
  }
  delete []map;
}

void calc_force_hs(double (*x)[dim],double (*f)[dim],double *a,double *U,int (*list)[Nn],double L){
  double dx,dy,dr2,dr,dUr,aij,t;
  ini_array(f);
  *U=0;
  for(int i=0;i<Np;i++){
    for(int j=1;j<=list[i][0];j++){
      dx=x[i][0]-x[list[i][j]][0];
      dy=x[i][1]-x[list[i][j]][1];

      dx-=L*floor((dx+0.5*L)/L);
      dy-=L*floor((dy+0.5*L)/L);
      dr2=dx*dx+dy*dy;
      aij = (a[i] + a[list[i][j]]) / 2.0;
      if(dr2<aij*aij){
        dr = sqrt(dr2);
        t=dr/aij;
        dUr=-(1.0-t)/aij;
        f[i][0]-=dUr*dx/dr;
        f[list[i][j]][0]+=dUr*dx/dr;
        f[i][1]-=dUr*dy/dr;
        f[list[i][j]][1]+=dUr*dy/dr;
      }
    }
  }
}

void eom_langevin_hs(double (*v)[dim],double (*x)[dim],double (*f)[dim],double *a,double *U,double dt,int (*list)[Nn],double *kine,double L){
  double zeta=1.0;
  *kine=0.0;
  calc_force_hs(x,f,a,&(*U),list,L);
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++){
      v[i][j] = f[i][j]/zeta;
      x[i][j]+=v[i][j]*dt;
      //*kine +=v[i][j]*v[i][j]/3.0/Np;                                                                                                                                                                                                                                         
    }
  p_boundary(x,L);
}

void update(double (*x_update)[dim],double (*x)[dim]){
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++)
      x_update[i][j]=x[i][j];
}

void calc_disp_max(double *disp_max,double (*x)[dim],double (*x_update)[dim],double L)
{
  double dx,dy;
  double disp;
  for(int i=0;i<Np;i++){
    dx=x[i][0]-x_update[i][0];
    dy=x[i][1]-x_update[i][1];
    dx-=L*floor((dx+0.5*L)/L);
    dy-=L*floor((dy+0.5*L)/L);
    disp = dx*dx+dy*dy;
    if(disp > *disp_max)
      *disp_max =disp;
  }
}

void auto_list_update(double *disp_max,double (*x)[dim],double (*x_update)[dim],int (*list)[Nn],int M,double Rcell,double L){
  static int count=0;
  count++;
  calc_disp_max(&(*disp_max),x,x_update,L);
  if(*disp_max > skin*skin*0.25){
    cell_list(list,x,M,Rcell,L);
    update(x_update,x);
    *disp_max=0.0;
    count=0;
  }
}

void sizechange(double (*x)[dim],double *L,double d_phi,double alla,double *phi){
  double dL;
  int j,k;
  dL = sqrt(M_PI*alla)/2.0*(1.0/sqrt(*phi+d_phi)-1.0/sqrt(*phi));
  for(j=0;j<Np;j++){
    for(k=0;k<dim;k++){
      x[j][k] *= (1.0-dL/(*L));
    }
  }
  *L += dL;
  *phi += d_phi;
}

void steepest_descent(double (*x)[dim],double (*f)[dim],int (*list)[Nn],double *a,int M,double *U,double Rcell,double L){
  double dx,dy,dy_temp,dt0=0.0001,zeta=0.0,sum_force =0.0,dxmax=0.0,dymax=0.0,dummy=0.0;
  double x0[Np][dim],v[Np][dim];
  cell_list(list,x,M,Rcell,L);

  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++){
      x0[i][j] = x[i][j];
      v[i][j] = 0.0;
    }

  for(;;){
    calc_force_hs(x,f,a,&(*U),list,L);
    sum_force=0.0;
    for(int i=0;i<Np;i++){
      for(int j=0;j<dim;j++){
        v[i][j] = f[i][j];
        x[i][j] += v[i][j]*dt0;
      }
      sum_force += sqrt(f[i][0]*f[i][0]+f[i][1]*f[i][1])/Np;
      dy = x[i][1]-x0[i][1];
      dx = x[i][0]-x0[i][0];
      dy -= L*floor((dy+0.5*L)/L);
      dx -= L*floor((dx+0.5*L)/L);
    }

    //    printf("SD: sum_force=%.16f,x=%f\n",sum_force,x[0][0]);                                                                                                                                                                                                               

    if(dx*dx+dy*dy > dxmax*dxmax+dymax*dymax){
      dxmax = dx;
      dymax = dy;
    }
    if(dxmax*dxmax+dymax*dymax > 0.5*skin*skin*0.25){
      //printf("cell update: SD %f,%f,%f\n",dxmax,dymax,dzmax);                                                                                                                                                                                                                 
      cell_list(list,x,M,Rcell,L);
      for(int i=0;i<Np;i++){
        for(int j=0;j<dim;j++){
          x0[i][j]=x[i][j];
          dxmax=0.0;
          dymax=0.0;
        }
      }
    }
    if(sum_force<5.e-2)
      break;
  }
}

void copy_array(double (*x)[dim],double (*x0)[dim]){
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++)
      x0[i][j]=x[i][j];
}

void norm_array(double *d,double (*x)[dim]){
  *d=0.0;
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++)
      *d += x[i][j]*x[i][j];
  *d = sqrt(*d);
}

int FIRE(double (*x)[dim],double (*f)[dim],int (*list)[Nn],double *a,int M,double *U,double Rcell,double L)
{
  int i,j,imax;
  double alpha = 0.1,P;
  double dt0=0.0001;
  double v[Np][dim],x0[Np][dim],dx,dy,dy_temp,dxmax,dymax,v_nor,f_nor;
  int count=0,count0=0;
  int conv=0;
  double sum_force=0.0;
  double finc=1.1,falpha=0.99,fdec=0.5;

  cell_list(list,x,M,Rcell,L);

  for(i=0;i<Np;i++){
    for(j=0;j<dim;j++){
      v[i][j]  = 0.0;
      x0[i][j] = x[i][j];
    }
  }
  P=0.0;
  dxmax = 0.0;
  dymax = 0.0;

  calc_force_hs(x,f,a,&(*U),list,L);

    for(i=0;i<Np;i++){
      for(j=0;j<dim;j++){
        v[i][j] += 0.5*f[i][j]*dt0;
      }
    }

    norm_array(&f_nor,f);
    norm_array(&v_nor,v);

    for(i=0;i<Np;i++){
      for(j=0;j<dim;j++){
        v[i][j] = (1.0-alpha)*v[i][j]+alpha*f[i][j]/(f_nor+DBL_EPSILON)*v_nor;
      }

      P +=v[i][0]*f[i][0]+v[i][1]*f[i][1]+v[i][2]*f[i][2];
      sum_force += sqrt(f[i][0]*f[i][0]+f[i][1]*f[i][1]+f[i][2]*f[i][2])/Np;

      dy = x[i][1]-x0[i][1];
      dx = x[i][0]-x0[i][0];
      dy -= L*floor((dy+0.5*L)/L);
      dx -= L*floor((dx+0.5*L)/L);

      if(dx*dx+dy*dy > dxmax*dxmax+dymax*dymax){
        dxmax = dx;
        dymax = dy;
        imax = i;
      }
    }

    if(dxmax*dxmax+dymax*dymax > 0.5*skin*skin*0.25){
      //printf("cell update: FIRE %f,%f,%f,imax=%d,U=%f\n",dxmax,dymax,dzmax,imax,*U/Np);                                                                                                                                                                                       
      cell_list(list,x,M,Rcell,L);

      for(i=0;i<Np;i++){
        for(j=0;j<dim;j++){
          x0[i][j]=x[i][j];
        }
      }

      dxmax=0.0;
      dymax=0.0;
    }

    //    printf("FIRE: sum_force=%.16f,dt=%f,x=%f, alpha=%f,P=%f,gamma=%f \n",sum_force,dt0,x[0][0],alpha,P,gamma);                                                                                                                                                            
    if(P>=0){
      conv++;
      if(conv>5){
        dt0*=finc;
        if(dt0>0.002)
          dt0=0.002;
        alpha*=falpha;
        conv=0;
      }
    }

    if(P<0){
      alpha=0.1;
      dt0*=fdec;
      conv=0.0;
      for(i=0;i<Np;i++){
        for(j=0;j<dim;j++){
          v[i][j]= 0.0;
        }
      }
    }
    if(sum_force <1.e-12||count>1e9){
      p_boundary(x,L);
      //std::cout<<"gamma="<<gamma<<"\t"<<"Iteration times ="<< count <<"\t"<<"energy="<<*U/(double)Np<<"\t"<<"stress ="<<*rfxy<<std::endl;                                                                                                                                     
      break;
    }
  }
  return 0;
}

void output(double phi,double U)
{
  char filename[128];
  ofstream file;
  sprintf(filename,"energy.dat");
  file.open(filename,ios::app);
  file << fixed << setprecision(10) << phi << " " << U << endl;
  file.close();
}

int main()
{
  double x[Np][dim],x_update[Np][dim],v[Np][dim],f[Np][dim],kine,a[Np];
  int list[Np][Nn];
  double t,U,disp_max=0.0,d_phi=1.e-5,dt=0.01;
  double r1 = 1.0,r2 = 1.4;
  double phi=0.835;
  double alla = 0.5*(Np * r1 * r1) + 0.5*(Np * r2 * r2);
  double L = sqrt(M_PI / 4.0 / phi * alla);
  double RCHK=3.0;
  int M = (int)(L / RCHK);
  double Rcell = L/M;

  set_diameter(a,r1,r2);
  ini_coord_rand(x,L);
  ini_array(v);
  cell_list(list,x,M,Rcell,L);

  for(t=0.0;t<100.0;t+=dt){
    auto_list_update(&disp_max,x,x_update,list,M,Rcell,L);
    eom_langevin_hs(v,x,f,a,&U,dt,list,&kine,L);
  }

  auto_list_update(&disp_max,x,x_update,list,M,Rcell,L);

  steepest_descent(x,f,list,a,M,&U,Rcell,L);
  FIRE(x,f,list,a,M,&U,Rcell,L);

  while(1){
    sizechange(x,&L,d_phi,alla,&phi);
    M = (int)(L / RCHK);
    Rcell = L/M;
    steepest_descent(x,f,list,a,M,&U,Rcell,L);
    FIRE(x,f,list,a,M,&U,Rcell,L);
    output(phi,U);
    if(phi>0.86){
      break;
    }
  }

  return 0;
}
