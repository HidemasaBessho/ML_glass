#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <iomanip>
using namespace std;
#define Npm 4096
#define Pm 500

int copy(double *x_update,double *y_update,double *x,double *y,int Np){
  int i;
  for(i=0;i<Np;i++){
    x_update[i]=x[i];
    y_update[i]=y[i];
  }
  return 0;
}

int calc_disp_max(double *disp_max,double *x,double *y,double *x_update,double *y_update,int Np,int L){
  int i;
  double dx,dy;
  double disp;

  for(i=0;i<Np;i++){
    dx=x[i]-x_update[i];
    if(dx>0.5*L) dx-=L;
    else if(dx<-0.5*L)dx+=L;
    
    dy=y[i]-y_update[i];
    if(dy>0.5*L) dy-=L;
    else if(dy<-0.5*L)dy+=L;
    disp = dx*dx+dy*dy;

    if(disp > *disp_max)
      *disp_max =disp;    
  }
 
  
  return 0;
}


int com_correction(double *x,double *y,double *x_corr,double *y_corr,int Np,double L){
  int i;
  double dx,dy;
  static double x0[Npm],y0[Npm];
  static bool IsFirst = true;
  if(IsFirst){
    for(i=0;i<Np;i++){
      x0[i]=x[i];
      y0[i]=y[i];
    }
    IsFirst = false;
  }
 
  for(i=0;i<Np;i++){
    dx=x[i]-x0[i];
    if(dx>0.5*L) dx-=L;
    else if(dx<-0.5*L)dx+=L;
   
    dy=y[i]-y0[i];
    if(dy>0.5*L) dy-=L;
    else if(dy<-0.5*L)dy+=L;
   
    *x_corr+=dx/Np; //center of mass displacement.x
    *y_corr+=dy/Np;
   
    x0[i]=x[i];
    y0[i]=y[i];
  }
  return 0;
}

int output_disp1(double* x,double* y,double* a,double* dx, double* dy,int Np,double T,char* argv[]){
  int i;
  char filename[128];
  ofstream file;
  sprintf(filename,"./%s/disp1_N%d_T%.2f.dat",argv[1],Np,T);
  file.open(filename);
  for(i=0;i<Np;i++){
    file << fixed << setprecision(20) << x[i] << " " << y[i] << " " << a[i] << " " << dx[i] << " " << dy[i] << " " << sqrt(dx[i]*dx[i]+dy[i]*dy[i]) << endl;
  }
  file.close();
  return 0;
}

int output_disp2(double* x,double* y,double* a,double* dr,int Np,double T,int count){
  int i;
  char filename[128];
  ofstream file;
  sprintf(filename,"disp1_N%d_T%.2d.dat",Np,T);
  file.open(filename);
  for(i=0;i<Np;i++){
    file << fixed << setprecision(20) << x[i] << " " << y[i] <<" " << a[i] << " " << dr[i] << endl;
  }
  file.close();
  return 0;
}

int output_t(double *x,double *y,double avU, double avK,int Np,double t,double time_stamp,double x_corr,double y_corr,double T,char* argv[])
{
  int i;
  char filename[64];
  char filename2[64];
  FILE *fp;
  FILE *fp2;
  
  sprintf(filename,"./%s/time_coord_N%d_T%.2f.dat",argv[1],Np,T);
  fp=fopen(filename,"a+");
  
  for(i=0;i<Np;i++)
    fprintf(fp,"%f\t%f\t%f\n",t-time_stamp,x[i]-x_corr,y[i]-y_corr);
  fclose(fp);
 
  sprintf(filename2,"./%s/time_energy_N%d_T%.2f.dat",argv[1],Np,T);
  fp2=fopen(filename2,"a+");
  fprintf(fp2,"%f\t%f\t%f\t%f\t%f\n",t-time_stamp, avK,avU,x_corr,y_corr);
  fclose(fp2);
  
  return 0;
}

double unif_rand(double left, double right)
{
  return left + (right - left) * rand() / RAND_MAX;
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
      rsq = v1 * v1 + v2 * v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0 * log(rsq) / rsq);
    gset = v1 * fac;
    iset = 0.50;
    return v2 * fac;
  }
  else {
    iset = 0;
    return gset;
  }
}

int ini_coord_rand(double* x, double* y, double* a, int Np, double L, double r1, double r2){
  int k, p;
  p = 0;
  for (k = 0; k < Np; k++) {
    x[k] = unif_rand(0, 1) * L;
    y[k] = unif_rand(0, 1) * L;
    if (p == 0) {
      a[k] = r1;
      p = p + 1;
    }
    else {
      a[k] = r2;
      p = p - 1;
    }
  }
  return 0;
}

int ini(double* vx, double* vy, int Np) {
  int j;
  for (j = 0; j < Np; j++) {
    vx[j] = 0.0;
    vy[j] = 0.0;
  }
  return 0;
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

void update(double L,int Np,double *x,double *y,int M,double RCHK,int (*list)[Pm])
{
  int i,j,k;
  int nx,ny;
  int l,m;
  double dx,dy,r;
  
  int (*map)[Npm]=new int[M*M][Npm];
  
  for(i=0;i<M;i++)
    for(j=0;j<M;j++)
      map[i+M*j][0]=0;
  
  for(i=0;i<Np;i++){
    nx=f((int)(x[i]*M/L),M);
    ny=f((int)(y[i]*M/L),M);
    
    for(m=ny-1;m<=ny+1;m++){
      for(l=nx-1;l<=nx+1;l++){
        map[f(l,M)+M*f(m,M)][map[f(l,M)+M*f(m,M)][0] +1]=i;
        map[f(l,M)+M*f(m,M)][0]++;
        //printf("%d\n", map[f(l,M)+M*f(m,M)][0]);
      }
    }
  }
  for(i=0;i<Np;i++){
    list[i][0]=0;
    nx = f((int)(x[i]*M/L),M);
    ny = f((int)(y[i]*M/L),M);
    
    for (k=1; k<=(map[nx+M*ny][0]); k++){
      j = map[nx+M*ny][k];
      if(j>i){
        dx =x[i] - x[j];
        dy =y[i] - y[j];
        
        if(dx<-L/2.0)
          dx+=L;
        else if(dx> L/2.0)
          dx-=L;

        if(dy<-L/2.0){
          dy+=L;
        }
        else if(dy> L/2.0)
          dy-=L;

        r = dx*dx + dy*dy;
        
        if(r<RCHK*RCHK){
          list[i][0]++;
          list[i][list[i][0]]=j;
        }
      }
    }
  }
  delete []map;
}

int calc_force_hs(double* x, double* y, double L, int Np, double* a, double* kx, double* ky, double* avU,int (*list)[Pm]) {
  int i, j, k;
  *avU = 0.0;
  double r;
  double t, drU;
  double dx, dy;
  double aij;
  double cut;
  for (k = 0; k < Np; k++) {
    kx[k] = 0.0;
    ky[k] = 0.0;
  }
  for (i = 0; i < Np; i++){
    for (j = 1; j <=list[i][0]; j++){
      dx = x[list[i][j]] - x[i];
      dy = y[list[i][j]] - y[i];
      if (dx > (0.5 * L))
        dx -= L;
      if (dx < -(0.5 * L))
        dx += L; 
      if (dy > (0.5 * L))
        dy -= L;
      if (dy < -(0.5 * L))
        dy += L;
      aij = (a[i] + a[list[i][j]]) / 2.0;
      r = sqrt(dx * dx + dy * dy); 
      t = r / aij;
      cut = aij;
      if (r < cut) {
        drU = -(1 - t) / aij; //analytical calculation of the 1'st derivative
      }
      else {
        drU = 0.0;
        continue;
      } 
      kx[list[i][j]] -= drU * dx / r;
      kx[i] += drU * dx / r;
      ky[list[i][j]] -= drU * dy / r;
      ky[i] += drU * dy / r;
      *avU += (1 - t) * (1 - t);      
    }
  }
  *avU /= double(Np);
  return 0;
}

int calc_force(double* x, double* y, double L, int Np, double* a, double* kx, double* ky, double* avU,int (*list)[Pm]) {
  int i, j, k;
  *avU = 0.0;
  double r2;
  double w2,w4,w12,drU;
  double dx, dy;
  double aij;
  double cut;
  cut = 3.0;
  for (k = 0; k < Np; k++) {
    kx[k] = 0.0;
    ky[k] = 0.0;
  }
  for (i = 0; i < Np; i++)
    {
      for (j = 1; j <=list[i][0]; j++)
        {
          dx = x[list[i][j]] - x[i];
          dy = y[list[i][j]] - y[i];
          if (dx > (0.5 * L))
            dx -= L;
          if (dx < -(0.5 * L))
            dx += L;
          if (dy > (0.5 * L))
            dy -= L;
          if (dy < -(0.5 * L))
            dy += L;
          aij = (a[list[i][j]] + a[i]) / 2.0;
          r2 = dx * dx + dy * dy; //avoid using sqrt()  
          w2 = aij*aij / r2;
          w4=w2*w2;
          w12=w4*w4*w4;
          if (r2 < cut*cut) {
            drU = (-12.0) * w12 / r2; //analytical calculation of the 1'st derivative / r
          }
          else {
            drU = 0.0;
            continue;
          }
          kx[list[i][j]] -= drU * dx;
          kx[i] += drU * dx;
          ky[list[i][j]] -= drU * dy;
          ky[i] += drU * dy;
          *avU += w12;
        }
    }
  *avU /= double(Np);
  return 0; 
}

int eq_motion(double* x, double* y, double* vx, double* vy, double dt, double* kx, double* ky, int Np, double* avK, double T,double* a){
  double zeta,alpha;
  zeta = 1.;
  int k;
  for (k = 0; k < Np; k++) {
    if(a[k]==1.0){
      alpha=1.0;
    }
    else if(a[k]==1.4){
      alpha=1.96;
    }
    vx[k] += (-vx[k] * zeta * dt + kx[k] * dt + sqrt(2. * zeta * T * dt) * gaussian_rand())/alpha;
    vy[k] += (-vy[k] * zeta * dt + ky[k] * dt + sqrt(2. * zeta * T * dt) * gaussian_rand())/alpha;
    x[k] += vx[k] * dt;
    y[k] += vy[k] * dt;
    *avK += vx[k] * vx[k] + vy[k] * vy[k];
  }
  *avK = *avK / Np / 2.0;
  return 0;
}

int p_bound(double* x, double* y, int Np, double L) {
  int k;
  for (k = 0; k < Np; k++) {
    if (x[k] < 0.0) {
      x[k] = x[k] + L;
    }
    if (x[k] >  L) {
      x[k] = x[k] - L;
    }
    if (y[k] < 0.0) {
      y[k] = y[k] + L;
    }
    if (y[k] >  L) {
      y[k] = y[k] - L;
    }
  }
  return 0;
}

int main(int argc, char *argv[])
{
  double t, avU=0.0, avK=0.0;
  double x_corr=0.0,y_corr=0.0;
  int i,count=0,count_noupdate=0;
  double sampling_time,time_stamp=0.;
  double disp_max=0.0;
  int Np = 4096;
  double r1=1.0, r2=1.4;
  double dt = 0.01, time_max = 1.e+8; //parameters
  double sampling_time_max=10000.;
  double time_coord=1.e+5;
  double time_stable = 100.;
  double T;
  double RCHK=4.5;
  double L = sqrt(double(Np) / 0.8);
  int    M=(int)(L/RCHK);
  //  cout << "L=" << L <<" "<< "M="<<M <<endl;
 
  double* x, * y, * vx, * vy, * a, * kx, * ky,*x_update,*y_update,*dx,*dy,*x0,*y0;
  int (*list)[Pm]=new int[Npm][Pm];
 
  x = new double[Np];
  y = new double[Np];
  x_update = new double[Np];
  y_update = new double[Np];
  vx = new double[Np];
  vy = new double[Np];
  a = new double[Np];
  kx = new double[Np];
  ky = new double[Np];
  dx = new double[Np];
  dy = new double[Np];
  x0 = new double[Np];
  y0 = new double[Np];
  char filename[128];

  for(i=1;i<argc;i++)
    T = atof(argv[i]); // reading temp from shell script  atof : floot  atoi : int
  // change ensemble number : srand(ensemble number)

  ini_coord_rand(x, y, a, Np, L, r1, r2);
  ini(vx, vy, Np);  
 
  //HP

  for (t = 0.; t < time_stable; t += dt) {
    update(L,Np,x,y,M,RCHK,list);
    calc_force_hs(x, y, L, Np, a, kx, ky, &avU,list);
    eq_motion(x, y, vx, vy, dt, kx, ky, Np, &avK, 0.0, a);
    p_bound(x, y, Np, L);
    //    cout << t <<endl;
  }
 
  update(L,Np,x,y,M,RCHK,list);

  //BHHP

  sampling_time=5.*dt;
  copy(x_update,y_update,x,y,Np);
  for (t = dt; t < time_max; t += dt) {
    count++;
    calc_force(x, y, L, Np, a, kx, ky, &avU,list);
    eq_motion(x, y, vx, vy, dt, kx, ky, Np, &avK, T, a);
    com_correction(x,y,&x_corr,&y_corr,Np, L);
    p_bound(x, y, Np, L);

    //logarithmic sampling
    if(int(t/dt) == int((sampling_time + time_stamp)/dt)){
      output_t(x,y,avU,avK,Np,t,time_stamp,x_corr,y_corr,T,argv);
      sampling_time*=pow(10.,0.1);
      sampling_time=int(sampling_time/dt)*dt;
      if(sampling_time > sampling_time_max/pow(10.,0.1)){
        time_stamp=t;
        sampling_time=5.*dt;
      }
    }

    //auto update
    calc_disp_max(&disp_max,x,y,x_update,y_update,Np,L);
    //  cout <<disp_max <<endl;
    count_noupdate++;
    if(disp_max>0.8*0.8){
      update(L,Np,x,y,M,RCHK,list);
      disp_max=0.0;
      copy(x_update,y_update,x,y,Np);
      // cout <<count_noupdate <<endl;
      count_noupdate=0;
    }

    if(t==time_coord){
      for(i=0;i<Np;i++){
	x0[i] = x[i];
	y0[i] = y[i];
      }
    }

    if(t==(time_coord+dt)){
      for(i=0;i<Np;i++){
	dx[i] = x[i] - x0[i];
	dy[i] = y[i] - y0[i];
      }
      output_disp1(x,y,a,dx,dy,Np,T,argv);
    }
  }
  delete[] x;
  delete[] y;
  delete[] x_update;
  delete[] y_update;
  delete[] vx;
  delete[] vy;
  delete[] a;
  delete[] kx;
  delete[] ky;
  delete[] list;
  delete[] dx;
  delete[] dy;
  delete[] x0;
  delete[] y0;
  return 0;
}
