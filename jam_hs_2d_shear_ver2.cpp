#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cfloat>
using namespace std;

#define Np 2048 //# of the particles
#define phi 0.9 //packing fraction
#define Nn 100 //# of the neigbour lists
#define a1 1.0 //diameter1
#define a2 1.4 //diameter2
#define L sqrt((a1*a1+a2*a2)*Np/2.0/4.0/phi)
#define dtbdhs 0.005 //dt for harmonic sphere potential
#define skin 1.0// skin size for list update

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

void ini_coord_rand(double* x, double* y, double* a){
  int k,p;
  p=0;
  for(k=0;k<Np;k++){
    x[k] = unif_rand(0, 1)*L;
    y[k] = unif_rand(0, 1)*L;
    if(p==0){
      a[k]=a1;
      p=p+1;
    }
    else{
      a[k]=a2;
      p=p-1;
    }
  }
}

void ini_vec(*x){
  for(int i=0;i<Np;i++)
    x[i] = 0.0;
}

void p_boundary(double *x,double *y){
  for(int i=0;i<Np;i++){
    x[i]-=L*floor(x[i]/L);
    y[i]-=L*floor(y[i]/L);
  }
}

void LE_boundary(double *x,double *y,double gamma){
  for(int i=0;i<Np;i++){
    if(x[i]<gamma*y[i])
      x[i]+=L;
    if(x[i]>L+gamma*y[i])
      x[i]-=L;

    if(y[i]<0.0){
      y[i]+=L;
      x[i]+=gamma*L;
    }
    if(y[i]>L){
      y[i]-=L;
      x[i]-=gamma*L;
    }
  }
}

void ini_matrix(double (*x)[dim]){
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

void cell_list(double* x, double* y, int M, double gamma,int(*list)[Nn]){
  int i,j,k;
  int nx,ny;
  int l,m;
  double dx,dy,r2,dy_temp;
  double thresh=cut*7./6.+skin;

  int (*map)[Np]=new int[M*M][Np];

  for(i=0;i<M;i++)
    for(j=0;j<M;j++)
      map[i+M*j][0]=0;

  for(i=0;i<Np;i++){
    nx=f((int)((x[i]-gamma*y[i])*M/L),M);
    ny=f((int)(y[i]*M/L),M);

    for(m=ny-1;m<=ny+1;m++){
      for(l=nx-1;l<=nx+1;l++){
	      map[f(l,M)+M*f(m,M)][map[f(l,M)+M*f(m,M)][0] +1]=i;
	      map[f(l,M)+M*f(m,M)][0]++;
      }
    }
  }
  
  for(i=0;i<Np;i++){
    list[i][0]=0;
    nx = f((int)((x[i]-gamma*y[i])*M/L),M);
    ny = f((int)(y[i]*M/L),M);

    for (k=1; k<=(map[nx+M*ny][0]); k++){
      j = map[nx+M*ny][k];
      if(j>i){
	      dx =x[i]-x[j];
	      dy =y[i]-y[j];
        dy_temp=dy;
	      dy -= L*floor((dy+0.5*L)/L);
	      dx -= gamma*L*floor((dy_temp+0.5*L)/L);
	      dx -= L*floor((dx+0.5*L)/L);
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

void calc_force_hs(double* x, double* y, double* a, double* fx, double* fy, double *U, int(*list)[Nn], double gamma, double *txy){
  int i,j;
  double r,t,f,dx,dy,aij,cut,dy_temp;
  ini_vec(fx);
  ini_vec(fy);
  *txy=0.0;
  *U = 0.0;
  for(i=0;i<Np;i++){
    for(j=1;j<=list[i][0];j++){
      dx = x[i] - x[list[i][j]];
      dy = y[i] - y[list[i][j]];
      dy_temp=dy;
	    dy -= L*floor((dy+0.5*L)/L);
	    dx -= gamma*L*floor((dy_temp+0.5*L)/L);
	    dx -= L*floor((dx+0.5*L)/L);
      aij = (a[i] + a[list[i][j]])/2.0;
      r = sqrt(dx*dx+dy*dy);
      cut = aij;
      if(r < cut){
        t = r/aij;
        f = -(1-t)/aij;
        *txy += dx*dy/r*f/L/L; //shear stress
      }
      else{
        f = 0.0;
        continue;
      }
      fx[i] -= 1.0*f*dx/r;
      fy[i] -= 1.0*f*dy/r;
      fx[list[i][j]] += 1.0*f*dx/r;
      fy[list[i][j]] += 1.0*f*dy/r;
      *U += (1 - t)*(1 - t)/2.0/double(Np);
    }
  }
}

void eq_motion_sd(double* x, double* y, double* vx, double* vy, double dt, double* fx, double* fy){
  int i;
  for(i=0;i<Np;i++){
    vx[k] = fx[k];
    vy[k] = fy[k];
    x[k] += vx[k] * dt;
    y[k] += vy[k] * dt;
  }
  p_boundary(x,y);
}

void copy_array(double *x,double *x0){
  for(int i=0;i<Np;i++)
    x0[i]=x[i];
}

void update(double *x_update,double *y_update,double *x,double *y)
{
  for(int i=0;i<Np;i++)
    x_update[i]=x[i];
    y_update[i]=y[i];
}

void calc_disp_max(double *disp_max,double *x,double *y,double *x_update,double *y_update,double gamma)
{
  double dx,dy,dy_temp;
  double disp;
  for(int i=0;i<Np;i++){
    dx=x[i]-x_update[i];
    dy=y[i]-y_update[i];
    dy_temp=dy;
    dy -= L*floor((dy+0.5*L)/L);
    dx -= gamma*L*floor((dy_temp+0.5*L)/L);
    dx -= L*floor((dx+0.5*L)/L);
    disp = dx*dx+dy*dy;
    if(disp > *disp_max)
      *disp_max =disp;
  }
}

void auto_list_update(double *disp_max,double *x,double *y,double *x_update,double *y_update,double gamma,int (*list)[Nn],int M){
  static int count=0;
  count++;
  calc_disp_max(&(*disp_max),x,y,x_update,y_update,gamma);
  if(*disp_max > 0.5*skin*skin*0.25){
    cell_list(x,y,M,gamma,list);
    update(x_update,y_update,x,y);
    *disp_max=0.0;
    count=0;
  }
}

int FIRE(int (*list)[Nn],double *x,double *y,double *x_update,double *y_update,double *fx,double *fy,double *a,double *U,double *txy,double *disp_max){
  double ux[Np],uy[Np],P,u_norm,f_tot;
  double dt=dt0,alpha=0.1;
  int count=0;
  
  // initialization
  ini_vec(vx);
  ini_vec(vy);
  cell_list(x,y,M,gamma,list);
  calc_force_hs(x,y,a,fx,fy,&(*U),list,gamma,&(*txy));
  
  for(;;){
    // initialization
    P=0.0;
    f_tot=0.0;
    
    // 1. MD integration with Velocity Verlet
    for(i=0;i<Np;i++){
      x[i]+=ux[i]*dt+0.5*fx[i]*dt*dt;
      y[i]+=uy[i]*dt+0.5*fy[i]*dt*dt;
      ux[i]+=0.5*fx[i]*dt;
      uy[i]+=0.5*fy[i]*dt;
    }
    calc_force_hs(x,y,a,fx,fy,&(*U),list,gamma,&(*txy));
    for(i=0;i<Np;i++){
      ux[i]+=0.5*fx[i]*dt;
      uy[i]+=0.5*fy[i]*dt;
    }
    LE_boundary(x,y,gamma);
    auto_list_update(&(*disp_max),x,y,x_update,y_update,gamma,list,int M);
    
    // 2. calculate power P
    for(i=0;i<Np;i++){
      u_norm=sqrt(ux[i]*ux[i]+uy[i]*uy[i]);
      f_norm=sqrt(fx[i]*fx[i]+fy[i]*fy[i]);
      f_tot += f_norm/Np;
      P+=fx[i]*vx[i]+fy[i]*vy[i];
      
      // 3. velocity modification
      ux[i] = (1.0-alpha)*ux[i]+A*fx[i]*u_norm/(f_norm+DBL_EPSILON);
      uy[i] = (1.0-alpha)*uy[i]+A*fy[i]*u_norm/(f_norm+DBL_EPSILON);
    }
    
    // 4. converge criterion
    if(f_tot<1.0e-12){
      break;
    }

    // 5. evaluate power P
    if(P>=0.0){
      count++;
      if(count>5){
        count=0;
        dt=min(1.1*dt,dt_max);
        alpha*=0.99;
      }
    }
    else{
      count=0;
      ini_vec(ux);
      ini_vec(uy);
      dt*=0.5;
      alpha=0.1;
    }
  }
}

void output_coord(double *x,double *y,double *a,double gamma,double *x0,double *y0){
  int i;
  double dx,dy,dy_temp;
  char filename[128];
  ofstream file;
  sprintf(filename,"coord_disp_gamma%.3f.csv",gamma);
  file.open(filename,ios::app);
  for(i=0;i<Np;i++){
    dx=x[i]-x0[i]-gamma*y0[i];
    dy=y[i]-y0[i];
    dy_temp=dy;
    dy -= L*floor((dy+0.5*L)/L);
    dx -= gamma*L*floor((dy_temp+0.5*L)/L);
    dx -= L*floor((dx+0.5*L)/L);
    file << fixed << setprecision(10) << dx << "," << dy << "," <<  a[i] << endl;
  }
  file.close();
}

void output_shear(double gamma,double txy,double U){
  char filename[128];
  ofstream file;
  sprintf(filename,"shear.csv");
  file.open(filename,ios::app);
  file << fixed << setprecision(10) << gamma << "," << txy << "," << U << endl;
  file.close();
}

int main()
{
  double *x,*y,*vx,*vy,*fx,*fy,*a,*x_update,*y_update,*x0,*y0;
  x = new double[Np];
  y = new double[Np];
  x_update = new double[Np];
  y_update = new double[Np];
  vx = new double[Np];
  vy = new double[Np];
  a = new double[Np];
  fx = new double[Np];
  fy = new double[Np];
  x0 = new double[Np];
  y0 = new double[Np];
  int (*list)[Nn] = new int[Np][Nn];
  double U,txy,gamma=0.0,disp_max=0.0,d_gamma=1.e-3;
  int M=(int)(L/(cut*sqrt(2.0)+skin));
  double t,time_stable=1.e+3;
  int i,j;
  char filename[128];
  ifstream file;
  char filename2[128];
  ofstream file2;
  
  ini_vec(vx);
  ini_vec(vy);
  ini_coord_rand(x,y,a);
  copy_array(x,x_update);
  copy_array(y,y_update);
  
  for(t=0.0;t<time_stable;t+=dths){
    auto_list_update(&disp_max,x,y,x_update,y_update,gamma,list,M);
    calc_force_hs(x,y,a,fx,fy,&U,list,gamma,&txy);
    eq_motion_sd(x,y,vx,vy,dths,fx,fy);
  }
  
  cell_list(x,y,M,gamma,list);
  copy_array(x,x_update);
  copy_array(y,y_update);
  
  FIRE(list,x,y,x_update,y_update,fx,fy,a,&U,&txy,&disp_max);
  
  copy_array(x,x0);
  copy_array(y,y0);
  while(gamma<0.5){
    gamma += d_gamma;
    for(i=0;i<Np;i++)
      x[i] += d_gamma*y[i];
    FIRE(list,x,y,x_update,y_update,fx,fy,a,&U,&txy,&disp_max);
    output_coord(x,y,a,gamma,x0,y0);
    output_shear(gamma,txy,U);
  }
  
  delete[] x;
  delete[] y;
  delete[] vx;
  delete[] vy;
  delete[] a;
  delete[] fx;
  delete[] fy;
  delete[] list;
  delete[] x0;
  delete[] y0;
  return 0;
}
