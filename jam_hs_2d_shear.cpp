#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cfloat>

#define Np 4096
#define phi 0.9
#define Nn 1000
#define r1 1.0
#define r2 1.4
#define L sqrt(M_PI*Np*(a1*a1+a2*a2)/8.0/phi)
#define dtbd 0.1
#define dim 2
#define cut 2.5
#define skin 1.0
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

void ini_coord_rand(double (*x)[dim]){
  for(int i = 0 ; i < Np; i++){
    x[i][0] = L*unif_rand(0.,1.);
    x[i][1] = L*unif_rand(0.,1.);
    x[i][2] = L*unif_rand(0.,1.);
  }
}

void set_diameter(double *a){
  int p=0,i;
  for(i=0;i<Np;i++){
    if(p==0){
      a[i] = a1;
      p += 1;
    }
    else if(p==1){
      a[i] = a2;
      p -= 1;
    }
  }
}

void p_boundary(double (*x)[dim]){
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++)
      x[i][j]-=L*floor(x[i][j]/L);
}

void scale_sliding_blick(double (*x)[dim],double gamma)
{
  for(int i=0;i<Np;i++){
    if(x[i][0]<gamma*(x[i][1]))
      x[i][0]+=L;
    if(x[i][0]>L+gamma*(x[i][1]))
      x[i][0]-=L;

    if(x[i][1]<0.0){
      x[i][1]+=L;
      x[i][0]+=gamma*L;
    }
    if(x[i][1]>L){
      x[i][1]-=L;
      x[i][0]-=gamma*L;
    }
  }
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

void cell_list(int (*list)[Nn],double (*x)[dim],int M,double gamma)
{
  int i,j,k;
  int nx,ny;
  int l,m;
  double dx,dy,dy_temp,r2;
  double thresh=cut+skin;
  
  int (*map)[Np]=new int[M*M][Np];
  
  for(i=0;i<M;i++)
    for(j=0;j<M;j++)
	map[i+M*j][0]=0;
  
  for(i=0;i<Np;i++){
    nx=f((int)((x[i][0]-gamma*x[i][1])*M/L),M);
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
    nx = f((int)((x[i][0]-gamma*x[i][1])*M/L),M);
    ny = f((int)(x[i][1]*M/L),M);
    // std::cout<<nx<<std::endl;
    // printf("map=%d,%d\n",nx+M*ny+M*M*nz,map[nx+M*ny+M*M*nz][0]);

    for (k=1;k<=(map[nx+M*ny][0]);k++){
      j = map[nx+M*ny][k];
      if(j>i){
        dx =x[i][0] - x[j][0];
        dy =x[i][1] - x[j][1];
        dy_temp=dy;
        dy -= L*floor((dy+0.5*L)/L);
        dx -= gamma*L*floor((dy_temp+0.5*L)/L);
        dx -= L*floor((dx+0.5*L)/L);
        r2 = dx*dx + dy*dy + dz*dz;
        if(r2<thresh*thresh){
          list[i][0]++;
          list[i][list[i][0]]=j;
        }
      }
    }
  }
  delete []map;
}

void calc_force_hs(double (*x)[dim],double (*f)[dim],double *a,double *U,int (*list)[Nn],double *txy,double *stress,double gamma){
  double dx,dy,dy_temp,dr2,dUr,aij,eij,dUrcut,Ucut,dr,t;
  ini_matrix(f);
  ini_array(stress);
  *U=0.0;
  *txy=0.0;
	
  for(int i=0;i<Np;i++){
    for(int j=1;j<=list[i][0];j++){
      dx=x[i][0]-x[list[i][j]][0];
      dy=x[i][1]-x[list[i][j]][1];
      dy_temp=dy;
      dy-=L*floor((dy+0.5*L)/L);
      dx-=gamma*L*floor((dy_temp+0.5*L)/L);
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
        *txy += dx*dy/dr*dUr/L/L; //shear stress
        stress[i] += 0.5*dx*dy/dr*dUr/L/L; //local shear stress
      }
    }
  }
}

void eom_langevin_hs(double (*v)[dim],double (*x)[dim],double (*f)[dim],double *a,double *U,double dt,int (*list)[Nn],double *kine,double *txy,double *stress,double gamma){
  double zeta=1.0;
  *kine=0.0;
  calc_force_hs(x,f,a,&(*U),list,&(*txy),stress,gamma);
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++){
      v[i][j]+=-zeta*v[i][j]*dt+f[i][j]*dt;
      x[i][j]+=v[i][j]*dt;
      *kine +=v[i][j]*v[i][j]/3.0/Np;
    }
  p_boundary(x);
}

void steepest_descent(double (*x)[dim],double (*f)[dim],double gamma,int (*list)[Nn],double *a,int M,double *U,double *txy,double *stress){
  double dx,dy,dy_temp,dz,dt0=0.001,zeta=0.0,sum_force =0.0,dxmax=0.0,dymax=0.0,dzmax=0.0;
  double x0[Np][dim],v[Np][dim];
  cell_list(list,x,M,gamma);

  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++){
      x0[i][j] = x[i][j];
      v[i][j] = 0.0;
    }

  for(;;){
    calc_force_hs(x,f,a,&(*U),list,&(*txy),stress,gamma);
    sum_force=0.0;
    for(int i=0;i<Np;i++){
      for(int j=0;j<dim;j++){
        v[i][j] = f[i][j];
        x[i][j] += v[i][j]*dt0;
      }
      sum_force += sqrt(f[i][0]*f[i][0]+f[i][1]*f[i][1])/Np;
      dy = x[i][1]-x0[i][1];
      dx = x[i][0]-x0[i][0];
      dy_temp=dy;
      dy -= L*floor((dy+0.5*L)/L);
      dx -= gamma*L*floor((dy_temp+0.5*L)/L);
      dx -= L*floor((dx+0.5*L)/L);
    } 
    
    if(dx*dx+dy*dy > dxmax*dxmax+dymax*dymax){
      dxmax = dx;
      dymax = dy;
    }
    if(dxmax*dxmax+dymax*dymax > 0.8*skin*skin*0.25){
      //printf("cell update %f,%f,%f\n",dxmax,dymax,dzmax);
      cell_list(list,x,M,gamma);
      for(int i=0;i<Np;i++)
        for(int j=0;j<dim;j++)
          x0[i][j]=x[i][j];
      dxmax=0.0;
      dymax=0.0;
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

int FIRE(double (*x)[dim],double (*f)[dim],double gamma,int (*list)[Nn],double *a,int M,double *U,double *txy,double *stress)
{
  int i,j,imax;
  double alpha = 0.1,P;
  double dt0=0.0001;
  double v[Np][dim],x0[Np][dim],dx,dy,dy_temp,dz,dxmax,dymax,dzmax,v_nor,f_nor;
  int count=0,count0=0;
  int conv=0;
  double sum_force=0.0;
  double finc=1.1,falpha=0.99,fdec=0.5;
  
  cell_list(list,x,M,gamma); 

  for(i=0;i<Np;i++){
       for(j=0;j<dim;j++){
      v[i][j]  = 0.0;
      x0[i][j] = x[i][j];
    }
  }
  P=0.0;
  dxmax = 0.0;
  dymax = 0.0; //call list update
  
  calc_force_hs(x,f,a,&(*U),list,&(*txy),stress,gamma);
   
  for(;;){
    P=0.0;
    sum_force=0.0;
    count++;
    count0++;
  
    for(i=0;i<Np;i++)
      for(j=0;j<dim;j++){
        x[i][j]+=v[i][j]*dt0+0.5*f[i][j]*dt0*dt0;
        v[i][j]+=0.5*f[i][j]*dt0;
      }
    calc_force_hs(x,f,a,&(*U),list,&(*txy),stress,gamma);
  
    for(i=0;i<Np;i++)
      for(j=0;j<dim;j++)
        v[i][j] += 0.5*f[i][j]*dt0;

    norm_array(&f_nor,f);
    norm_array(&v_nor,v);     

    for(i=0;i<Np;i++){
      for(j=0;j<dim;j++)
        v[i][j] = (1.0-alpha)*v[i][j]+alpha*f[i][j]/(f_nor+DBL_EPSILON)*v_nor;
    
      P +=v[i][0]*f[i][0]+v[i][1]*f[i][1]+v[i][2]*f[i][2];      
      sum_force += sqrt(f[i][0]*f[i][0]+f[i][1]*f[i][1]+f[i][2]*f[i][2])/Np;
      
      dy = x[i][1]-x0[i][1]; 
      dx = x[i][0]-x0[i][0];
      dy_temp=dy;
      dy -= L*floor((dy+0.5*L)/L);
      dx -= gamma*L*floor((dy_temp+0.5*L)/L);
      dx -= L*floor((dx+0.5*L)/L);
      
      if(dx*dx+dy*dy > dxmax*dxmax+dymax*dymax){
        dxmax = dx;
        dymax = dy;
        imax = i;
      }
    }
    
    if(dxmax*dxmax+dymax*dymax > 0.8*skin*skin*0.25){
      //printf("cell update %f,%f,%f,imax=%d\n",dxmax,dymax,dzmax,imax);
      cell_list(list,x,M,gamma);
     
      for(i=0;i<Np;i++)
        for(j=0;j<dim;j++) 
          x0[i][j]=x[i][j];
      
      dxmax=0.0;
      dymax=0.0;
    }
    
    // printf("FIRE: sum_force=%.16f,dt=%f,x=%f, alpha=%f,P=%f,gamma=%f \n",sum_force,dt0,x[0][0],alpha,P,gamma);
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
      for(i=0;i<Np;i++)
        for(j=0;j<dim;j++)
          v[i][j]= 0.0; 
    }
    if(sum_force <1e-8||count>1e9){
      scale_sliding_blick(x,gamma);
      //cout<<"gamma="<<gamma<< " " <<"Iteration times ="<< count << endl;
      break;  
    }
  }
  return 0;
}

void output_coord_NAD(double (*x)[dim],double(*x0)[dim],double *a,double gamma){
  double dx,dy,dy_temp,dz; 
  char filename[128];
  ofstream file;
  sprintf(filename,"coord_disp_gamma%.3f.csv",gamma);
  file.open(filename);
  file << setprecision(6)<<"# type,x,y,z,dx_na,dy_na,dz_na"<< endl;
  for(int i=0;i<Np;i++){
    dy = x[i][1]-x0[i][1];
    dx = x[i][0]-x0[i][0];
    dz = x[i][2]-x0[i][2];
    dy_temp=dy;
    dy -= L*floor((dy+0.5*L)/L);
    dx -= gamma*L*floor((dy_temp+0.5*L)/L);
    dx -= L*floor((dx+0.5*L)/L);
    dz -= L*floor((dz+0.5*L)/L);
    // std::cout<<std::setprecision(15)<<x[i][0]<<" "<<x0[i][0]+gamma*x[i][1]<<" "<<x[i][2]<<" "<<x0[i][2]<<std::endl;
    file<< setprecision(10)<<a[i]<<","<<x[i][0]<<","<<x[i][1]<<","<<dx-gamma*x[i][1]<<","<<dy<<","<< endl;
  }

  file.close();
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

void auto_list_update(double *disp_max,double (*x)[dim],double (*x_update)[dim],int (*list)[Nn],int M){
  static int count=0;
  count++;
  calc_disp_max(&(*disp_max),x,x_update);
  if(*disp_max > skin*skin*0.25){
    cell_list(list,x,M,0.0);
    update(x_update,x);
    *disp_max=0.0;
    count=0;
  }
}

void output_shear(double gamma,double txy,double txy0,double U){
  char filename[128];
  ofstream file;
  sprintf(filename,"shear_stress_energy.csv");
  file.open(filename,ios::app);
  file<< setprecision(10) << gamma << "," << txy-txy0 << "," << U << endl;
  file.close();
}

int main(){
  double x[Np][dim],x_update[Np][dim],v[Np][dim],f[Np][dim],kine,a[Np],stress[Np],x0[Np][dim];
  int list[Np][Nn];
  double tout=0.0,U,disp_max=0.0,temp_anneal,gamma=0.0,txy=0.0,txy0=0.0,t=0.0,teq=100.0,d_gamma=1.e-7;
  int M=(int)(L/(cut+skin));
  set_diameter(a);
  ini_coord_rand(x);
 
  ini_array(v);
  cell_list(list,x,M,0.0);
  
  for(t=0.0;t<teq;t+=dths){
    auto_list_update(&disp_max,x,x_update,list,M);
    eom_langevin_hs(v,x,f,a,&U,dths,list,&kine,&txy,stress,gamma);
  }
  
  steepest_descent(x,f,gamma,list,a,M,&U,&txy,stress);
  FIRE(x,f,gamma,list,a,M,&U,&txy,stress);
  copy_array(x,x0);
  txy0 = txy;
  
  while(gamma<1.e-3){
    gamma += d_gamma;
    for(int i=0;i<Np;i++)
      x[i][0] += d_gamma*x[i][1];
    steepest_descent(x,f,gamma,list,a,M,&U,&txy,stress);
    FIRE(x,f,gamma,list,a,M,&U,&txy,stress);
    d_gamma *= pow(10.0,0.1);
    output_coord_NAD(x,x0,a,gamma);
    output_shear(gamma,txy,txy0,U);
  }
  
  d_gamma=1.e-3;
  
  while(gamma<0.5){
    gamma += d_gamma;
    for(int i=0;i<Np;i++)
      x[i][0] += d_gamma*x[i][1];
    steepest_descent(x,f,gamma,list,a,M,&U,&txy,stress);
    FIRE(x,f,gamma,list,a,M,&U,&txy,stress);
    d_gamma *= pow(10.0,0.1);
    output_coord_NAD(x,x0,a,gamma);
    output_shear(gamma,txy,txy0,U);
  }
  
  return 0;
}
