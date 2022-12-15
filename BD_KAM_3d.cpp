#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cfloat>
using namespace utd;

#define Np 4000 //# of the particles
#define rho 1.2 //# density
#define Nn 1000  //# of the neigbour lists
#define L pow(Np/rho,1./3.)
#define teq 10000 //equilibration time
//#define tmax 1000 //production run time
#define dtbdhs 0.1
#define dtmd 0.001 //dt for molecular dynamics
#define dtbd 0.01 //dt for brownian dynamics
#define temp 1.0 // temperature
#define dim 3 //spatial dimension
#define cut 2.5 //potential cut off
#define skin 1.0// skin size for list update
#define d_gamma 0.001

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


void set_diameter(int *a){
  for(int i=0;i<0.2*Np;i++){
    a[2*i]=2;
    a[2*i+1]=1;
  }
  for(int i=0.4*Np;i<Np;i++)
    a[i]=1;
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
    if(x[i][2]<0.0)
      x[i][2]+=L;
    if(x[i][2]>L)
      x[i][2]-=L;
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

int com_correction(double (*x)[dim],double *x_corr,double *y_corr){
  int i,j;
  double dx,dy;
  static double X0[Np][dim];
  static bool IsFirst = true;
  if(IsFirst){
    for(i=0;i<Np;i++){
      for(j=0;j<dim;j++){
	      X0[i][j] = x[i][j];
      }
    }
    IsFirst = false;
  }

  for(i=0;i<Np;i++){
    dx = x[i][0]-X0[i][0];
    dy = x[i][1]-X0[i][1];
    dz = x[i][2]-X0[i][1];
    dx-=L*floor((dx+0.5*L)/L);
    dy-=L*floor((dy+0.5*L)/L);
    dz-=L*floor((dz+0.5*L)/L);

    *x_corr+=dx/double(Np); //center of mass displacement.x
    *y_corr+=dy/double(Np);
    *z_corr+=dz/double(Np);

    for(j=0;j<dim;j++){
      X0[i][j]=x[i][j];
    }
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

void cell_list(int (*list)[Nn],double (*x)[dim],int M,double gamma)
{
  int i,j,k;
  int nx,ny,nz;
  int l,m,n;
  double dx,dy,dy_temp,dz,r2;
  double thresh=cut+skin;
  
  int (*map)[Np+1]=new int[M*M*M][Np+1];
  
  for(i=0;i<M;i++)
    for(j=0;j<M;j++)
      for(k=0;k<M;k++)
	map[i+M*j+M*M*k][0]=0;
  
  for(i=0;i<Np;i++){
    nx=f((int)((x[i][0]-gamma*x[i][1])*M/L),M);
    ny=f((int)(x[i][1]*M/L),M);
    nz=f((int)(x[i][2]*M/L),M);
    for(n=nz-1;n<=nz+1;n++)
      for(m=ny-1;m<=ny+1;m++)
	for(l=nx-1;l<=nx+1;l++){
	  map[f(l,M)+M*f(m,M)+M*M*f(n,M)][map[f(l,M)+M*f(m,M)+M*M*f(n,M)][0] +1]=i;
	  map[f(l,M)+M*f(m,M)+M*M*f(n,M)][0]++;
	}  
  }
  
  for(i=0;i<Np;i++){
    list[i][0]=0;
    nx = f((int)((x[i][0]-gamma*x[i][1])*M/L),M);
    ny = f((int)(x[i][1]*M/L),M);
    nz = f((int)(x[i][2]*M/L),M);
    // std::cout<<nx<<std::endl;
    // printf("map=%d,%d\n",nx+M*ny+M*M*nz,map[nx+M*ny+M*M*nz][0]);

    for (k=1; k<=(map[nx+M*ny+M*M*nz][0]); k++){
      j = map[nx+M*ny+M*M*nz][k];
     
      if(j>i){
	dx =x[i][0] - x[j][0];
	dy =x[i][1] - x[j][1];
	dz =x[i][2] - x[j][2];
	dy_temp=dy;
	dy -= L*floor((dy+0.5*L)/L);
	dx -= gamma*L*floor((dy_temp+0.5*L)/L);
	dx -= L*floor((dx+0.5*L)/L);
	dz -= L*floor((dz+0.5*L)/L);

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

void calc_force_hs(double (*x)[dim],double (*f)[dim],int *a,double *U,int (*list)[Nn]){
  double dx,dy,dz,dr2,dUr,w2,w6,w12,w2cut,w6cut,w12cut,aij,eij,dUrcut,Ucut,dr,t;
  ini_matrix(f);
  *U=0;
  for(int i=0;i<Np;i++)
    for(int j=1;j<=list[i][0];j++){
      dx=x[i][0]-x[list[i][j]][0];
      dy=x[i][1]-x[list[i][j]][1];
      dz=x[i][2]-x[list[i][j]][2];
     
      dx-=L*floor((dx+0.5*L)/L);
      dy-=L*floor((dy+0.5*L)/L);
      dz-=L*floor((dz+0.5*L)/L);
      dr2=dx*dx+dy*dy+dz*dz;
      if(a[i]+a[list[i][j]] == 2){
	aij=1.0;
	eij=1.0;
  //      1.0 0.8 0.88
  //	  1.0 1.5 0.5
      }
      if(a[i]+a[list[i][j]] == 3){
	aij=0.8;
	eij=1.5;
      }
      if(a[i]+a[list[i][j]] == 4){
	aij=0.88;
	eij=0.5;
      }

      if(dr2<aij*aij){
	//	printf("%f\n",t);
	t=sqrt(dr2/aij*aij);
	dr=sqrt(dr2);
	dUr=-(1.-t)/aij;
	f[i][0]-=dUr*dx/dr;
	f[list[i][j]][0]+=dUr*dx/dr;
	f[i][1]-=dUr*dy/dr;
	f[list[i][j]][1]+=dUr*dy/dr;
	f[i][2]-=dUr*dz/dr;
	f[list[i][j]][2]+=dUr*dz/dr;
      }
    }
}

void calc_force(double (*x)[dim],double (*f)[dim],int *a,double *U,double *rfxy,int (*list)[Nn],double gamma,double *stress){
  double dx,dy,dy_temp,dz,dr2,dUr,w2,w6,w12,w2cut,w6cut,w12cut,aij,eij,dUrcut,Ucut,dr;
  double V = L*L*L;
 
  ini_matrix(f);
  ini_array(stress);
  
  *U=0,*rfxy=0.0;
  for(int i=0;i<Np;i++)
    for(int j=1;j<=list[i][0];j++){
      dx=x[i][0]-x[list[i][j]][0];
      dy=x[i][1]-x[list[i][j]][1];
      dz=x[i][2]-x[list[i][j]][2];
      dy_temp=dy;
      dy-=L*floor((dy+0.5*L)/L);
      dx-=gamma*L*floor((dy_temp+0.5*L)/L);
      dx-=L*floor((dx+0.5*L)/L);
      dz-=L*floor((dz+0.5*L)/L);
      
      dr2=dx*dx+dy*dy+dz*dz;
      
      if(a[i]+a[list[i][j]] == 2){
	aij=1.0;
	eij=1.0;
	//      1.0 0.8 0.88
	//	  1.0 0.5 1.5
      }
      if(a[i]+a[list[i][j]] == 3){
	aij=0.8;
	eij=1.5;
      }
      if(a[i]+a[list[i][j]] == 4){
	aij=0.88;
	eij=0.5;
      }
      
      if(dr2<cut*cut*aij*aij){
	dr=sqrt(dr2);
	w2=aij*aij/dr2;
	w6=w2*w2*w2;
	w12=w6*w6;
	
	w2cut=1./cut/cut;
	w6cut=w2cut*w2cut*w2cut;
	w12cut=w6cut*w6cut;
	dUrcut=-48.*eij*w12cut/(cut*aij)+24.*eij*w6cut/(cut*aij);
	Ucut=4.*eij*w12cut-4.*eij*w6cut;
	
       	dUr=(-48.*eij*w12+24*eij*w6)/dr2-dUrcut/dr;
	//	dUr=(-48.*eij*w12+24*eij*w6)/dr2;
	f[i][0] -=dUr*dx;
	f[list[i][j]][0] +=dUr*dx;
	f[i][1] -=dUr*dy;
	f[list[i][j]][1] +=dUr*dy;
	f[i][2] -=dUr*dz;
	f[list[i][j]][2] +=dUr*dz;
       	*U +=4.*eij*(w12-w6)-Ucut-dUrcut*(dr-cut*aij);
	//	*U+=4.*eij*(w12-w6)-Ucut;
	*rfxy += dUr*dx*dy/V;
	stress[i] += 0.5*dUr*dx*dy/V;
	stress[list[i][j]] += 0.5*dUr*dx*dy/V;
	//std::cout<<"f_ij="<<dUr*dx<<std::endl;
      }
    }
}

void eom_langevin_hs(double (*v)[dim],double (*x)[dim],double (*f)[dim],int *a,double *U,double dt,double temp0,int (*list)[Nn],double *kine){
  double zeta=1.0;
  double fluc=sqrt(2.*zeta*temp0*dt);
  *kine=0.0;
  calc_force_hs(x,f,a,&(*U),list);
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++){
      v[i][j]+=-zeta*v[i][j]*dt+f[i][j]*dt+fluc*gaussian_rand();
      x[i][j]+=v[i][j]*dt;
      *kine +=v[i][j]*v[i][j]/3.0/Np;
    }
  p_boundary(x);
}


void eom_langevin(double (*v)[dim],double (*x)[dim],double (*f)[dim],int *a,double *U,double dt,double temp0,int (*list)[Nn],double *kine,double *stress){
  double zeta=1.0;
  double fluc=sqrt(2.*zeta*temp0*dt);
  double dummy;
   *kine=0.0;
  calc_force(x,f,a,U,&dummy,list,0.0,stress);
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++){
      v[i][j]+=-zeta*v[i][j]*dt+f[i][j]*dt+fluc*gaussian_rand();
      x[i][j]+=v[i][j]*dt;
      *kine +=v[i][j]*v[i][j]/3.0/Np;
    }
  p_boundary(x);
}

void eom_md(double (*v)[dim],double (*x)[dim],double (*f)[dim],int *a,double *U,double dt,int (*list)[Nn],double *stress){
  double dummy=0.0;
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++){
      x[i][j]+=v[i][j]*dt+0.5*f[i][j]*dt*dt;
      v[i][j]+=0.5*f[i][j]*dt;
    }
  calc_force(x,f,a,U,&dummy,list,0.0,stress);

  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++){
      v[i][j]+=0.5*f[i][j]*dt;
    }
  p_boundary(x);
}

void output_coord(double (*x)[dim],double(*x0)[dim],int *a,double gamma,double *stress){
  double dx,dy,dy_temp,dz; 
  char filename[128];
  ofstream file;
  sprintf(filename,"coord_disp_gamma%.3f.csv",gamma);
  file.open(filename);
  file<< setprecision(6)<<"# type[i],x[i],y[i],z[i],stress_xy[i]"<<endl;
  for(int i=0;i<Np;i++){
    /*  dy = x[i][1]-x0[i][1];
    dx = x[i][0]-x0[i][0];
    dz = x[i][2]-x0[i][2];
    dy_temp=dy;
    dy -= L*floor((dy+0.5*L)/L);
    dx -= gamma*L*floor((dy_temp+0.5*L)/L);
    dx -= L*floor((dx+0.5*L)/L);
    dz -= L*floor((dz+0.5*L)/L);*/
    if(a[i]==1)
      file<< setprecision(10)<<(int)0<<","<<x[i][0]<<","<<x[i][1]<<","<<x[i][2]<<","<<stress[i]*Np<<endl;
    if(a[i]==2)
      file<< setprecision(10)<<(int)1<<","<<x[i][0]<<","<<x[i][1]<<","<<x[i][2]<<","<<stress[i]*Np<<endl;
  }
  file.close();
}

void output(int k,double (*v)[dim],double U){
  char filename[128];
  double K=0.0;
  
  ofstream file;
  sprintf(filename,"energy.dat");
  file.open(filename,ios::app); //append
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++)
      K+=0.5*v[i][j]*v[i][j];

  cout<< setprecision(6)<<k*dtmd<<"\t"<<K/Np*2./3.<<"\t"<<U/Np<<"\t"<<(K+U)/Np<<endl;
  file<< setprecision(6)<<k*dtmd<<"\t"<<K/Np*2./3.<<"\t"<<U/Np<<"\t"<<(K+U)/Np<<endl;
  file.close();
}

void output_t(double (*x)[dim],double t,double time_stamp,double x_corr,double y_corr)
{
  int i;
  char filename[128];
  ofstream file;
  sprintf(filename,"time_coord_N%d_T%.2f.dat",Np,temp);
  file.open(filename,ios::app);
  for(i=0;i<Np;i++){
    file << setprecision(6) << t-time_stamp << " " << x[i][0]-x_corr << " " << x[i][1]-y_corr << endl;
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
  double dx,dy,dz;
  double disp;
  for(int i=0;i<Np;i++){
    dx=x[i][0]-x_update[i][0];
    dy=x[i][1]-x_update[i][1];
    dz=x[i][2]-x_update[i][2];
    dx-=L*floor((dx+0.5*L)/L);
    dy-=L*floor((dy+0.5*L)/L);
    dz-=L*floor((dz+0.5*L)/L);
    disp = dx*dx+dy*dy+dz*dz;
    if(disp > *disp_max)
      *disp_max =disp;
  }
}

void auto_list_update(double *disp_max,double (*x)[dim],double (*x_update)[dim],int (*list)[Nn],int M){
  static int count=0;
  count++;
  calc_disp_max(&(*disp_max),x,x_update);
  if(*disp_max > 0.5*skin*skin*0.25){
    cell_list(list,x,M,0.0);
    update(x_update,x);
    *disp_max=0.0;
    count=0;
  }
}

int main(){
  double x[Np][dim],x0[Np][dim],x_update[Np][dim],v[Np][dim],f[Np][dim],stress[Np],xf[Np][dim];
  int list[Np][Nn],a[Np];
  double tout=0.0,U,rfxy,kine,disp_max=0.0,temp_anneal,gamma=0.0,fs=0.0,t=0.0;
  double sampling_time,time_stamp=0.0,sampling_time_max=1.e+4,t0=0.0,dT=8.33*1.e-5;
  int j=0,ens_count=0;
  int M = (int)(L/(cut*sqrt(2.0)+skin));
  // std::cout<<L/(cut*sqrt(2.0)+skin)<<std::endl;
  set_diameter(a);
  ini_coord_rand(x);
 
  ini_matrix(v);
  cell_list(list,x,M,0.0);
  //std::cout<<"L="<<L<<" "<<"M="<<M<<std::endl;
  
  while(j*dtbdhs < 100.){
    j++;
    auto_list_update(&disp_max,x,x_update,list,M);
    eom_langevin_hs(v,x,f,a,&U,dtbdhs,0.0,list,&kine);
    // std::cout<<x[0][2]<<" "<<list[0][0]<<" "<<kine<<std::endl;
  }
  
  j=0;
  while(t=0;t<tanneal;t+=dtbd){
    j++;
    temp_anneal=-dT*t+4.0;
    auto_list_update(&disp_max,x,x_update,list,M);
    eom_langevin(v,x,f,a,&U,dtbd,temp_anneal,list,&kine,stress);
    // std::cout<<kine<<" "<<U/Np<<std::endl;
    // std::cout<<x[0][2]<<" "<<f[0][0]<<" "<<kine<<std::endl;
  }
  
  j=0;
  while(j*dtbd < teq){
    j++;
    t=j*dtbd;
    auto_list_update(&disp_max,x,x_update,list,M);
    eom_langevin(v,x,f,a,&U,dtbd,temp,list,&kine,stress);
  }
  
  sampling_time=5.0*dtbd;
  for(t=dtbd;t<tmax;t+=dtbd){
    auto_list_update(&disp_max,x,x_update,list,M);
    eom_langevin(v,x,f,a,&U,dtbd,temp,list,&kine);
    com_correction(x,&x_corr,&y_corr);
    if(int(t/dtbd) == int((sampling_time + time_stamp)/dtbd)){
      output_t(x,t,time_stamp,x_corr,y_corr);
      sampling_time*=pow(10.,0.1);
      sampling_time=int(sampling_time/dtbd)*dtbd;
      if(sampling_time > sampling_time_max/pow(10.,0.1)){
        time_stamp=t;
        sampling_time=5.0*dtbd;
	      ens_count++;
      }
    }
    if(ens_count>=150){
     cout << ens_count << endl;
     break;
    }
  }
  
  return 0;
}
