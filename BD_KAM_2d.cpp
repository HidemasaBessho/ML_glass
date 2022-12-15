#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cfloat>
using namespace std;

#define Np 4000 //# of the particles
#define rho 1.09 //# density
#define Nn 100  //# of the neigbour lists
#define L sqrt(Np/rho)
#define teq 4000 //equilibration time
#define tmax 100000 //production run time
#define dtmd 0.001 //dt for molecular dynamics
#define dtbd 0.01 //dt for brownian dynamics
#define temp 1.0 // temperature
#define dim 2 //spatial dimension
#define cut 2.5 //potential cut off
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

void ini_coord_square(double (*x)[dim]){
  int num_x = (int)sqrt(Np)+1;
  int num_y = (int)sqrt(Np)+1;
  int i,j,k=0;
  for(j=0;j<num_y;j++){
    for(i=0;i<num_x;i++){
      x[i+num_x*j][0] = i*L/(double)num_x;
      x[i+num_x*j][1] = j*L/(double)num_y;
      k++;
      if(k>=Np)
        break;
    }
    if(k>=Np)
      break;
  }
}

void set_diameter(int *a){
  for(int i=0;i<0.35*Np;i++){
    a[2*i]=2;
    a[2*i+1]=1;
  }
  for(int i=0.7*Np;i<Np;i++)
    a[i]=1;
}

void p_boundary(double (*x)[dim]){
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++)
      x[i][j]-=L*floor(x[i][j]/L);
}

void ini_array(double (*x)[dim]){
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++)
      x[i][j]=0.0;
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
    dx-=L*floor((dx+0.5*L)/L);
    dy-=L*floor((dy+0.5*L)/L);

    *x_corr+=dx/double(Np); //center of mass displacement.x
    *y_corr+=dy/double(Np);

    for(j=0;j<dim;j++){
      X0[i][j]=x[i][j];
    }
  }

  return 0;
}

void list_verlet(int (*list)[Nn],double (*x)[dim]){
  double dx,dy,dr2;
  double thresh=cut+skin;
  for(int i=0;i<Np;i++)
    for(int j=0;j<Nn;j++)
      list[i][j]=0;

  for(int i=0;i<Np;i++)
    for(int j=0;j<Np;j++){
      if(j>i){
	dx=x[i][0]-x[j][0];
	dy=x[i][1]-x[j][1];
	dx-=L*floor((dx+0.5*L)/L);
	dy-=L*floor((dy+0.5*L)/L);
	dr2=dx*dx+dy*dy;
	if(dr2<thresh*thresh){
	  list[i][0]++;
	  list[i][(int)list[i][0]]=j;
	}
      }
    }
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


void calc_force(double (*x)[dim],double (*f)[dim],int *a,double *U,int (*list)[Nn]){
  double dx,dy,dr2,dUr,w2,w6,w12,w2cut,w6cut,w12cut,aij,eij,dUrcut,Ucut,dr;
  ini_array(f);
  *U=0;
  for(int i=0;i<Np;i++)
    for(int j=1;j<=list[i][0];j++){
      dx=x[i][0]-x[list[i][j]][0];
      dy=x[i][1]-x[list[i][j]][1];
      dx-=L*floor((dx+0.5*L)/L);
      dy-=L*floor((dy+0.5*L)/L);
      dr2=dx*dx+dy*dy;
      if(a[i]+a[list[i][j]] == 2){
	aij=1.0;
	eij=1.0;
	//      1:0 0:8 0:88
	//  1:0 0:5 1:5
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
	f[i][0]-=dUr*dx;
	f[list[i][j]][0]+=dUr*dx;
	f[i][1]-=dUr*dy;
	f[list[i][j]][1]+=dUr*dy;
	*U+=4.*eij*w12-4.*eij*w6-Ucut-dUrcut*(dr-cut*aij);
      }
    }
}

void eom_langevin(double (*v)[dim],double (*x)[dim],double (*f)[dim],int *a,double *U,double dt,double temp0,int (*list)[Nn],double *kine){
  double zeta=1.0;
  double fluc=sqrt(2.*zeta*temp0*dt);
  *kine=0.0;
  calc_force(x,f,a,&(*U),list);
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++){
      v[i][j]+=-zeta*v[i][j]*dt+f[i][j]*dt+fluc*gaussian_rand();
      x[i][j]+=v[i][j]*dt;
      *kine +=v[i][j]*v[i][j]/2.0/Np;
    }
  p_boundary(x);
}

void eom_md(double (*v)[dim],double (*x)[dim],double (*f)[dim],int *a,double *U,double dt,int (*list)[Nn]){
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++){
      x[i][j]+=v[i][j]*dt+0.5*f[i][j]*dt*dt;
      v[i][j]+=0.5*f[i][j]*dt;
    }
  calc_force(x,f,a,&(*U),list);
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++){
      v[i][j]+=0.5*f[i][j]*dt;
    }
  p_boundary(x);
}

void output(int k,double (*v)[dim],double U,double T){
  char filename[128];
  double K=0.0;

  ofstream file;
  sprintf(filename,"energy.dat");
  file.open(filename,ios::app); //append
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++)
      K+=0.5*v[i][j]*v[i][j];

  //  std::cout<< std::setprecision(6)<<k*dtbd<<"\t"<<T<<"\t"<<K/Np<<"\t"<<U/Np<<"\t"<<(K+U)/Np<<std::endl;
  file<< setprecision(6) << k*dtbd << " " << T << " " << K/Np << " " << U/Np << " " << (K+U)/Np << endl;
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

void auto_list_update(double *disp_max,double (*x)[dim],double (*x_update)[dim],int (*list)[Nn],int M){
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

void copy_array(double (*x)[dim],double (*x0)[dim]){
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++)
      x0[i][j]=x[i][j];
}

int main(){
  double x[Np][dim],x_update[Np][dim],v[Np][dim],f[Np][dim],kine,t=0.0,fs=0.0,xf[Np][dim];
  int list[Np][Nn],a[Np];
  double tout=0.0,U,disp_max=0.0,temp_anneal,sampling_time,time_stamp=0.0,sampling_time_max=1.e+4,t0=0.0,dT=8.33*1.e-5;
  double tanneal=(4.0-temp)/dT;
  int i,j=0,count=0,tcount=0,ens_count=0;
  int M=(int)(L/(cut+skin));
  double x_corr=0.0,y_corr=0.0;
  set_diameter(a);
  ini_coord_square(x);
  ini_array(v);
  cell_list(list,x,M);
  std::cout<<"L="<<L<<"M="<<M<<std::endl;

  j=0;
  while(t=0;t<tanneal;t+=dtbd){
    j++;
    tcount++;
    temp_anneal=-dT*t+4.0;
    auto_list_update(&disp_max,x,x_update,list,M);
    eom_langevin(v,x,f,a,&U,dtbd,temp_anneal,list,&kine);
    com_correction(x,&x_corr,&y_corr);
    //  std::cout<<f[0][0]<<" "<<kine<<std::endl;
  }

  j=0;
  tcount=0;

  while(j*dtbd < teq){
    j++;
    tcount++;
    auto_list_update(&disp_max,x,x_update,list,M);
    eom_langevin(v,x,f,a,&U,dtbd,temp,list,&kine);
    com_correction(x,&x_corr,&y_corr);
    //  std::cout<<f[0][0]<<" "<<kine<<std::endl;
  }

  j=0;
  t=0.0;

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
