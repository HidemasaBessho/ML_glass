#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cfloat>


#define Np 4096 //# of the particles
#define rho 1.2 //# density
#define Nn 100  //# of the neigbour lists
# define L sqrt(Np/rho)
#define teq 4000 //equilibration time
#define tmax 1000000 //production run time
#define dtmd 0.001 //dt for molecular dynamics
#define dtbd 0.005 //dt for brownian dynamics
#define temp 0.4 // temperature
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

void set_diameter(double *a){
  for(int i=0;i<Np/2;i++){
    a[2*i]=5./6.;
    a[2*i+1]=7./6.;
  }
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


void calc_force(double (*x)[dim],double (*f)[dim],double *a,double *U,int (*list)[Nn]){
  double dx,dy,dr2,dUr,w2,w6,w12,w2cut,w6cut,w12cut,aij,dUrcut,Ucut,dr;
  ini_array(f);
  *U=0;
  for(int i=0;i<Np;i++)
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

void eom_langevin(double (*v)[dim],double (*x)[dim],double (*f)[dim],double *a,double *U,double dt,double temp0,int (*list)[Nn],double *kine){
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

void eom_md(double (*v)[dim],double (*x)[dim],double (*f)[dim],double *a,double *U,double dt,int (*list)[Nn]){
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

  std::ofstream file;
  sprintf(filename,"energy.dat");
  file.open(filename,std::ios::app); //append
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++)
      K+=0.5*v[i][j]*v[i][j];

  //  std::cout<< std::setprecision(6)<<k*dtbd<<"\t"<<K/Np<<"\t"<<U/Np<<"\t"<<(K+U)/Np<<std::endl;
  file<< std::setprecision(6)<<k*dtbd<<"\t"<<T<<"\t"<<K/Np<<"\t"<<U/Np<<"\t"<<(K+U)/Np<<std::endl;
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
    //    std::cout<<"update"<<*disp_max<<" "<<count<<std::endl;
    *disp_max=0.0;
    count=0;
  }
}

void Fs(double (*x)[dim],double (*xf)[dim],double *fs,double *x_corr,double *y_corr){
  double q=2.0*M_PI/1.0;
  int i,j;
  double dx,dy;
  *fs=0.0;
  for(i=0;i<Np;i++){
    dx = (x[i][0]-(*x_corr))-xf[i][0];
    dy = (x[i][1]-(*y_corr))-xf[i][1];
    dx-=L*floor((dx+0.5*L)/L);
    dy-=L*floor((dy+0.5*L)/L);
    *fs += (cos(-q*dx)+cos(-q*dy))/Np/2.0;
  }
}

void output_Fs(double t,double fs,double time_stamp,double t0){
  char filename[128];
  std::ofstream file;
  sprintf(filename,"fs.dat");
  file.open(filename,std::ios::app);
  file<< std::setprecision(6)<<t0-(t-time_stamp)<<"\t"<<fs<< std::endl;
  file.close();
}

void copy_array(double (*x)[dim],double (*x0)[dim]){
  for(int i=0;i<Np;i++)
    for(int j=0;j<dim;j++)
      x0[i][j]=x[i][j];
}

int main(){
  double x[Np][dim],x_update[Np][dim],v[Np][dim],f[Np][dim],a[Np],kine,xf[Np][dim];
  int list[Np][Nn];
  double tout=0.0,U,disp_max=0.0,temp_anneal,t=0.0,sampling_time,time_stamp=0.0,sampling_time_max=1.e+4,fs=0.0,t0=0.0;
  int i,j=0,count=0,tcount=0;
  int M=(int)(L/(cut*7./6.+skin));
  double x_corr=0.0,y_corr=0.0;
  set_diameter(a);
  ini_coord_square(x);
  ini_array(v);
  cell_list(list,x,M);
  std::cout<<"L="<<L<<"M="<<M<<std::endl;

  j=0;
  while(j*dtbd < teq){
    j++;
    tcount++;
    temp_anneal=4.0-j*dtbd*(4.0-temp)/teq;
    auto_list_update(&disp_max,x,x_update,list,M);
    eom_langevin(v,x,f,a,&U,dtbd,temp_anneal,list,&kine);
    com_correction(x,&x_corr,&y_corr);
    if(tcount==int(teq/1000.0/dtbd)){
      output(j,v,U,temp_anneal);
      tcount=0;
    }
    //  std::cout<<f[0][0]<<" "<<kine<<std::endl;
  }

  j=0;
  tcount++;

  while(j*dtbd < teq){
    j++;
    tcount++;
    auto_list_update(&disp_max,x,x_update,list,M);
    eom_langevin(v,x,f,a,&U,dtbd,temp,list,&kine);
    com_correction(x,&x_corr,&y_corr);
    if(tcount==int(teq/1000.0/dtbd)){
      output(j,v,U,temp);
      tcount=0;
    }
  }

  j=0;
  sampling_time=10.0*dtbd;
  for(t=0.0;t<tmax;t+=dtbd){
    j++;
    auto_list_update(&disp_max,x,x_update,list,M);
    eom_langevin(v,x,f,a,&U,dtbd,temp,list,&kine);
    com_correction(x,&x_corr,&y_corr);
    if(int(t/dtbd) == int((sampling_time + time_stamp)/dtbd)){
      if(count==0){
	t0=t-time_stamp;
	for(i=0;i<Np;i++){
	  xf[i][0] = x[i][0]-x_corr;
	  xf[i][1] = x[i][1]-y_corr;
	}
	sampling_time*=pow(10.,0.1);
	sampling_time=int(sampling_time/dtbd)*dtbd;
	count++;
      }
      else{
	Fs(x,xf,&fs,&x_corr,&y_corr);
	output_Fs(t,fs,time_stamp,t0);
	sampling_time*=pow(10.,0.1);
	sampling_time=int(sampling_time/dtbd)*dtbd;
      }
      if(sampling_time > sampling_time_max/pow(10.,0.1)){
        time_stamp=t;
        sampling_time=10.0*dtbd;
	count=0;
      }
    }
  }
  return 0;
}
