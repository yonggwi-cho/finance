#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mm_malloc.h>

#define FileName "BS_FDM.txt"

/* Option payoff (for call options) */
double CallPayoff(double x, double k){
  double y;
  y=exp(0.5*(k+1)*x)-exp(0.5*(k-1)*x);
  if (y<0) y=0.0;
  return y;
}

/* Boundary condition at S=0 */
double CallBC1(double x, double tau){
  return 0.0;
}

/* Boundary condition for S -> infinity */
double CallBC2(double x, double tau, double alpha, double beta, double k){
  double z;
  z=exp((1.0-alpha)*x-beta*tau)-exp(-(k+beta)*tau);
  return z;
}

/* Solution of the Black-Scholes equation by explicit FDM */
void ExplicitFDM(double *V, double xmax, int M, int N, double K, double r, double sigma, double T){
  double *u0, *u;
  double k, alpha, beta, a, dt, dx, tau;
  int i, j;
  
  u0=(double *)malloc(sizeof(double)*(2*N+1));
  u =(double *)malloc(sizeof(double)*(2*N+1));

  k=2.0*r/(sigma*sigma);
  alpha=-1.0/2.0*(k-1);
  beta=-1.0/4.0*(k+1)*(k+1); // ?
  
  dt=T*(1.0/2.0)*sigma*sigma/M;
  dx=xmax/N;
  a=dt/(dx*dx);
  
  /* Set the terminal condition. */
  for (j=-N; j<=N; j++)
    u0[j+N]=CallPayoff(j*dx,k);
  
  for (i=1; i<=M; i++)
    {
      tau=i*dt;
      
      /* Set the boundary conditions. */
      u[0]  =CallBC1(-N*dx,tau);
      u[2*N]=CallBC2( N*dx,tau,alpha,beta,k);
      
      for (j=-N+1; j<N; j++)
	u[j+N]=u0[j+N]+a*(u0[j+N-1]-2*u0[j+N]+u0[j+N+1]); 
      
      for (j=-N; j<=N; j++)
	u0[j+N]=u[j+N]; // ?
    }
  
  for (j=-N; j<=N; j++)
    V[j+N]=K*exp(alpha*j*dx+beta*0.5*sigma*sigma*T)*u0[j+N]; // U -> V
  free(u0);
  free(u);
  
  return;
}

void TridiagonalSolve(int np, double *alpha, double *beta, double *gamma, double *dX, double *dY, int iopt)
{
  // Solve linear simultaneous equations with a tridiagonal coefficient matrix.
  //    np: dimension of the matrix.
  //    alpha[0] - alpha[np-1]: diagonal elements of the matrix.
  //    beta[1] - beta[np-1]:   subdiagonal elements of the matrix.
  //    gamma[0] - gamma[np-2]: superdiagonal elements of the matrix.
  //    dX[np]: solution vector
  //    dY[np]: right-hand-side vector 
  //    iopt: LU decomposition and forward/backward substitution when iopt=1.
  //          LU decomposition only when iopt=2.
  //          forward/backward substitution only when iopt=3.
  
  int i;
  double piv;
  
  if (iopt==1 || iopt==2){
    for (i=0; i<np-1; i++){
      piv = beta[i+1] / alpha[i];
      alpha[i+1] -= gamma[i]*piv;
      beta[i+1] = piv;
    }
  }
  
  if (iopt==1 || iopt==3){
    for (i=0; i<np-1; i++)
      dY[i+1] -= dY[i]*beta[i+1];
    dX[np-1] = dY[np-1] / alpha[np-1];
    for (i=np-2; i>=0; i--)
      dX[i] = (dY[i] - gamma[i]*dX[i+1]) / alpha[i];
  }
  
  return;
}

/* Solution of the Black-Scholes equation by implicit FDM */
void ImplicitFDM(double *V, double xmax, int M, int N, double K, double r, double sigma, double T)
{
  double *u0, *u;
  double k, alpha, beta, a, dt, dx, tau;
  double *e, *f, *g, *x, *y;
  int i, j, iopt;
  
  u0=(double *)malloc(sizeof(double)*(2*N+1));
  u =(double *)malloc(sizeof(double)*(2*N+1));
  e =(double *)malloc(sizeof(double)*(2*N-1)); // diagonal elements of the coefficient matrix.
  f =(double *)malloc(sizeof(double)*(2*N-1)); // subdiagonal elements of the coefficient matrix.
  g =(double *)malloc(sizeof(double)*(2*N-1)); // superdiagonal elements of the coefficient matrix.
  x =(double *)malloc(sizeof(double)*(2*N-1)); // solution vector.
  y =(double *)malloc(sizeof(double)*(2*N-1)); // right-hand-side vector.
  
  k=2.0*r/(sigma*sigma);
  alpha=-1.0/2.0*(k-1);
  beta=-1.0/4.0*(k+1)*(k+1);
  
  dt=T*(1.0/2.0)*sigma*sigma/M;
  dx=xmax/N;
  a=dt/(dx*dx);
  
  /* Set the coefficient matrix. */
  //e[0]=1.0+2.0*a;
  //g[0]=-a;
  for (j=0; j<2*N-2; j++){
    e[j]=1.0+2.0*a;
    f[j]=-a;
    g[j]=-a;
  }
  e[2*N-2]=1.0+2.0*a;
  f[2*N-2]=-a;
  
  /* LU decomposition of the coefficient matrix. */
  iopt=2;
  TridiagonalSolve(2*N-1, e, f, g, x, y, iopt);
  
  /* Set the terminal condition. */
  for (j=-N; j<=N; j++)
    u0[j+N]=CallPayoff(j*dx,k);
  
  for (i=1; i<=M; i++)
    {
      tau=i*dt;
      
      /* Set the boundary conditions. */
      u[0]  =CallBC1(-N*dx,tau);
      u[2*N]=CallBC2( N*dx,tau,alpha,beta,k);
      
      /* Set the right-hand-side vector. */
      for (j=-N+1; j<N; j++)
	y[j+N-1]=u0[j+N];
      y[0]    +=a*u[0];
      y[2*N-2]+=a*u[2*N];

      /* Solve for new u. */
      iopt=3;
      TridiagonalSolve(2*N-1, e, f, g, x, y, iopt);

      for (j=-N+1; j<N; j++){
	u[j+N]=x[j+N-1];
      }
      for (j=-N; j<=N; j++){
	u0[j+N]=u[j+N];
	
      }
    }
  
  for (j=-N; j<=N; j++)
    V[j+N]=K*exp(alpha*j*dx+beta*0.5*sigma*sigma*T)*u0[j+N];
  
  free(u0);
  free(u);
  free(e);
  free(f);
  free(g);
  free(x);
  free(y);
  
  return;
}

/* Solution of the Black-Scholes equation by the Crank-Nicolson method */
void CrankNicolson(double *V, double xmax, int M, int N, double K, double r, double sigma, double T)
{
	double *u0, *u;
	double k, alpha, beta, a, dt, dx, tau;
	double *e, *f, *g, *x, *y;
	int i, j, iopt;
	
	u0=(double *)malloc(sizeof(double)*(2*N+1));
	u =(double *)malloc(sizeof(double)*(2*N+1));
	e =(double *)malloc(sizeof(double)*(2*N-1)); // diagonal elements of the coefficient matrix.
	f =(double *)malloc(sizeof(double)*(2*N-1)); // subdiagonal elements of the coefficient matrix.
	g =(double *)malloc(sizeof(double)*(2*N-1)); // superdiagonal elements of the coefficient matrix.
	x =(double *)malloc(sizeof(double)*(2*N-1)); // solution vector.
	y =(double *)malloc(sizeof(double)*(2*N-1)); // right-hand-side vector.

	k=2.0*r/(sigma*sigma);
	alpha=-1.0/2.0*(k-1);
	beta=-1.0/4.0*(k+1)*(k+1);
	
	dt=T*(1.0/2.0)*sigma*sigma/M;
	dx=xmax/N;
	a=dt/(dx*dx);

	/* Set the coefficient matrix. */
	e[0]=1.0+a;
	g[0]=-0.5*a;
	for (j=0; j<2*N-2; j++){
		e[j]=1.0+a;
		f[j]=-0.5*a;
		g[j]=-0.5*a;
	}
	e[2*N-2]=1.0+a;
	f[2*N-2]=-0.5*a;

	/* LU decomposition of the coefficient matrix. */
	iopt=2;
	TridiagonalSolve(2*N-1, e, f, g, x, y, iopt);

	/* Set the terminal condition. */
	for (j=-N; j<=N; j++)
		u0[j+N]=CallPayoff(j*dx,k);
	
	for (i=1; i<=M; i++)
	{
		tau=i*dt;

		/* Set the boundary conditions. */
		u[0]  =CallBC1(-N*dx,tau);
		u[2*N]=CallBC2( N*dx,tau,alpha,beta,k);

		/* Set the right-hand-side vector. */
		for (j=-N+1; j<N; j++)
			y[j+N-1]=0.5*a*u0[j+N-1] + (1.0-a)*u0[j+N] + 0.5*a*u0[j+N+1];
		y[0]    +=0.5*a*u[0];
		y[2*N-2]+=0.5*a*u[2*N];
		
		/* Solve for new u. */
		iopt=3;
		TridiagonalSolve(2*N-1, e, f, g, x, y, iopt);
		for (j=-N+1; j<N; j++)
			u[j+N]=x[j+N-1];

		for (j=-N; j<=N; j++)
			u0[j+N]=u[j+N];
	}

	for (j=-N; j<=N; j++)
		V[j+N]=K*exp(alpha*j*dx+beta*0.5*sigma*sigma*T)*u0[j+N];

	free(u0);
	free(u);
	free(e);
	free(f);
	free(g);
	free(x);
	free(y);
	
	return;
}

// main 
int main(){
  double K, r, sigma, T, xmax;
  int N, M, l;
  double *V;
  FILE *OutputFile;
  
  OutputFile=fopen(FileName,"w");
  if (!OutputFile) printf("Output file cannot be opened.\n");
  
  K=100;
  r=0.1;
  sigma=0.3;
  T=1.0;

  /* T, xmax, N and M must satisfy dt <= (1/2)dx*dx, where */
  /* dt = (1/2)*T*sigma*sigma/M and dx = xmax/N.           */
  xmax=5.0;
  N=200;
  M=150;
  //printf("dt/(dx*dx) = %lf\n", 0.5*T*sigma*sigma/M/((xmax/N)*(xmax/N)));
  fprintf(OutputFile, "dt/(dx*dx) = %lf\n", 0.5*T*sigma*sigma/M/((xmax/N)*(xmax/N)));

  V=(double*)malloc(sizeof(double)*(2*N+1));
  //ExplicitFDM(V,xmax,M,N,K,r,sigma,T);
  ImplicitFDM(V,xmax,M,N,K,r,sigma,T);
  //	CrankNicolson(V,xmax,M,N,K,r,sigma,T);
  
  for(l=0; l<=2*N; l++) printf("%f  %lf\n", K*exp(xmax/N*(l-N)),V[l]);
  for(l=0; l<=2*N; l++) fprintf(OutputFile, "%f  %lf\n", K*exp(xmax/N*(l-N)),V[l]);

  free(V);
  fclose(OutputFile);
  
  return 0;	
}
