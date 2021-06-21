#include <stdio.h>
#include <math.h>
#include "complex.h"
#include "gsl_rng.h"

#define N 500
#define nciclos 125
#define lambda 0.4
#define pi 3.14159265359
#define TMAX 10000
#define nD 2000

gsl_rng *tau;

double V[N];
fcomplex f[N];
fcomplex chi[N];
fcomplex A[N];
fcomplex a[N-1];
fcomplex b[N-1];

int CalcV(double k);
int finicial(double k);
int calca();
int calcb(double s);
int calcchi();
double calcf(int n);
double PD(double Norm);
double PI(double Norm);

int main()
{
	int j, n, m, mt, l;
	int semilla=1465124;
	extern gsl_rng *tau;
	double k, s, x, Norm, coef;
	
	m = mt = 0;
	
	tau=gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(tau,semilla);
	
	for(l=0;l<1000;l++)
	{
		chi[0] = chi[N-1] = Complex(0.0,0.0);
		
		k = (2*pi*nciclos)/(N);
		s = (1)/(4*pow(k,2));
		
		CalcV(k);
		
		for(j=0;j<N;j++)
		{
			A[j] = Complex(-2-V[j],2/s);
		}
			
		finicial(k);
		
		calca();
		
		for(n=1;n<TMAX;n++) //Este es el bucle que hace correr el tiempo, debe pararse cada nD pasos
		{
			calcb(s);
			
			calcchi();
			
			Norm = calcf(n); //Se calcula la función y a la vez la norma
			
			if(nD%n==0)
			{					
				x = gsl_rng_uniform(tau); //Se calcula un número aleatorio
				
				if(x<PD(Norm))
				{
					mt = mt + 1;
					m = m + 1;
					
					n = TMAX; //De esta forma el bucle no se volverá a ejecutar y acaba el programa
				}
				else
				{
					for(j=(4*N/5)-1;j<N;j++)
					{
						f[j] = Complex(0.0,0.0);
					}
					
					for(j=0;j<N;j++) //Se calcula la nueva norma
					{
						Norm = Norm + pow(Cabs(f[j]),2);
					}
		
					Norm = sqrt(Norm);
					
					for(j=0;j<N;j++) //Se normaliza la función
					{
						f[j].r = f[j].r/Norm;
						f[j].i = f[j].i/Norm;
					}
					
					for(j=0;j<N;j++) //Se calcula la norma tras normalizar (debería ser 1)
					{
						Norm = Norm + pow(Cabs(f[j]),2);
					}
		
					Norm = sqrt(Norm);
					
					x = gsl_rng_uniform(tau);
					
					if(x<PI(Norm))
					{
						m = m + 1;
						
						n = TMAX; //De esta forma el bucle no se volverá a ejecutar y acaba el programa
					}
					else
					{
						for(j=0;j<N/5;j++)
						{
							f[j] = Complex(0.0,0.0);
						}
						
						for(j=0;j<N;j++) //Se calcula la nueva norma
						{
							Norm = Norm + pow(Cabs(f[j]),2);
						}
		
						Norm = sqrt(Norm);
					
						for(j=0;j<N;j++) //Se normaliza la función
						{
							f[j].r = f[j].r/Norm;
							f[j].i = f[j].i/Norm;
						}
					
						for(j=0;j<N;j++) //Se calcula la norma tras normalizar (debería ser 1)
						{
							Norm = Norm + pow(Cabs(f[j]),2);
						}
		
						Norm = sqrt(Norm);
					}
				}
			}
		}
	}
	
	coef = 1.0*mt/m;
	
	printf("Se ha detectado la partícula %i veces. %i veces ha sido a la derecha.\n El coeficiente de transmisión es de %lf",m,mt,coef);
	
	return 0;
}

int CalcV(double k)
{
	int j;
	
	for(j=0;j<N;j++)
	{
		if(j>=((2*N/5)-1) && j<=((3*N/5)-1))
		{
			V[j]=lambda*pow(k,2);
		}
		else
		{
			V[j]=0;
		}
	}
	
	return 0;
}

int finicial(double k)
{	
	int j;
	
	f[0] = f[N-1] = Complex(0.0,0.0);
	
	for(j=1;j<N-1;j++)
	{	
		f[j] = Cgauss(k*j,exp(-8*pow(4*j-N,2)/pow(N,2)));
	}
	
	return 0;
}

int calca()
{
	int j;
	
	a[N-2] = Complex(0.0,0.0);
	
	for(j=N-2;j>0;j--)
	{
		a[j-1] = Cdiv(Complex(-1.0,0.0),Cadd(A[j],a[j]));
	}
	
	return 0;
}

int calcb(double s)
{
	int j;
	
	b[N-2] = Complex(0.0,0.0);
	
	for(j=N-2;j>0;j--)
	{
		b[j-1] = Cdiv(Csub(Cmul(Complex(0.0,4.0/s),f[j]),b[j]),Cadd(A[j],a[j]));
	}
	
	return 0;
}

int calcchi()
{
	int j;
	
	for(j=1;j<N-1;j++)
	{
		chi[j] = Cadd(Cmul(a[j-1],chi[j-1]),b[j-1]);
	}
	
	return 0;
}

double calcf(int n)
{
	int j;
	double Norm;
	
	for(j=1;j<N-1;j++)
	{
		f[j] = Csub(chi[j],f[j]);
	}
	
	Norm = 0.0;
	
	for(j=0;j<N;j++)
	{
		Norm = Norm + pow(Cabs(f[j]),2);
	}
	
	Norm = sqrt(Norm);

	return Norm;
}

double PD(double Norm)
{
	int j;
	double PD;
	
	PD = 0.0;
	
	for(j=(4*N/5)-1;j<N;j++)
	{
		PD = PD + pow(Cabs(f[j]),2);
	}
	
	PD = PD/Norm; //Se divide la probabilidad por la norma para normalizarla
	
	return PD;
}

double PI(double Norm)
{
	int j;
	double PI;
	
	PI = 0.0;
	
	for(j=0;j<N/5;j++)
	{
		PI = PI + pow(Cabs(f[j]),2);
	}
	
	PI = PI/Norm; //Se divide la probabilidad por la norma para normalizarla
	
	return PI;
}
