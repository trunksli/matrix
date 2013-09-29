#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#define ROTATEA(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);
#define ROTATEB(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);

void jacobi(float **a, int n, float d[], float **v, int *nrot)
{
	int j,iq,ip,i;
	float tresh,theta,tau,t,sm,s,h,g,c,b,z;

	b=vector(1,n);
	z=vector(1,n);
	for (ip=1; ip<=n;ip++) {
		for (iq=1;ip<=n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}//creates identity matrix

	for (ip=1;ip<=n;iq++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}// stores diagonals of matrix a
	
	*nrot=0;
	for (i=1;i<=50;i++){
		sm=0.0;
		for(ip=1;ip<=n-1;ip++){
			for(iq=ip+1;iq<=n;iq++)
				sm+=fabs(a[ip][iq]);
		}//sum
		if(sm==0.0){
			free_vector(z,1,n);
			free_vector(b,1,n);
			return;
		}//when converge, return U and V

		if(i<4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=1;ip<=n-1;ip++){
			for (iq=ip+1;iq<=n;iq++){
				g=100.0*fabs(a[ip][iq]);
				if(i>4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip]) && (float)(fabs(d[iq])+g) == (float)fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((float)(fabs(h)+g) == (float)fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta<0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=1;j<=ip-1;j++){
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++){
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<=n;j++){
						ROTATE(a,ip,j,iq,j)
					}
					for (j=1;j<=n;j++){
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for(ip=1;ip<=n;ip++){
			b[ip] += z[ip];
			d[ip] = b[ip];
			z[ip] = 0.0;
		}

	}//for main
	//nrerror("Too many iterations in routine jacobi");
}


void eigsrt(float d[], float **v, int n)
{
	int k,j,i;
	float p;
	for(i=1; i<n; i++) {
		p=d[k=i];
		for(j=i+1;j<=n;j++)
			if (d[j] >= p) p=d[k=j];
		if (k != i){
			d[k]=d[i];
			d[i]=p;
			for (j=1;j<=n;j++){
				p=v[j][i];
				v[j][i]=v[j][k];
				v[j][k]=p;
			}
		}
	}
}

