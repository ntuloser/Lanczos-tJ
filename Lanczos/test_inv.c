#include<stdio.h> 
#include<stdlib.h>
#include <math.h>
#include "nrutil.h"
#define TINY 1.0e-20

void ludcmp(float **a, int n, int *indx, float *d){
	int i,imax,j,k;
	float big,dum,sum,temp;
	float *vv;

	vv=vector(1,n);
	*d=1.0;

	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++){
			if ((temp=fabs(a[i][j])) > big) big=temp;
		}
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}

	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;


		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}

		if (j != imax) {
			for (k=1;k<=n;k++)
			{dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}

		indx[j]=imax;

		if (a[j][j] == 0.0){
			a[j][j]=TINY;

		}

		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum; }
	}

	free_vector(vv,1,n); 
}

double determinant(float** a,int N){
	float d;
	int j,*indx;
	indx=ivector(1, N);
	ludcmp(a,N,indx,&d);
	for(j=1;j<=N;j++) {d *= a[j][j];}
	return d;

}


int main(){
	printf("hello\n");
	//double** a=malloc(4* sizeof(double*));
	//for(int i=0; i<4;i++){
	//	a[i]=malloc(4*sizeof(double));
	//}
	float **a; a=matrix(1,3,1,3);
	//for(int i=1; i<=4;i++){
	//	for(int j=1;j<=4;j++){
	//		a[i][j]=(i+j)*j;
	//	}
	//}
	a[1][1]=3;
	a[2][2]=2;
	a[3][3]=1;
	a[1][3]=1;
	a[3][1]=1;
	printf("%f\n",determinant(a,3));
	printf("could not possilby reach here\n");	
	return 0;
}
