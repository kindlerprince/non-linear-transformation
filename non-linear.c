#include<stdio.h>
#include<stdlib.h>
#include<math.h>
float *gaussian_elimination(float **arr, int n);
void print_mat(float **arr, int n,int m);
float **transpose(float **arr,int n,int m);
float **matrix_multiply(float **a, float **b,int n,int m,int p);
float **input_data(FILE *fp,int N,int d,float *y);
void error(float *w,int N,int d);

int main()
{
	int n,i,j,N=160,d=3;
	FILE *fp;
	fp=fopen("data_non-linear.csv","r");
	float **x,*y,r,*w,**AB,**x_T,*B,sum;
	y=(float *)malloc(N*sizeof(float));
	x=input_data(fp,N,d,y);
	fclose(fp);
	//print_mat(x,N,d+1);
	x_T=transpose(x,N,d+1);
	//print_mat(x_T,d+1,N);

	/*for(i=0; i<N; i++)
	{
		printf("%f\n",y[i]);
	}*/

	AB=matrix_multiply(x_T,x,d+1,N,d+1);

	B=(float *)malloc((d+1)*sizeof(float));

	for(i=0; i<d+1; i++)
	{
		sum=0;
		for(j=0; j<N; j++)
		{
			sum+=x_T[i][j]*y[j];
		}
		B[i]=sum;
	}
	/*for(i=0; i<d+1; i++)
	{
		printf("%f\n",B[i]);
	}*/
	for(i=0; i<d+1; i++)
	{
		AB[i]=(float *)realloc(AB[i],(d+2)*sizeof(float));
	}
	for(i=0; i<d+1; i++)
	{
		AB[i][d+1]=B[i];
	}
	w=gaussian_elimination(AB,d+1);
	printf("+------------------------+\n");
	for(i=0; i<d+1; i++)
	{
		printf("%f\n",w[i]);
	}
	printf("+------------------------+\n");
	error(w,40,d);
	return 0;
}
float *gaussian_elimination(float **arr, int n)
{
	int i,j,k;
	float r,subs,sum,*x;
	x=(float *)malloc(n*sizeof(float));
	for(i=0; i<n; i++)
	{
		for(j=i+1; j<n; j++)
		{
			r=arr[j][i]/arr[i][i];
			for(k=0; k<n+1; k++)
			{
				arr[j][k]=arr[j][k]-r*arr[i][k];
			}
		}
	}
	//print_mat(arr,n,n+1);
	for(i=n-1; i>=0; i--)
	{
		subs=0;
		/*  The below loop won't run for the first time */
		for(j=n-1; j>i; j--)
		{
			subs+=arr[i][j]*x[j];
		}
		x[j]=(arr[i][n]-subs)/arr[i][j];

	}
	return x;
}

float **matrix_multiply(float **a,float **b,int n,int m,int p)
{
	int i,j,k;
	float **c,sum;
	c=(float **)malloc(n*sizeof(float *));
	for(i=0; i<n; i++)
	{
		*(c+i)=(float *)malloc(p*sizeof(float));
	}
	for(i=0; i<n; i++)
	{
		for(k=0; k<p; k++)
		{
			sum=0;
			for(j=0; j<m; j++)
			{
				sum+=a[i][j]*b[j][k];
			}
			c[i][k]=sum;
		}
	}
	return c;
}

float **transpose(float **arr,int n,int m)
{
	int i,j;
	float **tarr;
	tarr=(float **)malloc(m*sizeof(float *));
	for(i=0; i<m; i++)
	{
		*(tarr+i)=(float *)malloc(n*sizeof(float));
	}
	for(i=0; i<n; i++)
	{
		for(j=0; j<m; j++)
		{
			tarr[j][i]=arr[i][j];
		}
	}
	return tarr;
}
void error(float *w,int N,int d)
{
	int i,j;
	float **x,*y,e,E,sum;
	FILE *fp;
	fp=fopen("data_non-linear_check.csv","r");
	y=(float *)malloc(N*sizeof(float));
	x=input_data(fp,N,d,y);
	for(i=0; i<N; i++)
	{
		sum=0;
		for(j=0; j<d+1; j++)
		{
			sum+=w[j]*x[i][j];
		}
		e+=(y[i]-sum)*(y[i]-sum);
	}
	E=e/N;
	printf("Square Error = %f \n",E);
	
}
void print_mat(float **arr, int n,int m)
{
	int i,j;
	printf("+--------------------------------------------------+\n");
	for(i=0; i<n; i++)
	{
		for(j=0; j<m; j++)
		{
			printf("%f ",arr[i][j]);
		}
		printf("\n");
	}
	printf("+--------------------------------------------------+\n");
}
float **input_data(FILE *fp,int N,int d,float *y)
{
	int i,j;
	float **x,z;
	x=(float **)malloc(N*sizeof(float *));

	for(i=0; i<N; i++)
	{
		*(x+i)=(float *)malloc((d+1)*sizeof(float));
	}

	for(i=0; i<N; i++)
	{
		x[i][0]=1;
		fscanf(fp,"%f,",&z);
		x[i][1]=z*z;
		x[i][2]=1/z;
		x[i][3]=exp(3*z);
		fscanf(fp,"%f\n",&y[i]);
	}
	return x;
}
