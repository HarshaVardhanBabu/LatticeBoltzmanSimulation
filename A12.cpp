	#include <iostream>
	#include <string.h>
	#include <fstream> 

	using namespace std;

	//A.1.2 The LBM Code (2DQ4)
	const int n = 100;
	const int m = 100;

	int main()
	{

		double f1[n+1][m+1]={0},f2[n+1][m+1]={0},f3[n+1][m+1]={0},f4[n+1][m+1]={0};
		double rho[n+1][m+1]= {0},feq=0,x[n+1]={0},y[m+1]={0};
		int i=0;
		ofstream results("/home/harshavardhanbabu/A12results.dat");
		ofstream midplbm("/home/harshavardhanbabu/A12midplbm.dat");
		ofstream r("/home/harshavardhanbabu/A12_100.dat");
		double dx = 1.0;
		double dy = dx;
		double dt = 1.0;
		x[0] = 0.0;
		y[0] = 0.0;
		for(i=1;i<= n; i++)
		{
			x[i] = x[i-1] + dx;
		}
		for(int j=1;j<= m; j++)
		{
			y[j] = y[j-1] + dy;
		}
		
			// for(int j=0;j<= m; j++)
			// {
			// 	cout<<j<<"\t"<<x[j]<<"\t"<<y[j]<<endl;
			// }
		
		double csq = (dx*dx)/(dt*dt);
		double alpha = 0.25;
		double omega = 1.0/(2.0*alpha/(dt*csq)+0.5);
		double mstep = 400;

		for(int j=0;j<=m;j++)
		{
			for(int i=0;i<=n;i++)
			{
				rho[i][j] = 0.0;
				
			}
		}

		for(int j=0;j<=m;j++)
		{
			for(int i=0;i<=n;i++)
			{
				f1[i][j] = 0.25*rho[i][j];
				f2[i][j] = 0.25*rho[i][j];
				f3[i][j] = 0.25*rho[i][j];
				f4[i][j] = 0.25*rho[i][j];
			}
		}
		for(int k =1;k<=mstep;k++)
		{
			for(int j=0;j<=m;j++)
				{
					for(int i=0;i<=n;i++)
						{
							feq = 0.25*rho[i][j];
							f1[i][j] = omega*feq+(1.0-omega)*f1[i][j];
							f2[i][j] = omega*feq+(1.0-omega)*f2[i][j];
							f3[i][j] = omega*feq+(1.0-omega)*f3[i][j];
							f4[i][j] = omega*feq+(1.0-omega)*f4[i][j];
						}
				}
		// Streaming
		for(int j=0;j<=m;j++)
		{
			for(int i=0;i<=n;i++)
			{
				f1[n-i][j] = f1[n-i-1][j];
				f2[i-1][j] = f2[i][j];
			}
		}
		for(int i=0;i<=n;i++)
		{
			for(int j=0;j<=m;j++)
			{
				f3[i][m-j] = f1[i][m-j-1];
				f4[i][j-1] = f4[i][j];
			}
		}
		for(int j=1;j<=m;j++)
		{
			f1[0][j] = 0.5 - f2[0][j];
			f3[0][j] = 0.5 - f4[0][j];
			f1[n][j] = 0.0;
			f2[n][j] = 0.0;
			f3[n][j] = 0.0;
			f4[n][j] = 0.0;

		}
	for(int i=1;i<=m;i++)
		{
			f1[i][m] = 0.0;
			f2[i][m] = 0.0;
			f3[i][m] = 0.0;
			f4[i][m] = 0.0;
			f1[i][0] = f1[i][1];
			f2[i][0] = f2[i][1];
			f3[i][0] = f3[i][1];
			f4[i][0] = f4[i][1];
		}
		for(int i=0;i<=n;i++)
		{
			for(int j=0;j<=m;j++)
			{
				rho[i][j] = f1[i][j]+f2[i][j]+f3[i][j]+f4[i][j];
				// if(k == 100)
				// {
					
				// 	r<<i<<"\t"<<j<<"\t"<<rho[i][j]<<endl;

					
				// }
				
			}
		}
		// write to file the rho[i][j] with time step;


	}	
	for (int i = 0; i <= n; ++i)
	{
		midplbm<<x[i]<<"\t"<<rho[i][m/2]<<"\r\n";
	}
	midplbm.close();
	results.close();
	r.close();
		return 0;
}

		


