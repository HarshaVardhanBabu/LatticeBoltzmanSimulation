#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>

using namespace std;

const double pi = 3.141592653589793;
const int n = 101;
const int m = 101;
const int d = 9;

class D2Q9
{
	void init_Basic();
	vector<int> cy,cx;
	vector<double> w;

	double f[d][n][m];
	double feq;
	double rho[n][m],x[n],y[m];

	double u[n][m];
	double v[n][m];
	double strf[n][m];

	double tw ;
	double csq ;
	double rhoo;
	double dx,dy;
	double dt ;
	double alpha;
	double Re;
	double omega;
	double mstep;
public:
	D2Q9();
	~D2Q9();

	void collision();
	void streaming();
	void boundaryConditions();
	void result();
	void startSimulation();

	/* data */
};
void D2Q9 :: result()
{
	double sum =0.0;
	for (int j = 0; j < m; ++j)
	{
		for (int i = 0; i < n; ++i)
		{
			sum = 0.0;
			for(int k = 0; k < d; ++k)
			{
				sum += f[k][i][j];
			}
			rho[i][j] = sum;
		}
	}
	ofstream results("/home/harshavardhanbabu/Qresu.dat");
	results<<"Variables = X\tY\tT\r\n";
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			results<<x[i]<<"\t"<<y[j]<<"\t"<<rho[i][j]<<"\r\n";
		}
	}
	results.close();
	ofstream results1("/home/harshavardhanbabu/midtlbmtc.dat");
	for (int i = 0; i < n; ++i)
	{
		results1<<x[i]<<"\t"<<rho[i][m/2]<<"\r\n";
	}
	results1.close();


}
void D2Q9 :: boundaryConditions()
{
	for (int j = 0; j < m; ++j)
	{
		f[1][0][j] = w[1]*tw + w[3]*tw - f[3][0][j];
		f[5][0][j] = w[5]*tw + w[7]*tw - f[7][0][j];
		f[8][0][j] = w[8]*tw + w[6]*tw - f[6][0][j];
		f[3][n-1][j] = -f[1][n-1][j];
		f[6][n-1][j] = -f[8][n-1][j];
		f[7][n-1][j] = -f[5][n-1][j];
	}
	for (int i = 0; i < n; ++i)
	{
		f[4][i][m-1] = -f[2][i][m-1];
		f[7][i][m-1] = -f[5][i][m-1];
		f[8][i][n-1] = -f[6][i][m-1];
		for (int k = 1; k < d; ++k)
		{
			f[k][i][0] = f[k][i][1];
		}
	}
}

void D2Q9 :: streaming()
{
	for (int j = m-1; j >= 0; --j)
	{
		for (int i = 0; i < n; ++i)
		{
			f[2][i][j] = f[2][i][j-1];
			f[6][i][j] = f[6][i+1][j-1];
		}

	}
	for (int j = m-1; j >= 0; --j)
	{
		for (int i = n-1; i >= 0; --i)
		{
		    for(int k = 0; k < d; ++k)
			{
			f[1][i][j] = f[1][i-1][j];
			f[5][i][j] = f[5][i-1][j-1];
		}

	}
	for (int j = 0; j < m; ++j)
	{
		for (int i = n-1; i >= 0 ; --i)
		{
			f[4][i][j] = f[4][i][j+1];
			f[8][i][j] = f[8][i-1][j+1];
		}
	}
	for (int j = 0; j < m; ++j)
	{
		for (int i = 0; i < n; ++i)
		{
			f[3][i][j] = f[3][i+1][j];
			f[7][i][j] = f[7][i+1][j+1];
		}
	}
}

}
void D2Q9 :: collision()
{
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			for(int k = 0; k < d; ++k)
			{
				feq = w[k]*rho[i][j];
				f[k][i][j] = omega*feq + (1.0-omega)*f[k][i][j];
			}
		}
	}

}
void D2Q9 :: init_Basic()
{
	cx.resize(9);
	cy.resize(9);
	w.resize(9);
	cx[0] = 0;
	cy[0] = 0;
	w[0]  = 4.0/9;

	for(int i = 1; i <= 4; i++)
	{
		cx[i] = cos((pi/2)*(i-1));
		cy[i] = sin((pi/2)*(i-1));
		w[i]  = 1.0/9;
	}

	for (int i = 5; i <= 8; ++i)
	{
		cx[i] = sqrt(2)*cos((pi/2)*(i-4-0.5));
		cy[i] = sqrt(2)*sin((pi/2)*(i-4-0.5));
		w[i]  = 1.0/36;
	}
	x[0] = 0;
	y[0] = 0;
	for (int i = 1; i < n; ++i)
	{
		x[i] = x[i-1] + dx;
	}
	for (int i = 1; i < m; ++i)
	{
		y[i] = y[i-1] + dy;
	}
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			rho[i][j] = 0.0;
			for(int k = 0; k < d; ++k)
			{
				f[k][i][j] = w[k]*rho[i][j];
				if(i == 0)
				{
					f[k][i][j] = w[k]*tw;
				}
			}
		}

	}
	// Code to print the directions in the lattice..
	// cout<<"The directions assigned are as follows:"<<endl;
	// for(int i=0;i<=8;++i)
	// {
	// 	cout<<"e"<<i<<" = "<<"("<<cx[i]<<","<<cy[i]<<");"<<endl;
	// }
}
D2Q9 :: D2Q9()
{
		cout<<"Start of main"<<endl;
		tw = 1.0;
		dx = 1.0;
		dy = dx;
		dt =1.0;
		alpha = 0.25;
		csq = (dx*dx)/(dt*dt);
		//Re = u0*m/alpha;
		//cout<<"Reynolds Number is :"<<Re<<endl;
		omega = 1.0/(3.0*alpha/(csq*dt)+0.5);
		mstep = 400;
		init_Basic();
}
void D2Q9 :: startSimulation()
{
	double sum;
	for (int i = 1; i <= mstep; ++i)
	{
		for (int i = 0; i < m; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				sum = 0.0;
				for(int k = 0; k < d; ++k)
				{
					sum += f[k][i][j];
				}
				rho[i][j] = sum;
			}
		}
		collision();
		streaming();
		boundaryConditions();
	}

}
D2Q9 :: ~D2Q9()
{
	cout<<"End"<<endl;
}
int main(int argc, char const *argv[])
{
	D2Q9 d;
	d.startSimulation();
	d.result();
	return 0;
}
