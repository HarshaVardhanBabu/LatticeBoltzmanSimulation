#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>

using namespace std;

const double pi = 3.141592653589793;
const int n = 101;
const int m = 101;
const int d = 9;

//Computer Code for lid-driven cavity...

class LidDrivenCavity
{
	vector<int> cy,cx;
	vector<double> w;

	 double f[d][n][m];//={0};
	 double feq[d][n][m];//={0};
	 double rho[n][m];//={0};
	 double u[n][m];//={0};
	 double v[n][m];//={0};
	 double strf[n][m];//={0};

	double u0 ;
	double sumvelo ;
	double rhoo;
	double dx,dy;
	double dt ;
	double alpha;
	double Re;
	double omega;
	double mstep;
public:
	LidDrivenCavity();
	~LidDrivenCavity();
	void init_Basic();
	void collision();
	void streaming();
	void sfbound();
	void rhouv();
	void result();

	/* data */
};
void LidDrivenCavity :: result()
{
	strf[0][0]  =  0;
	double rhoav = 0.0,rhom = 0.0;
	for (int i = 0; i < n; ++i)
	{
		rhoav = 0.5*(rho[i-1][0] + rho[i][0]);
		if(i != 0)
		{
			strf[i][0] -= rhoav*0.5*(v[i-1][0]+v[i][0]);
		}
		for(int j=1;j<m;j++)
		{
			rhom = 0.5*(rho[i][j] + rho[i][j-1]);
			strf[i][j] += rhom*0.5*(u[i][j-1]+u[i][j]);
		}
	}
	ofstream datfile("/home/harshavardhanbabu/stream.dat");
	if (datfile.is_open())
		{
			datfile.close();
		}
  	else
  		{
  			cout<<"File Handling failed"<<endl;
  		}
    		

}
void LidDrivenCavity :: rhouv()
{
	double ssum;
	for (int j = 0; j < m; ++j)
	{
		for (int i = 0; i < n; ++i)
		{
			ssum = 0.0;
			for (int k = 0; k < d; ++k)
			{
				ssum += f[k][i][j];
			}
			rho[i][j] = ssum;
		}
	}
	for (int i = 0; i < n; ++i)
	{
		rho[i][m-1] = f[0][i][m-1] + f[1][i][m-1] + f[3][i][m-1] + 2.0*(f[2][i][m-1] + f[6][i][m-1] + f[5][i][m-1]);
	}
	double usum,vsum;
	for (int i = 0; i < n; ++i)
	{
		for (int j = 1; j < m-1; ++j)
		{
			usum = 0;
			vsum = 0;
			for(int k = 0; k < d; ++k)
			{
				usum += f[k][i][j]*cx[k];
				vsum += f[k][i][j]*cy[k];
			}
			u[i][j] = usum/rho[i][j];
			v[i][j] = vsum/rho[i][j];
			//cout<<i<<"\t"<<j<<"\t"<<u[i][j]<<v[i][j]<<rho[i][j]<<endl;
		}
	}
}
void LidDrivenCavity :: sfbound()
{
	for (int j = 0; j < m; ++j)
		{
			// Bounce Back on West Boundary ..
			f[1][0][j] = f[3][0][j];
			f[5][0][j] = f[5][0][j];
			f[8][0][j] = f[8][0][j];
			// Bounce Back on East Boundary..
			f[3][n-1][j] = f[1][n-1][j];
			f[7][n-1][j] = f[5][n-1][j];
			f[6][n-1][j] = f[8][n-1][j];
		}
		//  Bounce Back on South Boundary ..
		for (int i = 0; i < n; ++i)
		{
			f[2][i][0] = f[4][i][0];
			f[5][i][0] = f[7][i][0];
			f[6][i][0] = f[8][i][0];
		}
		// moving lid,north boundary
		double rhon;
		for (int i = 1; i < n-1; ++i)
		{
			rhon = f[0][i][m-1] + f[1][i][m-1] + f[3][i][m-1] + 2*(f[2][i][m-1] + f[6][i][m-1] + f[5][i][m-1]);
			f[4][i][m-1] = f[2][i][m-1];
			f[8][i][m-1] = f[6][i][m-1] + rhon*u0/6.0;
			f[7][i][m-1] = f[5][i][m-1] - rhon*u0/6.0;
		}

}
void LidDrivenCavity :: collision()
{
	double t1,t2;
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			t1 = u[i][j]*u[i][j] + v[i][j]*v[i][j] ;
			for (int k = 0; k < d; ++k)
			{
				t2 = u[i][j]*cx[k] + v[i][j]*cy[k];
				feq[k][i][j] = rho[i][j]*w[k]*(1.0 + 3.0*t2 + 4.5*t2*t2 -1.5*t1);
				f[k][i][j] = omega*feq[k][i][j] + (1.0 - omega)*f[k][i][j];
			}
		}
	}
}
void LidDrivenCavity :: streaming()
{
	for (int j = 0; j < m; ++j)
	{
		for (int i = n-1; i > 0; --i)//right to left
		{
			f[1][i][j] = f[1][i-1][j];
		}
		for (int i = 0; i < n-1; ++i)//left to right
		{
			f[3][i][j] = f[3][i+1][j];
		}
	}
	for (int j = m-1; j > 0; --j)
	{
		for (int i = 0; i < n; ++i)
		{
			f[2][i][j] = f[2][i][j-1];
		}
		for (int i = n-1; i > 0; --i)
		{
			f[5][i][j] = f[5][i-1][j-1];
		}
		for (int i = 0; i < n-1; ++i)
		{
			f[6][i][j] = f[6][i+1][j-1];
		}
	}
	for (int j = 0; j < m-1; ++j)
	{
		for (int i = 0; i < n; ++i)
		{
			f[4][i][j] = f[4][i][j+1];
		}
		for(int i=0;i<n-1;i++)
		{
			f[7][i][j] = f[7][i+1][j+1];
		}
		for (int i = n-1; i > 0; --i)
		{
			f[8][i][j] = f[8][i-1][j+1];
		}
	}

}
LidDrivenCavity :: LidDrivenCavity()
{
		cout<<"Start of main"<<endl;
		u0 = 0.10;
		sumvelo = 0.0;
		rhoo = 5.0;
		dx = 1.0;
		dy = dx;
		dt =1.0;
		alpha = 0.01;
		Re = u0*m/alpha;
		cout<<"Reynolds Number is :"<<Re<<endl;
		omega = 1.0/(3*alpha+0.5);
		mstep = 40000;
		init_Basic();
		for (int j = 0; j < m; ++j)
		{
			for (int i = 0; i < n; ++i)
			{
				rho[i][j] = rhoo;
				u[i][j] = 0.0;
				v[i][j] = 0.0;
			}
		}
		// Top Plate Velocity condition
		for (int i = 1; i < n-1; ++i)
		{
			u[i][m] = u0;
			v[i][m] = 0.0;
		}
		init_Basic();
		ofstream myfile ("/home/harshavardhanbabu/timeu.dat");
		myfile<<"Time\tu\tv\r\n";
		for (int i = 0; i < 1000; ++i)
		{
			cout<<i<<endl;
			collision();
			streaming();
			sfbound();
			rhouv();
			
						myfile<<i<<"\t"<<u[n/2][m/2]<<"\t"<<v[n/2][m/2]<<"\r\n";
						// if(i == 100)
			// {
			// 	ofstream myfile ("/home/harshavardhanbabu/uv.dat");
  	// 			if (myfile.is_open())
  	// 				{
  	// 					for (int i = 0; i < n; ++i)
  	// 					{
	  // 						for (int j = 0; j < m; ++j)
	  // 						{
	  // 							myfile<<u[i][j]<<"\t"<<v[i][j]<<endl;
	  // 						}
  	// 					}
    					
   //  					myfile.close();
    					
   //  				}
   //  			else
   //  			{
   //  				cout << "Unable to open file";
					
   //  			} 
   //  			break;	
			// }
		}
		myfile.close();
		cout<<"Done!"<<endl;

}
LidDrivenCavity :: ~LidDrivenCavity()
{
	cout<<"That's it.."<<endl;	
}
void LidDrivenCavity :: init_Basic()
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
	// cout<<"The directions assigned are as follows:"<<endl;
	// for(int i=0;i<=8;++i)
	// {
	// 	cout<<"e"<<i<<" = "<<"("<<cx[i]<<","<<cy[i]<<");"<<endl;
	// }
	/*fill_n( f, sizeof f / sizeof **f, 0 );
	fill_n( feq, sizeof feq / sizeof **feq, 0 );
	fill_n( rho, sizeof rho / sizeof **rho, 0 );
	fill_n( u, sizeof u / sizeof **u, 0 );
	fill_n( v, sizeof v / sizeof **v, 0 );
	fill_n( strf, sizeof strf / sizeof **strf, 0 );*/
	double f[d][n][m];//={0};
	 double feq[d][n][m];//={0};
	 double rho[n][m];//={0};
	 double u[n][m];//={0};
	 double v[n][m];//={0};
	 double strf[n][m];//={0};
	 
	 	for (int i = 0; i < n; ++i)
	 		{
	 			for (int j = 0; j < m; ++j)
	 			{
	 				rho[i][j]=0;
	  u[i][j]=0;
	 v[i][j]=0;
	  strf[i][j]=0;
	 				for (int k = 0; k < d; ++k		)
						 {
						 	 f[d][n][m]=0;
	 feq[k][i][j]=0;
	 

	 			}
	  		}
	  	}
	
}

int main()
{
	cout<<"In main"<<endl;
	LidDrivenCavity l;

	return 0;
}