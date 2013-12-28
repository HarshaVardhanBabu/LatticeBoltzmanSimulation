#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>

using namespace std;

// Code is taken from A.A. Mohamad Lattice Boltzmann method Fundamentals and Engineering Applicaitons with Computer codes..

// Check for the correct path 

const double pi = 3.141592653589793;
const int m = 101; // No.of.Lattice nodes + 1;
class Diffusion
{
	void init();
	void collision();
	void streaming();
	void boundaryCond();
	double fo[m],f1[m],f2[m],rho[m],feq[m],x[m];
	double dx,dt,csq,omega,alpha,mstep,twall;
public:
	Diffusion();	
	Diffusion(int dx,int dt,double alpha,int mstep,double twall);
	void startSimulation();
	void writeResults();
	void generateGnuPlotFile();
	~Diffusion();

	/* data */
};
void Diffusion :: collision()
{
	for (int i = 0; i < m; ++i)
	{
		  rho[i] = f1[i] + f2[i];
          feq[i] = 0.5*rho[i];
          f1[i] = (1.0-omega)*f1[i]+omega*feq[i];
          f2[i] = (1.0-omega)*f2[i]+omega*feq[i];
	}
}
void Diffusion :: streaming()
{
		for(int i =1;i<m-1;i++)
        {
          f1[m-i-1] = f1[m-i-2];
          f2[i-1] = f2[i];
        }
}
void Diffusion :: boundaryCond()
{
	f1[0] = twall - f2[0];
	f1[m-1] = f1[m-2];
	f2[m-1] = f2[m-2];
}
void Diffusion :: startSimulation()
{
	for (int i = 1; i < mstep; ++i)
		{
			collision();
			streaming();
			boundaryCond();
		}
		cout<<"Simulation Completed Successfully!!"<<endl;
	
}
void Diffusion :: writeResults()
{
	ofstream results("/home/harshavardhanbabu/A11.dat");
	if (results.is_open())
	{
		for (int i = 0; i < m; ++i)
		{
			results<<x[i]<<"\t"<<rho[i]<<"\r\n";
		}
		cout<<"File Written Successfully.."<<endl;
	}
	else
	{
		cout<<"Writing failed!! Check for the path.."<<endl;
	}
}
void Diffusion :: init()
{
	x[0] = 0;
	for (int i = 1; i < m; ++i)
	{
		x[i] = x[i-1] + dx;
	}
	csq = (dx*dx)/(dt*dt);
	omega = 1.0/(alpha/(dt*csq)+0.5);
	for (int i = 0; i < m; ++i)
	{
		rho[i] = 0.0;
		f1[i] = 0.5*rho[i];
		f2[i] = 0.5*rho[i];
	}
}

Diffusion :: Diffusion(int dx,int dt,double alpha,int mstep,double twall)
{
	this->dx = dx;
	this->dt = dt;
	this->alpha = alpha;
	this->mstep = mstep;
	this->twall = twall;
	init();
}

Diffusion :: ~Diffusion()
{

}
Diffusion :: Diffusion()
{

}
void Diffusion :: generateGnuPlotFile()
{
	ofstream gnuplot("/home/harshavardhanbabu/A11.plt");
	gnuplot<<"set terminal svg enhanced size 1000 1000 fname \" 1D Diffusion \" fsize 36\r\n";
	gnuplot<<"set output \"Diffusion.svg\" \r\n ";
	gnuplot<<"set title \"Plot of T Vs x\" \r\n";
	gnuplot<<"set xlabel \"x\" \r\n ";
	gnuplot<<"set ylabel \"T\" \r\n";
	gnuplot<<"plot \"A11.dat\" using 1:2 title \"\" \r\n";
	gnuplot.close();
}
int main(int argc, char const *argv[])
{
	Diffusion d(1.0,1.0,0.25,200,1.0);
	d.startSimulation();
	d.writeResults();
	d.generateGnuPlotFile();
	return 0;
}
