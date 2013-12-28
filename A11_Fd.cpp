//The finitie difference method;

#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>

using namespace std;

const int m = 101;
const int mstep = 400;


class Temp
{
	double dens,fo[m],f[m];
	double alpha;
	double mstep;
	double dx;
	double dt;
	double x;
public:
	void init();
	Temp();
	~Temp();

	/* data */
};
Temp :: ~Temp()
{
	cout<<"That's it.."<<endl;
}
void Temp :: init()
{
	dx = 1.0;
	dt = 0.5;
	alpha = 0.25;
	for (int i = 0; i < m; ++i)
	{
		fo[i] = 0.0;
		f[i] = 0.0;
	}
	fo[0] = 1.0;
	f[0] = 1.0;
	fo[m-1] = fo[m-2];
	f[m-1] = f[m-2];	
}
Temp :: Temp()
{
	init();
	for(int i=1;i <= mstep;i++)
	{
		for (int j = 1; j < m-1; ++i)
		{
			f[j] = fo[j] + dt*alpha*(fo[j+1] - 2*fo[j] + fo[j-1])/(dx*dx);
		}
		for (int j = 1; j < m-1; ++i)
		{
			fo[i] = f[i];
		}
		fo[m-1] = f[m-2];
	}
	x = 0.0;
	ofstream datfile("/home/harshavardhanbabu/temperature.dat");
	for (int i = 0; i < m; ++i)
	{
		datfile<<x<<"\t"<<f[i]<<"\r\n";
		x = x+dx;
	}
	datfile.close();

}
int main()
{
	Temp t;
	return 0;
}