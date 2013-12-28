/*-----------------------------------------------
 Poiseuille Flow due to Pressure Gradient alone
---------------------------------------------------*/

#include <stdio.h>
#include <math.h>
  //9-speed square lattice
#define np 8	

// Input variable

#define tau 1.1//2.155
#define rho_in 2.0      //
#define rho_out 1.0      //
#define rho0 1.5      //
#define uxo 0
#define uyo 0
#define ref 0.7  

#define taut 1.4//
#define T_in 1.0      //
#define T_out 1.0      //
#define T0 1.0          //  
#define tref 0.65 //
#define dkt 10.0 
#define delH 0.0 

#define tauc 1.4//
#define C_in 0.2      //
#define C_out 0.001      //
#define C0 0.001          //
#define cref 0.65 //
#define dka 0.00 //0.1

//dely=width/(ny+1); delx=dely
//nx=(int)length/delx;
//x_step=nx/4;

#define nx 339    // nx=no of x-grids -1  //
#define ny 33  // ny=no of y-grids -1   //
#define time_step 250 //write up
#define time 5000
#define time_step_vel 500 //for vel
#define time_vel 10000
#define x_step 25 // print out
#define y_step 4
#define length 12.0//length of the channel(in meter)  //
#define width  1.2 //Width of the channel(in m)   //
#define esp 1.0 //Speed of the sound
/*----------------------------------------------------------------------------
			Declaring global variables
------------------------------------------------------------------------------*/
double f[np+3][nx+3][ny+3]={0}; //Local distribution function
double feq[np+3][nx+3][ny+3]={0}; //Equilibrium distribution function
double fn[np+3][nx+3][ny+3]={0}; // A step ahead euilibrium function
double rho[nx+3][ny+3]={0}; // Local density
double ux[nx+3][ny+3]={0}; //Local velocity in x-direction
double uy[nx+3][ny+3]={0}; //Local velocity in Y-direction


int ex[]={0,1,0,-1,0,1,-1,-1,1}; //ei's are direction vecors of 9-speed square lattice
int ey[]={0,0,1,0,-1,1,1,-1,-1};

double rho_grad=0.0,length_correct=0.0,delx=0.0,dely=0.0,delt=0.0,dknud[nx+3]={0.0}, dknud0=0.0, Betat=0.0,T_grad=0.0, Dm = 0.0, Betacc = 0.0, Betatt = 0.0;

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
 			 Initializing macroscopic properties
------------------------------------------------------------------------------------*/			
 			 void initialization(){
 			 	int i,j;
 			 	for(i=0;i<=nx+1;i++){
 			 		for(j=0;j<=ny+1;j++){
 			 			rho[i][j]=rho0;
 			 			ux[i][j]=uxo;
 			 			uy[i][j]=uyo;
 			 		}
 			 	}
 			 }
 			 double pLinear(double x)
 			 {
 			 	return ((rho_in) + (x/length)*(rho_out-rho_in));
 			 }
/*----------------------------------------------------------------------------------


------------------------------------------------------------------------------------
	Defining Equilibrium Distribution function for 9 speed square lattice
-----------------------------------------------------------------------------------*/
	void equilibrium(){
		int i,j,k;

	// Equilibrium distribution function at i=0
		for(i=0;i<=nx+1;i++){
			for(j=0;j<=ny+1;j++){
				feq[0][i][j]=4.0*rho[i][j]*(1.0-3.0*(ux[i][j]*ux[i][j]+uy[i][j]*uy[i][j])/2.0)/9.0;
			}
		}
	// Equilibrium distribution function for i=1,2,3,4

		for(i=0;i<=nx+1;i++){
			for(j=0;j<=ny+1;j++){
				for(k=1;k<=4;k++){
					feq[k][i][j]=rho[i][j]*(1.0+3.0*(ux[i][j]*ex[k]+uy[i][j]*ey[k])-3.0*(ux[i][j]*ux[i][j]+uy[i][j]*uy[i][j])/2.0+9.0*(ex[k]*ux[i][j]+ey[k]*uy[i][j])*(ex[k]*ux[i][j]+ey[k]*uy[i][j])/2.0)/9.0;
				}
			}
		}

	// Equilibrium distribution function for i=5,6,7,8

		for(i=0;i<=nx+1;i++){
			for(j=0;j<=ny+1;j++){
				for(k=5;k<=np;k++){
					feq[k][i][j]=rho[i][j]*(1.0+3.0*(ux[i][j]*ex[k]+uy[i][j]*ey[k])-3.0*(ux[i][j]*ux[i][j]+uy[i][j]*uy[i][j])/2.0+9.0*(ex[k]*ux[i][j]+ey[k]*uy[i][j])*(ex[k]*ux[i][j]+ey[k]*uy[i][j])/2.0)/36.0;
				}
			}
		}

	}
/*----------------------------------------------------------------------------------
 *
 -----------------------------------------------------------------------------------
 			Initializing distribution function
----------------------------------------------------------------------------------*/

 			void initialize_distribution_function(){
 				int i,j,k;
 				for(i=0;i<=nx+1;i++){
 					for(j=0;j<=ny+1;j++){
 						for(k=0;k<=np;k++){
 							f[k][i][j]=feq[k][i][j];
 						}
 					}
 				}
 			}

/*----------------------------------------------------------------------------------

  ----------------------------------------------------------------------------------
  		Boundary conditions and Streaming
----------------------------------------------------------------------------------*/

  		void boundary_conditions_and_streaming(){
  			int i,j,k,p,q;
	//First of all calculate unknown fi's
	//Wall below
  			for(i=1;i<=nx;i++){
	                               // f(r+dr,t+dt)=f(r,t) (Streaming)
  				f[0][i][0]=fn[0][i][0];
  				f[1][i][0]=fn[1][i-1][0];
  				f[3][i][0]=fn[3][i+1][0];
  				f[4][i][0]=fn[4][i][1];
  				f[7][i][0]=fn[7][i+1][1];
  				f[8][i][0]=fn[8][i-1][1];
					// Above f's are known

		f[2][i][0]=f[4][i][0];//f[4][i][0];
		f[5][i][0]=ref*f[7][i][0]+(1.0-ref)*f[8][i][0];//f[7][i][0]-(f[1][i][0]-f[3][i][0])/2;
		f[6][i][0]=ref*f[8][i][0]+(1.0-ref)*f[7][i][0];//f[8][i][0]+(f[1][i][0]-f[3][i][0])/2;
	}

	//Wall above
	for(i=1;i<=nx;i++){
	                   	// f(r+dr,t+dt)=f(r,t) (Streaming)
		f[0][i][ny+1]=fn[0][i][ny+1];
		f[1][i][ny+1]=fn[1][i-1][ny+1];
		f[2][i][ny+1]=fn[2][i][ny];
		f[3][i][ny+1]=fn[3][i+1][ny+1];
		f[5][i][ny+1]=fn[5][i-1][ny];
		f[6][i][ny+1]=fn[6][i+1][ny];
		
		f[4][i][ny+1]=f[2][i][ny+1];//f[2][i][ny+1];
		f[8][i][ny+1]=ref*f[6][i][ny+1]+(1.0-ref)*f[5][i][ny+1];//f[6][i][ny+1]-(f[1][i][ny+1]-f[3][i][ny+1])/2;
		f[7][i][ny+1]=ref*f[5][i][ny+1]+(1.0-ref)*f[6][i][ny+1];//f[5][i][ny+1]+(f[1][i][ny+1]-f[3][i][ny+1])/2;
	}

	//Inlet
	for(j=1;j<=ny;j++){


		f[0][0][j]=fn[0][0][j];
		f[2][0][j]=fn[2][0][j-1];
		f[3][0][j]=fn[3][1][j];
		f[4][0][j]=fn[4][0][j+1];
		f[6][0][j]=fn[6][1][j-1];
		f[7][0][j]=fn[7][1][j+1];

		ux[0][j]=1-(f[0][0][j]+f[2][0][j]+f[4][0][j] +2*(f[3][0][j]+f[6][0][j]+f[7][0][j]))/rho_in;

		f[1][0][j]=f[3][0][j]+(2*rho_in*ux[0][j])/3;
		f[8][0][j]=rho_in*ux[0][j]/6+f[6][0][j]+(f[2][0][j]-f[4][0][j])/2;
		f[5][0][j]=rho_in-(f[0][0][j]+f[1][0][j]+f[2][0][j]+f[3][0][j]+f[4][0][j]+f[6][0][j]+f[7][0][j]+f[8][0][j]);
	}

	//outlet
	for(j=1;j<=ny;j++){


		f[0][nx+1][j]=fn[0][nx+1][j];
		f[1][nx+1][j]=fn[1][nx][j];
		f[2][nx+1][j]=fn[2][nx+1][j-1];
		f[4][nx+1][j]=fn[4][nx+1][j+1];
		f[5][nx+1][j]=fn[5][nx][j-1];
		f[8][nx+1][j]=fn[8][nx][j+1];

		ux[nx+1][j]=-1+(f[0][nx+1][j]+f[2][nx+1][j]+f[4][nx+1][j]+2*(f[1][nx+1][j]+f[5][nx+1][j]+f[8][nx+1][j]))/rho_out;


		f[3][nx+1][j]=f[1][nx+1][j]-2*rho_out*ux[nx+1][j]/3;
		f[7][nx+1][j]=f[5][nx+1][j]+(f[2][nx+1][j]-f[4][nx+1][j])/2-rho_out*ux[nx+1][j]/6;
		f[6][nx+1][j]=rho_out-(f[0][nx+1][j]+f[1][nx+1][j]+f[2][nx+1][j]+f[3][nx+1][j]+f[4][nx+1][j]+f[5][nx+1][j]+f[7][nx+1][j]+f[8][nx+1][j]);
	}

	//Bottom inlet corner


	f[0][0][0]=fn[0][0][0];
	f[3][0][0]=fn[3][1][0];
	f[4][0][0]=fn[4][0][1];
	f[7][0][0]=fn[7][1][1];

	f[1][0][0]=f[3][0][0];
	f[2][0][0]=f[4][0][0];
	f[5][0][0]=f[7][0][0];
	f[8][0][0]=(rho_in-(f[0][0][0]+2*(f[3][0][0]+f[4][0][0]+f[7][0][0])))/2;
	f[6][0][0]=f[8][0][0];


	// Top Inlet corner


	f[0][0][ny+1]=fn[0][0][ny+1];
	f[2][0][ny+1]=fn[2][0][ny];
	f[3][0][ny+1]=fn[3][1][ny+1];
	f[6][0][ny+1]=fn[6][1][ny];


	f[1][0][ny+1]=f[3][0][ny+1];
	f[4][0][ny+1]=f[2][0][ny+1];
	f[5][0][ny+1]=(rho_in-(f[0][0][ny+1]+2*(f[2][0][ny+1]+f[3][0][ny+1]+f[6][0][ny+1])))/2;
	f[7][0][ny+1]=f[5][0][ny+1];
	f[8][0][ny+1]=f[6][0][ny+1];


	//Bottom outlet corner

	f[0][nx+1][0]=fn[0][nx+1][0];
	f[1][nx+1][0]=fn[1][nx][0];
	f[4][nx+1][0]=fn[4][nx+1][1];
	f[8][nx+1][0]=fn[8][nx][1];

	f[3][nx+1][0]=f[1][nx+1][0];
	f[2][nx+1][0]=f[4][nx+1][0];
	f[5][nx+1][0]=(rho_out-(f[0][nx+1][0]+2*(f[1][nx+1][0]+f[4][nx+1][0]+f[8][nx+1][0])))/2;
	f[7][nx+1][0]=f[5][nx+1][0];
	f[6][nx+1][0]=f[8][nx+1][0];


	//Top outlet corner

	f[0][nx+1][ny+1]=fn[0][nx+1][ny+1];
	f[1][nx+1][ny+1]=fn[1][nx][ny+1];
	f[2][nx+1][ny+1]=fn[2][nx+1][ny];
	f[5][nx+1][ny+1]=fn[5][nx][ny];

	f[3][nx+1][ny+1]=f[1][nx+1][ny+1];
	f[4][nx+1][ny+1]=f[2][nx+1][ny+1];
	f[7][nx+1][ny+1]=f[5][nx+1][ny+1];	
	f[6][nx+1][ny+1]=(rho_out-(f[0][nx+1][ny+1]+2*(f[1][nx+1][ny+1]+f[2][nx+1][ny+1]+f[5][nx+1][ny+1])))/2;
	f[8][nx+1][ny+1]=f[6][nx+1][ny+1];


}

 /* --------------------------------------------------------------------------------------------------
 				Streaming of particle inside channel (Excluding boundaries)
 -----------------------------------------------------------------------------------------------------*/

 				void streaming(){
       				  // f(r+dr,t+dt)=f(r,t)   For all points lying inside channel except boundary points
 					int i,j,k;
 					double sum;
 					for(i=1;i<=nx;i++){
 						for(j=1;j<=ny;j++){
 							f[0][i][j]=fn[0][i][j];
 							f[1][i][j]=fn[1][i-1][j];
 							f[2][i][j]=fn[2][i][j-1];
 							f[3][i][j]=fn[3][i+1][j];
 							f[4][i][j]=fn[4][i][j+1];
 							f[5][i][j]=fn[5][i-1][j-1];
 							f[6][i][j]=fn[6][i+1][j-1];
 							f[7][i][j]=fn[7][i+1][j+1];
 							f[8][i][j]=fn[8][i-1][j+1];
 						}
 					}
 				}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------	
void calculate_velocity_density(){ 	// It calculates velocity and density by updated distribution function i.e at t=t+dt
	int i,j,k;
	double sum,sumux,sumuy;
	// CALCULATION OF DENSITY
	for(i=0;i<=nx+1;i++){
		for(j=0;j<=ny+1;j++){
			rho[i][j]=0;
		}
	}
	
	for(i=0;i<=nx+1;i++){
		for(j=0;j<=ny+1;j++){
			sum=0;
			for(k=0;k<=np;k++){
				sum+=f[k][i][j];
			}
			rho[i][j]=sum;
			
		}
	}

	for(j=0;j<=ny+1;j++){
		rho[0][j]=rho_in;
		rho[nx+1][j]=rho_out;
	}
	

	// CALCULATION OF VELOCITIES
	for(i=0;i<=nx+1;i++){
		for(j=0;j<=ny+1;j++){
			sumux=0;
			sumuy=0;
			for(k=0;k<=np;k++){
				sumux+=ex[k]*f[k][i][j];
				sumuy+=ey[k]*f[k][i][j];
			}	
			ux[i][j]=sumux/rho[i][j];
			uy[i][j]=sumuy/rho[i][j];
		}
	}
}


/*-----------------------------------------------------------------------------

  -----------------------------------------------------------------------------
  	 Calculation of distribution function from boltzmann transport eqn
------------------------------------------------------------------------------*/

  	 void collision(){	
	//f'(r,t)=f(r,t)-1/tau(f(r,t)-feq(r,t)) this is due to collision only
	//f(r+dr,t+dt)=f(r,t) Particle undergoes streaming just after collision
  	 	int i,j,k;
  	 	double tauv;
  	 	for(k=0;k<=np;k++){
  	 		for(i=0;i<=nx+1;i++){
  	 			for(j=0;j<=ny+1;j++){
  	 				tauv=0.5 + rho0/rho[i][j]*(tau-0.5); 
  	 				fn[k][i][j]=f[k][i][j]*(1.0-1.0/tauv)+(1.0/tauv)*feq[k][i][j];
  	 			}
  	 		}
  	 	}
  	 }



/*-------------------------------------------------------------------MAIN FUNCTION----------------------------------------------------------------------------------------------------*/
  	 main()
  	 {

  	 	int i,j,k,itr,step,num_file=0,num_file1, cnum_file, tnum_file;
  	 	double sum,ux_avg,rhox,delp,dneu,uxmax,Pr,y,x,uavg[nx+3]={0.0},up[nx+3][ny+3]={0.0},term[nx+3]={0.0},Nusselt[nx+3]={0.0},Betav=0.0, Re=0.0, zstar=0.0,Beta=0.0, Betac = 0.0, Sc = 0.0, Sh[nx+3]={0.0};
	double max[nx+3]; //Maximum velocity at each x-point (Numerical)
	double uavgcal[nx+3]={0};
	double Tm[nx+3]={0}, Cm[nx+3]={0};
	char *file_name[]={"ut0.dat","ut1000.dat","ut2000.dat","ut3000.dat","ut4000.dat","ut5000.dat","ut6000.dat","ut7000.dat","ut8000.dat","ut9000.dat","ut10000.dat","ut11000.dat","ut12000.dat","ut13000.dat","ut14000.dat","ut15000.dat","ut16000.dat","ut17000.dat","ut18000.dat","ut19000.dat","ut20000.dat","ut21000"};
	char *file_name1[]={"y_uy_0","y_uy_1","y_uy_2","y_uy_3","y_uy_4","y_uy_5"};
	char *file_name2[]={"y_ux_0","y_ux_1","y_ux_2","y_ux_3","y_ux_4","y_ux_5"};
	char *file_name3[]={"rhox_0","rhox_1","rhox_2","rhox_3","rhox_4","rhox_5","rhox_6"};	
	char *file_name4[]={"rhoy_0","rhoy_1","rhoy_2","rhoy_3","rhoy_4","rhoy_5","rhoy_6"};	
	
	FILE *fp1,*fp2,*fp3, *confp, *Tonfp;
//--------------Calculation of delx,dely,delt,delp,dneu,---------------------------------------------------	
	dely=width/(ny+1);
	delx=dely;
       //delx = length/(nx+1);
       //dely=delx;
	length_correct=delx*(nx+1);
	delt=delx/esp;
        delp=(rho_in-rho_out)*esp*esp/3.0;      //Pressure drop across channel
        dneu=(2.0*tau-1.0)/6.0*(esp*esp*delt);     //viscocity
        uxmax=(width*width/4.0)*(delp/length_correct)/(2.0*rho0*dneu); //Maximum velocity
		dknud0 = sqrt(22.0/7.0/24.0)*(2*tau-1.0)/ny*rho0/rho_out; //Tang et al. use rho0 = rho_out in simul to match Zhang's Kn
	//	dknud = sqrt(8.0*7.0/3.0/22.0)*(tau-0.5)/ny; // Zhang et al.
		uavg[nx+1] = (1.0/6.0 + (1.0-ref)/ref*1.1466*dknud0+2.0*0.31*dknud0*dknud0)*delp/length_correct*width*width/(2*dneu*rho0); //Zhang et al.
		Dm = (2.0*tauc-1.0)/6.0*(esp*esp*delt); 
		Betat = 1.66*(1.0-tref)/tref;
		Betav = (1.0-ref)/ref;
		Betac = (1.0-cref)/cref;
		Betacc = Betac * dknud0 * dka / Dm;
		Re = uavg[nx+1]*(2*width)/dneu;
		zstar = Re * 0.7 * 4 * width;
		Beta = Betat/Betav;
		Pr = (2*tau-1.0)/(2*taut-1.0);
		Sc = (2*tau-1.0)/(2*tauc-1.0);
		Betatt = -dka * delH/ dkt ;  

		printf("length_correct = %lf dknud0 = %lf Re = %lf  uavg = %lf   Betat = %lf   Betatt =  %lf   Betav = %lf  Beta = %lf  Betac =  %lf   Betacc = %lf  Bv*Kn = %lf   Pr = %lf   Sc = %lf   Dm =  %lf\n",length_correct, dknud0, Re, uavg[nx+1],Betat, Betatt, Betav, Beta, Betac, Betacc, Betav*dknud0, Pr, Sc, Dm);

//-----------------------------------------------------------------------------------------------------------	
	initialization(); // It initialises all macroscopic variables
 	equilibrium(); // It calculates equilibrium distribution function
	initialize_distribution_function(); //It initializes all distribution function
	
	
	num_file=0;
	//Put constraint at inlet and outlet
	for(j=0;j<=ny+1;j++)
	{
		rho[0][j]=rho_in;
		rho[nx+1][j]=rho_out;
	}
	 //Time loop starts from here..
	for(itr=0;itr<=time_vel;itr++){
	         	equilibrium();// It calculates equilibrium distribution function at each point	
				collision();  // f'(r+dr,t+dt)=f(r,t)-1/tau(f(r,t)-feq(r,t)
		       boundary_conditions_and_streaming(); // f(r+dr,d+dt)=f'(r+dr,d+dt) at boundary*/  
		       streaming();  //f(r+dr,d+dt)=f'(r+dr,d+dt) inside the channel excluding boundary
		       calculate_velocity_density(); // It calculates density and velocity at (r+dr,t+dt)

		       for(i=0;i<nx+1;i++)
		       {
		       	uavgcal[i] = 0.0;
		       	y=0.0;
		       	for(j=0;j<=ny;j++)
		       	{
		       		uavgcal[i] = uavgcal[i] + 0.5*dely/width*(ux[i][j] + ux[i][j+1]); 
		       		y+=dely;
		       	}
		       }

		       for(i=0;i<nx+1;i++)
		       {
		       	dknud[i]=dknud0/rho[i][0]*rho_out; 
		       	term[i] = Betav*1.1466*dknud[i] + 2*0.31*dknud[i]*dknud[i];
		       	uavg[i] = 1.0/6.0 + term[i];
		       	y=0;
		       	for(j=0;j<=ny+1;j++)
		       	{    
		       		up[i][j] = -(y/width)*(y/width) + y/width + term[i];
		       		y+=dely;
		       	}
		       }     
		       if(itr%time_step_vel == 0){
		       	printf("iteration is going for Time=%d to Time=%d\n",itr,(itr+time_step_vel));
		       	fp1=fopen(file_name[num_file],"w");
		       	num_file++;

		       	y=0;
		       	for(j=0;j<=ny+1;j++){                                   
		       		//fprintf(fp1,"%lf %lf %lf %lf  %lf   %lf   %lf    %lf   %lf  %lf  %lf\n",ux[2][j]/uavgcal[2],ux[nx/4][j]/uavgcal[nx/4],ux[nx/2][j]/uavgcal[nx/2],ux[3*nx/4][j]/uavgcal[3*nx/4],ux[nx-1][j]/uavgcal[nx-1],up[2][j]/uavg[2],up[nx/4][j]/uavg[nx/4],up[nx/2][j]/uavg[nx/2],up[3*nx/4][j]/uavg[3*nx/4],up[nx-1][j]/uavg[nx-1],y/width);
		       		fprintf(fp1, "%lf\t%lf\n",(y/width),ux[nx-1][j]/uavgcal[nx-1]);
		       		y+=dely;
		       	}
		       	fclose(fp1);

		       }

		   }
		   // Writing SlipVelocity Vs Length
			x=0.0;
			i=0;
			fp1=fopen("SlipVelocity.dat","w");
			for (i = 0; i <= nx+1; ++i)
			{
				fprintf(fp1,"%lf\t%lf\r\n",x,ux[i][0]);
				x=x+delx;
			}
				
			fclose(fp1);	     

  //------------------------------------------------------------------------------------------------------------	
	// Writing density and x in file
			fp2=fopen("rho_x.dat","w");
			x=0.0;
			for(i=0;i<=nx+1;i++){
				fprintf(fp2,"%lf %lf\n",x,rho[i][(ny+1)/2]);
				x=x+delx;
			}
			fclose(fp2);
			// Writing Non-Dimensional Pressure 
			fp3 = fopen("P_x.dat","w");
			x = 0.0;
			double p_linear = 0.0;
			double p_nond  = 0.0;
			for(i=0;i<=nx+1;i++)
			{
				p_linear = pLinear(x);
				p_nond = (rho[i][(ny+1)/2] - p_linear)/(rho_out);
				fprintf(fp3,"%lf\t%lf\n",(x/length),p_nond);
				x=x+delx;
			}
			fclose(fp3);
/*
//-----------------writting ux and y at different x---------------------------
			num_file=0;
			x=0.0;
			i=0;
			do{
				fp1=fopen(file_name1[num_file],"w");

				y=0.0;
				for(j=0;j<=ny+1;j++){
					fprintf(fp1,"%lf %lf\n",ux[i][j],y);
					y=y+dely;
				}
				i=i+(int)x_step;
				x=x+x_step*delx;
				fclose(fp1);
				num_file++;
			}while(x<=length);

//--------------------- Writting uy and y at different x in file---------------------------------
			num_file=0;
			x=0.0;
			i=0;
			do{
				fp1=fopen(file_name2[num_file],"w");
				y=0.0;
				for(j=0;j<=ny+1;j++){
					fprintf(fp1,"%lf ux[nx-1][j]/uavgcal[nx-1]%30.28lf\n",y,uy[i][j]);
					y=y+dely;
				}
				x=x+x_step*delx;
				i=i+(int)x_step;
				fclose(fp1);
				num_file++;
			}while(x<=length);

//------------------Writing rho and x ux[nx-1][j]/uavgcal[nx-1]at different y-ordinate-----------------------------------------
			num_file=0;
			y=0.0;
			j=0;
			do{
				fp1=fopen(file_name3[num_file],"w");
				x=0;
				for(i=0;i<=nx+1;i++){
					fprintf(fp1,"%lf %lf\n",x,rho[i][j]);
					x=x+delx;
				}
				j=j+(int)y_step;
				y=y+y_step*dely;
				fclose(fp1);
				num_file++;
			}while(y<=width);

//-------------- Writting rho and y at different x---------------------------------------------------

			num_file=0;
			x=0.0;
			i=0;
			do{
				fp1=fopen(file_name4[num_file],"w");
				y=0;
				for(j=0;j<=ny+1;j++){
					fprintf(fp1,"%lf %lf\n",y,rho[i][j]);
					y=y+dely;
				}
				i=i+(int)x_step;
				x=x+x_step*delx;
				fclose(fp1);
				num_file++;
			}while(x<=length);

*/
		}
