/*---------------------------------------------------
 Poiseuille Flow due to Pressure Gradient and Gravity
-----------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define np 8	 //9-speed square lattice
// Input variable
#define tau 2.0

#define flag 1

// flag = 1 Only Pressure driven.
// flag = 2 Only Gravity driven.
// flag = 3 Both Pressure and Gravity. 

//Make rho_in = rho_out = rho0 to make pressure gradient zero 

#define rho_in 2.01
#define rho_out 1.99
#define rho0 2.00

#define uxo 0
#define uyo 0

 // Gravity value. Zero for no gravity.
#define g 0.01 
//dely=width/(ny+1); delx=dely
//nx=(int)length/delx;
//x_step=nx/4;

#define nx 199    // nx=no of x-grids -1
#define ny 9   // ny=no of y-grids -1

 // Make the Changes in timestep depending up on the flag
// Make sure that time is ten times timestep or else increase the files in the array in main
#define time_step 10 //write up
#define time 100

#define x_step 131 // print out
#define y_step 25
#define length 200.0//length of the channel(in meter)
#define width  10.0 //Width of the channel(in m)
#define esp 1.0 //Speed of the sound
#define snapshot_time 1000 // for 4 points in the tube
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
// Weights and direction vectors for gravity inclusion
double w[] = {0,1.0/3,1.0/3,1.0/3,1.0/3,1.0/12,1.0/12,1.0/12,1.0/12};// 
double grav[]={0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0 };
double delt;
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
/*----------------------------------------------------------------------------------


------------------------------------------------------------------------------------
	Defining Equilibrium Distribution function for 9 speed square lattice
-----------------------------------------------------------------------------------*/
	void equilibrium(){
		int i,j,k;

	// Equilibrium distribution function at i=0
		for(i=0;i<=nx+1;i++){
			for(j=0;j<=ny+1;j++){
				feq[0][i][j]=4.0*rho[i][j]*(1.0-3.0*(ux[i][j]*ux[i][j]+uy[i][j]*uy[i][j])/2.0)/9.;
			}
		}
	// Equilibrium distribution function for i=1,2.0,3,4

		for(i=0;i<=nx+1;i++){
			for(j=0;j<=ny+1;j++){
				for(k=1;k<=4;k++){
					feq[k][i][j]=rho[i][j]*(1.0+3.0*(ux[i][j]*ex[k]+uy[i][j]*ey[k])-3.0*(ux[i][j]*ux[i][j]+uy[i][j]*uy[i][j])/2.0+9*(ex[k]*ux[i][j]+ey[k]*uy[i][j])*(ex[k]*ux[i][j]+ey[k]*uy[i][j])/2.0)/9.0;
				}
			}
		}

	// Equilibrium distribution function for i=5,6,7,8

		for(i=0;i<=nx+1;i++){
			for(j=0;j<=ny+1;j++){
				for(k=5;k<=np;k++){
					feq[k][i][j]=rho[i][j]*(1.0+3.0*(ux[i][j]*ex[k]+uy[i][j]*ey[k])-3.0*(ux[i][j]*ux[i][j]+uy[i][j]*uy[i][j])/2.0+9*(ex[k]*ux[i][j]+ey[k]*uy[i][j])*(ex[k]*ux[i][j]+ey[k]*uy[i][j])/2.0)/36.0;
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

  				f[2][i][0]=f[4][i][0];
  				f[5][i][0]=f[7][i][0]-(f[1][i][0]-f[3][i][0])/2;
  				f[6][i][0]=f[8][i][0]+(f[1][i][0]-f[3][i][0])/2;
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

  				f[4][i][ny+1]=f[2][i][ny+1];
  				f[8][i][ny+1]=f[6][i][ny+1]-(f[1][i][ny+1]-f[3][i][ny+1])/2;
  				f[7][i][ny+1]=f[5][i][ny+1]+(f[1][i][ny+1]-f[3][i][ny+1])/2;
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
  	 	for(k=0;k<=np;k++)
  	 	{
  	 		for(i=0;i<=nx+1;i++)
  	 		{
  	 			for(j=0;j<=ny+1;j++)
  	 			{
  	 				if (flag == 1 )
  	 				{
  	 					fn[k][i][j]=f[k][i][j]*(1.0-1.0/tau)+(1.0/tau)*feq[k][i][j];
  	 				}
  	 				if(flag == 2 || flag ==3)
  	 				{
  	 					
  	 					fn[k][i][j]=f[k][i][j]*(1.0-1.0/tau)+(1.0/tau)*feq[k][i][j] + g*rho[i][j]*grav[k]*delt*w[k];
  	 				}

  	 				
  	 			}
  	 		}
  	 	}
  	 }	



/*-------------------------------------------------------------------MAIN FUNCTION----------------------------------------------------------------------------------------------------*/
  	 main()
  	 {
  	 	int i,j,k,itr,step,num_file=0,num_file1;
  	 	double sum,ux_avg,rhox,delp,dneu,delx,dely,uxmax,y,x;
	double max[nx+3]; //Maximum velocity at each x-point (Numerical)
	double shear_stress[nx+3][ny+3]={0};
	double dux[nx+3][ny+3]={0};
	char *file_nameue[]={"eut0.dat","eut100.dat","eut200.dat","eut300.dat","eut400.dat","eut500.dat","eut600.dat","eut700.dat","eut800.dat","eut900.dat","eut1000.dat"};
	char *file_nameum[]={"mut0.dat","mut100.dat","mut200.dat","mut300.dat","mut400.dat","mut500.dat","mut600.dat","mut700.dat","mut800.dat","mut900.dat","mut1000.dat"};
	char *file_nameus[]={"sut0.dat","sut100.dat","sut200.dat","sut300.dat","sut400.dat","sut500.dat","sut600.dat","sut700.dat","sut800.dat","sut900.dat","sut1000.dat"};
	char *file_namerho[]={"rhot0.dat","rhot100.dat","rhot200.dat","rhot300.dat","rhot400.dat","rhot500.dat","rhot600.dat","rhot700.dat","rhot800.dat","rhot900.dat","rhot1000.dat"};
	// flag = 1 Only Pressure driven.
// flag = 2 Only Gravity driven.
// flag = 3 Both Pressure and Gravity. 
	char *title[] = {"Only Pressure driven.","Only Gravity driven.","Both Pressure and Gravity."};
	FILE *fp1,*fp2,*fp3,*fp4,*fp5;
	

	if ( (flag == 3 || flag == 1) && (rho_in-rho_out) == 0)
	{
		printf("You set flag to %d but haven't changed the rho_in and rho_out values..\n",flag);
		printf("Terminating...\n");
		exit(0);
	}
	if(flag == 2 && (rho_in-rho_out) != 0)
	{
		printf("You set flag to 2 but have rho_in and rho_out values different ...\n");
		printf("Terminating...\n");
		exit(0);
	}
//--------------Calculation of delx,dely,delt,delp,dneu,---------------------------------------------------	
	dely = width/(ny+1);
	printf("dely is %lf\n Width is%lf\n ", dely,width);
	delx=dely;
	delt=delx/esp;
        delp=(rho_in-rho_out)*esp*esp/3.0;      //Pressure drop across channel
        dneu=(2.0*tau-1.0)/6.0*(esp*esp*delt);     //viscocity
        if(flag == 1)
        {
        	uxmax=(width*width/4.0)*(delp/length)/(2.0*rho0*dneu); //Maximum velocity Pressure Driven alone
        }
        if(flag == 2)
        {
        	 uxmax = (g*width*width)/(8.0*dneu); //Maximum velocity with gravity alone
        }
        if(flag == 3)
        {

        	uxmax = ( (width*width) / (8.0*rho0*dneu) )*fabs((delp/length) - (rho0*g)) ;//Combined
        	printf("Maximum Velocity of combined is %lf\n",uxmax);
        	printf("Delp is %lf=\t, rhog is %lf\n",(delp/length),(rho0*g));
        	printf("Value of absolute term is %lf\n",fabs((delp/length) - (rho0*g)));
        }
        
//-----------------------------------------------------------------------------------------------------------	
	initialization(); // It initialises all macroscopic variables
 	equilibrium(); // It calculates equilibrium distribution function
	initialize_distribution_function(); //It initializes all distribution function
	
	    //Time loop starts from here..
	num_file=0;
	for(itr=0;itr<=time;itr++)
	{
               		//Put constraint at inlet and outlet
		for(j=0;j<=ny+1;j++)
		{
			rho[0][j]=rho_in;
			rho[nx+1][j]=rho_out;
		}

			equilibrium();// It calculates equilibrium distribution function at each point	
			collision();  // f'(r+dr,t+dt)=f(r,t)-1/tau(f(r,t)-feq(r,t)
				      // here f'(r+dr,d+dt) is distribution function just after collision
	       boundary_conditions_and_streaming(); // f(r+dr,d+dt)=f'(r+dr,d+dt) at boundary
	       streaming();  //f(r+dr,d+dt)=f'(r+dr,d+dt) inside the channel excluding boundary
	       calculate_velocity_density(); // It calculates density and velocity at (r+dr,t+dt)
		//Writting data in file for plotting the graph between ux vs y at x=L,x=L/2,x=0

	       if(itr%time_step == 0)
	       {
	       	printf("iteration is going for Time=%d to Time=%d\n",itr,(itr+time_step));
	       	fp1=fopen(file_nameue[num_file],"w");
	       	fprintf(fp1,"#\"For Time=%d to Time=%d\"\n",itr,(itr+time_step));
	 		fprintf(fp1, "#Mode:\t%s\n",title[flag-1]);      	
	       	fp3=fopen(file_nameum[num_file],"w");
	       	fprintf(fp3,"#\"For Time=%d to Time=%d\"\n",itr,(itr+time_step));
	       	fprintf(fp3, "#Mode:\t%s\n",title[flag-1]);
	       	fp4=fopen(file_nameus[num_file],"w");
	       	fprintf(fp4,"#\"For Time=%d to Time=%d\"\n",itr,(itr+time_step));
	       	fprintf(fp4, "#Mode:\t%s\n",title[flag-1]);
	       	fp2=fopen(file_namerho[num_file],"w");
	       	fprintf(fp2,"#\"For Time=%d to Time=%d\"\n",itr,(itr+time_step));
	       	fprintf(fp2, "#Mode:\t%s\n",title[flag-1]);
	       	num_file++;

	       	y=0;
	       	for(j=0;j<=ny+1;j++){
	       		fprintf(fp1,"%lf %lf\n",(ux[nx+1][j]/uxmax),y);
	       		y+=dely;
				//printf("%lf\n",dely );
	       	}
	       	fclose(fp1);
	       	y=0;
	       	for(j=0;j<=ny+1;j++){
	       		fprintf(fp3,"%lf %lf\n",(ux[nx/2][j]/uxmax),y);
	       		y+=dely;
	       	}
	       	fclose(fp3);
	       	y=0;
	       	for(j=0;j<=ny+1;j++)
	       	{
	       		fprintf(fp4,"%lf %lf\n",(ux[0][j]/uxmax),y);
	       		y+=dely;
	       	}
	       	fclose(fp4);

	       	x=0.0;
	       	for(i=0;i<=nx+1;i++){
	       		fprintf(fp2,"%lf %lf\n",x,rho[i][(ny+1)/2]);
	       		x=x+delx;
	       	}
	       	fclose(fp2);
	       }
	   }
	}
