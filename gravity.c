/* THIS PROGRAM GIVES VELOCITY PROFILE FOR COMPRESSIBLE FLUID AND INCOMPRESSIBLE FLUID IN THE CHANNEL
NOTE: SINCE WE ARE USING ARRAY IN THIS PROGRAM.SO DUE TO INSUFFICIENT MEMORY  WE CAN NOT MAKE ANY CHOICE FOR LENGTH AND DIAMETER.THIS PROGRAM WORKS BEST FOR L/D<6
*/

#include<stdio.h>
#define np 8	 //9-speed square lattice
		// Input variable
#define tau 1.0
#define rho_in 2.01
#define rho_out 1.99
#define rho0 2
#define uxo 0
#define uyo 0
#define g 3.0

//dely=width/(ny+1); delx=dely
//nx=(int)length/delx;
//x_step=nx/4;

#define nx 199   // nx=no of x-grids -1
#define ny 9 // ny=no of y-grids -1
// #define time_step 1//write up
// #define time 1
#define time_step 10000 //write up
#define time 50000
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
double w[] = {4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};
double grav[]={0.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0, 1.0};
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
	double delx = width/(ny+1);
	// Equilibrium distribution function at i=0
	for(i=0;i<=nx+1;i++){
		for(j=0;j<=ny+1;j++){
				feq[0][i][j]=4*rho[i][j]*(1-3*(ux[i][j]*ux[i][j]+uy[i][j]*uy[i][j])/2)/9;
		}
	}
	//printf("%lf\n",feq[0][0][0]);
	// Equilibrium distribution function for i=1,2,3,4

	for(i=0;i<=nx+1;i++){
		for(j=0;j<=ny+1;j++){
			for(k=1;k<=4;k++){
	feq[k][i][j]=rho[i][j]*(1+3*(ux[i][j]*ex[k]+uy[i][j]*ey[k])-3*(ux[i][j]*ux[i][j]+uy[i][j]*uy[i][j])/2+9*(ex[k]*ux[i][j]+ey[k]*uy[i][j])*(ex[k]*ux[i][j]+ey[k]*uy[i][j])/2)/9;
			}
		}
	}

	// Equilibrium distribution function for i=5,6,7,8

	for(i=0;i<=nx+1;i++){
		for(j=0;j<=ny+1;j++){
			for(k=5;k<=np;k++){
	feq[k][i][j]=rho[i][j]*(1+3*(ux[i][j]*ex[k]+uy[i][j]*ey[k])-3*(ux[i][j]*ux[i][j]+uy[i][j]*uy[i][j])/2+9*(ex[k]*ux[i][j]+ey[k]*uy[i][j])*(ex[k]*ux[i][j]+ey[k]*uy[i][j])/2)/36;
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
				//printf("%lf is \n",f[k][i][j]);
				}
			//	printf("The velocity @ n/2,m/2 is%lf\n",rho[i][j]);
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
	//double force = 0.0,g = 10.0;
	for(k=0;k<=np;k++){
		for(i=0;i<=nx+1;i++){
			for(j=0;j<=ny+1;j++){
				//force = ex[k]*g*w[k];
		fn[k][i][j]=f[k][i][j]*(1-1/tau)+(1/tau)*feq[k][i][j] + g*rho[i][j]*grav[k]*delt*w[k];
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
	char *file_nameue[]={"geut0.dat","geut100.dat","geut200.dat","geut300.dat","geut400.dat","geut500.dat","geut600.dat","geut700.dat","geut800.dat","geut900.dat","geut1000.dat"};
	char *file_nameum[]={"gmut0.dat","gmut100.dat","gmut200.dat","gmut300.dat","gmut400.dat","gmut500.dat","gmut600.dat","gmut700.dat","gmut800.dat","gmut900.dat","gmut1000.dat"};
	char *file_nameus[]={"gsut0.dat","gsut100.dat","gsut200.dat","gsut300.dat","gsut400.dat","gsut500.dat","gsut600.dat","gsut700.dat","gsut800.dat","gsut900.dat","gsut1000.dat"};
	char *file_namerho[]={"rhot0.dat","rhot100.dat","rhot200.dat","rhot300.dat","rhot400.dat","rhot500.dat","rhot600.dat","rhot700.dat","rhot800.dat","rhot900.dat","rhot1000.dat"};
	//char *gnuplot[] = {"gt0.plt","gt100.plt","gt200.plt","gt300.plt","gt400.plt","gt500.plt","gt600.plt","gt700.plt","gt800.plt","gt900.plt","gt1000.plt"};
	FILE *fp1,*fp2,*fp3,*fp4,*fp5;

//--------------Calculation of delx,dely,delt,delp,dneu,---------------------------------------------------
        dely = width/(ny+1);
        printf("dely is %lf\n Width is%lf\n ", dely,width);
        delx=dely;
        delt=delx/esp;
       // delp=(rho_in-rho_out)*esp*esp/3.0;      //Pressure drop across channel
        dneu=(2.0*tau-1.0)/6.0*(esp*esp*delt);     //viscocity
       // uxmax=(width*width/4.0)*(delp/length)/(2.0*rho0*dneu); //Maximum velocity
//-----------------------------------------------------------------------------------------------------------
	initialization(); // It initialises all macroscopic variables
 	equilibrium(); // It calculates equilibrium distribution function
	initialize_distribution_function(); //It initializes all distribution function

	    //Time loop starts from here..
	 num_file=0;
	for(itr=0;itr<=time;itr++)
	{
               		//Put constraint at inlet and outlet
               	for(j=0;j<=ny+1;j++){
               		rho[0][j]=rho_in;
               		rho[nx+1][j]=rho_out;
               		}

			equilibrium();// It calculates equilibrium distribution function at each point
			collision();  // f'(r+dr,t+dt)=f(r,t)-1/tau(f(r,t)-feq(r,t)
				      // here f'(r+dr,d+dt) is distribution function just after collision

	       boundary_conditions_and_streaming(); // f(r+dr,d+dt)=f'(r+dr,d+dt) at boundary
	       streaming();  //f(r+dr,d+dt)=f'(r+dr,d+dt) inside the channel excluding boundary
	       calculate_velocity_density(); // It calculates density and velocity at (r+dr,t+dt)

	       //Writting data in file for plotting the graph between ux vs y at exit
	       printf("%d Time is\n",itr );
		if(itr%time_step == 0)
		{
		printf("iteration is going for Time=%d to Time=%d\n",itr,(itr+time_step));
	       	fp1=fopen(file_nameue[num_file],"w");
	       	fp3=fopen(file_nameum[num_file],"w");
	       	fp4=fopen(file_nameus[num_file],"w");
	       	fp2=fopen(file_namerho[num_file],"w");
	       	num_file++;
		y=0;
		for(j=0;j<=ny+1;j++){
				fprintf(fp1,"%lf %lf\n",(ux[nx+1][j]),y);
				y+=dely;
				//printf("%lf\n",dely );
			}
			fclose(fp1);
			y=0;
		for(j=0;j<=ny+1;j++){
				fprintf(fp3,"%lf %lf\n",(ux[nx/2][j]),y);
				y+=dely;
			}
			fclose(fp3);
			y=0.0;
		for(j=0;j<=ny+1;j++)
		{
				fprintf(fp4,"%lf %lf\n",(ux[0][j]),y);
				y+=dely;
			}
			fclose(fp4);

		 	x=0.0;
		       for(i=0;i<=nx+1;i++){
			fprintf(fp2,"%lf %lf\n",x,rho[i][(ny+1)/2]);
				x=x+delx;
				}
		      fclose(fp2);
		// }
	   // Writting ux and y at different x at time=snapshot
	 //   if(itr==snapshot_time){
		// num_file1=0;
  // 		x=0.0;
  //   		i=0;
  // 	     do{
	 //       fp1=fopen(file_name5[num_file1],"w");

	 //       y=0.0;
	 //       for(j=0;j<=ny+1;j++){
	 //       fprintf(fp1,"%lf %lf\n",ux[i][j],y);
	 //       y=y+dely;
	 //       }
	 //       i=i+(int)x_step;
	 //       x=x+x_step;
	 //       fclose(fp1);
	 //       num_file1++;
  //      }while(x<=length);
  // }
}
	// printf("/n Generating GnuPlot .plt files for plotting\n");
	// for (int i = 0; i <= 10; ++i)
	// {
	// 	    fp5 =fopen(,"w");
	// fprintf(fp5,"set terminal svg enhanced size 1000 1000 fname \" Poiseuille Flow \" fsize 36\r\n");
	// fprintf(fp5,"set output \"poiseuilleFlow.svg\" \r\n ");
	// fprintf(fp5,"set title \"Plot of T Vs x\" \r\n");
	// fprintf(fp5,"set xlabel \"x\" \r\n ");
	// fprintf(fp5,"set ylabel \"T\" \r\n");
	// fprintf(fp5,"plot \"A11.dat\" using 1:2 title \"\" \r\n");
	// }



//End of time loop

// //------------------------------------------------------------------------------------------------------------
// 	// Writing density and x in file
// 	fp2=fopen("rho_x","w");
// 	x=0.0;
//         for(i=0;i<=nx+1;i++){
// 		fprintf(fp2,"%lf %lf\n",x,rho[i][(ny+1)/2]);
// 		x=x+delx;
//  		}
//        fclose(fp2);

// //-----------------writting ux and y at different x---------------------------
//        num_file=0;
//        x=0.0;
//        i=0;
//        do{
// 	       fp1=fopen(file_name1[num_file],"w");

// 	       y=0.0;
// 	       for(j=0;j<=ny+1;j++){
// 	       fprintf(fp1,"%lf %lf\n",ux[i][j],y);
// 	       y=y+dely;
// 	       }
// 	       i=i+(int)x_step;
// 	       x=x+x_step*delx;
// 	       fclose(fp1);
// 	       num_file++;
//        }while(x<=length);

// //--------------------- Writting uy and y at different x in file---------------------------------
//        num_file=0;
//        x=0.0;
//        i=0;
//        do{
//                fp1=fopen(file_name2[num_file],"w");
//                y=0.0;
// 	       for(j=0;j<=ny+1;j++){
// 	       fprintf(fp1,"%lf %30.28lf\n",y,uy[i][j]);
// 	       y=y+dely;
// 	       }
// 	       x=x+x_step*delx;
// 	       i=i+(int)x_step;
// 	       fclose(fp1);
// 	       num_file++;
//        }while(x<=length);

// //------------------Writing rho and x at different y-ordinate-----------------------------------------
//        num_file=0;
//        y=0.0;
//        j=0;
//        do{
//                fp1=fopen(file_name3[num_file],"w");
//                x=0;
// 	       for(i=0;i<=nx+1;i++){
// 	       fprintf(fp1,"%lf %lf\n",x,rho[i][j]);
// 	       x=x+delx;
// 	       }
// 	       j=j+(int)y_step;
// 	       y=y+y_step*dely;
// 	       fclose(fp1);
// 	       num_file++;
//        }while(y<=width);

// //-------------- Writting rho and y at different x---------------------------------------------------

//        num_file=0;
//        x=0.0;
//        i=0;
//        do{
//                fp1=fopen(file_name4[num_file],"w");
//                y=0;
// 	       for(j=0;j<=ny+1;j++){
// 	       fprintf(fp1,"%lf %30.28lf\n",y,rho[i][j]);
// 	       y=y+dely;
// 	       }
// 	       i=i+(int)x_step;
// 	       x=x+x_step*delx;
// 	       fclose(fp1);
// 	       num_file++;
//        }while(x<=length);

// //---------------------Writting stress and y at exit-----------------------------------------
// 	fp1=fopen("stress_y","w");
// 	fp2=fopen("ux_y_at_exit","w");
// 	y=0.0;

// 	for(j=0;j<=ny;j++){
// 		fprintf(fp2,"%lf %20.18lf\n",y,ux[nx+1][j]);
// 		dux[nx+1][j+1]=ux[nx+1][j]-ux[nx+1][j+1];
// 		shear_stress[nx+1][j+1]=rho[nx+1][j+1]*dneu*(dux[nx+1][j+1]/dely);
// 		y=y+dely;
// 		fprintf(fp1,"%20.18lf %lf\n",shear_stress[nx+1][j+1],y);
// 		}
// 		fclose(fp1);
// 		fclose(fp2);
   }
}
