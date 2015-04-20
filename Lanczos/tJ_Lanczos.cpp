#include<iostream>
#include <stdlib.h>
#include <ctime>
#include <cmath>
#include <stdio.h>
#include"config.h"
#include <vector>
//#include <omp.h>

using namespace std;

#define PI 3.1415926535897932

#define SIZE 4
#define DELTA 6//6


const double t=1.;
const double J=2.;//0.33;
const int num_e=SIZE*SIZE-DELTA;

const int Variational_steps=1;
const double del_t = 0.015;//0.03
const double differStep=0.008; //little larger then the variation of energy per site.


////Variational parameter////
double D=0.880081;//0.0200;        // the first  parameter
double Mu=-1.2273;//0.520803;     // the second parameter
double g=1;             // the third  parameter
/////////////////////////////
std::vector<double> Energy_log;
std::vector<double> D_log;
std::vector<double> Mu_log;




void update_aij(double ** a_ij, double kk[SIZE*SIZE][2]){
	for (int i=0; i< (SIZE*2-1); i++) {
		for (int j=0; j <(SIZE*2 -1 ); j++) {

			a_ij[i][j]=0;

			for (int m=0; m<SIZE*SIZE; m++) {
				double ek,Dk;
				ek = (-2.)*t*(cos(kk[m][0])+cos(kk[m][1])) - Mu;
				Dk = D*(cos(kk[m][0])-cos(kk[m][1])) ;

				if (D>0.0001) {
					//If D is not vanishing
					//sine term is cancel when summing over the BZ zone.
					a_ij[i][j] += Dk/(ek+pow(pow(ek,2)+pow(Dk,2),0.5) )*cos(kk[m][0]*(i-(SIZE-1))+kk[m][1]*(j-(SIZE-1)) )  / (SIZE*SIZE);
				}
				else{
					//cout<<"error in vanishing D"<<endl;

					a_ij[i][j]+=0.5*(1-std::abs(ek)/ek)*cos(kk[m][0]*(i-(SIZE-1))+kk[m][1]*(j-(SIZE-1)) )/(SIZE*SIZE);}

			}
			//Checked
			//printf("a[%i][%i]=%f\n",i,j,a_ij[i][j]);
		}
	}
}


double determinant(double** a)
{
	double x[num_e/2][num_e/2];
	for (int i=0;i< num_e/2;i=i+1)
		for (int j=0;j< num_e/2;j=j+1)
			x[i][j]=a[i][j];
	int k=0,n=0;
	while (k< num_e/2)
	{
		int l=k;
		while (std::abs(x[l][k])<0.00001)
		{
			l=l+1;
			if (l== num_e/2)
				return 0;
		}
		if (l!=k)
		{
			n=n+1;
			for (int i=0;i< num_e/2;i=i+1)
			{
				double b;
				b=x[k][i];
				x[k][i]=x[l][i];
				x[l][i]=b;
			}
		}

		for (int i=k+1;i< num_e/2;i=i+1)
		{
			double r=x[i][k]/x[k][k];
			for (int j=k;j< num_e/2;j=j+1)
			{
				x[i][j]=x[i][j]-r*x[k][j];
			}
		}
		k=k+1;
	}
	double det=1;
	for (int i=0;i< num_e/2;i=i+1)
	{det=det*x[i][i];
	}
	return det*pow(-1.0,n);
}

void Traversal(int level, int bound, config* config_level, double ** a_ij, double *** slater,const double deter_origin,const double deriv_D, const double deriv_Mu, double* ptr_tot_E, double* ptr_tot_O_DtimesE, double* ptr_tot_O_MutimesE, double* ptr_temp_EperSample);




int main(){
	srand((unsigned)time(0));
    
	double E_tmp[4];




	// The Brillouin Zone with periodic in x, antiperiodic in y
	double kk[SIZE*SIZE][2];    // 1-dim vector representation of kx,ky
	for (int idx =0; idx<SIZE*SIZE; idx++) {
		int i=      idx % SIZE;
		int j=      idx / SIZE;
		kk[idx][0]= PI*( 1.-2.0*i/SIZE);
		kk[idx][1]= PI*( (SIZE-1.0)/SIZE-2.0*j/SIZE);
		//Checked
		//printf("%f, %f\n",kk[idx][0],kk[idx][1]);
	}




	//// all possible a_ij
	//// consider all possible r_i - r_j
	//// This generate a 2S-1 X 2S-1 matrix
	//// We shall call the matrix by giving the displacement vector

	double ** a_ij = new double*[2*SIZE-1];
	for (int k=0; k < (2*SIZE-1); k++) {
		a_ij[k] =new double[2*SIZE-1];
	}
	//we set(update) a_ij in each iteration of the variation for-loop.

    
	//config config_level[0](SIZE,DELTA);
    /////////Construct a intermidiate config called config_level[1] !!
	//config config_level[1](SIZE,DELTA);
	/////////config_level[0].copy_config_to( &config_level[1] );
	//config config_level[2](SIZE,DELTA);
    //config config_level[3](SIZE,DELTA);
    config config_level[4]={config(SIZE,DELTA),config(SIZE,DELTA),config(SIZE,DELTA),config(SIZE,DELTA)};
    //config config1[4];
    
    int lvl=4;
    double *** slater;
    slater = new double**[lvl];
    for (int i=0; i<lvl; i++) {
        slater[i] = new double*[num_e/2];
        for (int k=0; k<num_e/2; k++) {
            slater[i][k] =new double[num_e/2];
        }
    }
    /*
	double ** slater[0];
	slater[0] = new double*[num_e/2];
	for (int k=0; k<num_e/2; k++) {
		slater[0][k] =new double[num_e/2];
	}
	double ** slater[1];
	slater[1] = new double*[num_e/2];
	for (int k=0; k<num_e/2; k++) {
		slater[1][k] =new double[num_e/2];
	}
	double ** slater[2];
	slater[2] = new double*[num_e/2];
	for (int k=0; k<num_e/2; k++) {
		slater[2][k] =new double[num_e/2];
	}
    double ** slater[3];
    slater[3] = new double*[num_e/2];
    for (int k=0; k<num_e/2; k++) {
        slater[3][k] =new double[num_e/2];
    }


	double ** inv_slater[0];
	inv_slater[0] = new double*[num_e/2];
	for (int k=0; k<num_e/2; k++) {
		inv_slater[0][k] =new double[num_e/2];
	}

	double ** inv_slater[1];
	inv_slater[1] = new double*[num_e/2];
	for (int k=0; k<num_e/2; k++) {
		inv_slater[1][k] =new double[num_e/2];
	}*/







	/////////////////////////////
	////Variational Procedure////
	/////////////////////////////

	for (int stp=0;stp<Variational_steps ; stp++) {

		// SET a_ij matrix
        update_aij(a_ij,kk);

		///Reset all variable///
		int     tot_doub    =0;
		double  tot_E       =0;
		double  tot_Esquare =0;
		double  tot_accept  =0;

        double  tot_O[3]={0,0,0};
		double  tot_O_DtimesE =0;
		double  tot_O_MutimesE =0;
		double  tot_O_gtimesE =0;
        double  tot_OtimesO[3][3];
        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++) {
                tot_OtimesO[i][j]=0;
            }
        }

        
		double  tot_E_sqr      =0;
		double  tot_doub2   =0;

		/////////////////////////////////
		//**the monte carlo procedure**//
		/////////////////////////////////
        
        config_level[0].rand_init_no_d();
        config_level[0].printconfig();
        
		int sample=50000;
		int interval=10;
		int warmup=3000;
		//		int takeInv=500;
		int totsteps=sample*interval+warmup;

		for (int steps=0;steps<=totsteps;steps++){

			//Random
			//Generating a new config from config_level[0] !
			//
			int flipped=0;
			int idx =   rand()%(SIZE*SIZE);     //choose electron
			int move =  rand()%4;
			int overbound=0;

            config_level[0].swap(&config_level[1], int(idx/SIZE),idx%SIZE, move,&flipped,&overbound);
            //int tJ = config_level[0].swap(&config_level[1], x,y,move,&flipped,&overbound);

			/*
			   cout<<"\n config_level[0]:\n";
			   config_level[0].printconfig();
			   cout<<"config_level[1]:\n";
			   config_level[1].printconfig();
			   cout<<"\n edited:"<<flipped<<"\n";
			 */



			///////////////////////////////////////
			///   Metropolis algorithm      ///////
			///////////////////////////////////////

			/// Given the random probability 0<p<1 ///

			double p=0;
			while (p==0){
				p=(rand()+1)/(double)(RAND_MAX);
			}

			config_level[0].set_slater(a_ij,slater[0]);
			config_level[1].set_slater(a_ij,slater[1]);

			if (flipped==1) {
				//cout<<"prob:"<<pow(determinant(slater[1],num_e/2)/determinant(slater[0],num_e/2),2)<<"\n";
				int num_d_a = config_level[0].num_doublon();
				int num_d_b = config_level[1].num_doublon();

				if (p<pow(determinant(slater[1])/determinant(slater[0]),2) *pow(g,2*(num_d_b-num_d_a))|| abs(determinant(slater[0]))<0.00001){

					//updated config.
					config_level[1].copy_config_to( &config_level[0] );
					tot_accept+=1;
				}
			}





			////////////////////////////////////////////////////////////////
			///////////  guys, Time to Sampling  ~~  /////////////
			//////////////////////////////////////////

			if (steps%interval==0 && steps>warmup){

				config_level[0].set_slater(a_ij,slater[0]);
				//if (std::abs(determinant(slater[0]))<0.000001)   cout<<"small_deter!!!\n";
				//createInverse(slater[0],inv_slater[0],num_e/2);
                
				/////optimization//
				double  deter_origin = determinant(slater[0])*config_level[0].SX;
                
				D+=differStep;
				update_aij(a_ij,kk);
				config_level[0].set_slater(a_ij,slater[0]);
				double  deter_D     =   determinant(slater[0])*config_level[0].SX;
                
				D-=differStep;
				Mu+=differStep;
				update_aij(a_ij,kk);
				config_level[0].set_slater(a_ij,slater[0]);
				double  deter_Mu    =   determinant(slater[0])*config_level[0].SX;
                
				Mu-=differStep;
				update_aij(a_ij,kk);
				config_level[0].set_slater(a_ij,slater[0]);
                
				double  deriv_D     =   (deter_D-deter_origin)/(differStep*deter_origin);
				double  deriv_Mu    =   (deter_Mu-deter_origin)/(differStep*deter_origin);
				//double  deriv_g     =   ( pow(g,num_d_a) - pow((g+0.001),num_d_a) )/0.001;
				tot_O[0] +=  deriv_D;
				tot_O[1] +=  deriv_Mu;
                tot_OtimesO[0][1]   += (deriv_D*deriv_Mu);
                tot_OtimesO[1][0]   += (deriv_D*deriv_Mu);
                tot_OtimesO[0][0]   += (deriv_D*deriv_D);
                tot_OtimesO[1][1]   += (deriv_Mu*deriv_Mu);
				//tot_O_g +=  deriv_g;
				//////////


                //////////////////
				//    t-J       //
                //////////////////
                
                //  config_level[0], config_level[1], a_ij, slater[1],deter_origin,
                //  tot_E_sqr, tot_E, tot_O_DtimesE, tot_O_MutimesE
                //  deriv_D, deriv_Mu
                double temp_EperSample=0;
                int N=1;

                Traversal(0,1, config_level, a_ij, slater, deter_origin, deriv_D,deriv_Mu, &tot_E, &tot_O_DtimesE, &tot_O_MutimesE, &temp_EperSample);
                
                tot_E_sqr += pow(temp_EperSample,2);
                
			}//End of Sampling.


		}//end of monte carlo loop;

        
        ///////////////////////////////////////
        /// Optimization is hard T______T  ////
        ///////////////////////////////////////

		cout<<" D = "<<D<<" Mu = "<<Mu<<" g = "<<g<<"\n";
		double avg_E= tot_E/sample;
        Energy_log.push_back(avg_E);
        D_log.push_back(D);
        Mu_log.push_back(Mu);
		double num_doub = double(tot_doub)/double(sample);

        double err_E=std::pow( double((tot_E_sqr/sample-pow(avg_E,2))/sample) , 0.5)/pow(double(SIZE),2);
		double err_doub=pow( (tot_doub2/sample-pow(num_doub,2))/sample,0.5)/pow(double(SIZE),2);


		//cout<<"tot_E = "<<tot_E<<", tot_doub = "<<tot_doub<<endl;
		cout<<"sample = "<<sample<<"\n";
		cout<<"avg_E = "<<avg_E /SIZE/SIZE<<" err: "<<err_E<<endl;
		//cout<<"Esquare ="<<tot_Esquare /pow(SIZE,4)/sample<<endl;
		//cout<<"E_variance = "<<pow((std::abs(pow(avg_E,2)-tot_Esquare/sample))/sample,0.5)/pow(SIZE,2)<<endl;
		cout<<"doublon number="<<num_doub/SIZE/SIZE<<" err: "<<err_doub<<endl;
		cout<<"acceptance ratio:" <<tot_accept/totsteps<<endl;

        // adjust the order, avoiding the round off error.//
		double grad_D   =   2*(tot_O_DtimesE -  (tot_O[0]*tot_E/sample)) /sample;
		double grad_Mu  =   2*(tot_O_MutimesE - (tot_O[1]*tot_E/sample)) /sample;
		//double grad_g   =   2*(tot_O_gtimesE/sample - (tot_O_g/sample)*(tot_E/sample));
        double Smatrix[2][2];
        for (int i=0; i<2; i++) {
            for (int j=0; j<2; j++) {
                Smatrix[i][j] = (tot_OtimesO[i][j]-tot_O[i]*tot_O[j]/sample)/sample;
                cout<<"i,j:"<<i<<" "<<j<<" --> "<<Smatrix[i][j]<<endl;
            }
        }
        double inverse_Smatrix[2][2];
        double deter_Smatrix = (Smatrix[0][0]*Smatrix[1][1]-Smatrix[0][1]*Smatrix[1][0]);
        inverse_Smatrix[0][0]= Smatrix[1][1]/deter_Smatrix;
        inverse_Smatrix[0][1]= -Smatrix[1][0]/deter_Smatrix;
        inverse_Smatrix[1][0]= -Smatrix[0][1]/deter_Smatrix;
        inverse_Smatrix[1][1]= Smatrix[0][0]/deter_Smatrix;
                                                      

		cout<<"\n";
		cout<<"D = "<<D<<" , grad_D = "<<grad_D<<endl;
		cout<<"Mu = "<<Mu<<" , grad_Mu = "<<grad_Mu<<endl;
		//cout<<"g  = "<<g <<" , grad_g  = "<<grad_g<<endl;
		cout<<"\n";

        double dD0= inverse_Smatrix[0][0]*grad_D+inverse_Smatrix[0][1]*grad_Mu;
        double dD1= inverse_Smatrix[1][0]*grad_D+inverse_Smatrix[1][1]*grad_Mu;
        
        cout<<"\n";
        cout<<"dD0 = "<<dD0<<endl;
        cout<<"dD1 = "<<dD1<<endl;
        //cout<<"g  = "<<g <<" , grad_g  = "<<grad_g<<endl;
        cout<<"\n";

        
        
        // The Variational method //
        
        //del_t define in the beginning
        ////////////////////////////////////
        // The Steepest Decend(SD) method //
        ////////////////////////////////////
		//D   -= del_t*grad_D*5;
		//Mu  -= del_t*grad_Mu*5;
		//g   -= del_t*grad_g;
        
        ////////////////////////////////////
        // Th SR method                 //
        /////////////////////////////////////

        D   -= del_t*dD0;
        Mu  -= del_t*dD1;
        //g   -= del_t*grad_g;


	}//end of variational loop;
    
    
    
    // Output the log //
    cout<<"log:"<<endl;
    for (int i=0; i<Energy_log.size(); i++) {
        cout<<"E: "<<Energy_log[i]<<",D: "<<D_log[i]<<",Mu: "<<Mu_log[i]<<endl;
    }

	return 0;
}




void Traversal(int lvl_now, int bound, config* config_level, double ** a_ij, double *** slater, const double deter_origin, const double deriv_D, const double deriv_Mu,double* ptr_tot_E, double* ptr_tot_O_DtimesE, double* ptr_tot_O_MutimesE, double* ptr_temp_EperSample){
    
    if (lvl_now==bound) {
        return;
    }
    
    for (int x=0; x<SIZE; x++) {
        for (int y=0; y<SIZE; y++) {
            for (int move=0; move<3; move++) {
                if (move==1)    continue;
                //cout<<"move:"<<move<<"x:"<<x<<"y:"<<y<<endl;
                
                double E_term[3]={0,0,0};
                double ratio=0;
                double Ek_square=0;
                int flipped=0;
                int overbound=0 ;
                
                int tJ = config_level[lvl_now].swap(&config_level[lvl_now+1], x,y,move,&flipped,&overbound);
                
                if (flipped==1) {
                    config_level[lvl_now+1].set_slater(a_ij,slater[1]);
                    
                    //for (int i=0; i<num_e/2; i++) {
                    //	ratio += slater[1][idx_1][i]*inv_slater[0][i][idx_1];
                    //}
                    
                    //Ek=-t*ratio *pow(g,num_d_b-num_d_a);
                    if (tJ==1) {// ele -- empty
                        if (not overbound) {
                            
                            E_term[0]  =  -t*determinant(slater[1])/deter_origin*config_level[lvl_now+1].SX;
                            {
                                //Traversal.
                                Traversal(lvl_now+1,bound, config_level, a_ij, slater, deter_origin, deriv_D,deriv_Mu, ptr_tot_E, ptr_tot_O_DtimesE, ptr_tot_O_MutimesE, ptr_temp_EperSample);

                            }
                            //function --> H_square
                        }
                        else{//overbound
                            E_term[0]  =  t*determinant(slater[1])/deter_origin*config_level[lvl_now+1].SX;
                            
                        }
                    }
                    else if(tJ==2){// eleup -- eledown
                        E_term[1]  = +J/2*determinant(slater[1])/deter_origin*config_level[lvl_now+1].SX;
                        {
                            //Traversal.
                        }
                        E_term[2]  = (-J/4*2);  // contribution from the S1Z S2Z term
                        // contribution from the n_up*n_down term
                        {
                            config_level[lvl_now].copy_config_to(&config_level[lvl_now+1]);
                            //Traversal.
                        }
                        //function --> H_square
                        
                    }
                    else cout<<"GG\n";
                    
                    //Ek_square=-t  ;
                }//end if flipped==1
                else{//if flipped==0
                    if (tJ==0) {
                        ;   //emtpy -- empty
                    }
                    else if(tJ==3){
                        ;   //ele -- ele
                    }
                    
                    //E=0;
                }
                
                
                
                if (lvl_now==0) {
                    
                    for (int i=0; i<3; i++) {
                        (*ptr_tot_E)       +=  E_term[i];
                        (*ptr_tot_O_DtimesE) += E_term[i]*deriv_D;
                        (*ptr_tot_O_MutimesE) += E_term[i]*deriv_Mu;
                        //tot_O_gtimesE += E*deriv_g;
                    }
                    
                    for (int i=0; i<3; i++) {
                        (*ptr_temp_EperSample)       +=  E_term[i];
                    }

                }
            }// endfor move
        }//endfor y
    }//endfor x
}


//  better data structure --> memory cost rather than speed cost
//  parallel computing in openMP

