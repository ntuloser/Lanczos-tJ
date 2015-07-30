#include<iostream>
#include <stdlib.h>
#include <ctime>
#include <cmath>
#include <stdio.h>
#include"config.h"
#include <vector>

using namespace std;

#define PI 3.1415926535897932

#define SIZE 4
#define DELTA 6//6


const double t=1.;
const double J=2.;//0.33;
const int num_e=SIZE*SIZE-DELTA;

const int Variational_steps=50;
const double del_t = 0.015;//0.03
const double differStep=0.008; //little larger then the variation of energy per site.


////Variational parameter////
double D=1.23848;//0.0200;        // the first  parameter
double Mu=-0.348159;//0.520803;     // the second parameter
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

    
	config alpha(SIZE,DELTA);
	//Construct a intermidiate config called beta !!
	config beta(SIZE,DELTA);
	//alpha.copy_config_to( &beta );
	config gamma(SIZE,DELTA);
	//...

	double ** slater_a;
	slater_a = new double*[num_e/2];
	for (int k=0; k<num_e/2; k++) {
		slater_a[k] =new double[num_e/2];
	}

	double ** slater_b;
	slater_b = new double*[num_e/2];
	for (int k=0; k<num_e/2; k++) {
		slater_b[k] =new double[num_e/2];
	}


	double ** slater_g;
	slater_g = new double*[num_e/2];
	for (int k=0; k<num_e/2; k++) {
		slater_g[k] =new double[num_e/2];
	}


	double ** inv_slater_a;
	inv_slater_a = new double*[num_e/2];
	for (int k=0; k<num_e/2; k++) {
		inv_slater_a[k] =new double[num_e/2];
	}

	double ** inv_slater_b;
	inv_slater_b = new double*[num_e/2];
	for (int k=0; k<num_e/2; k++) {
		inv_slater_b[k] =new double[num_e/2];
	}







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

        
		double  tot_E2      =0;
		double  tot_doub2   =0;

		/////////////////////////////////
		//**the monte carlo procedure**//
		/////////////////////////////////
        
        alpha.rand_init_no_d();
        alpha.printconfig();
        
		int sample=50000;
		int interval=10;
		int warmup=30000;
		//		int takeInv=500;
		int totsteps=sample*interval+warmup;

		for (int steps=0;steps<=totsteps;steps++){

			//Random
			//Generating a new config from alpha !
			//
			int flipped=0;
			int idx =   rand()%(SIZE*SIZE);     //choose electron
			int move =  rand()%4;
			int overbound=0;

			 alpha.swap(&beta, int(idx/SIZE),idx%SIZE, move,&flipped,&overbound);
            //int tJ = alpha.swap(&beta, x,y,move,&flipped,&overbound);

			/*
			   cout<<"\n alpha:\n";
			   alpha.printconfig();
			   cout<<"beta:\n";
			   beta.printconfig();
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

			alpha.set_slater(a_ij,slater_a);
			beta.set_slater(a_ij,slater_b);

			if (flipped==1) {
				//cout<<"prob:"<<pow(determinant(slater_b,num_e/2)/determinant(slater_a,num_e/2),2)<<"\n";
				int num_d_a = alpha.num_doublon();
				int num_d_b = beta.num_doublon();

				if (p<pow(determinant(slater_b)/determinant(slater_a),2) *pow(g,2*(num_d_b-num_d_a))|| abs(determinant(slater_a))<0.00001){

					//updated config.
					beta.copy_config_to( &alpha );
					tot_accept+=1;
				}
			}





			////////////////////////////////////////////////////////////////
			///////////  guys, Time to Sampling  ~~  /////////////
			//////////////////////////////////////////

			if (steps%interval==0 && steps>warmup){

				alpha.set_slater(a_ij,slater_a);
				//if (std::abs(determinant(slater_a))<0.000001)   cout<<"small_deter!!!\n";
				//createInverse(slater_a,inv_slater_a,num_e/2);
                
				/////optimization//
				double  deter_origin = determinant(slater_a);
                
				D+=differStep;
				update_aij(a_ij,kk);
				alpha.set_slater(a_ij,slater_a);
				double  deter_D     =   determinant(slater_a);
                
				D-=differStep;
				Mu+=differStep;
				update_aij(a_ij,kk);
				alpha.set_slater(a_ij,slater_a);
				double  deter_Mu    =   determinant(slater_a);
                
				Mu-=differStep;
				update_aij(a_ij,kk);
				alpha.set_slater(a_ij,slater_a);
                
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
                double temp_EperSample=0;

				for (int x=0; x<SIZE; x++) {
					for (int y=0; y<SIZE; y++) {
						for (int move=0; move<3; move++) {
                            if (move==1)    continue;
                            //cout<<"move:"<<move<<"x:"<<x<<"y:"<<y<<endl;
                            

							double E_term[4]={0,0,0,0};
							double ratio=0;
							double Ek_square=0;
                            flipped=0;
                            overbound=0 ;


							//int site_idx=(x*SIZE+y);

							int tJ = alpha.swap(&beta, x,y,move,&flipped,&overbound);
							
                            if (flipped==1) {
								beta.set_slater(a_ij,slater_b);

								//for (int i=0; i<num_e/2; i++) {
								//	ratio += slater_b[idx_1][i]*inv_slater_a[i][idx_1];
								//}

								if (not overbound) {
									//Ek=-t*ratio *pow(g,num_d_b-num_d_a);
                                    if (tJ==1) {
                                        E_term[0]  =  -t*determinant(slater_b)/deter_origin;
                                    }
                                    else if(tJ==2){
                                        E_term[1]  = -J/2*determinant(slater_b)/deter_origin;// exchange
                                        E_term[2]  = (-J/4*1);  // contribution from the S1Z S2Z term
                                        E_term[3]  = (-J/4*1);  // contribution from the n_up*n_down term
                                    }
                                    else cout<<"GG\n";


									//Ek_square=-t  ;
								}
								else{
									//Ek=t*ratio *pow(g,num_d_b-num_d_a);
                                    if (tJ==1) {
                                        E_term[0] =  t*determinant(slater_b)/deter_origin ;
                                    }
                                    else if(tJ==2){
                                        E_term[1] =  -J/2*determinant(slater_b)/deter_origin ;// exchange
                                        E_term[2] =  (-J/4*1);// contri from the S1Z S2Z term
                                        E_term[3] =  (-J/4*1);  // contribution from the n_up*n_down term

                                    }
                                    else cout<<"GG\n";

									//Ek_square=t  ;
								}

							}//end if flipped==1
							else{//if flipped==0
                                if (tJ==0) {
                                    ;
                                }
                                else if(tJ==3){
                                    ;
                                }
                                
								//E=0;
							}


							//tot_Esquare +=  Ek*U*beta.num_doublon(); //contri ~ -6
							//tot_Esquare +=  U*alpha.num_doublon()*Ek; //contri ~ -6.2

                            for (int i=0; i<4; i++) {
                                tot_E       +=  E_term[i];
                                tot_O_DtimesE += E_term[i]*deriv_D;
                                tot_O_MutimesE += E_term[i]*deriv_Mu;
                                //tot_O_gtimesE += E*deriv_g;
                            }
                            
                            for (int i=0; i<4; i++) {
                                temp_EperSample       +=  E_term[i];
                            }
                            
                            
						}// endfor move

						//tot_Esquare += pow(U*alpha.num_doublon(),2);//contribution ~ 1


					}//endfor y



				}//endfor x
                tot_E2 += pow(temp_EperSample,2);

			}


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

        double err_E=std::pow( double((tot_E2/sample-pow(avg_E,2))/sample) , 0.5)/pow(double(SIZE),2);
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
		D   -= del_t*grad_D*5;
		Mu  -= del_t*grad_Mu*5;
		//g   -= del_t*grad_g;
        
        ////////////////////////////////////
        // Th SR method                 //
        /////////////////////////////////////

        //D   -= del_t*dD0;
        //Mu  -= del_t*dD1;
        //g   -= del_t*grad_g;
        
        
        



	}//end of variational loop;

    cout<<"log:"<<endl;
    for (int i=0; i<Energy_log.size(); i++) {
        cout<<"E: "<<Energy_log[i]<<",D: "<<D_log[i]<<",Mu: "<<Mu_log[i]<<endl;
        
    }

	return 0;
}
//term 1 ~ -0.7, term2 ~ 0.15, term3 ~ -0.35 , total ~-0.9

