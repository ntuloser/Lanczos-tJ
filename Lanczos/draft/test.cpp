#include<iostream>
#include <stdlib.h>
#include <ctime>
#include <cmath>
#include <stdio.h>
#include"config.h"

using namespace std;

#define PI 3.1415926535897932

#define SIZE 4
#define DELTA 6


const double t=1.;
const double J=0;//0.33;
const double U=1.;
const int num_e=SIZE*SIZE-DELTA;

////Variational parameter////
double D=0.119767;
double Mu=0.;
double g=0.754434;
/////////////////////////////



void update_aij(double ** a_ij, double kk[SIZE*SIZE][2]){
    for (int i=0; i< (SIZE*2-1); i++) {
        for (int j=0; j <(SIZE*2 -1 ); j++) {
            
            a_ij[i][j]=0;
            
            for (int m=0; m<SIZE*SIZE; m++) {
                double ek,Dk;
                ek = (-2.)*t*(cos(kk[m][0])+cos(kk[m][1]))-Mu;
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


void createInverse(double** c,double** b,int n)
{
	double a[n][n];
	for (int i=0;i<n;i=i+1)
		for (int j=0;j<n;j=j+1)
			a[i][j]=c[i][j];
	double Inv[n][n];
	for (int i=0;i<n;i=i+1){
		for (int j=0;j<n;j=j+1){
			if(i==j)
				Inv[i][j]=1;
			else
				Inv[i][j]=0;
		}
    }
    
    
	for (int i=0;i<n;i=i+1)
	{
		int k=i;
		while (abs(a[k][i])<0.000001&&k<n)
			k=k+1;
		if (k==n){
            cout<<"gg in deter: "<<determinant(c);
            cout<<"stop"<<endl;system("pause");
        }
		if (k!=i){ //a[i][i] is small
			for (int j=0;j<n;j=j+1)
			{
				double b=a[i][j];
				a[i][j]=a[k][j];
				a[k][j]=b;
				b=Inv[i][j];
				Inv[i][j]=Inv[k][j];
				Inv[k][j]=b;
			}
        }
        
        
		double diag=a[i][i];
		for (int j=0;j<n;j=j+1)
		{
			Inv[i][j]=Inv[i][j]/diag;
			a[i][j]=a[i][j]/diag;
		}
        
		for (int j=0;j<n;j=j+1)
		{
			double b=a[j][i];
			for (int l=0;l<n;l=l+1)
			{
				if (j!=i)
				{
					Inv[j][l]=Inv[j][l]-Inv[i][l]*b;
					a[j][l]=a[j][l]-a[i][l]*b;
				}
			}
		}
	}
	//cout<<determinant(c)*determinant(Inv)<<endl;
	for (int i=0;i<n;i=i+1)
		for (int j=0;j<n;j=j+1)
			b[i][j]=Inv[i][j];
}




int main(){
    srand((unsigned)time(0));
    
    int Variational_steps=1;
    double step_size=0.001;
    double E_tmp[4];
    
    double del_t = 0.001;
    double deriv;
    
    
    
    
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

        double  tot_O_D     =0;
        double  tot_O_Mu    =0;
        double  tot_O_g     =0;
        double  tot_O_DtimesE =0;
        double  tot_O_MutimesE =0;
        double  tot_O_gtimesE =0;
        
        
        double  tot_E2      =0;
        double  tot_doub2   =0;
        
        alpha.rand_init();
        alpha.printconfig();
        
        
        
        
        
		/////////////////////////////////
		//**the monte carlo procedure**//
		/////////////////////////////////
		int sample=30000;
		int interval=50;
		int warmup=30000;
        //		int takeInv=500;
		int totsteps=sample*interval+warmup;
        
		for (int steps=0;steps<=totsteps;steps++){

            //Random
            //Generating a new config from alpha !
            //
            int flipped=0;
            int idx =   rand()%num_e;     //choose electron
            int move =  rand()%4;
            int overbound=0;
            
            if (idx<num_e/2) {
                //printf("idx:%d move:%d up\n",idx,move);
                alpha.hop_up(&beta, idx,move,&flipped,&overbound);
            }
            else{
                //printf("idx:%d move:%d down\n",idx-num_e/2,move);
                alpha.hop_down(&beta, idx-num_e/2,move,&flipped,&overbound);
            }
            
            
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
                int num_d_a = alpha.num_doublon();
                
                if (std::abs(determinant(slater_a))<0.001) {
                    cout<<"small_deter!!!\n";
                }
                createInverse(slater_a,inv_slater_a,num_e/2);
                
                /////optimization//
                double  deter_origin = determinant(slater_a);
                
                D+=0.001;
                update_aij(a_ij,kk);
                alpha.set_slater(a_ij,slater_a);
                double  deter_D     =   determinant(slater_a);
                
                D-=0.001;
                Mu+=0.001;
                update_aij(a_ij,kk);
                alpha.set_slater(a_ij,slater_a);
                double  deter_Mu    =   determinant(slater_a);
                
                Mu-=0.001;
                update_aij(a_ij,kk);
                alpha.set_slater(a_ij,slater_a);
                
                double  deriv_D     =   (deter_D-deter_origin)/(0.001*deter_origin);
                double  deriv_Mu    =   (deter_Mu-deter_origin)/(0.001*deter_origin);
                double  deriv_g     =   ( pow(g,num_d_a) - pow((g+0.001),num_d_a) )/0.001;
                tot_O_D +=  deriv_D;
                tot_O_Mu+=  deriv_Mu;
                tot_O_g +=  deriv_g;
                //////////
                
                //////////
                //  U   //
                //////////
                tot_doub    += alpha.num_doublon();
                tot_E       += U * alpha.num_doublon();
                tot_O_DtimesE   +=  U *alpha.num_doublon()*deriv_D;
                tot_O_MutimesE  +=  U *alpha.num_doublon()*deriv_Mu;
                tot_O_gtimesE   +=  U *alpha.num_doublon()*deriv_g;
                
                tot_doub2   += pow(alpha.num_doublon(),2);
                double temp_E2      = U * alpha.num_doublon();
                
                
                //t
                for (int idx_1=0; idx_1<num_e; idx_1++) {
                    for (int move=0; move<4; move++) {
                        
                        double Ek=0;
                        overbound=0 ;
                        double ratio=0;
                        double Ek_square=0;
                        
                        if (idx_1<num_e/2) {//electron up
                            
                            alpha.hop_up(&beta, idx_1,move,&flipped,&overbound);
                            if (flipped==1) {
                                
                                beta.set_slater(a_ij,slater_b);
                                int num_d_b = beta.num_doublon();
                                
                                for (int i=0; i<num_e/2; i++) {
                                    ratio += slater_b[idx_1][i]*inv_slater_a[i][idx_1];
                                }
                                
                                if (not overbound) {
                                    Ek=-t*ratio *pow(g,num_d_b-num_d_a);
                                    //Ek=-t*determinant(slater_b)/determinant(slater_a) *pow(g,num_d_b-num_d_a);
                                    
                                    
                                    Ek_square=-t  ;
                                }
                                else{
                                    Ek=t*ratio *pow(g,num_d_b-num_d_a);
                                    //Ek=t*determinant(slater_b)/determinant(slater_a) *pow(g,num_d_b-num_d_a);
                                    
                                    
                                    Ek_square=t  ;
                                }
                                
                            }//end if flipped==1
                            else{
                                Ek=0;
                            }
                            
                        }//end electron up
                        
                        else{//idx_1 >=  num_e/2
                            //electron down
                            
                            
                            alpha.hop_down(&beta, (idx_1-num_e/2),move,&flipped,&overbound);
                            if (flipped==1) {
                                int num_d_b = beta.num_doublon();
                                beta.set_slater(a_ij,slater_b);
                                
                                for (int i=0; i<num_e/2; i++) {
                                    ratio += inv_slater_a[idx_1-num_e/2][i]*slater_b[i][idx_1-num_e/2];
                                }
                                
                                if (not overbound) {
                                    Ek=-t*ratio *pow(g,num_d_b-num_d_a);
                                    //Ek=-t*determinant(slater_b)/determinant(slater_a)  *pow(g,num_d_b-num_d_a);
                                    
                                    Ek_square=-t;
                                }
                                else{
                                    Ek=t*ratio *pow(g,num_d_b-num_d_a);
                                    //Ek=t*determinant(slater_b)/determinant(slater_a)  *pow(g,num_d_b-num_d_a);
                                    
                                    Ek_square=t;
                                }
                                
                                
                            }
                            else
                                Ek=0;
                            
                        }//end electron down
                        
                        tot_Esquare +=  Ek*U*beta.num_doublon(); //contri ~ -6
                        tot_Esquare +=  U*alpha.num_doublon()*Ek; //contri ~ -6.2
                        
                        tot_E       +=  Ek;
                        tot_O_MutimesE += Ek*deriv_Mu;
                        tot_O_DtimesE += Ek*deriv_D;
                        tot_O_gtimesE += Ek*deriv_g;
                        
                        temp_E2     +=  Ek;
                    }// end move forloop
                }//end idx forloop
                
                tot_E2 += pow(temp_E2,2);
                tot_Esquare += pow(U*alpha.num_doublon(),2);//contribution ~ 1
                
            }//end sampling-if.
            
            
            
            
            
        }//end of monte carlo loop;
        
        cout<<" D = "<<D<<" Mu = "<<Mu<<" g = "<<g<<"\n";
        double avg_E= tot_E/sample;
        double num_doub = double(tot_doub)/double(sample);
        
        double err_E=pow( (tot_E2/sample-pow(avg_E,2))/sample , 0.5)/pow(SIZE,2);
        double err_doub=pow( (tot_doub2/sample-pow(num_doub,2))/sample,0.5)/pow(SIZE,2);
        
        
        //cout<<"tot_E = "<<tot_E<<", tot_doub = "<<tot_doub<<endl;
        cout<<"sample = "<<sample<<"\n";
        cout<<"avg_E = "<<avg_E /SIZE/SIZE<<" err: "<<err_E<<endl;
        //cout<<"Esquare ="<<tot_Esquare /pow(SIZE,4)/sample<<endl;
        cout<<"E_variance = "<<pow((std::abs(pow(avg_E,2)-tot_Esquare/sample))/sample,0.5)/pow(SIZE,2)<<endl;
        cout<<"doublon number="<<num_doub/SIZE/SIZE<<" err: "<<err_doub<<endl;
        cout<<"acceptance ratio:" <<tot_accept/totsteps<<endl;
        
        
        double grad_D   =   2*(tot_O_DtimesE/sample - (tot_O_D/sample)*(tot_E/sample));
        double grad_Mu  =   2*(tot_O_MutimesE/sample - (tot_O_Mu/sample)*(tot_E/sample));
        double grad_g   =   2*(tot_O_gtimesE/sample - (tot_O_g/sample)*(tot_E/sample));

        cout<<"\n";
        cout<<"D = "<<D<<" , grad_D = "<<grad_D<<endl;
        cout<<"Mu = "<<Mu<<" , grad_Mu = "<<grad_Mu<<endl;
        cout<<"g  = "<<g <<" , grad_g  = "<<grad_g<<endl;
        cout<<"\n";
        
        
        step_size=0.01;
        D   -= step_size*grad_D;
        Mu  -= step_size*grad_Mu;
        g   -= step_size*grad_g;
        
            
            
        
        
        
    }//end of variational loop;
    
    
    
    return 0;
}


//tested: constuctor(SIZE,dope) , ran_init()
//tested: printconfig(), swap(site_idx,move)
//tested: hop_up(idx, move), hop down(idx, move)


/*
 printf("size of config: %lu\n",sizeof(config));
 config alpha(4,0);
 alpha.printconfig();
 printf("doublon number: %d\n",alpha.num_doublon());
 
 
 alpha.rand_init();
 alpha.printconfig();
 printf("doublon number: %d\n",alpha.num_doublon());
 
 
 config beta(4,0);
 beta=alpha.swap(1,1);
 beta.printconfig();
 printf("doublon number: %d\n",alpha.num_doublon());
 */



/*
 stp_size 0.0003 too small to move away.
 stp_size 0.003  too big to stay stable.
 sample 20000 --> err 0.0006
 sample 100000--> err 0.0003
 sample 10000 --> err 0.001
 
 */
