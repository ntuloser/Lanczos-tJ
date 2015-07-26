#include<iostream>
#include <stdlib.h>
#include <ctime>
#include <cmath>
#include <stdio.h>
#include <vector>
//#include <omp.h>
#include <mpi.h>
#include <functional>
//#include <unordered_map>
#include "khash.h"
KHASH_MAP_INIT_INT64(kMap, double)
//sturct kh_kMap_t

#define SIZE 4
#define DELTA 6//6
#define LEV 3// calculation given to H^lev
#include"config.h"

using namespace std;

#define PI 3.1415926535897932
//int sss=4;
//int ddd=6;


const double t=1.;
const double J=2.;//0.33;
const int num_e=SIZE*SIZE-DELTA;

const int Variational_steps=1;
const double del_t = 0.015;//0.03
const double differStep=0.008; //little larger then the variation of energy per site.


////Variational parameter////
double D=  0.862516;//0.880081;        // the first  parameter
double Mu= -1.29518;//-1.2273;     // the second parameter
double g=1;             // the third  parameter
/////////////////////////////
std::vector<double> Energy_log;
std::vector<double> D_log;
std::vector<double> Mu_log;

int countaccept=0;
int counttotal=0;


void createInverse(double c[SIZE*SIZE/2][SIZE*SIZE/2],double b[SIZE*SIZE/2][SIZE*SIZE/2],int n){
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
    for (int i=0;i<n;i=i+1){
        int k=i;
        while (abs(a[k][i])<0.00001&&k<n)
            k=k+1;
        if (k==n)
        {cout<<"stop"<<endl;system("pause");}
        if (k!=i)
            for (int j=0;j<n;j=j+1)
            {
                double b=a[i][j];
                a[i][j]=a[k][j];
                a[k][j]=b;
                b=Inv[i][j];
                Inv[i][j]=Inv[k][j];
                Inv[k][j]=b;
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

void Traversal(int level, int bound, config* config_level, double ** a_ij, double *** slater,const double& deter_origin,const double& deriv_D, const double& deriv_Mu, double* ptr_tot_E, double* ptr_tot_O_DtimesE, double* ptr_tot_O_MutimesE, double* ptr_temp_EperSample,double *Energy_level, double * tot_E_power,  kh_kMap_t* khashMap);


long long tonumber(const config& alpha){
    long long number=0;
/*    for (int idx=0; idx<alpha.num_ele/2 ; idx++) {
        double tmp = 4.0*idx;
        number += alpha.electronsup[idx][0] * pow(4,(tmp+0.));
        number += alpha.electronsup[idx][1] * pow(4,(tmp+1.));
        number += alpha.electronsdown[idx][0]*pow(4,(tmp+2.));
        number += alpha.electronsdown[idx][1]*pow(4,(tmp+3.));
        
    }*/
    //cout<<" The number transform from the config :"<<number<<endl;
    
        int power[16]={1,4,16,64,256,1024,4096,16384,65536,262144,1048576,4194304,16777216,67108864,268435456,1073741824  };
        number += alpha.electronsup[0][0] * power[0];
        number += alpha.electronsup[0][1] * power[1];
        number += alpha.electronsdown[0][0]*power[2];
        number += alpha.electronsdown[0][1]*power[3];

        number += alpha.electronsup[1][0] * power[4];
        number += alpha.electronsup[1][1] * power[5];
        number += alpha.electronsdown[1][0]*power[6];
        number += alpha.electronsdown[1][1]*power[7];

        number += alpha.electronsup[2][0] * power[8];
        number += alpha.electronsup[2][1] * power[9];
        number += alpha.electronsdown[2][0]*power[10];
        number += alpha.electronsdown[2][1]*power[11];

        number += alpha.electronsup[3][0] * power[12];
        number += alpha.electronsup[3][1] * power[13];
        number += alpha.electronsdown[3][0]*power[14];
        number += (long long)alpha.electronsdown[3][1]* power[15];

        number += (long long)alpha.electronsup[4][0] * power[15]*power[1];
        number += (long long)alpha.electronsup[4][1] * power[15]*power[2];
        number += (long long)alpha.electronsdown[4][0]*power[15]*power[3];
        number += (long long)alpha.electronsdown[4][1]*power[15]*power[4];


	return number;
}



int main(int argc, char** argv){
    
    double buffer[3]={D,Mu,0};//D,Mu,g
    double ksum, Nsum;
    ///////////  Using MPI ////////
    MPI::Status status;
    MPI_Request request;
    MPI::Init(argc,argv);
    int numnodes = MPI::COMM_WORLD.Get_size();
    int mynode = MPI::COMM_WORLD.Get_rank();
    std::cout << "Process " << mynode<< " of " << numnodes << std::endl;
    
    //srand(  pow((unsigned)time(0),(mynode+1))  );
    srand((mynode+1)*(unsigned)time(0));
    
    
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
    
    
    int ret;
    khiter_t k;
    khash_t(kMap) *khashMap = kh_init(kMap);
    //unordered_map<long long,double,hash<long long> > dMap;
    
    /////////////////////////////
    ////Variational Procedure////
    /////////////////////////////
    
    for (int stp=0;stp<Variational_steps ; stp++) {
        
        
        ///Reset all variable///
        int     tot_doub    =0;
        double  tot_E       =0;
        double  tot_E_power[LEV];
        for (int i=0; i<LEV; i++) {
            tot_E_power[i]=0;
        }
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
        int sample=20000/numnodes;
        int interval=30;
        int warmup=3000;
        int totsteps=sample*interval+warmup;
        
        ///////////  Using OpenMP ////////
        //int num_threads=omp_get_max_threads() ;
        //////////////////////////////////
        
        //#pragma omp parallel for
        //for (int prl_seed=0; prl_seed<numnodes; prl_seed++) {
        
        //// all possible a_ij
        //// consider all possible r_i - r_j
        //// This generate a 2S-1 X 2S-1 matrix
        //// We shall call the matrix by giving the displacement vector
        
        double ** a_ij = new double*[2*SIZE-1];
        for (int k=0; k < (2*SIZE-1); k++) {
            a_ij[k] =new double[2*SIZE-1];
        }
        //we set(update) a_ij in each iteration of the variation for-loop.
        
        
        
        
        /*try {
         config* config_level = new config[LEV];
         }catch (bad_alloc xa) {
         cout << "Allocation Failure\n";
         return 1;
         }*/
        
        
        int lvl=LEV+1;
        config* config_level = new config[lvl];
        for (int i=0; i<lvl; i++) {
            ;
            //config_level[i]=config(SIZE,DELTA);
        }
        //config config_level[3]={config(SIZE,DELTA),config(SIZE,DELTA),config(SIZE,DELTA)};
        
        double *** slater;
        slater = new double**[lvl];
        for (int i=0; i<lvl; i++) {
            slater[i] = new double*[num_e/2];
            for (int k=0; k<num_e/2; k++) {
                slater[i][k] =new double[num_e/2];
            }
        }
        double *** inv_slater;
        inv_slater = new double**[lvl];
        for (int i=0; i<lvl; i++) {
            inv_slater[i] = new double*[num_e/2];
            for (int k=0; k<num_e/2; k++) {
                inv_slater[i][k] =new double[num_e/2];
            }
        }
        
        
        // SET a_ij matrix
        update_aij(a_ij,kk);
        
        
        config_level[0].rand_init_no_d();
        config_level[0].printconfig();
        
        
        //		int takeInv=500;
        
        for (int steps=0;steps<=totsteps;steps++){
            
            //Random
            //Generating a new config from config_level[0] !
            //
            int flipped=0;
            int idx =   rand()%(SIZE*SIZE);     //choose electron
            int move =  rand()%4;
            int overbound=0;
            
            config_level[0].swap(&config_level[1], int(idx/SIZE),idx%SIZE, move,&flipped,&overbound);
            
            ///////////////////////////////////////
            ///   Metropolis algorithm      ///////
            ///////////////////////////////////////
            
            /// Given the random probability 0<p<1 ///
            
            double p=0;
            while (p==0){
                p=(rand()+1)/(double)(RAND_MAX);
            }
            
            
            if (flipped==1) {
                int num_d_a = config_level[0].num_doublon();
                int num_d_b = config_level[1].num_doublon();
                
                counttotal+=2;
                
                double deter_0;
                long long number=tonumber(config_level[0]);
                k = kh_get(kMap, khashMap, number);
                if ( k == kh_end(khashMap) ){
                    config_level[0].set_slater(a_ij,slater[0]);
                    deter_0 = determinant(slater[0]);
                    k = kh_put(kMap, khashMap, number, &ret);
                    kh_value(khashMap, k) = deter_0;
                }
                else{
                    countaccept+=1;
                    deter_0 = kh_value(khashMap, k);
                }
                
                double deter_1;
                number =tonumber(config_level[1]);
                k = kh_get(kMap, khashMap, number);
                if ( k == kh_end(khashMap) ){
                    config_level[1].set_slater(a_ij,slater[1]);
                    deter_1 = determinant(slater[1]);
                    k = kh_put(kMap, khashMap, number, &ret);
                    kh_value(khashMap, k) = deter_1;
                }
                else{
                    countaccept+=1;
                    deter_1 = kh_value(khashMap, k);
                }
                
                if (p<pow(deter_1/deter_0,2) *pow(g,2*(num_d_b-num_d_a))|| abs(deter_0)<0.00001){
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
                double  deter_origin = determinant(slater[0])*config_level[0].Sign;
                
                D+=differStep;
                update_aij(a_ij,kk);
                config_level[0].set_slater(a_ij,slater[0]);
                double  deter_D     =   determinant(slater[0])*config_level[0].Sign;
                
                D-=differStep;
                Mu+=differStep;
                update_aij(a_ij,kk);
                config_level[0].set_slater(a_ij,slater[0]);
                double  deter_Mu    =   determinant(slater[0])*config_level[0].Sign;
                
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
                /////////// ///////
                
                //  config_level[0], config_level[1], a_ij, slater[1],deter_origin,
                //  tot_E_sqr, tot_E, tot_O_DtimesE, tot_O_MutimesE
                //  deriv_D, deriv_Mu
                double temp_EperSample=0;
                int N=1;
                double Energy_level[LEV];
                for (int i=0; i<LEV; i++) {
                    Energy_level[i]=0;
                }
                
                Traversal(0,LEV, config_level, a_ij, slater, deter_origin, deriv_D, deriv_Mu, &tot_E, &tot_O_DtimesE, &tot_O_MutimesE, &temp_EperSample,Energy_level,tot_E_power, khashMap);
                //get the value of tot_E, tot_O_DtimesE, tot_O_MutimesE, temp_EperSample, Energy_level, tot_E_power
                
                tot_E_sqr += pow(temp_EperSample,2);
                
            }//End of Sampling.
            
        }//end of monte carlo loop
        
        //};  //end of parallel #openmp
        //end of parallel
        
        //
        
        
        
        double SUM_tot_E=0;
        double SUM_tot_O_DtimesE, SUM_tot_O_MutimesE, SUM_tot_E_sqr;
        double SUM_tot_E_power[LEV];
        for (int i=0; i <LEV; i++) {
            SUM_tot_E_power[i]=0;
        }
        double SUM_tot_O[3]={0,0,0};
        double SUM_tot_OO[4]={0,0,0,0};
        double local_tot_OO[4]={tot_OtimesO[0][0],tot_OtimesO[0][1],tot_OtimesO[1][0],tot_OtimesO[1][1]};
        MPI_Reduce(&tot_E,&SUM_tot_E,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Reduce(&tot_O_DtimesE,&SUM_tot_O_DtimesE,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Reduce(&tot_O_MutimesE,&SUM_tot_O_MutimesE,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Reduce(&tot_E_sqr,&SUM_tot_E_sqr,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Reduce(tot_E_power,SUM_tot_E_power,LEV,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Reduce(tot_O,SUM_tot_O,3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Reduce(local_tot_OO,SUM_tot_OO,4,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        
        
        
        if (mynode==0) {//This is the master branch
            
            sample=sample*numnodes;
            
            tot_E           =SUM_tot_E;
            tot_O_DtimesE   =SUM_tot_O_DtimesE;
            tot_O_MutimesE  =SUM_tot_O_MutimesE;
            tot_E_sqr       =SUM_tot_E_sqr;
            for (int i=0; i<LEV; i++) {
                tot_E_power[i]     =SUM_tot_E_power[i];
            }
            for (int i=0; i<3; i++) {
                tot_O[i]     =SUM_tot_O[i];
            }
            tot_OtimesO[0][0]=SUM_tot_OO[0];
            tot_OtimesO[1][0]=SUM_tot_OO[1];
            tot_OtimesO[0][1]=SUM_tot_OO[2];
            tot_OtimesO[1][1]=SUM_tot_OO[3];
            
            
            
            
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
            for (int i=0; i<LEV; i++) {
                cout<<"Epower"<<i<<" = "<<tot_E_power[i]/pow(SIZE,double((i+1)*2))/sample  <<endl;
            }
            
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
            ///////////////////////////////////////
            /// Optimization is hard T______T  ////
            ///////////////////////////////////////
            
            
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
            
            ///////////////////////////
            //  Preparing Broadcast ///
            ///////////////////////////
            buffer[0]=D;
            buffer[1]=Mu;
        }
        
        //////////////////
        //  Broadcast ///
        ////////////////
        MPI_Bcast(buffer,3,MPI_DOUBLE,0,MPI_COMM_WORLD);
        //int MPI_Bcast( void *buffer, int count, MPI_Datatype datatype, int root,  MPI_Comm comm )
        
        
        D  = buffer[0];
        Mu = buffer[1];
        
        
    }//end of variational loop;
    
    
    if (mynode==0) {
        // Output the log //
        cout<<"log:"<<endl;
        for (int i=0; i<Energy_log.size(); i++) {
            cout<<"E: "<<Energy_log[i]<<",D: "<<D_log[i]<<",Mu: "<<Mu_log[i]<<endl;
        }
        cout<<"hash_map_efficiency:"<<endl;
        cout<<"total = "<<counttotal<< "  accept = "<<countaccept<<"ratio"<<double(countaccept)/counttotal<<endl;
    }
    MPI::Finalize();
    
    return 0;
}





void Traversal(int lvl_now, int bound, config* config_level, double ** a_ij, double *** slater, const double& deter_origin, const double& deriv_D, const double& deriv_Mu,double* ptr_tot_E, double* ptr_tot_O_DtimesE, double* ptr_tot_O_MutimesE, double* ptr_temp_EperSample, double *Energy_level, double * tot_E_power, kh_kMap_t* khashMap){
    
    
    double deter[ (bound+1) ];
    khiter_t k;
    int ret;
    
    
    if (lvl_now==bound) {
        return;
    }
    
    for (int x=0; x<SIZE; x++) {
        for (int y=0; y<SIZE; y++) {
            for (int move=0; move<3; move++) {
                if (move==1)    continue;
                //cout<<"move:"<<move<<"x:"<<x<<"y:"<<y<<endl;
                
                
                double ratio=0;
                int flipped=0;
                int overbound=0 ;
                
                int tJ = config_level[lvl_now].swap(&config_level[lvl_now+1], x,y,move,&flipped,&overbound);
                
                if (flipped==1) {
                    counttotal+=1;
                    
                    long long number =tonumber(config_level[lvl_now+1]);
                    k = kh_get(kMap, khashMap, number);
                    if ( k == kh_end(khashMap) ){
                        config_level[lvl_now+1].set_slater(a_ij,slater[lvl_now+1]);
                        deter[lvl_now+1] = determinant(slater[lvl_now+1]);
                        k = kh_put(kMap, khashMap, number, &ret);
                        kh_value(khashMap, k) = deter[lvl_now+1];
                    }
                    else{
                        countaccept+=1;
                        deter[lvl_now+1] = kh_value(khashMap, k);
                    }
                    

                    
                    
                    //for (int i=0; i<num_e/2; i++) {
                    //	ratio += slater[1][idx_1][i]*inv_slater[0][i][idx_1];
                    //}
                    
                    //Ek=-t*ratio *pow(g,num_d_b-num_d_a);
                    if (tJ==1) {// ele -- empty
                        if (not overbound) {
                            //Traversal.
                            //New config generate
                            //Energy get.
                            Energy_level[lvl_now]  =  -t;
                            
                            
                            Traversal(lvl_now+1,bound, config_level, a_ij, slater, deter_origin, deriv_D,deriv_Mu, ptr_tot_E, ptr_tot_O_DtimesE, ptr_tot_O_MutimesE, ptr_temp_EperSample,Energy_level,tot_E_power, khashMap);
                            
                            if (lvl_now==0) {

                                Energy_level[lvl_now]*=(deter[lvl_now+1]/deter_origin*config_level[lvl_now+1].Sign);
                                //Energy_level[lvl_now]*=(determinant(slater[lvl_now+1])/deter_origin*config_level[lvl_now+1].Sign);
                                (*ptr_tot_E)       +=  Energy_level[0];
                                (*ptr_tot_O_DtimesE) += Energy_level[0]*deriv_D;
                                (*ptr_tot_O_MutimesE) += Energy_level[0]*deriv_Mu;
                                //tot_O_gtimesE += E*deriv_g;
                                (*ptr_temp_EperSample)       +=  Energy_level[0];
                            }
                            else{
                                double temp=1.0;
                                for (int i=0; i<=lvl_now; i++) {
                                    temp*=Energy_level[i];
                                }
                                temp*=(deter[lvl_now+1]/deter_origin*config_level[lvl_now+1].Sign);
                                //temp*=(determinant(slater[lvl_now+1])/deter_origin*config_level[lvl_now+1].Sign);
                                tot_E_power[lvl_now] += temp;
                            }
                            
                        }
                        else{
                            //In the case of overbound
                            //Energy_level[lvl_now]  =  t*determinant(slater[lvl_now+1])/deter_origin*config_level[lvl_now+1].Sign;
                            Energy_level[lvl_now]  =  t;
                            
                            
                            Traversal(lvl_now+1,bound, config_level, a_ij, slater, deter_origin, deriv_D,deriv_Mu, ptr_tot_E, ptr_tot_O_DtimesE, ptr_tot_O_MutimesE, ptr_temp_EperSample,Energy_level,tot_E_power, khashMap);
                            
                            if (lvl_now==0) {
                                Energy_level[lvl_now]*=(deter[lvl_now+1]/deter_origin*config_level[lvl_now+1].Sign);
                                //Energy_level[lvl_now]*=(determinant(slater[lvl_now+1])/deter_origin*config_level[lvl_now+1].Sign);
                                
                                (*ptr_tot_E)       +=  Energy_level[0];
                                (*ptr_tot_O_DtimesE) += Energy_level[0]*deriv_D;
                                (*ptr_tot_O_MutimesE) += Energy_level[0]*deriv_Mu;
                                //tot_O_gtimesE += E*deriv_g;
                                (*ptr_temp_EperSample)       +=  Energy_level[0];
                            }
                            else{
                                double temp=1.0;
                                for (int i=0; i<=lvl_now; i++) {
                                    temp*=Energy_level[i];
                                }
                                temp*=(deter[lvl_now+1]/deter_origin*config_level[lvl_now+1].Sign);
                                //temp*=(determinant(slater[lvl_now+1])/deter_origin*config_level[lvl_now+1].Sign);
                                
                                tot_E_power[lvl_now] += temp;
                            }
                            
                            
                            
                        }
                    }
                    
                    
                    
                    else if(tJ==2){// eleup -- eledown
                        
                        {// contri from the superexchange.
                            //Energy_level[lvl_now]  = +J/2*determinant(slater[lvl_now+1])/deter_origin*config_level[lvl_now+1].Sign;
                            Energy_level[lvl_now]  = +J/2;
                            
                            Traversal(lvl_now+1,bound, config_level, a_ij, slater, deter_origin, deriv_D,deriv_Mu, ptr_tot_E, ptr_tot_O_DtimesE, ptr_tot_O_MutimesE, ptr_temp_EperSample,Energy_level,tot_E_power, khashMap);
                            
                            
                            if (lvl_now==0) {
                                Energy_level[lvl_now]*=(deter[lvl_now+1]/deter_origin*config_level[lvl_now+1].Sign);
                                //Energy_level[lvl_now]*=(determinant(slater[lvl_now+1])/deter_origin*config_level[lvl_now+1].Sign);
                                
                                (*ptr_tot_E)       +=  Energy_level[0];
                                (*ptr_tot_O_DtimesE) += Energy_level[0]*deriv_D;
                                (*ptr_tot_O_MutimesE) += Energy_level[0]*deriv_Mu;
                                //tot_O_gtimesE += E*deriv_g;
                                (*ptr_temp_EperSample)       +=  Energy_level[0];
                            }
                            else{
                                double temp=1.0;
                                for (int i=0; i<=lvl_now; i++) {
                                    temp*=Energy_level[i];
                                }
                                temp*=(deter[lvl_now+1]/deter_origin*config_level[lvl_now+1].Sign);
                                //temp*=(determinant(slater[lvl_now+1])/deter_origin*config_level[lvl_now+1].Sign);
                                
                                tot_E_power[lvl_now] += temp;
                            }
                            
                            
                        }
                        
                        // contribution from the n_up*n_down term
                        {// giving the identical configuration in the next level.
                            Energy_level[lvl_now]  = (-J/4*2);  // contribution from the S1Z S2Z term
                            
                            config_level[lvl_now].copy_config_to(&config_level[lvl_now+1]);
                            config_level[lvl_now+1].set_slater(a_ij,slater[lvl_now+1]);
                            
                            counttotal+=1;
                            number= tonumber(config_level[lvl_now+1]);
                            k = kh_get(kMap, khashMap, number);
                            if ( k == kh_end(khashMap) ){
                                config_level[lvl_now+1].set_slater(a_ij,slater[lvl_now+1]);
                                deter[lvl_now+1] = determinant(slater[lvl_now+1]);
                                k = kh_put(kMap, khashMap, number, &ret);
                                kh_value(khashMap, k) = deter[lvl_now+1];
                            }
                            else{
                                countaccept+=1;
                                deter[lvl_now+1] = kh_value(khashMap, k);
                            }
                            
                            

                            

                            
                            
                            Traversal(lvl_now+1,bound, config_level, a_ij, slater, deter_origin, deriv_D,deriv_Mu, ptr_tot_E, ptr_tot_O_DtimesE, ptr_tot_O_MutimesE, ptr_temp_EperSample,Energy_level,tot_E_power, khashMap);
                            
                            if (lvl_now==0) {
                                
                                (*ptr_tot_E)       +=  Energy_level[0];
                                (*ptr_tot_O_DtimesE) += Energy_level[0]*deriv_D;
                                (*ptr_tot_O_MutimesE) += Energy_level[0]*deriv_Mu;
                                //tot_O_gtimesE += E*deriv_g;
                                (*ptr_temp_EperSample)       +=  Energy_level[0];
                            }
                            
                            else{
                                double temp=1.0;
                                for (int i=0; i<=lvl_now; i++) {
                                    temp*=Energy_level[i];
                                }
                                temp*=(deter[lvl_now+1]/deter_origin*config_level[lvl_now+1].Sign);
                                //temp*=(determinant(slater[lvl_now+1])/deter_origin*config_level[lvl_now+1].Sign);
                                
                                tot_E_power[lvl_now] += temp;
                            }
                            
                            
                            
                            
                        }
                    }
                    
                    else cout<<"GG\n";
                    
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
                
                
            }// endfor move
        }//endfor y
    }//endfor x
}


//  better data structure --> memory cost rather than speed cost
//  parallel computing in openMP
//
//  complexity is about
// ( Interval(30) + 32^LEV(5) )*SAMPLE( 30000 )*order(N^3)

