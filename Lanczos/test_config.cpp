#include<iostream>
#include<cmath>
#include<ctime>
#include"config.h"
#include"function.h"


using namespace std;

#define SIZE 4
#define PI 3.1415926535897932
#define DELTA 0

const double t=1.0;
const double U=1.0;



int main(){
    srand((unsigned)time(0));
    
	////Variational parameter////
	double D=1.;
	double Mu=0.;
    double g=1.;
    
    ///////
    
    int num_e=SIZE*SIZE-DELTA;
    double tot_E    = 0;
    double tot_doub = 0;
    double tot_accept = 0;
    double tot_det=0;
    
    
    
    double kk[SIZE*SIZE][2];    // 1-dim vector representation of kx,ky
    for (int idx =0; idx<SIZE*SIZE; idx++) {
        int i=      idx % SIZE;
        int j=      idx / SIZE;
        kk[idx][0]= PI*( 1.-2.0*i/SIZE);
        kk[idx][1]= PI*( (SIZE-1.0)/SIZE-2.0*j/SIZE);
    }
    
    double** a_ij = new double*[2*SIZE-1];
    for (int i=0; i<2*SIZE-1; i++) {
        a_ij[i] = new double[2*SIZE-1];
    }
    
    for (int idx=0; idx<(2*SIZE-1)*(2*SIZE-1); idx++) {
        int i=      idx % (2*SIZE-1);
        int j=      idx / (2*SIZE-1);
        
        a_ij[i][j]=0;
        
        for (int m=0; m<SIZE*SIZE; m++) {
            double ek,Dk;
            ek = (-2.)*t*(cos(kk[m][0])+cos(kk[m][1]))-Mu;
            Dk = D*(cos(kk[m][0])-cos(kk[m][1])) ;
            
            if (D>0.0001) {
                //sine term is cancel when summing over the BZ zone.
                a_ij[i][j] += Dk/(ek+pow(pow(ek,2)+pow(Dk,2),0.5) )*cos(kk[m][0]*(i-(SIZE-1))+kk[m][1]*(j-(SIZE-1)) )  / (SIZE*SIZE);
            }
            else
                cout<<"error in vanishing D"<<endl;
        }
    }
    
    ////////////////////////
    ///////////////////////
    //////////////////////
    /////////////////////
    ////////////////////
    ///////////////////
    
    config alpha(SIZE,DELTA);
    alpha.rand_init();
    alpha.printconfig();
    
    config beta(SIZE,DELTA);
    
    double** slater_a = new double*[num_e/2];
    for (int i=0; i<num_e/2; i++) {
        slater_a[i] = new double[num_e/2];
    }
    double** slater_b = new double*[num_e/2];
    for (int i=0; i<num_e/2; i++) {
        slater_b[i] = new double[num_e/2];
    }
    
    
    
    int interval = 200;
    int sample   = 20000;
    int warmup   = 20000;
    int tot_steps= warmup+ sample*interval;
    for (int steps=0; steps<tot_steps ; steps++) {
        
        alpha.copy_config_to(&beta);
        
        int idx=rand() % num_e;
        rand();
        int move=rand() % 4;
        int flipped=0;
        
        if (idx<num_e/2) {
            alpha.hop_up(&beta, idx,move, &flipped);
        }
        else{
            alpha.hop_down(&beta, (idx-num_e/2),move, &flipped);
        }
        
        
        alpha.set_slater( a_ij,slater_a);
        beta.set_slater( a_ij,slater_b);
        
        
        if (steps>warmup && steps%interval==0) {
            //cout<<"det: "<<pow( determinant(slater_b,num_e/2)/determinant(slater_a,num_e/2) ,2)<<endl;
            
            //tot_det+= pow( determinant(slater_b,num_e/2)/determinant(slater_a,num_e/2) ,2);
        }
        
        if (flipped==1) {
            
            //alpha.set_slater( a_ij,slater_a);
            //beta.set_slater( a_ij,slater_b);
            
            double p=0;
            while (p==0){
                p=(rand()+1)/(double)(RAND_MAX);
            }
            
            
            if (p < pow( determinant(slater_b,num_e/2)/determinant(slater_a,num_e/2) ,2) || abs(determinant(slater_a,num_e/2))<0.00001  ){
                
                beta.copy_config_to( &alpha);
                tot_accept+=1;
                
            }
        }
        
        
        if (steps>warmup && steps%interval==0) {
            
            int num_doub =    alpha.num_doublon();
            tot_doub+= num_doub;
            //cout<<tot_doub<<endl;
            
        }
        
    }
    
    cout<<"D = "<<D<<", Mu = "<<Mu<<endl;
    cout<<"tot_doub = "<<tot_doub<<endl;
    cout<<"num_doub = "<<tot_doub/sample/SIZE/SIZE<<endl;
    cout<<"acceptance ratio:" <<tot_accept/tot_steps<<endl;
    cout<<"tot_det = "<<tot_det/sample<<endl;
    
    return 0;
}
