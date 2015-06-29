#include<iostream>
#include<ctime>
#include <stdlib.h>
#include <unordered_map>
#include <cmath>


#define SIZE 4
#define DELTA 6//6
#define LEV 2// calculation given to H^lev

#include "config.h"

#define PI 3.1415926535897932


using namespace std;

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

int lvl = LEV+1;

long long tonumber(const config& alpha){
    long long number=0;
    for (int idx=0; idx<alpha.num_ele/2 ; idx++) {
        double tmp = 4.0*idx;
        number += alpha.electronsup[idx][0] * pow(4,(tmp+0.));
        number += alpha.electronsup[idx][1] * pow(4,(tmp+1.));
        number += alpha.electronsdown[idx][0]*pow(4,(tmp+2.));
        number += alpha.electronsdown[idx][1]*pow(4,(tmp+3.));
        
    }
    //cout<<" The number transform from the config :"<<number<<endl;
    return number;
}

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

int main(){
    srand((unsigned)time(0));
    
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

    double ** a_ij = new double*[2*SIZE-1];
    for (int k=0; k < (2*SIZE-1); k++) {
        a_ij[k] =new double[2*SIZE-1];
    }

    
    
    double *** slater;
    slater = new double**[lvl];
    for (int i=0; i<lvl; i++) {
        slater[i] = new double*[num_e/2];
        for (int k=0; k<num_e/2; k++) {
            slater[i][k] =new double[num_e/2];
        }
    }
    unordered_map<long long,double,hash<long long> > dMap;

    
    

    config alpha;
    config beta;
    alpha.rand_init_no_d();
    //alpha.printconfig();
    alpha.set_slater(a_ij,slater[0]);

    std::clock_t start;

    double duration_rand_init_no_d;
    double duration_set_slater;
    double duration_determinant;
    double dur_swap;
    double dur_copy;
    double dur_mapFindInsert;
    
    
    
    start = std::clock();
    for (int i=0; i<1000000; i++) {
        alpha.rand_init_no_d();
    //    alpha.set_slater(a_ij,slater[0]);
//        determinant(slater[0]);
    }
    duration_rand_init_no_d = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    
    
    
    start = std::clock();
    for (int i=0; i<1000000; i++) {
        alpha.rand_init_no_d();
        alpha.set_slater(a_ij,slater[0]);
        //        determinant(slater[0]);
    }
    duration_set_slater = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    duration_set_slater -= duration_rand_init_no_d;
    
    
    start = std::clock();
    for (int i=0; i<1000000; i++) {
        alpha.rand_init_no_d();
        alpha.set_slater(a_ij,slater[0]);
        determinant(slater[0]);
    }
    duration_determinant = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    duration_determinant-= duration_rand_init_no_d;
    duration_determinant-= duration_set_slater;
    
    
    
    
    
    
    start = std::clock();
    for (int i=0; i<1000000; i++) {
        
        int flipped=0;
        int idx =   rand()%(SIZE*SIZE);     //choose electron
        int move =  rand()%4;
        int overbound=0;
        
        alpha.swap(&beta, int(idx/SIZE),idx%SIZE, move,&flipped,&overbound);
        
        //beta.copy_config_to(&alpha);
        
    }
    dur_swap = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    
    
    
    start = std::clock();
    for (int i=0; i<1000000; i++) {
        
        int flipped=0;
        int idx =   rand()%(SIZE*SIZE);     //choose electron
        int move =  rand()%4;
        int overbound=0;
        
        alpha.swap(&beta, int(idx/SIZE),idx%SIZE, move,&flipped,&overbound);
        
        beta.copy_config_to(&alpha);
        
    }
    dur_copy = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    dur_copy -= dur_swap;
    
    
    double dur_tonumber;
    start = std::clock();
    for (int i=0; i<1000000; i++) {
        
        int flipped=0;
        int idx =   rand()%(SIZE*SIZE);     //choose electron
        int move =  rand()%4;
        int overbound=0;
        
        alpha.swap(&beta, int(idx/SIZE),idx%SIZE, move,&flipped,&overbound);
        beta.copy_config_to(&alpha);
        tonumber(alpha);
    }
    dur_tonumber = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    dur_tonumber -= dur_swap;
    dur_tonumber -= dur_copy;
    
    
    
    int countaccept=0;
    start = std::clock();
    for (int i=0; i<1000000; i++) {
        
        int flipped=0;
        int idx =   rand()%(SIZE*SIZE);     //choose electron
        int move =  rand()%4;
        int overbound=0;
        
        alpha.swap(&beta, int(idx/SIZE),idx%SIZE, move,&flipped,&overbound);
        beta.copy_config_to(&alpha);
        
        double deter_0;
        long long number = tonumber(alpha);
        unordered_map<long long,double>::const_iterator got = dMap.find (number);
        if ( got == dMap.end() ){
            deter_0 = determinant(slater[0]);
            pair<long long,double> tmp (number,deter_0);
            dMap.insert( tmp);
            //dMap.insert ( {tonumber(config_level[0]), deter_0} );
        }
        else{
            countaccept+=1;
            deter_0 = got->second;
        }

    }
    dur_mapFindInsert = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    dur_mapFindInsert -= dur_tonumber;
    dur_mapFindInsert -= dur_swap;
    dur_mapFindInsert -= dur_copy;
    
    
    double dur_find;
    start = std::clock();
    for (int i=0; i<1000000; i++) {
        
        
        double deter_0;
        long long number = tonumber(alpha);
        unordered_map<long long,double>::const_iterator got = dMap.find (number);
        if ( got == dMap.end() ){
            deter_0 = determinant(slater[0]);
            pair<long long,double> tmp (number,deter_0);
            dMap.insert( tmp);
            //dMap.insert ( {tonumber(config_level[0]), deter_0} );
        }
        else{
            //countaccept+=1;
            deter_0 = got->second;
        }
        
    }
    dur_find = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    dur_find -= dur_tonumber;
    
    
    
    
    cout<<"test performance with running 100,0000 times"<<endl;
    cout<<"duration_rand_init_no_d : "<<duration_rand_init_no_d<<endl;
    cout<<"duration_set_slater : "<<duration_set_slater<<endl;
    cout<<"duration_determinant : "<<duration_determinant<<endl;
    cout<<"duration_swap : "<<dur_swap<<endl;
    cout<<"duration_copy : "<<dur_copy<<endl;
    cout<<"duration_tonumber : "<<dur_tonumber<<endl;
    cout<<endl;
    cout<<"dur_find : "<<dur_find<<endl;
    cout<<"dur_mapFindInsert : "<<dur_mapFindInsert<<endl;
    cout<<"hash_map_efficiency:"<<endl;
    cout<<"total = "<<1000000<< "  accept = "<<countaccept<<"ratio"<<double(countaccept)/1000000<<endl;


    

return 0;
}
