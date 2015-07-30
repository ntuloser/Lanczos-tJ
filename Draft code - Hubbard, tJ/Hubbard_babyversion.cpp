#include<iostream>
#include<cstdlib>
#include<ctime>
#include<cmath>
#include<stdio.h>
#include<stdlib.h>
//#include<config.h>

#define Size 4
#define Pi 3.1415926535897932
#define delta 0
#define avoid_doublon 0 //Hubbard = 0   t-J = 1

using namespace std;


const int t=1;
const double J=0;//0.33;
const double U=1;
const int num_ele=Size*Size-delta;


class onsite{
	private:
	public:
		int idx_up;
		int idx_down;
		int state;
};

void swap(onsite*** lattice,int x,int y,int x_new,int y_new,int up_new[num_ele/2][2], int down_new[num_ele/2][2],double slater_new[num_ele/2][num_ele/2],const double a_ij[2*Size-1][2*Size-1] ){

	// Simultaneously
	// Updating the electron config.
	// and

	// construct the new slater determinant by
	/////////////////////////////////////
	// UPDATING the slater determinant
	/////////////////////////////////////
	// case 0 hopping
	// 0-1 hopping the up electron
	// 0-2 hopping the down electron
	// case 1 superexchange
	//

	if (lattice[x][y]->state==0) {

		//empty<-->up , hopping up elec
		if (lattice[x_new][y_new]->state==1) {

			int idx_up=lattice[x_new][y_new]->idx_up;

			//empty
			//empty
			up_new[ idx_up ][0]=x;
			up_new[ idx_up ][1]=y;

			for (int i=0;i<num_ele/2;i++){
				int x=up_new[idx_up][0]-down_new[i][0];
				int y=up_new[idx_up][1]-down_new[i][1];
				slater_new[idx_up][i]=a_ij[x+(Size-1)][y+(Size-1)];
			}
		}

		//empty<-->down , hopping down elec
		else if (lattice[x_new][y_new]->state==2) {

			int idx_down = lattice[x_new][y_new]->idx_down;

			//empty
			//empty
			down_new[ idx_down ][0] = x;
			down_new[ idx_down ][1] = y;

			for (int i=0;i<num_ele/2;i++)
			{
				int x=up_new[i][0]-down_new[idx_down][0];
				int y=up_new[i][1]-down_new[idx_down][1];
				slater_new[i][idx_down]=a_ij[x+(Size-1)][y+(Size-1)];
			}
		}

	}
	else if (lattice[x][y]->state==1) {
		//up<-->empty, hopping up elec
		if (lattice[x_new][y_new]->state==0) {

			int idx_up = lattice[x][y]->idx_up;

			up_new[ idx_up ][0] =x_new;
			up_new[ idx_up ][1] =y_new;
			//empty
			//empty

			for (int i=0;i<num_ele/2;i++){
				int x=up_new[idx_up][0]-down_new[i][0];
				int y=up_new[idx_up][1]-down_new[i][1];
				slater_new[idx_up][i]=a_ij[x+(Size-1)][y+(Size-1)];
			}

		}

		//up<-->down, superexchange !!
		else if (lattice[x_new][y_new]->state==2) {

			int idx_up = lattice[x][y]->idx_up;
			int idx_down=lattice[x_new][y_new]->idx_down;

			up_new[ lattice[x][y]->idx_up ][0] =x_new;
			up_new[ lattice[x][y]->idx_up ][1] =y_new;
			down_new[ lattice[x_new][y_new]->idx_down ][0] = x;
			down_new[ lattice[x_new][y_new]->idx_down ][1] = y;


			for (int i=0;i<num_ele/2;i++){
				int x=up_new[idx_up][0]-down_new[i][0];
				int y=up_new[idx_up][1]-down_new[i][1];
				slater_new[idx_up][i]=a_ij[x+(Size-1)][y+(Size-1)];
			}
			for (int i=0;i<num_ele/2;i++)
			{
				int x=up_new[i][0]-down_new[idx_down][0];
				int y=up_new[i][1]-down_new[idx_down][1];
				slater_new[i][idx_down]=a_ij[x+(Size-1)][y+(Size-1)];
			}

		}

	}

	else if (lattice[x][y]->state==2) {
		//down<-->empty, hopping down elec.
		if (lattice[x_new][y_new]->state==0) {

			int idx_down = lattice[x][y]->idx_down;

			down_new[ lattice[x][y]->idx_down ][0] =x_new;
			down_new[ lattice[x][y]->idx_down ][1] =y_new;
			//empty
			//empty

			for (int i=0;i<num_ele/2;i++)
			{
				int x=up_new[i][0]-down_new[idx_down][0];
				int y=up_new[i][1]-down_new[idx_down][1];
				slater_new[i][idx_down]=a_ij[x+(Size-1)][y+(Size-1)];
			}

		}

		//down<-->up, superexchange !!
		else if (lattice[x_new][y_new]->state==1) {

			int idx_down = lattice[x][y]->idx_down;
			int idx_up   = lattice[x_new][y_new]->idx_up;

			down_new[ lattice[x][y]->idx_down ][0] = x_new;
			down_new[ lattice[x][y]->idx_down ][1] = y_new;
			up_new[ lattice[x_new][y_new]->idx_up ][0] = x;
			up_new[ lattice[x_new][y_new]->idx_up ][1] = y;

			for (int i=0;i<num_ele/2;i++){
				int x=up_new[idx_up][0]-down_new[i][0];
				int y=up_new[idx_up][1]-down_new[i][1];
				slater_new[idx_up][i]=a_ij[x+(Size-1)][y+(Size-1)];
			}
			for (int i=0;i<num_ele/2;i++)
			{
				int x=up_new[i][0]-down_new[idx_down][0];
				int y=up_new[i][1]-down_new[idx_down][1];
				slater_new[i][idx_down]=a_ij[x+(Size-1)][y+(Size-1)];
			}

		}
	}
}//end function swap


void hop(int updown, onsite*** lattice,int x,int y,int x_new,int y_new,int up_new[num_ele/2][2], int down_new[num_ele/2][2],double slater_new[num_ele/2][num_ele/2],const double a_ij[2*Size-1][2*Size-1]  ){
    
    
	if (lattice[x][y]->state==0) {

        cout<<"idx_up = "<<lattice[x][y]->idx_up<<endl;
        cout<<"idx_down = "<<lattice[x][y]->idx_down<<endl;
        cout<<x<<", "<<y<<", "<<lattice[x][y]->state<<", GG at function hop"<<endl;
        cout<<"idx_up_new = "<<lattice[x_new][y_new]->idx_up<<endl;
        cout<<"idx_down_new = "<<lattice[x_new][y_new]->idx_down<<endl;
        cout<<endl;

	}
    else if (updown==0){//moving up electron
        
        int idx_up = lattice[x][y]->idx_up;
        if (idx_up<0) {
            cout<<"FUCK"<<endl;
        }
        
        up_new[ idx_up ][0] =x_new;
        up_new[ idx_up ][1] =y_new;
        //empty
        //empty
        
        for (int i=0;i<num_ele/2;i++){
            int x=up_new[idx_up][0]-down_new[i][0];
            int y=up_new[idx_up][1]-down_new[i][1];

            slater_new[idx_up][i]=a_ij[x+(Size-1)][y+(Size-1)];
        }
        

    }
    else if (updown==1){//moving down electron
        
        int idx_down = lattice[x][y]->idx_down;
        
        if (idx_down<0) {
            cout<<"FUCK"<<endl;
        }
        
        down_new[ lattice[x][y]->idx_down ][0] =x_new;
        down_new[ lattice[x][y]->idx_down ][1] =y_new;
        //empty
        //empty
        
        for (int i=0;i<num_ele/2;i++)
        {
            int x=up_new[i][0]-down_new[idx_down][0];
            int y=up_new[i][1]-down_new[idx_down][1];

            slater_new[i][idx_down]=a_ij[x+(Size-1)][y+(Size-1)];
        }
        

    }

}//End function hop








double determinant(double a[num_ele/2][num_ele/2])
{
	double x[num_ele/2][num_ele/2];
	for (int i=0;i<num_ele/2;i=i+1)
		for (int j=0;j<num_ele/2;j=j+1)
			x[i][j]=a[i][j];
	int k=0,n=0;
	while (k<num_ele/2)
	{
		int l=k;
		while (abs(x[l][k])<0.00001)
		{
			l=l+1;
			if (l==num_ele/2)
				return 0;
		}
		if (l!=k)
		{
			n=n+1;
			for (int i=0;i<num_ele/2;i=i+1)
			{
				double b;
				b=x[k][i];
				x[k][i]=x[l][i];
				x[l][i]=b;
			}
		}

		for (int i=k+1;i<num_ele/2;i=i+1)
		{
			double r=x[i][k]/x[k][k];
			for (int j=k;j<num_ele/2;j=j+1)
			{
				x[i][j]=x[i][j]-r*x[k][j];
			}
		}
		k=k+1;
	}
	double det=1;
	for (int i=0;i<num_ele/2;i=i+1)
	{det=det*x[i][i];
	}
	return det*pow(-1.0,n);
}




int main(){

	////Variational parameter////
	double D=1.;
	double Mu=0.;
    double g=1.;

    
    
	/////////////////////////////
    double step_size=0.003;
    double E_tmp[3];
    double num_doub=0;
    int tot_doub=0;
    
    double Etot=0;
    double errE=0;
    double avgE=0;
    double tot_accept=0;

    double tot_det=0;
    
    
    
    
    // The Brillouin Zone with periodic in x, antiperiodic in y
    double kk[Size*Size][2];    // 1-dim vector representation of kx,ky
    for (int idx =0; idx<Size*Size; idx++) {
        int i=      idx % Size;
        int j=      idx / Size;
        kk[idx][0]= Pi*( 1.-2.0*i/Size);
        kk[idx][1]= Pi*( (Size-1.0)/Size-2.0*j/Size);
        //Checked
        //printf("%f, %f\n",kk[idx][0],kk[idx][1]);
    }
    
    //// all possible a_ij
    //// consider all possible r_i - r_j
    //// This generate a 2S-1 X 2S-1 matrix
    //// We shall call the matrix by giving the displacement vector
    
    double a_ij[2*Size-1][2*Size-1];
    for (int idx=0; idx<(2*Size-1)*(2*Size-1); idx++) {
        int i=      idx % (2*Size-1);
        int j=      idx / (2*Size-1);
        
        a_ij[i][j]=0;
        
        for (int m=0; m<Size*Size; m++) {
            double ek,Dk;
            ek = (-2.)*t*(cos(kk[m][0])+cos(kk[m][1]))-Mu;
            Dk = D*(cos(kk[m][0])-cos(kk[m][1])) ;
            
            if (D>0.0001) {
                //If D is not vanishing
                //sine term is cancel when summing over the BZ zone.
                a_ij[i][j] += Dk/(ek+pow(pow(ek,2)+pow(Dk,2),0.5) )*cos(kk[m][0]*(i-(Size-1))+kk[m][1]*(j-(Size-1)) )  / (Size*Size);
            }
            else
                cout<<"error in vanishing D"<<endl;
            
            //element[i][j]+=0.5*(1-abs(ek)/ek)*cos(kk[m][0]*(i-(s-1))+kk[m][1]*(j-(s-1)))/(s*s);
            
        }
        //Checked
        //printf("a[%i][%i]=%f\n",i,j,a_ij[i][j]);
    }


	int Variational_steps=1;
	for (int stp=0;stp<Variational_steps ; stp++) {










		///////////////////////////
		///// initial config  /////     no doublon.
		///////////////////////////

		onsite*** squarelattice = new onsite**[Size];
		for (int i=0; i<Size; i++) {
			squarelattice[i]= new onsite*[Size];
			for (int j=0; j<Size; j++) {
				squarelattice[i][j]=new onsite;
			}
		}

		srand((unsigned)time(0));
		int electronsup[num_ele/2][2];
		int electronsdown[num_ele/2][2];
		//	int site[Size][Size];
		int count=0;
		while (count < num_ele/2 ) // Filling the up spin part, total s*s/2-delta/2 electrons.
		{
			int r=rand()%(Size*Size);
			int x_ran=r%Size;
			int y_ran=r/Size;
			int filled=0;
			for (int i=0; i<count; i++) {
				if (electronsup[i][0]==x_ran && electronsup[i][1]==y_ran) {
					filled=1;
				}
			}
			if (filled==0) {
				electronsup[count][0]=x_ran;
				electronsup[count][1]=y_ran;
				count+=1;
			}

		}
		count = 0;
		while (count < num_ele/2 ) // Filling the down spin part, total s*s/2-delta/2 electrons.
		{
			int r=rand()%(Size*Size);
			int x_ran=r%Size;
			int y_ran=r/Size;
			int filled=0;
			for (int i=0; i<count; i++) {
				if (electronsdown[i][0]==x_ran && electronsdown[i][1]==y_ran) {
					filled=1;
				}
			}
			//// avoid doubly occupied////
            if (avoid_doublon == 1) {
                
                for (int i=0; i<num_ele/2; i++) {
                    if (electronsup[i][0]==x_ran && electronsup[i][1]==y_ran) {
                            filled=1;
                    }
                }
            }
			//////////////////////////////
			if (filled==0) {
				electronsdown[count][0]=x_ran;
				electronsdown[count][1]=y_ran;
				count+=1;
			}
		}


		//////////////////////////////////////////
		////// Construct the site matrix /////////
		//////////////////////////////////////////

		// Filling in the occupied, empty->0, up->1, down->2, updown->3
		for (int i=0; i<Size; i++) {
			for (int j=0; j<Size; j++) {
				//			site[i][j]=0;
				squarelattice[i][j]->state=0;
				squarelattice[i][j]->idx_up=-1;
				squarelattice[i][j]->idx_down=-1;
			}
		}
        
		for (int i=0;i<num_ele/2;i++){
            
			int x=electronsup[i][0];
			int y=electronsup[i][1];
			//		site[x][y]+=1;
			squarelattice[x][y]->state+=1;
			squarelattice[x][y]->idx_up=i;

			x=electronsdown[i][0];
			y=electronsdown[i][1];
			//		site[x][y]+=2;
			squarelattice[x][y]->state+=2;
			squarelattice[x][y]->idx_down=i;
		}



		/////////////////////////////
		//    The wave function    //
		/////////////////////////////
		//determining the slater matrix of the initial state
		double slater[num_ele/2][num_ele/2];
		double slater_new[num_ele/2][num_ele/2];
		for (int i=0;i<num_ele/2;i++){
			for (int j=0;j<num_ele/2;j++){
				// The ith electronup, and the jth electrondown.
				// Call the matrix a_ij by displacement vector
				// between two electron.
				//
				// We shift the displacement so that it is positive
				// to call matrix element.

				int x=electronsup[i][0]-electronsdown[j][0];
				int y=electronsup[i][1]-electronsdown[j][1];
				slater[i][j]=a_ij[x+(Size-1)][y+(Size-1)];
			}
		}

		/////////////////////////////////
		//**the monte carlo procedure**//
		/////////////////////////////////
		int sample=100000;
		int interval=500;
		int warmup=100000;
		int takeInv=500;
		int totsteps=sample*interval+warmup;
		for (int steps=0;steps<=totsteps;steps++){



			/////////////////////////////////////////////
			//Randomly generating a new configuration////
			/////////////////////////////////////////////

			int electronsup_new[num_ele/2][2];         //The new configuration
			int electronsdown_new[num_ele/2][2];
			for (int i=0;i<(num_ele/2);i++)
			{
				electronsup_new[i][0]=electronsup[i][0];
				electronsup_new[i][1]=electronsup[i][1];
				electronsdown_new[i][0]=electronsdown[i][0];
				electronsdown_new[i][1]=electronsdown[i][1];
			}

			for (int i=0;i<num_ele/2;i++){
				for (int j=0;j<num_ele/2;j++){
					slater_new[i][j]=slater[i][j];
				}
			}

			int flipped=0;
			int move;
			int idx;
            
            
            move=rand()%4;              //choose direction
            int move_x,move_y;
            if (move==0) {
                move_x=1;
                move_y=0;
            }
            else if(move==1){
                move_x=-1;
                move_y=0;
            }
            else if(move==2){
                move_x=0;
                move_y=-1;
            }
            else if(move==3){
                move_x=0;
                move_y=1;
            }
            
            
            ////////t_J model////////
            if (avoid_doublon==1) {
                idx=rand()%(Size*Size);     //choose site
                int x=idx%Size;
                int y=idx/Size;
                int x_new = (x+move_x+Size)%Size;
                int y_new = (y+move_y+Size)%Size;
                
                
                ///////////////////////////
                ///////Super exchange /////
                if (squarelattice[x_new][y_new]->state== squarelattice[x][y]->state) continue;
                
                swap(squarelattice,x,y,x_new,y_new,electronsup_new,electronsdown_new,slater_new,a_ij);
                flipped=1;
                /////////////////////////////
            }
            
            
            ////Hubbard model////////
            else{  //avoid_doublon ==0; Hubbard!
                idx = rand()%num_ele;     //choose electron
                
                if (idx<num_ele/2) {
                    int x=electronsup[idx][0];
                    int y=electronsup[idx][1];
                    int x_new = (x+move_x+Size)%Size;
                    int y_new = (y+move_y+Size)%Size;
                    
                    if (squarelattice[x_new][y_new]->idx_up==-1) {
                        //electronsup_new[idx][0]=x_new;
                        //electronsup_new[idx][1]=y_new;
                        
                        hop(0,squarelattice,x,y,x_new,y_new,electronsup_new,electronsdown_new,slater_new,a_ij);
                        flipped=1;
                    }
                    
                }//end if
                else{//idx >= num_ele/2
                    //idx=idx-num_ele/2;
                    int x=electronsdown[idx-num_ele/2][0];
                    int y=electronsdown[idx-num_ele/2][1];
                    int x_new = (x+move_x+Size)%Size;
                    int y_new = (y+move_y+Size)%Size;
                    
                    if (squarelattice[x_new][y_new]->idx_down==-1) {
                        //electronsdown_new[idx][0]=x_new;
                        //electronsdown_new[idx][1]=y_new;
                        hop(1,squarelattice,x,y,x_new,y_new,electronsup_new,electronsdown_new,slater_new,a_ij);
                        flipped=1;
                    }
                    
                    
                }//end else (up or down electron)
            }//end else(tJ or Hubbard)
            

            
            if (steps>warmup && steps%interval==0) {
                
                tot_det+= pow(determinant(slater_new)/determinant(slater),2);
            }


			////////////////////////////////////////////////////////////////
			///////////   Warm up   //////////////////
			//////////////////////////////////////////


			///////////////////////////////////////
			///   Metropolis condition      ///////
			///////////////////////////////////////

			/// Given the random probability 0<p<1 ///
			double p=0;
			while (p==0){
				p=(rand()+1)/(double)(RAND_MAX);
			}

			if (1/*steps<warmup*/){
				if (p<pow(determinant(slater_new)/determinant(slater),2)/* *pow(g,2*(numd1-numd))*/|| abs(determinant(slater))<0.00001){


					///////////
					//SWAP  lattice[x][y] and lattice[x+move_x][y+move_y]
					///////////
					int move_x,move_y;
					if (move==0) {
						move_x=1;
						move_y=0;
					}
					else if(move==1){
						move_x=-1;
						move_y=0;
					}
					else if(move==2){
						move_x=0;
						move_y=-1;
					}
					else if(move==3){
						move_x=0;
						move_y=1;
					}


                    if (flipped==1) {
                        
                        tot_accept+=1;
                        
                        if (avoid_doublon==1) {
                            
                            
                            int x=idx%Size;
                            int y=idx/Size;
                            
                            int x_new = (x+move_x+Size)%Size;
                            int y_new = (y+move_y+Size)%Size;
                            
                            int temp_state      =squarelattice[x][y]->state;
                            int temp_idx_up     =squarelattice[x][y]->idx_up;
                            int temp_idx_down   =squarelattice[x][y]->idx_down;
                            
                            squarelattice[x][y]->state=squarelattice[x_new][y_new]->state;
                            squarelattice[x][y]->idx_up=squarelattice[x_new][y_new]->idx_up;
                            squarelattice[x][y]->idx_down=squarelattice[x_new][y_new]->idx_down;
                            
                            squarelattice[x_new][y_new]->state = temp_state;
                            squarelattice[x_new][y_new]->idx_up = temp_idx_up;
                            squarelattice[x_new][y_new]->idx_down=temp_idx_down;
                            
                        }
                        else{
                            
                            if (idx<num_ele/2) {
                                ///moving up electron from x,y to x_new, y_new
                                int x=electronsup[idx][0];
                                int y=electronsup[idx][1];
                                int x_new = (x+move_x+Size)%Size;
                                int y_new = (y+move_y+Size)%Size;
                                
                                
                                int temp_idx_up     =squarelattice[x][y]->idx_up;
                                
                                squarelattice[x][y]->state -= 1;
                                squarelattice[x][y]->idx_up = -1;
                                
                                squarelattice[x_new][y_new]->state += 1;
                                squarelattice[x_new][y_new]->idx_up = temp_idx_up;
                                
                            }
                            else{
                                ///moving down electron from x,y to x_new, y_new
                                int x=electronsdown[idx-num_ele/2][0];
                                int y=electronsdown[idx-num_ele/2][1];
                                int x_new = (x+move_x+Size)%Size;
                                int y_new = (y+move_y+Size)%Size;
                                
                                
                                int temp_idx_down   =squarelattice[x][y]->idx_down;
                                
                                squarelattice[x][y]->state -= 2;
                                squarelattice[x][y]->idx_down = -1;
                                
                                squarelattice[x_new][y_new]->state += 2;
                                squarelattice[x_new][y_new]->idx_down = temp_idx_down;
                                
                            }
                            
                        }
                    }
                    
                    
					////////////////////////////////////////////
					//////// Switch electron list       ////////
					//////// Switch slater determinant  ////////
                    ////////////////////////////////////////////
					for (int i=0;i<num_ele/2;i=i+1){
						for (int j=0;j<=1;j=j+1)
						{
							electronsup[i][j]=electronsup_new[i][j];
							electronsdown[i][j]=electronsdown_new[i][j];
						}
					}
					for (int i=0;i<num_ele/2;i=i+1){
						for (int j=0;j<num_ele/2;j=j+1){
							slater[i][j]=slater_new[i][j];
						}
					}
                        
					///////////////////////////////////////////////////////



				}//end if
			}

			////////////////////////////////////////////////////////////////////
			///////////  The Post Warm up  Era  //////////////
			//////////////////////////////////////////
			else {
			}

			////////////////////////////////////////////////////////////////
			///////////  guys, Time to Sampling  ~~  /////////////
			//////////////////////////////////////////

			//
			// f_\alpha = \sum_\beta  < \alpha | O | \beta >* <\beta | \Psi> / <\alpha | \Psi>
			//
			if (steps%interval==0 && steps>warmup){


				double Ek=0;
				double ratio=0;
				int electronsup_beta[num_ele/2][2],electronsdown_beta[num_ele/2][2];
				double slater_beta[num_ele/2][num_ele/2];

				for (int i=0;i<num_ele/2;i=i+1){
					for (int j=0;j<2;j=j+1)
					{
						electronsup_beta[i][j]=electronsup[i][j];
						electronsdown_beta[i][j]=electronsdown[i][j];
					}
				}

				for (int m=0;m<num_ele/2;m=m+1){
					for (int n=0;n<num_ele/2;n=n+1){
						slater_beta[m][n]=slater[m][n];
					}
				}

                
                // Counting doublon///
                num_doub=0;
                for (int i=0; i < Size ; i++) {
                    for (int j=0; j<Size; j++) {
                        if (squarelattice[i][j]->state==3) {
                            num_doub+=1;
                        }
                        
                    }
                }

                Etot +=num_doub*U;
                tot_doub+=num_doub;


                if (avoid_doublon==1) {
                    
                    //scaning updown and leftright direction n.n. neighbor
                    for (int i=0; i < Size; i++) {
                        for (int j=0; j<Size; j++) {
                            
                            for (int k=0; k<2; k++) {
                                int move_i, move_j;
                                if (k==0) {
                                    move_i=1;
                                    move_j=0;
                                }
                                else{
                                    move_i=0;
                                    move_j=1;
                                }
                                
                                int i_new = (i+move_i)%Size;
                                int j_new = (j+move_j)%Size;
                                
                                int overbound  = ( (i+move_i)/Size );  // || (j+move_j)/Size );
                                //periodic boundary condition only in x axis.
                                
                                
                                double E_beta=0;
                                int state = squarelattice[i][j]->state;
                                int state_new = squarelattice[i_new][j_new]->state;
                                
                                
                                if (state==0) {
                                    if(state_new==1){
                                        
                                        //switching
                                        //(i)(j)
                                        //emtpy
                                        //empty
                                        //(i_new)(j_new)
                                        int idx_up = squarelattice[i_new][j_new]->idx_up;
                                        electronsup_beta[idx_up][0]=i;
                                        electronsup_beta[idx_up][1]=j;
                                        
                                        
                                        for (int m=0; m<num_ele/2; m++) {
                                            int x=electronsup_beta[idx_up][0]-electronsdown_beta[m][0];
                                            int y=electronsup_beta[idx_up][1]-electronsdown_beta[m][1];
                                            slater_beta[idx_up][m]=a_ij[x+(Size-1)][y+(Size-1)];
                                        }
                                        
                                        if (not overbound) {
                                            E_beta=-t*determinant(slater_beta)/determinant(slater);
                                        }
                                        else{
                                            E_beta=t*determinant(slater_beta)/determinant(slater);
                                        }
                                        
                                        
                                        /// Resume the inital config ///
                                        electronsup_beta[idx_up][0] = i_new;
                                        electronsup_beta[idx_up][1] = j_new;
                                        
                                        
                                        for (int m=0; m<num_ele/2; m++) {
                                            slater_beta[idx_up][m]=slater[idx_up][m];
                                        }
                                        /////////////////////////////////
                                        
                                        
                                    }//end of case 0,1
                                    
                                    
                                    else if(state_new==2){
                                        
                                        //switching
                                        //(i)(j)
                                        //emtpy
                                        //empty
                                        //(i_new)(j_new)
                                        int idx_down = squarelattice[i_new][j_new]->idx_down;
                                        electronsdown_beta[idx_down][0]=i;
                                        electronsdown_beta[idx_down][1]=j;
                                        
                                        
                                        for (int m=0; m<num_ele/2; m++) {
                                            int x=electronsup_beta[m][0]-electronsdown_beta[idx_down][0];
                                            int y=electronsup_beta[m][1]-electronsdown_beta[idx_down][1];
                                            slater_beta[m][idx_down]=a_ij[x+(Size-1)][y+(Size-1)];
                                        }
                                        
                                        if (not overbound) {
                                            E_beta=-t*determinant(slater_beta)/determinant(slater);
                                        }
                                        else{
                                            E_beta=t*determinant(slater_beta)/determinant(slater);
                                        }
                                        
                                        
                                        
                                        /// Resume the inital config//////
                                        electronsdown_beta[idx_down][0]=i_new;
                                        electronsdown_beta[idx_down][1]=j_new;
                                        
                                        
                                        for (int m=0; m<num_ele/2; m++) {
                                            slater_beta[m][idx_down]=slater[m][idx_down];
                                        }
                                        //////////////////////////////////
                                        
                                        
                                    }//end of case 0,2
                                }
                                else if(state==1){
                                    if (state_new==0) {
                                        
                                        
                                        //switching
                                        //(i)(j)
                                        int idx_up = squarelattice[i][j]->idx_up;
                                        electronsup_beta[idx_up][0] = i_new;
                                        electronsup_beta[idx_up][1] = j_new;
                                        //(i+1)(j)
                                        //empty
                                        //empty
                                        
                                        for (int m=0; m<num_ele/2; m++) {
                                            int x=electronsup_beta[idx_up][0]-electronsdown_beta[m][0];
                                            int y=electronsup_beta[idx_up][1]-electronsdown_beta[m][1];
                                            slater_beta[idx_up][m]=a_ij[x+(Size-1)][y+(Size-1)];
                                        }
                                        
                                        
                                        if (not overbound) {
                                            E_beta=-t*determinant(slater_beta)/determinant(slater);
                                        }
                                        else{
                                            E_beta=t*determinant(slater_beta)/determinant(slater);
                                        }
                                        
                                        
                                        /// Resume the inital config ///
                                        electronsup_beta[idx_up][0] = i;
                                        electronsup_beta[idx_up][1] = j;
                                        
                                        for (int m=0; m<num_ele/2; m++) {
                                            slater_beta[idx_up][m]=slater[idx_up][m];
                                        }
                                        /////////////////////////////////
                                        
                                        
                                        
                                    }//end of case 1,0
                                    else if(state_new==2){
                                        
                                        
                                        //switching
                                        //(i)(j)
                                        int idx_up = squarelattice[i][j]->idx_up;
                                        electronsup_beta[idx_up][0] = i_new;
                                        electronsup_beta[idx_up][1] = j_new;
                                        //(i_new)(j_new)
                                        int idx_down = squarelattice[i_new][j_new]->idx_down;
                                        electronsdown_beta[idx_down][0] = i;
                                        electronsdown_beta[idx_down][1] = j;
                                        
                                        for (int m=0; m<num_ele/2; m++) {
                                            int x=electronsup_beta[idx_up][0]-electronsdown_beta[m][0];
                                            int y=electronsup_beta[idx_up][1]-electronsdown_beta[m][1];
                                            slater_beta[idx_up][m]=a_ij[x+(Size-1)][y+(Size-1)];
                                        }
                                        for (int m=0; m<num_ele/2; m++) {
                                            int x=electronsup_beta[m][0]-electronsdown_beta[idx_down][0];
                                            int y=electronsup_beta[m][1]-electronsdown_beta[idx_down][1];
                                            slater_beta[m][idx_down]=a_ij[x+(Size-1)][y+(Size-1)];
                                        }
                                        
                                        
                                        E_beta=J/2*determinant(slater_beta)/determinant(slater);
                                        E_beta+=(-J)/2;
                                        
                                        /// Resume the inital config ///
                                        electronsup_beta[idx_up][0] = i;
                                        electronsup_beta[idx_up][1] = j;
                                        for (int m=0; m<num_ele/2; m++) {
                                            slater_beta[idx_up][m]=slater[idx_up][m];
                                        }
                                        electronsdown_beta[idx_down][0] = i_new;
                                        electronsdown_beta[idx_down][1] = j_new;
                                        for (int m=0; m<num_ele/2; m++) {
                                            slater_beta[m][idx_down]=slater[m][idx_down];
                                        }
                                        //////////////////////////////////
                                        
                                        
                                    }//end of case 1,2
                                    
                                }
                                else if(state==2){
                                    if (state_new==0) {
                                        
                                        //switching
                                        //(i)(j)
                                        int idx_down = squarelattice[i][j]->idx_down;
                                        electronsdown_beta[idx_down][0] = i_new;
                                        electronsdown_beta[idx_down][1] = j_new;
                                        //(i_new)(j_new)
                                        //empty
                                        //empty
                                        
                                        for (int m=0; m<num_ele/2; m++) {
                                            int x=electronsup_beta[m][0]-electronsdown_beta[idx_down][0];
                                            int y=electronsup_beta[m][1]-electronsdown_beta[idx_down][1];
                                            slater_beta[m][idx_down]=a_ij[x+(Size-1)][y+(Size-1)];
                                        }
                                        
                                        
                                        if (not overbound) {
                                            E_beta=-t*determinant(slater_beta)/determinant(slater);
                                        }
                                        else{
                                            E_beta=t*determinant(slater_beta)/determinant(slater);
                                        }
                                        
                                        /// Resume the inital config//////
                                        electronsdown_beta[idx_down][0] = i;
                                        electronsdown_beta[idx_down][1] = j;
                                        for (int m=0; m<num_ele/2; m++) {
                                            slater_beta[m][idx_down]=slater[m][idx_down];
                                        }
                                        //////////////////////////////////
                                        
                                        
                                    }//end of case 2,0
                                    else if(state_new==1){
                                        
                                        
                                        //switching
                                        //(i)(j)
                                        int idx_down = squarelattice[i][j]->idx_down;
                                        electronsdown_beta[idx_down][0] = i_new;
                                        electronsdown_beta[idx_down][1] = j_new;
                                        //(i_new)(j_new)
                                        int idx_up = squarelattice[i_new][j_new]->idx_up;
                                        electronsup_beta[idx_up][0] = i;
                                        electronsup_beta[idx_up][1] = j;
                                        
                                        for (int m=0; m<num_ele/2; m++) {
                                            int x=electronsup_beta[idx_up][0]-electronsdown_beta[m][0];
                                            int y=electronsup_beta[idx_up][1]-electronsdown_beta[m][1];
                                            slater_beta[idx_up][m]=a_ij[x+(Size-1)][y+(Size-1)];
                                        }
                                        for (int m=0; m<num_ele/2; m++) {
                                            int x=electronsup_beta[m][0]-electronsdown_beta[idx_down][0];
                                            int y=electronsup_beta[m][1]-electronsdown_beta[idx_down][1];
                                            slater_beta[m][idx_down]=a_ij[x+(Size-1)][y+(Size-1)];
                                        }
                                        
                                        
                                        E_beta=J/2*determinant(slater_beta)/determinant(slater);
                                        E_beta+=(-J)/2;
                                        
                                        
                                        /// Resume the inital config ///
                                        electronsdown_beta[idx_down][0] = i;
                                        electronsdown_beta[idx_down][1] = j;
                                        for (int m=0; m<num_ele/2; m++) {
                                            slater_beta[m][idx_down]=slater[m][idx_down];
                                        }
                                        electronsup_beta[idx_up][0] = i_new;
                                        electronsup_beta[idx_up][1] = j_new;
                                        for (int m=0; m<num_ele/2; m++) {
                                            slater_beta[idx_up][m]=slater[idx_up][m];
                                        }
                                        //////////////////////////////////
                                        
                                        
                                    }//end of case 2,1
                                    
                                }
                                
                                Etot+=E_beta;
                            }
                        }//End for_j
                    }//End for_i
                    
                }
                
                else{ // avoid_doublon==0 , Hubbard model !!
                 
                    for (int idx=0; idx<num_ele; idx++) {
                        for (int move=0; move<4; move++) {
                            
                            int move_x,move_y;
                            if (move==0) {
                                move_x=1;
                                move_y=0;
                            }
                            else if(move==1){
                                move_x=-1;
                                move_y=0;
                            }
                            else if(move==2){
                                move_x=0;
                                move_y=-1;
                            }
                            else if(move==3){
                                move_x=0;
                                move_y=1;
                            }
                            
                            double E_beta=0;
                            int overbound=0 ;
                            
                            if (idx<num_ele/2) {
                                int x=electronsup[idx][0];
                                int y=electronsup[idx][1];
                                int x_new = (x+move_x+Size)%Size;
                                int y_new = (y+move_y+Size)%Size;
                                
                                if ((y+move_y)==Size || (y+move_y)==-1  ) {
                                    overbound=1;
                                }
                                else overbound=0;
                                
                                if (squarelattice[x_new][y_new]->idx_up==-1) {

                                    
                                    //switching
                                    //(i)(j)
                                    int idx_up = idx;
                                    electronsup_beta[idx_up][0] = x_new;
                                    electronsup_beta[idx_up][1] = y_new;
                                    //(i_new)(j_new)
                                    //empty
                                    //empty
                                    
                                    for (int m=0; m<num_ele/2; m++) {
                                        int del_x=electronsup_beta[idx_up][0]-electronsdown_beta[m][0];
                                        int del_y=electronsup_beta[idx_up][1]-electronsdown_beta[m][1];
                                        slater_beta[idx_up][m]=a_ij[del_x+(Size-1)][del_y+(Size-1)];
                                    }
                                    
                                    
                                    if (not overbound) {
                                        E_beta=-t*determinant(slater_beta)/determinant(slater);
                                    }
                                    else{
                                        E_beta=t*determinant(slater_beta)/determinant(slater);
                                    }
                                    
                                    /// Resume the inital config//////
                                    electronsup_beta[idx_up][0] = x;
                                    electronsup_beta[idx_up][1] = y;
                                    for (int m=0; m<num_ele/2; m++) {
                                        slater_beta[idx_up][m]=slater[idx_up][m];
                                    }
                                    //////////////////////////////////
                                }
                                else E_beta=0;

                                
                                
                            }//end if
                            else{//idx >= num_ele/2
                                //idx=idx-num_ele/2;
                                int x=electronsdown[idx-num_ele/2][0];
                                int y=electronsdown[idx-num_ele/2][1];
                                int x_new = (x+move_x+Size)%Size;
                                int y_new = (y+move_y+Size)%Size;
                                
                                if (  (y+move_y)==Size || (y+move_y)==-1 ) {
                                    overbound=1;
                                }
                                else overbound=0;
                                
                                
                                if (squarelattice[x_new][y_new]->idx_down==-1) {

                                    
                                    //switching
                                    //(i)(j)
                                    int idx_down = idx-num_ele/2;
                                    electronsdown_beta[idx_down][0] = x_new;
                                    electronsdown_beta[idx_down][1] = y_new;
                                    //(i_new)(j_new)
                                    //empty
                                    //empty
                                    
                                    for (int m=0; m<num_ele/2; m++) {
                                        int del_x=electronsup_beta[m][0]-electronsdown_beta[idx_down][0];
                                        int del_y=electronsup_beta[m][1]-electronsdown_beta[idx_down][1];
                                        slater_beta[m][idx_down]=a_ij[del_x+(Size-1)][del_y+(Size-1)];
                                    }
                                    
                                    
                                    if (not overbound) {
                                        E_beta=-t*determinant(slater_beta)/determinant(slater);
                                    }
                                    else{
                                        E_beta=t*determinant(slater_beta)/determinant(slater);
                                    }
                                    
                                    /// Resume the inital config//////
                                    electronsdown_beta[idx_down][0] = x;
                                    electronsdown_beta[idx_down][1] = y;
                                    for (int m=0; m<num_ele/2; m++) {
                                        slater_beta[m][idx_down]=slater[m][idx_down];
                                    }
                                    //////////////////////////////////
                                }
                                else E_beta=0;

                                
                                
                            }//end else (up or down electron)

                            Etot+=E_beta;
                            
                        }//end for (move)
                    }//end for (idx)
                    
                }//if else t-J || Hubbard.



			}//End if sampling;


			//
		}//End of all steps
        
        
        
        
        
        
        /////////////////////
        /////////
        //      Derivative
        //
        ////////////////
        
        avgE=Etot/sample;
        double num_doub=tot_doub/double(sample);
        
        

        if (stp%3==0) {
            D+=0.001;
            Mu+=0;
            E_tmp[0]=avgE;
            
            cout<<"D = "<<D-0.001<<endl;
            cout<<"Mu = "<<Mu<<endl;
            cout<<"avgE = "<<E_tmp[0]/Size/Size<<endl;
            cout<<"Ek = "<<E_tmp[0]-num_doub<<endl;
            cout<<"doublon number="<<num_doub/Size/Size<<endl;
            cout<<"acceptance ratio:" <<tot_accept/totsteps<<endl;
            cout<<"tot_det = "<<tot_det/sample<<endl;


            
        }
        else if(stp%3==1){
            
            //With D+0.001, Mu
            E_tmp[1]=avgE;
            //
            D-=0.001;
            Mu+=0.001;
            //
        }
        else if(stp%3==2){
            //With D, Mu+0.001
            E_tmp[2]=avgE;
            //
            Mu-=0.001;
            //
            double grad_D  = (E_tmp[1]-E_tmp[0])/0.001;
            double grad_Mu = (E_tmp[2]-E_tmp[0])/0.001;
            
            cout<<"D = "<<D<<" , grad_D = "<<grad_D<<endl;
            cout<<"Mu = "<<Mu<<" , grad_Mu = "<<grad_Mu<<endl;
            cout<<"avgE = "<<E_tmp[0]<<endl;
            
            D   -= step_size*grad_D;
            Mu  -= step_size*grad_Mu;
        }
        
        
        

	}
	/*

	   avgE=Etot/sample;
	   errE=pow((E2tot/sample-pow(avgE,2))/sample,0.5);
	   cout<<" doublon number="<<avgnumd<<" err="<<errnumd<<endl;
	   cout<<" energy="<<avgE<<" err="<<errE<<endl;
	   cout<<" binding number="<<avgnumbind<<" err="<<errnumbind<<endl;
	   end_time=clock();
	   tot_time=(float)(end_time-start_time)/CLOCKS_PER_SEC;
	   cout<<" time="<<tot_time<<"sec";
	   system("pause");
	   return 0;*/




	//Without variation
	return 0;
}




/*
    Should normalized the variational parameter.
 
 
    Test result:
    1m24.294s (~65s in the loop) for no operation but running monte carlo, site 4*4

 */
