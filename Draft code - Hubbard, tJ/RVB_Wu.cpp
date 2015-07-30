#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#define s 4

using namespace std;
double determinant(double a[s*s/2][s*s/2]);
void createInverse(double c[s*s/2][s*s/2],double b[s*s/2][s*s/2],int n);
int main(void)
{
	clock_t start_time,end_time;
	float tot_time=0;
	start_time=clock();
	double pi=3.1415926535897932;
	const int t=1,u=1;
	const double g=0.81898098,D=0.52431606;
	double Etot=0,E2tot=0,errE,avgE;
	double numdtot=0,numd2tot=0,errnumd,avgnumd;

	double numbindtot=0,numbind2tot=0,errnumbind,avgnumbind;

	double II[s*s/2][s*s/2],III[s*s/2][s*s/2];   
	double ratio1=0;

	//creating the bloch wave numbers
	double kk[s*s][2];
	double kx,ky;
	int j=0,i=0,a=0;

	while (j<s)
	{
		i=0;
		while (i<s)
		{
			kk[a][0]=pi*(1-(2.0*i/s));
			kk[a][1]=pi*((s-1.0)/s-(2.0*j/s));
			a=a+1;
			i=i+1;
		}
		j=j+1;
	}

    //////////??????????????????????????
	//recording any possible matrix element in the slater determinant

	double element[s*2-1][s*2-1];
	for (int i=0;i<s*2-1;i=i+1)
	{
		for (int j=0;j<s*2-1;j=j+1)
		{
			element[i][j]=0;
			for (int m=0;m<s*s;m=m+1)
			{
				double uk2,vk2,ek,Dk;
				ek=(-2)*t*(cos(kk[m][0])+cos(kk[m][1])); //-2t(coskx+cosky)
				Dk=D*(cos(kk[m][0])-cos(kk[m][1]));     //D(coskx-cosky)
				//    uk2=0.5*(1+ek/pow(pow(ek,2)+pow(Dk,2),0.5));
				//    vk2=1-uk2;
				//    element[i][j]+=pow(abs(vk2/uk2),0.5)*abs(Dk)/Dk*cos(kk[m][0]*(i-(s-1))+kk[m][1]*(j-(s-1)))/(s*s);
				if (D>0.0001)  
					element[i][j]+=Dk/(pow(pow(ek,2)+pow(Dk,2),0.5)+ek)*cos(kk[m][0]*(i-(s-1))+kk[m][1]*(j-(s-1)))/(s*s);
				else
					element[i][j]+=0.5*(1-abs(ek)/ek)*cos(kk[m][0]*(i-(s-1))+kk[m][1]*(j-(s-1)))/(s*s); 
			}
            //printf("a[%i][%i]=%f\n",i,j,element[i][j]);

		}
	}

    
	//giving random initial configuration
	srand((unsigned)time(0));
	int electronsup[s*s/2][2];
	int electronsdown[s*s/2][2];
	int site[s][s];
	int count=0;
	while (count<s*s/2) // Filling the up spin part, total s*s/2 electrons.
	{
		int x=rand()%(s*s);
		int x_random=x/s;
		int y_random=x%s;
		int filled=0;
		for (int i=0;i<count;i=i+1) // check whether the electron is filled or not
		{
			if(electronsup[i][0]==x_random  &&  electronsup[i][1]==y_random)
				filled=1;
		}
		if (filled==0) // if not
		{
			electronsup[count][0]=x_random;
			electronsup[count][1]=y_random;
			count=count+1;
		}
	}
	count=0;
	while (count<s*s/2) // Filling the down spin part, total s*s/2 electrons.
	{
		int x=rand()%(s*s);
		int x_random=x/s;
		int y_random=x%s;
		int filled=0;
		for (int i=0;i<count;i=i+1) // check whether the electron is filled or not
		{
			if(electronsdown[i][0]==x_random &&electronsdown[i][1]==y_random)
				filled=1;
		}
		if (filled==0) // if not
		{
			electronsdown[count][0]=x_random;
			electronsdown[count][1]=y_random;
			count=count+1;
		}
	}
    

    
    
	//determining the slater matrix of the initial state

	double slater[s*s/2][s*s/2],slater1[s*s/2][s*s/2];
	for (int i=0;i<s*s/2;i=i+1)
		for (int j=0;j<s*s/2;j=j+1)
		{
			int x=electronsup[i][0]-electronsdown[j][0];
			int y=electronsup[i][1]-electronsdown[j][1];
			slater[i][j]=element[x+(s-1)][y+(s-1)];
		}

    
    
    
    
    
    
    /////////////////////////////////
	//**the monte carlo procedure**//
    /////////////////////////////////
    
	int sample=100000,interval=250,warmup=100000,takeInv=250;
	int totsteps=sample*interval+warmup;
	for (int steps=0;steps<=totsteps;steps=steps+1)
	{
        //////////////////////////////////////////
        ////// Construct the site matrix /////////
        //////////////////////////////////////////
        
		for (int i=0;i<s;i=i+1)
			for (int j=0;j<s;j=j+1)
				site[i][j]=0;
        // Filling in the occupied, empty->0, up->1, down->2, updown->3
		for (int i=0;i<s*s/2;i=i+1)
		{
			int x=electronsup[i][0];
			int y=electronsup[i][1];
			site[x][y]+=1;
			x=electronsdown[i][0];
			y=electronsdown[i][1];
			site[x][y]+=2;
		}
        
        
        /////////////////////////////////////////////
        //Randomly generating a new configuration////
        /////////////////////////////////////////////
        
		int electronsup1[s*s/2][2];         //The new configuration
		int electronsdown1[s*s/2][2];
		for (int i=0;i<s*s/2;i=i+1)
		{
			electronsup1[i][0]=electronsup[i][0];
			electronsup1[i][1]=electronsup[i][1];
			electronsdown1[i][0]=electronsdown[i][0];
			electronsdown1[i][1]=electronsdown[i][1];
		}

		//choose an electron and its hopping direction
		int m=rand()%(s*s); // choose electron
		int move=rand()%4;  // choose hopping
		int z,z1,n;
		int judge=0;
		if (m<(s*s/2))      // moving electronup
		{
			n=0;
			z=electronsup1[m][move/2];  //0,1 -->0 move x  2,3-->1 move y
			z1=z+((move%2)-0.5)*2;      //0--> -1, 1-->1,  2--> -1, 3--> 1
 
			if (z1==-1)                 // periodic boundary condition
				z1=s-1;
			if (z1==s)
				z1=0;
			for (int i=0;i<s*s/2;i=i+1) //
			{
				if (electronsup1[i][move/2]==z1&&electronsup1[i][abs(1-move/2)]==electronsup1[m][abs(1-move/2)])
					judge=1;
                    //if already occupied
			}

			if (judge==0)
				electronsup1[m][move/2]=z1;
		}
		else                // moving electrondown
		{
			n=1;
			m=m-s*s/2;

			z=electronsdown1[m][move/2];
			z1=z+((move%2)-0.5)*2;
			if (z1==-1)
				z1=s-1;
			if (z1==s)
				z1=0;
			for (int i=0;i<s*s/2;i=i+1)
			{
				if (electronsdown1[i][move/2]==z1&&electronsdown1[i][abs(1-move/2)]==electronsdown1[m][abs(1-move/2)])
					judge=1;
			}
			if (judge==0)
				electronsdown1[m][move/2]=z1;
		}

        
		//counting the doublon number of the original/new configuration

		int numd=0;
		int numd1=0;
		for (int i=0;i<s*s/2;i=i+1)
		{
			int x=electronsup[i][0];
			int y=electronsup[i][1];
			int x1=electronsup1[i][0];
			int y1=electronsup1[i][1];
			for (int j=0;j<s*s/2;j=j+1)
			{
				if (electronsdown[j][0]==x&&electronsdown[j][1]==y)
					numd=numd+1;
				if (electronsdown1[j][0]==x1&&electronsdown1[j][1]==y1)
					numd1=numd1+1;
			}
		}


		//start the process on inverse
		if (steps==warmup||(steps%takeInv==0&&steps>warmup))//() added
		{

			createInverse(slater,II,s*s/2);
			for (int i=0;i<s*s/2;i=i+1)
				for (int j=0;j<s*s/2;j=j+1)
					III[i][j]=II[i][j];
		}

        //
        // IF judge == 0 // That is the movement could be made , transition site not occupied.
        // A new config is generated.
        //
		if (judge==0)
		{ 
			for (int i=0;i<s*s/2;i=i+1)
				for (int j=0;j<s*s/2;j=j+1)
					slater1[i][j]=slater[i][j];
			if (n==0)   // x-case
			{
				for (int i=0;i<s*s/2;i=i+1)
				{
					int x=electronsup1[m][0]-electronsdown1[i][0];
					int y=electronsup1[m][1]-electronsdown1[i][1];
					slater1[m][i]=element[x+(s-1)][y+(s-1)];
				}
			}
			else        // y-case
			{
				for (int i=0;i<s*s/2;i=i+1)
				{
					int x=electronsup1[i][0]-electronsdown1[m][0];
					int y=electronsup1[i][1]-electronsdown1[m][1];
					slater1[i][m]=element[x+(s-1)][y+(s-1)];
				}
			}
            
            /// Given the random probability 0<p<1 ///
			double p=0;
			while (p==0)
				p=(rand()+1)/(double)(RAND_MAX);
            
            
            //////////////////////////////////////////
            ///////////   Warm up   //////////////////
            //////////////////////////////////////////
            
			if (steps<warmup)
			{
                ///// Taking determinant /////
				if (p<(pow(determinant(slater1)/determinant(slater),2)*pow(g,2*(numd1-numd)))||abs(determinant(slater))<0.00001)
				{   ////////////////////////
                    //////// Switch ////////
					for (int i=0;i<s*s/2;i=i+1)
						for (int j=0;j<=1;j=j+1)
						{
							electronsup[i][j]=electronsup1[i][j];
							electronsdown[i][j]=electronsdown1[i][j];
						}
					for (int i=0;i<s*s/2;i=i+1)
						for (int j=0;j<s*s/2;j=j+1)
							slater[i][j]=slater1[i][j];
					numd=numd1;
				}
			}
            
            /////////////////////////////////////
            ///  After Warm up /////////////////
            /////////////////////////////////////
			else
			{
				if (n==0){// moving up elec.
					for(int i=0;i<s*s/2;i=i+1)
					{
						ratio1+=slater1[m][i]*II[i][m];
					}
				}
				else{ //n==1 moving down elec.
					for(int i=0;i<s*s/2;i=i+1)
					{
						ratio1+=slater1[i][m]*II[m][i];
					}
				}

				if (p<pow(ratio1,2)*pow(g,2*(numd1-numd)))
				{//cout<<determinant(slater)<<"  "<<determinant(slater1)<<"  "<<p<<endl;
					for (int i=0;i<s*s/2;i=i+1)
						for (int j=0;j<=1;j=j+1)
						{
							electronsup[i][j]=electronsup1[i][j];
							electronsdown[i][j]=electronsdown1[i][j];
						}
					//更新INVERSE 
					if (n==0)
					{
						for(int i=0;i<s*s/2;i=i+1)
						{
							for(int j=0;j<s*s/2;j=j+1)
							{
								if (j!=m)
								{
									for(int k=0;k<s*s/2;k=k+1)
										II[i][j]-=III[i][m]/ratio1*slater1[m][k]*III[k][j];
								}
								else
									II[i][j]=III[i][j]/ratio1;
							}
						}
					}
					else
					{
						for (int i=0;i<s*s/2;i=i+1)
						{
							for (int j=0;j<s*s/2;j=j+1)
							{
								if(i!=m)
								{
									for (int k=0;k<s*s/2;k=k+1)
										II[i][j]-=III[m][j]/ratio1*slater1[k][m]*III[i][k];
								}
								else
									II[i][j]=III[i][j]/ratio1;
							}
						}
					}

					for (int i=0;i<s*s/2;i=i+1)
						for (int j=0;j<s*s/2;j=j+1)
							III[i][j]=II[i][j];
                    /// updating slater determinant
					for (int i=0;i<s*s/2;i=i+1)
						for (int j=0;j<s*s/2;j=j+1)
							slater[i][j]=slater1[i][j];
					numd=numd1;
				}
				ratio1=0;
			}
		}


        ///////////////////////////////////////////////////////
        // Sampling and calculate the physical properties./////
        ///////////////////////////////////////////////////////
		//   計
		//     算
		//       期
		//         望
		//           值
		//             !!   

		if (steps%interval==0 && steps>warmup)
		{

			int numbind=0;
			for (int i=0;i<s;i=i+1)          
				for (int j=0;j<s;j=j+1)
					if (site[i][j]==3)
					{
						int x1=i+1;
						if (x1==s)
							x1=0;
						int x2=i-1;
						if (x2==-1)
							x2=s-1;
						int y1=j+1;
						if (y1==s)
							y1=0;
						int y2=j-1;
						if (y2==-1)
							y2=s-1;

						if (site[x1][j]==0||site[x2][j]==0||site[i][y1]==0||site[i][y2]==0)
							numbind=numbind+1;
					}
			double Ek=0;
			double ratio=0;
			int electronsup2[s*s/2][2],electronsdown2[s*s/2][2];

			numbindtot=numbindtot+numbind/(double)(s*s);
			numbind2tot=numbind2tot+pow(numbind/double(s*s),2);
			numdtot=numdtot+numd/double(s*s);
			numd2tot=numd2tot+pow(numd/double(s*s),2);

			double slater2[s*s/2][s*s/2];
			for (int m=0;m<s*s/2;m=m+1)
				for (int n=0;n<s*s/2;n=n+1)
					slater2[m][n]=slater[m][n];
			for (int i=0;i<s*s/2;i=i+1)
				for (int j=0;j<2;j=j+1)
				{
					electronsup2[i][j]=electronsup[i][j];
					electronsdown2[i][j]=electronsdown[i][j];
				}
			for (int i=0;i<s*s/2;i=i+1)
			{
				int x,y,x1,y1,x2,y2,judge=0;
				x=electronsup[i][0];
				y=electronsup[i][1];
				x1=x+1;
				y1=y+1;
				x2=x-1;
				y2=y-1;
				if (x1==s)
					x1=0;
				if (y1==s)
					y1=0;
				if (x2==-1)
					x2=s-1;
				if (y2==-1)
					y2=s-1;

				int numd2=numd;




				for (int j=0;j<s*s/2;j=j+1)
					if (electronsup[j][0]==x1&&electronsup[j][1]==y)
						judge=1;

				if (judge==0)
				{ 
					electronsup2[i][0]=x1;
					for (int j=0;j<s*s/2;j=j+1)
					{
						int x=electronsup2[i][0]-electronsdown2[j][0];
						int y=electronsup2[i][1]-electronsdown2[j][1];
						slater2[i][j]=element[x+(s-1)][y+(s-1)];
					}
					for (int j=0;j<s*s/2;j=j+1)
					{
						if (electronsdown[j][0]==x&&electronsdown[j][1]==y)
							numd2=numd2-1;
						if (electronsdown[j][0]==x1&&electronsdown[j][1]==y)
							numd2=numd2+1;
					}
					for (int j=0;j<s*s/2;j=j+1)
					{
						ratio+=slater2[i][j]*II[j][i];
					}
					Ek=Ek-t*ratio*pow(g,numd2-numd);
					numd2=numd;
					ratio=0;
					for (int k=0;k<s*s/2;k=k+1)
						slater2[i][k]=slater[i][k];
					electronsup2[i][0]=x;
				}




				judge=0;
				for (int j=0;j<s*s/2;j=j+1)
					if (electronsup[j][0]==x2&&electronsup[j][1]==y)
						judge=1;
				if (judge==0)
				{
					electronsup2[i][0]=x2;
					for (int j=0;j<s*s/2;j=j+1)
					{
						int x=electronsup2[i][0]-electronsdown2[j][0];
						int y=electronsup2[i][1]-electronsdown2[j][1];
						slater2[i][j]=element[x+(s-1)][y+(s-1)];
					}
					for (int j=0;j<s*s/2;j=j+1)
					{
						if (electronsdown[j][0]==x&&electronsdown[j][1]==y)
							numd2=numd2-1;
						if (electronsdown[j][0]==x2&&electronsdown[j][1]==y)
							numd2=numd2+1;
					}
					for (int j=0;j<s*s/2;j=j+1)
					{
						ratio+=slater2[i][j]*II[j][i];
					}                 
					Ek=Ek-t*ratio*pow(g,numd2-numd);
					numd2=numd;
					ratio=0;
					for (int k=0;k<s*s/2;k=k+1)
						slater2[i][k]=slater[i][k];
					electronsup2[i][0]=x;
				}              


				judge=0;
				for (int j=0;j<s*s/2;j=j+1)
					if (electronsup[j][0]==x&&electronsup[j][1]==y1)
						judge=1;

				if (judge==0)
				{
					electronsup2[i][1]=y1;
					for (int j=0;j<s*s/2;j=j+1)
					{
						int x=electronsup2[i][0]-electronsdown2[j][0];
						int y=electronsup2[i][1]-electronsdown2[j][1];
						slater2[i][j]=element[x+(s-1)][y+(s-1)];
					}
					for (int j=0;j<s*s/2;j=j+1)
					{
						if (electronsdown[j][0]==x&&electronsdown[j][1]==y)
							numd2=numd2-1;
						if (electronsdown[j][0]==x&&electronsdown[j][1]==y1)
							numd2=numd2+1;
					}
					for (int j=0;j<s*s/2;j=j+1)
					{
						ratio+=slater2[i][j]*II[j][i];
					}
					if (y1==0)
						Ek=Ek+t*ratio*pow(g,numd2-numd);
					else
						Ek=Ek-t*ratio*pow(g,numd2-numd);
					numd2=numd;
					ratio=0;
					for (int k=0;k<s*s/2;k=k+1)
						slater2[i][k]=slater[i][k];
					electronsup2[i][1]=y;
				}


				judge=0;
				for (int j=0;j<s*s/2;j=j+1)
					if (electronsup[j][0]==x&&electronsup[j][1]==y2)
						judge=1;

				if (judge==0)
				{
					electronsup2[i][1]=y2;
					for (int j=0;j<s*s/2;j=j+1)
					{
						int x=electronsup2[i][0]-electronsdown2[j][0];
						int y=electronsup2[i][1]-electronsdown2[j][1];
						slater2[i][j]=element[x+(s-1)][y+(s-1)];
					}
					for (int j=0;j<s*s/2;j=j+1)
					{
						if (electronsdown[j][0]==x&&electronsdown[j][1]==y)
							numd2=numd2-1;
						if (electronsdown[j][0]==x&&electronsdown[j][1]==y2)
							numd2=numd2+1;
					}
					for (int j=0;j<s*s/2;j=j+1)
					{
						ratio+=slater2[i][j]*II[j][i];
					}
					if (y2==s-1)
						Ek=Ek+t*ratio*pow(g,numd2-numd);
					else
						Ek=Ek-t*ratio*pow(g,numd2-numd);
					numd2=numd;
					ratio=0;
					for (int k=0;k<s*s/2;k=k+1)
						slater2[i][k]=slater[i][k];
					electronsup2[i][1]=y;
				}
			}
			judge=0;
			for (int i=0;i<s*s/2;i=i+1)
			{
				int x,y,x1,y1,x2,y2,judge=0;
				x=electronsdown[i][0];
				y=electronsdown[i][1];
				x1=x+1;
				y1=y+1;
				x2=x-1;
				y2=y-1;
				if (x1==s)
					x1=0;
				if (y1==s)
					y1=0;
				if (x2==-1)
					x2=s-1;
				if (y2==-1)
					y2=s-1;
				int numd2=numd;


				for (int j=0;j<s*s/2;j=j+1)
					if (electronsdown[j][0]==x1&&electronsdown[j][1]==y)
						judge=1;

				if (judge==0)
				{
					electronsdown2[i][0]=x1;
					for (int j=0;j<s*s/2;j=j+1)
					{
						int x=electronsup2[j][0]-electronsdown2[i][0];
						int y=electronsup2[j][1]-electronsdown2[i][1];
						slater2[j][i]=element[x+(s-1)][y+(s-1)];
					}
					for (int j=0;j<s*s/2;j=j+1)
					{
						if (electronsup[j][0]==x&&electronsup[j][1]==y)
							numd2=numd2-1;
						if (electronsup[j][0]==x1&&electronsup[j][1]==y)
							numd2=numd2+1;
					}
					for (int j=0;j<s*s/2;j=j+1)
					{
						ratio+=slater2[j][i]*II[i][j];
					}
					Ek=Ek-t*ratio*pow(g,numd2-numd);
					numd2=numd;
					ratio=0;
					for (int k=0;k<s*s/2;k=k+1)
						slater2[k][i]=slater[k][i];
					electronsdown2[i][0]=x;
				}




				judge=0;
				for (int j=0;j<s*s/2;j=j+1)
					if (electronsdown[j][0]==x2&&electronsdown[j][1]==y)
						judge=1;

				if (judge==0)
				{
					electronsdown2[i][0]=x2;
					for (int j=0;j<s*s/2;j=j+1)
					{
						int x=electronsup2[j][0]-electronsdown2[i][0];
						int y=electronsup2[j][1]-electronsdown2[i][1];
						slater2[j][i]=element[x+(s-1)][y+(s-1)];
					}
					for (int j=0;j<s*s/2;j=j+1)
					{
						if (electronsup[j][0]==x&&electronsup[j][1]==y)
							numd2=numd2-1;
						if (electronsup[j][0]==x2&&electronsup[j][1]==y)
							numd2=numd2+1;
					}
					for (int j=0;j<s*s/2;j=j+1)
					{
						ratio+=slater2[j][i]*II[i][j];
					}
					Ek=Ek-t*ratio*pow(g,numd2-numd);
					numd2=numd;
					ratio=0;
					for (int k=0;k<s*s/2;k=k+1)
						slater2[k][i]=slater[k][i];
					electronsdown2[i][0]=x;
				}




				judge=0;
				for (int j=0;j<s*s/2;j=j+1)
					if (electronsdown[j][0]==x&&electronsdown[j][1]==y1)
						judge=1;

				if (judge==0)
				{
					electronsdown2[i][1]=y1;
					for (int j=0;j<s*s/2;j=j+1)
					{
						int x=electronsup2[j][0]-electronsdown2[i][0];
						int y=electronsup2[j][1]-electronsdown2[i][1];
						slater2[j][i]=element[x+(s-1)][y+(s-1)];
					}
					for (int j=0;j<s*s/2;j=j+1)
					{
						if (electronsup[j][0]==x&&electronsup[j][1]==y)
							numd2=numd2-1;
						if (electronsup[j][0]==x&&electronsup[j][1]==y1)
							numd2=numd2+1;
					}
					for (int j=0;j<s*s/2;j=j+1)
					{
						ratio+=slater2[j][i]*II[i][j];
					}
					if (y1==0)
						Ek=Ek+t*ratio*pow(g,numd2-numd);
					else
						Ek=Ek-t*ratio*pow(g,numd2-numd);
					numd2=numd;
					ratio=0;
					for (int k=0;k<s*s/2;k=k+1)
						slater2[k][i]=slater[k][i];
					electronsdown2[i][1]=y;
				}
				judge=0;
				for (int j=0;j<s*s/2;j=j+1)
					if (electronsdown[j][0]==x&&electronsdown[j][1]==y2)
						judge=1;

				if (judge==0)
				{
					electronsdown2[i][1]=y2;
					for (int j=0;j<s*s/2;j=j+1)
					{
						int x=electronsup2[j][0]-electronsdown2[i][0];
						int y=electronsup2[j][1]-electronsdown2[i][1];
						slater2[j][i]=element[x+(s-1)][y+(s-1)];
					}
					for (int j=0;j<s*s/2;j=j+1)
					{
						if (electronsup[j][0]==x&&electronsup[j][1]==y)
							numd2=numd2-1;
						if (electronsup[j][0]==x&&electronsup[j][1]==y2)
							numd2=numd2+1;
					}
					for (int j=0;j<s*s/2;j=j+1)
					{
						ratio+=slater2[j][i]*II[i][j];
					}
					if (y2==s-1)
						Ek=Ek+t*ratio*pow(g,numd2-numd);
					else
						Ek=Ek-t*ratio*pow(g,numd2-numd);
					numd2=numd;
					ratio=0;
					for (int k=0;k<s*s/2;k=k+1)
						slater2[k][i]=slater[k][i];
					electronsdown2[i][1]=y;
				}
				judge=0;
			}
			Etot=Etot+(Ek+u*numd)/(s*s);
			E2tot=E2tot+pow((Ek+u*numd)/(s*s),2);
		}
	}

	avgnumbind=numbindtot/sample;
	errnumbind=pow((numbind2tot/sample-pow(avgnumbind,2))/sample,0.5);

	avgnumd=numdtot/sample;
	errnumd=pow((numd2tot/sample-pow(avgnumd,2))/sample,0.5);

	avgE=Etot/sample;
	errE=pow((E2tot/sample-pow(avgE,2))/sample,0.5);
	cout<<" doublon number="<<avgnumd<<" err="<<errnumd<<endl;
	cout<<" energy="<<avgE<<" err="<<errE<<endl;
	cout<<" binding number="<<avgnumbind<<" err="<<errnumbind<<endl;
	end_time=clock();
	tot_time=(float)(end_time-start_time)/CLOCKS_PER_SEC;
	cout<<" time="<<tot_time<<"sec";
	system("pause");
	return 0;
}


double determinant(double a[s*s/2][s*s/2])
{      
	double x[s*s/2][s*s/2];
	for (int i=0;i<s*s/2;i=i+1)
		for (int j=0;j<s*s/2;j=j+1)
			x[i][j]=a[i][j];
	int k=0,n=0;
	while (k<s*s/2)
	{
		int l=k;
		while (abs(x[l][k])<0.00001)
		{
			l=l+1;
			if (l==s*s/2)
				return 0;
		}
		if (l!=k)
		{
			n=n+1;
			for (int i=0;i<s*s/2;i=i+1)
			{
				double b;
				b=x[k][i];
				x[k][i]=x[l][i];
				x[l][i]=b;          
			}
		}

		for (int i=k+1;i<s*s/2;i=i+1)
		{
			double r=x[i][k]/x[k][k];
			for (int j=k;j<s*s/2;j=j+1)
			{
				x[i][j]=x[i][j]-r*x[k][j];
			}
		}
		k=k+1;
	}
	double det=1;
	for (int i=0;i<s*s/2;i=i+1)
	{det=det*x[i][i];
	}
	return det*pow(-1.0,n);
}


void createInverse(double c[s*s/2][s*s/2],double b[s*s/2][s*s/2],int n)
{
	double a[n][n]; 
	for (int i=0;i<n;i=i+1)
		for (int j=0;j<n;j=j+1)
			a[i][j]=c[i][j];



	double Inv[n][n];
	for (int i=0;i<n;i=i+1)
		for (int j=0;j<n;j=j+1)
		{
			if(i==j)
				Inv[i][j]=1;
			else
				Inv[i][j]=0;
		}
	for (int i=0;i<n;i=i+1)
	{
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
