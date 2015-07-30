//extern int sss;
//extern int ddd;

class site {
	private:
    
	public:
		int idx_up;
		int idx_down;
		int state;
};


class config{
private:
    
public:
    int Sign;
    int Size;
    int Dope;
    int num_ele;
    
    //vector< vector< site > > squarelattice_not_ptr;
    //vector< vector< site > >* squarelattice;

    site*** squarelattice;
    int ** electronsup;
    int ** electronsdown;
    
    
    ////////////
    //function//
    ////////////
public:
    ////default_constructor
    //config();
    
    ////Constructor
    config(int s=SIZE, int d=DELTA );
    ////Destructor
    ~config();
    config(const config& a);
    /////////////////
    ////Operator
    config& operator=(const config& a);
    

    

    
    
    
    
    //////////////////////////
    //  random  initial   ////
    //////////////////////////
    void rand_init(){
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
           /* if (avoid_doublon == 1) {
                
                for (int i=0; i<num_ele/2; i++) {
                    if (electronsup[i][0]==x_ran && electronsup[i][1]==y_ran) {
                        filled=1;
                    }
                }
            }*/
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
				squarelattice[i][j]->state=0;
				squarelattice[i][j]->idx_up=-1;
				squarelattice[i][j]->idx_down=-1;
			}
		}
        
		for (int i=0;i<num_ele/2;i++){
            
			int x=electronsup[i][0];
			int y=electronsup[i][1];
			squarelattice[x][y]->state+=1;
			squarelattice[x][y]->idx_up=i;
            
			x=electronsdown[i][0];
			y=electronsdown[i][1];
			squarelattice[x][y]->state+=2;
			squarelattice[x][y]->idx_down=i;
		}
        

    }
    
    
    void rand_init_no_d(){
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
            
            for (int i=0; i<num_ele/2; i++) {
                    if (electronsup[i][0]==x_ran && electronsup[i][1]==y_ran) {
                            filled=1;
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
				squarelattice[i][j]->state=0;
				squarelattice[i][j]->idx_up=-1;
				squarelattice[i][j]->idx_down=-1;
			}
		}

        
		for (int i=0;i<num_ele/2;i++){
            
			int x=electronsup[i][0];
			int y=electronsup[i][1];
			squarelattice[x][y]->state+=1;
			squarelattice[x][y]->idx_up=i;
            
			x=electronsdown[i][0];
			y=electronsdown[i][1];
			squarelattice[x][y]->state+=2;
			squarelattice[x][y]->idx_down=i;
		}
        
        
    }

    
    
    
    
    
    
    void hop_up(config* beta_ptr, int idx, int move, int* f_ptr, int* o_ptr){

        this->copy_config_to(beta_ptr);
        
        //config* beta_ptr = (config*)malloc(sizeof(config));
        //memcpy(beta_ptr, this, sizeof(config));
        
        int move_x;
        int move_y;
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
        
        int x = beta_ptr->electronsup[idx][0];
        int y = beta_ptr->electronsup[idx][1];
        int x_new = (x+move_x+Size)%Size;
        int y_new = (y+move_y+Size)%Size;
        
        if ((y+move_y)==Size || (y+move_y)==-1  ) {
            (*o_ptr)=1;
        }
        else {
            (*o_ptr)=0;
        }

        
        if (beta_ptr->squarelattice[x_new][y_new]->idx_up==-1) {
            
            beta_ptr->electronsup[idx][0]=x_new;
            beta_ptr->electronsup[idx][1]=y_new;
            
            beta_ptr->squarelattice[x][y]->state  -= 1;
            beta_ptr->squarelattice[x][y]->idx_up = -1;
            beta_ptr->squarelattice[x_new][y_new]->state += 1;
            beta_ptr->squarelattice[x_new][y_new]->idx_up = idx;
            (*f_ptr)=1;
        }
        
        else
            (*f_ptr)=0;
        
    };
    
    
    
    
    
    void hop_down(config* beta_ptr,int idx, int move,int* f_ptr, int* o_ptr){
        
        this->copy_config_to(beta_ptr);
        //config* beta_ptr = (config*)malloc(sizeof(config));
        //memcpy(beta_ptr,this,sizeof(config));
        
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
        
        int x = beta_ptr->electronsdown[idx][0];
        int y = beta_ptr->electronsdown[idx][1];
        int x_new = (x+move_x+Size)%Size;
        int y_new = (y+move_y+Size)%Size;
        
        
        if ((y+move_y)==Size || (y+move_y)==-1  ) {
            (*o_ptr)=1;
        }
        else {
            (*o_ptr)=0;
        }
        
        
        if (beta_ptr->squarelattice[x_new][y_new]->idx_down==-1) {
            
            beta_ptr->electronsdown[idx][0]=x_new;
            beta_ptr->electronsdown[idx][1]=y_new;
            
            beta_ptr->squarelattice[x][y]->state -=  2;
            beta_ptr->squarelattice[x][y]->idx_down = -1;
            beta_ptr->squarelattice[x_new][y_new]->state += 2;
            beta_ptr->squarelattice[x_new][y_new]->idx_down = idx;
             /**/
            (*f_ptr)=1;
        
        }
        else
            (*f_ptr)=0;
        
    };
    
    
    
    
    
    
    int swap(config* beta_ptr, const int& x,const int& y, const int& move,int* f_ptr, int* o_ptr){
        
        int tJ=0;
        this->copy_config_to(beta_ptr);
        //config* beta_ptr = (config*)malloc(sizeof(config));
        //memcpy(beta_ptr,this,sizeof(config));
        
        
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

        
        int x_new = (x+move_x+Size)%Size;
        int y_new = (y+move_y+Size)%Size;
        
        if ((y+move_y)==Size || (y+move_y)==-1  ) {
            (*o_ptr)=1;
        }
        else {
            (*o_ptr)=0;
        }
        
        
        int temp_state = beta_ptr->squarelattice[x][y] -> state;
        int temp_idx_up = beta_ptr->squarelattice[x][y] -> idx_up;
        int temp_idx_down = beta_ptr->squarelattice[x][y] -> idx_down;
        
        /////////////////////////
        //// PAULI EXCLUSION ////
        /////////////////////////
        if (beta_ptr->squarelattice[x_new][y_new]->state == temp_state ) {
            (*f_ptr)=0;
            if (temp_state==0) {
                tJ=0;
                return tJ;
            }
            else if(temp_state==1 || temp_state==2 ){
                tJ=3;
                return tJ;
            }
        }
        /////////////////////////
        /////////////////////////
        
        
        beta_ptr->squarelattice[x][y]->state = beta_ptr->squarelattice[x_new][y_new]->state;
        beta_ptr->squarelattice[x][y]->idx_up = beta_ptr->squarelattice[x_new][y_new]->idx_up;
        beta_ptr->squarelattice[x][y]->idx_down = beta_ptr->squarelattice[x_new][y_new]->idx_down;
        
        beta_ptr->squarelattice[x_new][y_new]->state = temp_state;
        beta_ptr->squarelattice[x_new][y_new]->idx_up = temp_idx_up;
        beta_ptr->squarelattice[x_new][y_new]->idx_down = temp_idx_down;
        
        
        if (temp_idx_up != -1) {
            beta_ptr->electronsup[temp_idx_up][0]=x_new;
            beta_ptr->electronsup[temp_idx_up][1]=y_new;
            tJ+=1;
        }
        else if (temp_idx_down != -1) {
            beta_ptr->electronsdown[temp_idx_down][0]=x_new;
            beta_ptr->electronsdown[temp_idx_down][1]=y_new;
            tJ+=1;
        }
        
        temp_idx_up = beta_ptr->squarelattice[x][y] -> idx_up;
        temp_idx_down = beta_ptr->squarelattice[x][y] -> idx_down;
        
        if (temp_idx_up != -1) {
            beta_ptr->electronsup[temp_idx_up][0]=x;
            beta_ptr->electronsup[temp_idx_up][1]=y;
            tJ+=1;
        }
        else if (temp_idx_down != -1) {
            beta_ptr->electronsdown[temp_idx_down][0]=x;
            beta_ptr->electronsdown[temp_idx_down][1]=y;
            tJ+=1;
        }
        
        if (tJ==2) {
            // in the case of superexchange
            beta_ptr->Sign *= (-1);
        }
        (*f_ptr)=1;
        return tJ;
        
    };
    
    
    
    //Counting the number of doublon in the config.
    int num_doublon(){
        
        int count=0;
        for (int i=0; i<Size; i++) {
            for (int j=0; j<Size; j++) {
                if(this->squarelattice[i][j]->state == 3){
                    count+=1;
                }
            }
        }
        
        return count;
    };
    
    
    
    
    
    void printconfig(){
        
        printf("electronsup:\n");
        for (int idx=0; idx<num_ele/2; idx++) {
            printf("(%d, %d) ",electronsup[idx][0],electronsup[idx][1]);
        }
        printf("\n");
        printf("electronsdown:\n");
        for (int idx=0; idx<num_ele/2; idx++) {
            printf("(%d, %d) ",electronsdown[idx][0],electronsdown[idx][1]);
        }
        
        
        std::cout<<"\n";
        for (int i=0; i<Size; i++) {
            for (int j=0; j<Size; j++) {
                std::cout<<squarelattice[i][j]->state<<", ";
            }
            std::cout<<"\n";
        }
        //std::cout<<"\n";
    };
    
    
    void copy_config_to( config* beta_ptr){
        
        //We require a beta_ptr which has been constructed( memory allocated);

        //beta_ptr->Size = this->Size;
        //beta_ptr->Dope = this->Dope;
        //beta_ptr->num_ele = this->num_ele;
        
        
        
        beta_ptr->Sign = this->Sign;
        site* b_site=beta_ptr->squarelattice[0][0];
        site* this_site=this->squarelattice[0][0];
        memcpy ( b_site, this_site, Size*Size*sizeof(site) );
        memcpy ( &(beta_ptr->electronsup[0][0]), &(this->electronsup[0][0]), num_ele*sizeof(int) );
        memcpy ( &(beta_ptr->electronsdown[0][0]), &(this->electronsdown[0][0]), num_ele*sizeof(int) );
        


        /*
        for (int i=0; i<Size; i++) {
            for (int j=0; j<Size; j++) {
                
                site* b_site=beta_ptr->squarelattice[i][j];
                site* this_site=this->squarelattice[i][j];
                
                b_site->state = this_site->state;
                b_site->idx_up = this_site->idx_up;
                b_site->idx_down = this_site->idx_down;
                
            }
        }
        for (int i=0; i<num_ele/2; i++) {
            for (int j=0; j<2; j++) {
                beta_ptr->electronsup[i][j] = this->electronsup[i][j];
                beta_ptr->electronsdown[i][j] = this->electronsdown[i][j];
            }
        }*/
        
    }
    
    
    void set_slater( double** a_ij,double** slater ){
        
        //slater:  [num_ele/2 ][num_ele/2 ]
        //a_ij  :  [(2*Size-1)][(2*Size-1)]
        for (int i=0;i<num_ele/2;i++){
			for (int j=0;j<num_ele/2;j++){
				// The ith electronsup, and the jth electronsdown.
				// Call the matrix a_ij by displacement vector
				// between two electron.
				//
				// We shift the displacement so that it is positive
				// to call matrix element.
				int x=electronsup[i][0]-electronsdown[j][0];
				int y=electronsup[i][1]-electronsdown[j][1];
                
                int t1=x+(Size-1);
                int t2=y+(Size-1);
				slater[i][j]=a_ij[t1][t2];
               
                //slater[i][j]=*(*(a_ij + (sizeof(double *)*t1) ) +sizeof(double)*t2  )   ;
			}
		}

    };

    
    
};





config::config(int s, int d ){
    Sign=1;
    Size=s;
    num_ele=s*s-d;
    
    

    /*
    int i;
    site ** siteptrptr;
    squarelattice= (site ***)malloc(Size *sizeof(site **));
    siteptrptr = (site  **)malloc(Size*Size*sizeof(site *));
    for (i=0; i<Size ; i++, siteptrptr+=Size )
        squarelattice[i]=siteptrptr;
    */
    
    
    int i, j;
    site ** siteptrptr;
    site *  siteptr;
    squarelattice = ( site *** )malloc(Size*sizeof(site**)+Size*Size*sizeof(site*)+Size*Size*sizeof(site) );
    for (i=0, siteptrptr=(site **)(squarelattice+Size), siteptr=(site *)(siteptrptr+Size*Size); i<Size; i++, siteptrptr += Size, siteptr+= Size ) {
        squarelattice[i]=siteptrptr;
        for (j=0 ; j<Size; j++) {
            squarelattice[i][j]=(site*)(siteptr+j);
        }
        
    }

    int *intptr;
    electronsup = (int**) malloc(num_ele/2 * sizeof(int*)+ num_ele*sizeof(int) );
    for(i=0, intptr=(int*)(electronsup+num_ele/2); i<num_ele/2; i++, intptr+= 2){
        electronsup[i]=intptr;
    }
    electronsdown = (int**) malloc(num_ele/2 * sizeof(int*)+ num_ele*sizeof(int) );
    for(i=0, intptr=(int*)(electronsdown+num_ele/2); i<num_ele/2; i++, intptr+= 2){
        electronsdown[i]=intptr;
    }

    
    /*
    squarelattice = new site**[Size];
    for (int i=0; i<Size; i++) {
        squarelattice[i]= new site*[Size];
        for (int j=0; j<Size; j++) {
            squarelattice[i][j]=new site;
        }
    }
    
    electronsup = new int*[num_ele/2];
    for (int i=0; i<num_ele/2; i++) {
        electronsup[i] = new int[2];
    }
    electronsdown = new int*[num_ele/2];
    for (int i=0; i<num_ele/2; i++) {
        electronsdown[i] = new int[2];
    }*/
    
}


//copy constructor
config::config(const config& a){
    
    Sign = a.Sign;
    Size = a.Size;
    Dope = a.Dope;
    num_ele = a.num_ele;
    
    squarelattice = new site**[Size];
    for (int i=0; i<Size; i++) {
        squarelattice[i]= new site*[Size];
        for (int j=0; j<Size; j++) {
            squarelattice[i][j]=new site;
        }
    }
    
    electronsup = new int*[num_ele/2];
    for (int i=0; i<num_ele/2; i++) {
        electronsup[i] = new int[2];
    }
    electronsdown = new int*[num_ele/2];
    for (int i=0; i<num_ele/2; i++) {
        electronsdown[i] = new int[2];
    }

    
    for (int i=0; i<Size; i++) {
        for (int j=0; j<Size; j++) {
            squarelattice[i][j]->state = (a.squarelattice[i][j]->state);
            squarelattice[i][j]->idx_up = (a.squarelattice[i][j]->idx_up);
            squarelattice[i][j]->idx_down = (a.squarelattice[i][j]->idx_down);
        }
    }
    for (int i=0; i<num_ele/2; i++) {
        for (int j=0; j<2; j++) {
            electronsup[i][j] = a.electronsup[i][j];
            electronsdown[i][j] = a.electronsdown[i][j];
        }
    }

}

//destructor
config::~config(){
    
    free(squarelattice);
    free(electronsup);
    free(electronsdown);

    
    /*
    for (int i=0; i<Size; i++) {
        for (int j=0; j<Size; j++) {
            delete squarelattice[i][j];
        }
        delete [] squarelattice[i];
    }
    delete [] squarelattice;
    
    for (int i=0; i<num_ele/2; i++) {
        delete [] electronsup[i];
    }
    delete [] electronsup;

    for (int i=0; i<num_ele/2; i++) {
        delete [] electronsdown[i];
    }
    delete [] electronsdown;*/
    
}

config& config::operator=(const config& a){
    if (this != &a) {//avoid self assignment
        //destruct the original data
        for (int i=0; i<Size; i++) {
            for (int j=0; j<Size; j++) {
                delete squarelattice[i][j];
            }
            delete [] squarelattice[i];
        }
        delete [] squarelattice;
        
        for (int i=0; i<num_ele/2; i++) {
            delete [] electronsup[i];
        }
        delete [] electronsup;
        
        for (int i=0; i<num_ele/2; i++) {
            delete [] electronsdown[i];
        }
        delete [] electronsdown;
        // set new data from a
        
        Sign = a.Sign;
        Size = a.Size;
        Dope = a.Dope;
        num_ele = a.num_ele;
        
        squarelattice = new site**[Size];
        for (int i=0; i<Size; i++) {
            squarelattice[i]= new site*[Size];
            for (int j=0; j<Size; j++) {
                squarelattice[i][j]=new site;
            }
        }
        
        electronsup = new int*[num_ele/2];
        for (int i=0; i<num_ele/2; i++) {
            electronsup[i] = new int[2];
        }
        electronsdown = new int*[num_ele/2];
        for (int i=0; i<num_ele/2; i++) {
            electronsdown[i] = new int[2];
        }
        
        
        for (int i=0; i<Size; i++) {
            for (int j=0; j<Size; j++) {
                squarelattice[i][j]->state = (a.squarelattice[i][j]->state);
                squarelattice[i][j]->idx_up = (a.squarelattice[i][j]->idx_up);
                squarelattice[i][j]->idx_down = (a.squarelattice[i][j]->idx_down);
            }
        }
        for (int i=0; i<num_ele/2; i++) {
            for (int j=0; j<2; j++) {
                electronsup[i][j] = a.electronsup[i][j];
                electronsdown[i][j] = a.electronsdown[i][j];
            }
        }

    }
    return *this;
}


