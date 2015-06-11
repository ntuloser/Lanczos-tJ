#include<iostream>
#include<mpi.h>
#include <unistd.h>


int main(int argc, char** argv ){
	int mynode, numnodes;
	int sum, startval, endval, accum;

	MPI::Status status;
	MPI::Init( argc, argv );
	numnodes = MPI::COMM_WORLD.Get_size();
	mynode   = MPI::COMM_WORLD.Get_rank();


	startval = 1000*mynode/numnodes+1;
	endval   = 1000*(mynode+1)/numnodes;

	sum=0;


	for( int i=startval; i<=endval; ++i ){
		sum+=i;
		usleep(10000); //pass in microseconds
	}
	if( mynode != 0 ){
		MPI::COMM_WORLD.Send( &sum, 1, MPI::INT, 0, 1 );
	}
	else{
		for( int j=1; j<numnodes; ++j ){
			MPI::COMM_WORLD.Recv( &accum, 1, MPI::INT, j, 1, status );
			sum = sum + accum;
		}
	}


	if( mynode == 0 )
	{std::cout << "The sum from 0 to 1000 is: " << sum << std::endl;}
	MPI::Finalize();

	return 0;
}
