/* Helloword.c */
#include <stdio.h>
#include <mpi.h>
int main(int argc, char* argv[]) {
	int my_rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	printf( "Hello from process %d\n", my_rank); 
	MPI_Finalize();
	return 0;
}
