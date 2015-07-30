#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

/************************************************************
  This is a simple broadcast program in MPI
 ************************************************************/

int main(argc,argv)
	int argc;
	char *argv[];
{
	printf("Hello world\n");
	int i,myid, numprocs;
	int source,count;
	int buffer[6]={0,0,0,0,0,0};
	MPI_Status status;
	MPI_Request request;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	source=0;
	count=5;
	if(myid == source){
		for(i=0;i<count;i++)
			buffer[i]=i;
	}

	//int MPI_Bcast( void *buffer, int count, MPI_Datatype datatype, int root,  MPI_Comm comm )
	//count number of entries in buffer (integer)
	//datatype data type of buffer (handle)
	//root rank of broadcast root (integer)
	//comm communicator (handle)


	MPI_Bcast(buffer,3,MPI_INT,source,MPI_COMM_WORLD);
	for(i=0;i<6;i++)
		printf("%d ",buffer[i]);
	printf("\n");
	MPI_Finalize();
}

