# include <stdio.h>
# include <mpi.h>
#include <stdlib.h>

int const name_size = MPI_MAX_PROCESSOR_NAME;

int main (int argc, char *argv[]) {
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm shmcomm;

    int rank, size, namelen, shm_size;
    int verbose = 0;
    char name[name_size];

    MPI_Init( &argc, &argv );

    if (getenv("MSPARSM_VERBOSE")) verbose = 1;

    MPI_Get_processor_name( name, &namelen );

    // MPI_COMM_TYPE_SHARED: This type splits the communicator into subcommunicators, each of which can create a shared memory region.
    MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
    MPI_Comm_size (shmcomm, &shm_size);

    MPI_Comm_rank( comm, &rank );
    MPI_Comm_size( comm, &size );

    if ( verbose ) printf("Proceso [%d]/[%d], shmcom [%d], running on [%s]\n", rank, size, shm_size, name);


    MPI_Finalize();

    return 0;
}

