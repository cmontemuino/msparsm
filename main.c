# include <stdio.h>
# include <mpi.h>
#include <stdlib.h>

int const name_size = MPI_MAX_PROCESSOR_NAME;

int main (int argc, char *argv[]) {
    int rank, size, namelen;
    int verbose = 0;
    char name[name_size];

    MPI_Init( &argc, &argv );

    if (getenv("MSPARSM_VERBOSE")) verbose = 1;

    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    MPI_Get_processor_name( name, &namelen );

    if ( verbose ) printf("Proceso %d de %d procesos, nombre=%s\n", rank, size, name);

    MPI_Finalize();

    return 0;
}

