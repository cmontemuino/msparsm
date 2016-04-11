# include <stdio.h>
# include <mpi.h>
#include <stdlib.h>

int const name_size = MPI_MAX_PROCESSOR_NAME;

const int LOCAL_RESULT_TAG = 100;
const int GLOBAL_RESULT_TAG = 200;

/* defines global rank  -> shmcomm rank mapping;
    output: partners_map is array of ranks in shmcomm  */
void translate_ranks(MPI_Comm shmcomm, int partners[], int partners_map[])
{
    MPI_Group world_group, shared_group;
    int world_size;

    /* create MPI groups for global communicator and shm communicator */
    MPI_Comm_group (MPI_COMM_WORLD, &world_group);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_group (shmcomm, &shared_group);

    MPI_Group_translate_ranks (world_group, world_size, partners, shared_group, partners_map);
}

int main (int argc, char *argv[]) {
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm shmcomm;

    int global_rank, global_size, namelen;
    int i;
    int shm_rank, shm_size;
    int verbose = 0;
    char name[name_size];

    int global_result;

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( comm, &global_rank );
    MPI_Comm_size( comm, &global_size );

    if (getenv("MSPARSM_VERBOSE")) verbose = 1;

    MPI_Get_processor_name( name, &namelen );

    // MPI_COMM_TYPE_SHARED: This type splits the communicator into subcommunicators, each of which can create a shared memory region.
    MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
    MPI_Comm_size (shmcomm, &shm_size);

    MPI_Comm_rank( shmcomm, &shm_rank );
    MPI_Comm_size( shmcomm, &shm_size );

    int partners[global_size];
    int partners_map[global_size];
    for (i = 0; i < global_size; i++)
    {
        partners[i] = i;
    }

    translate_ranks(shmcomm, partners, partners_map);

    MPI_Barrier (comm);

    if ( global_rank == 0) {
        printf("Global master: [%d/%d] running on [%s]\n\n", global_rank, global_size, name);

//        for (i = 0; i < global_size; i++)
//        {
//            MPI_Recv(&global_result, 1, MPI_INT, i, GLOBAL_RESULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//            printf("Global master: received [%d] from process [%d]\n", global_result, i);
//        }

        for (i = 0; i < global_size; i++)
        {
            printf("Translated: Global[%d] -> Local[%d]\n", partners[i], partners_map[i]);
        }
    }

    if ( shm_rank == 0) {
        printf("\nLocal master: [%d/%d] running on [%s]\n", shm_rank, shm_size, name);
//        MPI_Send(&shm_rank, 1, MPI_INT, global_rank, GLOBAL_RESULT_TAG, MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}

