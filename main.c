# include <stdio.h>
# include <mpi.h>
#include <stdlib.h>

int const name_size = MPI_MAX_PROCESSOR_NAME;

void translate_ranks(MPI_Comm shmcomm, int rank_size, int global_ranks[], int shm_ranks[])
{
    MPI_Group world_group, shared_group;

    /* create MPI groups for global communicator and shm communicator */
    MPI_Comm_group (MPI_COMM_WORLD, &world_group);
    MPI_Comm_group (shmcomm, &shared_group);

    MPI_Group_translate_ranks (world_group, rank_size, global_ranks, shared_group, shm_ranks);
}

void print_ranks (int rank, int rank_size, int partners[], int partners_map[],
                       int n_node_partners, int n_inter_partners)
{
    int j, partner;
    char tmp_str_intra[rank_size*16];
    char tmp_str_inter[rank_size*10];
    int pos_intra = 0, pos_inter = 0;

    for (j=0; j<rank_size; j++)
    {
        partner = partners[j]; /* partner is in the world notation */
        if (partners_map[j] != MPI_UNDEFINED) /* partner j is on the same node  */
            pos_intra += sprintf (&tmp_str_intra[pos_intra], ", %d (%d)", partner, partners_map[j]);
        else
            pos_inter += sprintf (&tmp_str_inter[pos_inter], "jode, %d", partner);
    }

    if (n_inter_partners)
        printf ("i'm rank %d with %d internode partner%c%s \n",
                rank, n_inter_partners, n_inter_partners >1?'s':' ', tmp_str_inter);

    if (n_node_partners)
        printf ("i'm rank %d with %d intranode partner%c%s\n",
                rank, n_node_partners, n_node_partners > 1?'s':' ', tmp_str_intra);
}

int main (int argc, char *argv[]) {
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm shmcomm;

    int global_rank, global_size, namelen;
    int shm_rank, shm_size;
    int verbose = 0;
    char name[name_size];

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( comm, &global_rank );
    MPI_Comm_size( comm, &shm_rank );

    if (getenv("MSPARSM_VERBOSE")) verbose = 1;

    MPI_Get_processor_name( name, &namelen );

    // MPI_COMM_TYPE_SHARED: This type splits the communicator into subcommunicators, each of which can create a shared memory region.
    MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
    MPI_Comm_size (shmcomm, &shm_size);

    MPI_Comm_rank( shmcomm, &shm_rank );
    MPI_Comm_size( shmcomm, &shm_size );

    MPI_Barrier (comm);

    if ( global_rank == 0) {
        printf("Global master: [%d/%d] running on [%s]\n", global_rank, global_size, name);
    }

    if ( shm_rank == 0) {
        printf("Local master: [%d/%d] running on [%s]\n", shm_rank, shm_size, name);
    }

    int local_ranks[shm_size];
    int *global_ranks = (int*)malloc(shm_size*sizeof(int));
    local_ranks[0] = 0;
    local_ranks[1] = 1;
    local_ranks[2] = 2;
    local_ranks[3] = 3;
    translate_ranks(shmcomm, shm_size, local_ranks, global_ranks);

    int local_master_rank=0, global_master_rank=0;
    if (verbose) print_ranks (global_rank, shm_size, local_ranks, global_ranks, 1, 1);

    //if ( verbose ) printf("Proceso [%d]/[%d], shmcom [%d], running on [%s]\n", shm_rank, local_size, shm_size, name);


    MPI_Finalize();

    return 0;
}

