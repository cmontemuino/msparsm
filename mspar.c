#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ms.h"
#include "mspar.h"
//#include <mpi.h> /* OpenMPI library */

const int SAMPLES_NUMBER_TAG = 200;
const int RESULTS_TAG = 300;
const int GO_TO_WORK_TAG = 400;
const int NODE_MASTER_ASSIGNMENT = 500; // Used to assign a number of samples to a node master
const int ACK_TAG = 600; // Used by workers in the master node

typedef struct {
    int shm_rank;
    int world_rank;
} Rank_Struct;

// Following variables are with global scope in order to facilitate its sharing among routines.
// They are going to be updated in the masterWorkerSetup routine only, which is called only one, therefore there is no
// risk of race conditions or whatever other concurrency related problem.
MPI_Comm shmcomm; // shm intra-communicator
int world_rank, shm_rank;
int world_size, shm_size;
int shm_mode = 0;  // indicates whether MPI-3 SHM is going to be applied

// **************************************  //
// MASTER
// **************************************  //

int
masterWorkerSetup(int argc, char *argv[], int howmany, struct params parameters, unsigned int maxsites, int *excludeFrom, int *node_rank)
{
    // goToWork         : used by workers to realize if there is more work to do.
    // seedMatrix       : matrix containing the RNG seeds to be distributed to working processes.
    // localSeedMatrix  : matrix used by workers to receive RNG seeds from master.
    unsigned short *seedMatrix;
    unsigned short localSeedMatrix[3];

    // MPI-3 SHM related
    MPI_Win win;      // shm window object
    char *shm; // the shared memory
    char *shm_results; // memory place where all MPI process from one node will going to share
    int nodeHowmany;

    // MPI Initialization
    MPI_Init(&argc, &argv );
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // MPI_COMM_TYPE_SHARED: This type splits the communicator into subcommunicators, each of which can create a shared memory region.
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
    MPI_Comm_size( shmcomm, &shm_size );
    MPI_Comm_rank( shmcomm, &shm_rank );

    if (shm_size != world_size) // there are MPI process in more than 1 computing node
    {
        shm_mode = 1;
    }

    *excludeFrom = shm_mode;
    *node_rank = shm_rank;

    if(world_rank == 0)
    {
        int i;
        // Only the master process prints out the application's parameters
        for(i=0; i<argc; i++) {
            fprintf(stdout, "%s ",argv[i]);
            fflush(stdout);
        }
        // If there are (not likely) more processes than samples, then the process pull
        // is cut up to the number of samples. */
        if(world_size > howmany)
            world_size = howmany + 1; // the extra 1 is due to the master

        int nseeds;
        doInitializeRng(argc, argv, &nseeds, parameters);
        int dimension = nseeds * world_size;
        seedMatrix = (unsigned short *) malloc(sizeof(unsigned short) * dimension);
        for(i=0; i<dimension;i++)
            seedMatrix[i] = (unsigned short) (ran1()*100000);
    }

    // Filter out workers with rank higher than howmany, meaning there are more workers than samples to be generated.
    if(world_rank < howmany) {
        MPI_Scatter(seedMatrix, 3, MPI_UNSIGNED_SHORT, localSeedMatrix, 3, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);

        /*
        MPI_Barrier(MPI_COMM_WORLD);

        if ( shm_mode == 0)
        {
            parallelSeed(localSeedMatrix);
            MPI_Aint sz;
            int dispunit = 1;

            if (shm_rank != 0)
            {

                char *sample = generateSample(parameters, maxsites);
                int length = strlen(sample);

                MPI_Send(&length, 1, MPI_INT, 0, 10001, shmcomm);

                MPI_Win_allocate_shared(length, 1, MPI_INFO_NULL, shmcomm, &shm, &win);
                MPI_Win_shared_query(win, MPI_PROC_NULL, &sz, &dispunit, &shm_results);

                memcpy(shm_results, sample, length);
                free(sample);
                MPI_Win_fence(0, win);
                //MPI_Barrier(shmcomm);

                MPI_Win_free(&win);

            }
            else
            {
                int length;

                MPI_Status status;
                MPI_Recv(&length, 1, MPI_INT, 1, 10001, shmcomm, &status);
                //printf("[%d] - shm from [%d] ready for reading\n", shm_rank, status.MPI_SOURCE);


                MPI_Win_allocate_shared(0, 1, MPI_INFO_NULL, shmcomm, &shm, &win);
                MPI_Win_shared_query(win, MPI_PROC_NULL, &sz, &dispunit, &shm_results);

                MPI_Win_fence(0, win);
                //MPI_Barrier(shmcomm);
                printf("%s\n", shm_rank, shm_results);

                MPI_Win_free(&win);

            }
        }
        */

        Rank_Struct rank;
        MPI_Aint struct_rank_size, struct_rank_size_lb;
        MPI_Datatype struct_rank_type;

        buildRankDataType(&struct_rank_type);

        MPI_Type_get_extent(struct_rank_type, &struct_rank_size_lb, &struct_rank_size);

        rank.shm_rank = shm_rank;
        rank.world_rank = world_rank;

        void *ranks = malloc(struct_rank_size * (world_size));

        MPI_Gather (&rank, 1, struct_rank_type, ranks, 1, struct_rank_type, 0, MPI_COMM_WORLD);

        if (world_rank == 0) { // Global master
            int i, pendingNodeMasters = 0;

            // Distribute remaining samples.
            if (shm_mode) { // there is more than one node
                int node_count = numberOfNodes(ranks, struct_rank_size);

                // calculate how many samples are going to be distributed among all nodes
                nodeHowmany = howmany / node_count;
                int remainder = howmany % node_count;

                // Delegate samples on node where the global master resides to a secondary master node (world_rank = 1).
                int samples = nodeHowmany + remainder;
                MPI_Send(&samples, 1, MPI_INT, 1, NODE_MASTER_ASSIGNMENT, MPI_COMM_WORLD); // process 1 will output the results

                // Distribute samples among remaining nodes
                for (i = 1; i < world_size; ++i) {// don't include node with global master now
                    Rank_Struct *r;
                    r = ranks + struct_rank_size * i;
                    if (r->shm_rank == 0) {
                        MPI_Send(&nodeHowmany, 1, MPI_INT, r->world_rank, NODE_MASTER_ASSIGNMENT, MPI_COMM_WORLD);
                        pendingNodeMasters += 1;
                    }
                }

                // Receive results from master nodes
                int source;
                while (pendingNodeMasters) {
                    shm_results = readResults(MPI_COMM_WORLD, &source);
                    fprintf(stdout, "%s", shm_results);
                    fflush(stdout);
                    free(shm_results);
                    pendingNodeMasters -= 1;
                }
            } else { // There is only one node, hence a secondary master is not needed
                // Note: global master visits this branch if and only if there is one single node
                masterProcessingLogic(howmany, 0, parameters, maxsites);
            }
        } else {
            if (shm_rank == 0 || (shm_mode && world_rank == 1)) {
                // Note: global master never visits this branch
                MPI_Recv(&nodeHowmany, 1, MPI_INT, 0, NODE_MASTER_ASSIGNMENT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                masterProcessingLogic(nodeHowmany, 0, parameters, maxsites);

            }
            else
            {
                // Workers initialize the RGN
                parallelSeed(localSeedMatrix);
            }
        }

        MPI_Type_free(&struct_rank_type);
    }

    return world_rank;
}

/**
 * Calculates how many nodes are there taking the SHM ranks from all MPI processes.
 * This function is useful to determine the number of node masters.
 */
int
numberOfNodes(void *ranks, MPI_Aint rank_size)
{
    int i, result = 0;
    Rank_Struct *rank;
    for ( i = 0; i < world_size; i ++) {
        rank = ranks + rank_size * i;
        if (rank->shm_rank == 0)
            result += 1;
    }

    return result;
}

void
buildRankDataType(MPI_Datatype* result) {
    Rank_Struct rank;
    {
        int blocklen[] = {1, 1};
        MPI_Aint addr[3];
        MPI_Get_address(&rank, &addr[0]);
        MPI_Get_address(&rank.shm_rank, &addr[1]);
        MPI_Get_address(&rank.world_rank, &addr[2]);
        MPI_Aint displacements[] = {addr[1] - addr[0], addr[2] - addr[0]};
        MPI_Datatype types[2] = {MPI_INT, MPI_INT};
        MPI_Type_create_struct(2, blocklen, displacements, types, result);
        MPI_Type_commit(result);
    }
}

void
masterWorkerTeardown() {
    MPI_Finalize();
}

/*
 * Logic implemented by the master process.
 *
 * @param howmany total number of replicas to be generated
 * @param lastAssignedProcess last processes that has been assigned som work
 */
void
masterProcessingLogic(int howmany, int lastAssignedProcess, struct params parameters, unsigned int maxsites)
{
    int *processActivity = (int*) malloc(shm_size * sizeof(int));
    int node_offset = 0;
    if (world_rank == 1) { // it is the secondary master located at the main node (where the global master resides)
        node_offset = 1;
        processActivity[1] = 1;
    }

    if (howmany > 0) {
        processActivity[0] = 1; // Master initially does not generate replicas

        int i;
        for (i = node_offset+1; i < shm_size; i++) // secondary master initially does not generate replicas
            processActivity[i] = 0;

        int pendingJobs = howmany; // number of jobs already assigned but pending to be finalized

        char *results;
        char *sample;
        size_t offset, length;
        results = malloc(sizeof(char) * 1000);

        int idleProcess = 1;
        int samples = 1; // number of samples to be generated by workers
        while (howmany > 0 && idleProcess > 0) { // Send sample generation requests to all available workers
            idleProcess = findIdleProcess(processActivity, lastAssignedProcess, node_offset);
            if (idleProcess > 0) {
                assignWork(processActivity, idleProcess, samples);
                lastAssignedProcess = idleProcess;
                howmany--;
            }
        }
        while (howmany > 0) { // Collect previously assigned sample generation jobs
            if (world_rank > 1) { // it is not the main node, hence samples are received using MPI point-to-point
                sample = readResultsFromWorkers(1, processActivity);
                offset = strlen(results);
                length = strlen(sample);
                results = realloc(results, offset + length + 1);
                memcpy(results + offset, sample, length);
                free(sample);
                howmany-= samples; // in the readResultsFromWorkers we sent the goToWork, which means another sample is going to be generated
            } else { // In main node all workers directly outputs the samples, but...
                // ...we need to receive some ACK to verify a sample was indeed outputted
                int goToWork = 0; // when 1, then it does mean another sample was requested to a worker
                int additional_samples = readAckFromLocalWorker(howmany, processActivity, parameters, maxsites, &goToWork);

                // need to update counters if applicable
                pendingJobs -= additional_samples;
                howmany -= additional_samples + (goToWork * samples);
            }
            pendingJobs-= samples; // in either branch (at least) one sample was generated
        }

        while (pendingJobs > 0) { // collect any non yet received sample
            if (world_rank != 0 && shm_rank == 0) { // whatever node master, but the secondary master located at main node (i.e.: world_rank=  1)
                sample = readResultsFromWorkers(0, processActivity);
                offset = strlen(results);
                length = strlen(sample);
                results = realloc(results, offset + length + 1);
                memcpy(results + offset, sample, length);
                free(sample);
            } else {
                // we need to receive some ACK to verify a sample was outputted by a local worker process
                int goToWork = 0;
                readAckFromLocalWorker(0, processActivity, parameters, maxsites, &goToWork);
            }

            pendingJobs--;
        }


        if (world_rank == 0) { // global master, but there is only one single node, otherwise global master does not visit this function at all
            fprintf(stdout, "%s", results);
            fflush(stdout);
        } else { // node master in a multi-node scenario
            if (shm_rank == 0) // exclude the secondary master located at the main node (i.e. world_rank = 1)
                // Note: we do not want to inlude the secondary master here, as it is already sending the samples to the standard output (see readAckFromLocalWorker function)
                MPI_Send(results, strlen(results) + 1, MPI_CHAR, 0, RESULTS_TAG, MPI_COMM_WORLD);
        }
        free(results); // be good citizen
    }
}


int readAckFromLocalWorker(int remaining, int *workersActivity, struct params parameters, unsigned int maxsites, int *goToWork)
{
    int source;
    int additional_samples = 0;
    MPI_Status status;

    *goToWork = 0;

    MPI_Probe(MPI_ANY_SOURCE, ACK_TAG, shmcomm, &status);
//    while (!msg_avail) {
//        // This function is called by secondary master in main node only, hence it is safe to generate an extra sample
//        // and send it to the standard output.
//        MPI_Iprobe(MPI_ANY_SOURCE, ACK_TAG, shmcomm, &msg_avail, &status);
//        if (remaining) { // then generate sample while others reply
//            char *results = generateSamples(1, parameters, maxsites);
//            fprintf(stdout, "################%s", results);
//            fflush(stdout);
//            free(results);
//            remaining--;
//            additional_samples++;
//        }
//    }

    source = status.MPI_SOURCE;

    int ack;
    //MPI_Recv(&ack, 1, MPI_INT, source, ACK_TAG, shmcomm, MPI_STATUS_IGNORE);
    MPI_Request requests[2];
    MPI_Irecv(&ack, 1, MPI_INT, source, ACK_TAG, shmcomm, &requests[0]);

    if (remaining) { // then generate 1 sample while others reply
        char *results = generateSamples(1, parameters, maxsites);
        fprintf(stdout, "%s", results);
        fflush(stdout);
        free(results);
        remaining--;
        additional_samples++;
    }

    // Now let the worker know whether there are (still) any pending samples
    if (remaining > 0)
        *goToWork = 1;

    //MPI_Send(&goToWork, 1, MPI_INT, source, GO_TO_WORK_TAG, shmcomm);
    MPI_Isend(goToWork, 1, MPI_INT, source, GO_TO_WORK_TAG, shmcomm, &requests[1]);

    workersActivity[source]=0;

    MPI_Waitall(2, requests, MPI_STATUS_IGNORE);

    return additional_samples;
}

/*
 *
 * Esta función realiza dos tareas: por un lado hace que el Master escuche los resultados enviados por los workers y por
 * otro lado, se envía al worker que se ha recibido la información un mensaje sobre si debe seguir esperando por
 * trabajos o si ha de finalizar su contribución al sistema.
 *
 * @param goToWork indica si el worker queda en espera de más trabajo (1) o si ya puede finalizar su ejecución (0)
 * @param workersActivity el vector con el estado de actividad de los workers
 *
 * @return the generated sample
 */
char* readResultsFromWorkers(int goToWork, int* workersActivity)
{
    char *readResults(MPI_Comm, int*);

    int source;
    char *results = readResults(shmcomm, &source);

    MPI_Send(&goToWork, 1, MPI_INT, source, GO_TO_WORK_TAG, shmcomm);

    workersActivity[source]=0;
    return results;
}

char
*readResults(MPI_Comm comm, int *source)
{
    MPI_Status status;
    int size;

    MPI_Probe(MPI_ANY_SOURCE, RESULTS_TAG, comm, &status);
    MPI_Get_count(&status, MPI_CHAR, &size);
    *source = status.MPI_SOURCE;
    char *results = (char *) malloc(size*sizeof(char));

    MPI_Recv(results, size, MPI_CHAR, *source, RESULTS_TAG, comm, MPI_STATUS_IGNORE);

    return results;
}

/*
 * Finds an idle process from a list of potential worker processes.
 *
 * @param workersActivity status of all processes that can generate some work (0=idle; 1=busy)
 * @lastAssignedProcess last process assigned with some work
 *
 * @return idle process index or -1 if all processes are busy.
 */
int findIdleProcess(int *processActivity, int lastAssignedProcess, int node_offset) {
    /*
     * Implementation note: lastAssignedProcess is used to implement a fairness policy in which every available process
     * can be assigned with some work.
     */
    int result = -1;
    int i= lastAssignedProcess+1; // master process does not generate replicas
    while(i < shm_size && processActivity[i] == 1){
        i++;
    };

    if(i >= shm_size){
        i = node_offset + 1; // master process does not generate replicas
        while(i < lastAssignedProcess && processActivity[i] == 1){
            i++;
        }

        if(i <= lastAssignedProcess && processActivity[i] == 0){
            result = i;
        }
    } else {
        result = i;
    }

    return result;
}

/*
 * Assigns samples to the workers. This implies to send the number of samples to be generated..
 *
 * @param workersActivity worker's state (0=idle; 1=busy)
 * @param worker worker's index to whom a sample is going to be assigned
 * @param samples samples the worker is going to generate
 */
void assignWork(int* workersActivity, int worker, int samples) {
    MPI_Send(&samples, 1, MPI_INT, worker, SAMPLES_NUMBER_TAG, shmcomm);
    workersActivity[worker]=1;
}

// **************************************  //
// WORKERS
// **************************************  //

int
workerProcess(struct params parameters, unsigned maxsites)
{
    int samples;
    char *results;
    int master;
    int flag = 1;

    samples = receiveWorkRequest(&master); // TODO: could just skip this and start generating samples? (or it could be used to set the "samples" number once
//    results = generateSamples(samples, parameters, maxsites);
//
//    sendResultsToMasterProcess(results, master);
//
//    free(results); // be good citizen

    while (flag) {
        results = generateSamples(samples, parameters, maxsites);

        sendResultsToMasterProcess(results, master);

        free(results); // be good citizen

        flag = isThereMoreWork(master);
    }

    //return isThereMoreWork(master);
    return flag;
}

char *generateSamples(int samples, struct params parameters, unsigned maxsites)
{
    char *results;
    char *sample;
    size_t offset, length;

    results = generateSample(parameters, maxsites);

    int i;
    for (i = 1; i < samples; ++i) {
        sample = generateSample(parameters, maxsites);

        offset = strlen(results);
        length = strlen(sample);

        results = realloc(results, offset + length + 1);

        memcpy(results+offset, sample, length);
        free(sample);
    }

    return results;
}

/*
 * Receives the sample's quantity the Master process asked to be generated.
 *
 * @return samples to be generated
 */
int
receiveWorkRequest(int *master){
    int samples;

    *master = 0;
    if ( shm_mode && world_rank == shm_rank ) {
        *master = 1;
    }

    MPI_Recv(&samples, 1, MPI_INT, *master, SAMPLES_NUMBER_TAG, shmcomm, MPI_STATUS_IGNORE);
    return samples;
}

int isThereMoreWork(int master) {
    int goToWork;

    MPI_Recv(&goToWork, 1, MPI_INT, master, GO_TO_WORK_TAG, shmcomm, MPI_STATUS_IGNORE);

    return goToWork;
}

/*
 * Logic to generate a sample.
 *
 * @param samples samples to be generated
 *
 * @return the sample generated by the worker
 */
char*
generateSample(struct params parameters, unsigned maxsites)
{
    int segsites;
    size_t positionStrLength, gametesStrLenght, offset;
    double probss, tmrca, ttot;
    char *results;
    char **gametes;
    struct gensam_result gensamResults;

    if( parameters.mp.segsitesin ==  0 )
        gametes = cmatrix(parameters.cp.nsam,maxsites+1);
    else
        gametes = cmatrix(parameters.cp.nsam, parameters.mp.segsitesin+1 );

    gensamResults = gensam(gametes, &probss, &tmrca, &ttot, parameters, &segsites);
    results = doPrintWorkerResultHeader(segsites, probss, parameters, gensamResults.tree);
    offset = strlen(results);

    if(segsites > 0)
    {
        char *positionsStr = doPrintWorkerResultPositions(segsites, parameters.output_precision, gensamResults.positions);
        positionStrLength = strlen(positionsStr);


        char *gametesStr = doPrintWorkerResultGametes(segsites, parameters.cp.nsam, gametes);
        gametesStrLenght = strlen(gametesStr);

        results = realloc(results, offset + positionStrLength + gametesStrLenght + 1);
        //sprintf(results, "%s%s", results, positionsStr);
        memcpy(results+offset, positionsStr, positionStrLength+1);

        offset += positionStrLength;
        memcpy(results+offset, gametesStr, gametesStrLenght+1);
        free(positionsStr);
        free(gametesStr);
        if( parameters.mp.timeflag ) {
            free(gensamResults.tree);
        }
    }

    return results;
}


/*
 * Prints the number of segregation sites:
 *    \n
 *    // xxx.x xx.xx x.xxxx x.xxxx
 *    segsites: xxx
 */
char *doPrintWorkerResultHeader(int segsites, double probss, struct params pars, char *treeOutput){
    char *results;

    int length = 3 + 1; // initially "\n//" and optionally a "\n" when there is no treeOutput;
    if( (segsites > 0 ) || ( pars.mp.theta > 0.0 ) )
    {
        length += 21; // "segsites: " + estimation of segsites digits + CR/LF
        if (pars.mp.treeflag)
        {
            length += strlen(treeOutput);
        }

        if( (pars.mp.segsitesin > 0 ) && ( pars.mp.theta > 0.0 ))
        {
            length += 17; // "prob: " + estimation of probss digits + CR/LF
        }
    }
    results = malloc(sizeof(char)*length);

    sprintf(results, "\n//");

    if( (segsites > 0 ) || ( pars.mp.theta > 0.0 ) ) {
        if( pars.mp.treeflag ) {
            sprintf(results, "%s%s", results, treeOutput);
        } else {
            sprintf(results, "%s%s", results, "\n");
        }
        if( (pars.mp.segsitesin > 0 ) && ( pars.mp.theta > 0.0 )) {
            sprintf(results, "%sprob: %g\n", results, probss);
        }
        sprintf(results, "%ssegsites: %d\n", results, segsites);
    }

    return results;
}

/*
 * Prints the segregation site positions:
 *      positions: 0.xxxxx 0.xxxxx .... etc.
 */
char *doPrintWorkerResultPositions(int segsites, int output_precision, double *positions){
    int i;
    size_t offset;

    int positionStrLength = 3+output_precision;
    int length = 12 + positionStrLength*segsites; // "positions: " + LF/CR + digit + decimal point + space
    char *results = malloc(sizeof(char) * length);

    sprintf(results, "positions: ");
    offset = 11;

    char *positionStr = malloc(sizeof(char) * positionStrLength);

    for(i=0; i<segsites; i++){
        sprintf(positionStr, "%6.*lf ", output_precision, positions[i]);
        memcpy(results+offset, positionStr, positionStrLength+1);
        offset += positionStrLength;
    }

    free(positionStr);

    return results;
}

/*
 * Print the gametes
 */
char *doPrintWorkerResultGametes(int segsites, int nsam, char **gametes){
    int i;
    size_t offset;

    int gameteStrLength = segsites+1;
    int resultsLength = 1 + gameteStrLength*nsam; // LF/CR + (segsites + LF/CR)
    char *results = malloc(sizeof(char) * resultsLength);
    sprintf(results, "\n");
    offset=1;

    char *gameteStr = malloc(sizeof(char) * gameteStrLength * nsam);

    for(i=0;i<nsam; i++) {
        //sprintf(results, "%s%s\n", results, gametes[i]);
        sprintf(gameteStr, "%s\n ", gametes[i]);
        memcpy(results+offset, gameteStr, gameteStrLength+2);
        offset += gameteStrLength;
    }

    return results;
}

/*
 * Sent Worker's results to the Master process.
 *
 * @param results results to be sent
 *
 */
void sendResultsToMasterProcess(char* results, int master)
{
    if (shm_mode && world_rank != shm_rank) {
        MPI_Send(results, strlen(results)+1, MPI_CHAR, master, RESULTS_TAG, shmcomm);
    } else { // there is either one single node, or the worker is located at the master node
        int ack = 1;
        MPI_Send(&ack, 1, MPI_INT, master, ACK_TAG, shmcomm);
        fprintf(stdout, "%s", results);
        fflush(stdout);
    }
}

// **************************************  //
// UTILS
// **************************************  //

/*--------------------------------------------------------------
 *
 *  DESCRIPTION: (Append strings)  CMS
 *
 *    Given two strings, lhs and rhs, the rhs string is appended
 *    to the lhs string, which can later on can be safely accessed
 *    by the caller of this function.
 *
 *  ARGUMENTS:
 *
 *    lhs - The left hand side string
 *    rhs - The right hand side string
 *
 *  RETURNS:
 *    A pointer to the new string (rhs appended to lhs)
 *
 *------------------------------------------------------------*/
char *append(char *lhs, const char *rhs)
{
    const size_t len1 = strlen(lhs);
    const size_t len2 = strlen(rhs);
    const size_t newSize = len1 + len2 + 1; //+1 because of the terminating null

    char *const buffer = malloc(newSize);

    strcpy(buffer, lhs);
    strcpy(buffer+len1, rhs);

    return buffer;
} /* append */

/* Initialization of the random generator. */
unsigned short * parallelSeed(unsigned short *seedv){
    unsigned short *seed48();

    return seed48(seedv);
}

/*
 * name: doInitializeRng
 * description: En caso de especificarse las semillas para inicializar el RGN,
 *              se llama a la función commandlineseed que se encuentra en el
 *              fichero del RNG.
 *
 * @param argc la cantidad de argumentos que se recibió por línea de comandos
 * @param argv el vector que tiene los valores de cada uno de los argumentos recibidos
 */
void
doInitializeRng(int argc, char *argv[], int *seeds, struct params parameters)
{
    int commandlineseed(char **);
    int arg = 0;

    while(arg < argc){
        switch(argv[arg++][1]){
            case 's':
                if(argv[arg-1][2] == 'e'){
                    // Tanto 'pars' como 'nseeds' son variables globales
                    parameters.commandlineseedflag = 1;
                    *seeds = commandlineseed(argv+arg);
                }
                break;
        }
    }
}