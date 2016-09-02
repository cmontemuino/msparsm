#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "ms.h"
#include "mspar.h"

const int RESULTS_TAG = 300;

int diagnose = 0; // Used for diagnosing the application.

// Following variables are with global scope in order to facilitate its sharing among routines.
// They are going to be updated in the masterWorkerSetup routine only, which is called only one, therefore there is no
// risk of race conditions or whatever other concurrency related problem.
MPI_Comm shmcomm; // shm intra-communicator
int world_rank, shm_rank;
int world_size, shm_size;

// **************************************  //
// MASTER
// **************************************  //
void singleNodeProcessing(int howmany, struct params parameters, unsigned int maxsites, int *bytes)
{
    // No master process is needed. Every MPI process can just output the generated samples
    int samples = howmany / world_size;
    int remainder = howmany % world_size;

    if (world_rank == 0) // let the "global master" to generate the remainder samples as well
        samples += remainder;

    if (diagnose)
        fprintf(stderr, "[%d] -> Vamos a generar [%d] samples.\n", world_rank, samples);

    char *results = generateSamples(samples, parameters, maxsites, bytes);
    printSamples(results, *bytes);
}

void printSamples(char *results, int bytes)
{
    fprintf(stdout, "%s", results);
    fflush(stdout);

    if (diagnose)
        fprintf(stderr, "[%d] -> Printed [%d] bytes.\n", world_rank, bytes);

    free(results); // be good citizen
}

void secondaryNodeProcessing(int remaining, struct params parameters, unsigned int maxsites)
{
    int bytes = 0;
    char *results = malloc(sizeof(char) * 1000);
    if (remaining > 0)
        results = generateSamples(remaining, parameters, maxsites, &bytes);

    // Receive samples from workers in same node.
    int i;
    char *shm_results;
    for (i = 1; i < shm_size; i++){
        int source, length;

        shm_results = readResults(shmcomm, &source, &length);
        results = realloc(results, bytes + length + 1);
        memcpy(results + bytes, shm_results, length);
        bytes += length + 1;
        free(shm_results);
    }

    // Send gathered results to master in master-node
    sendResultsToMaster(results, bytes, MPI_COMM_WORLD);
}

void sendResultsToMaster(char *results, int bytes, MPI_Comm comm)
{
    MPI_Send(results, bytes + 1, MPI_CHAR, 0, RESULTS_TAG, comm);

    if (diagnose) {
        char *communicator = "MPI_COMM_WORLD";
        if (comm != MPI_COMM_WORLD)
            communicator = "SHM_COMM";

        fprintf(stderr, "[%d] -> Sent [%d] bytes to master in %s.\n", world_rank, bytes + 1, communicator);
    }

    free(results);
}

void principalMasterProcessing(int remaining, int nodes, struct params parameters, unsigned int maxsites)
{
    int bytes = 0;
    if (remaining > 0) {
        char *results = generateSamples(remaining, parameters, maxsites, &bytes);
        printSamples(results, bytes);
    }

    // Receive samples from other node masters, each one sending a consolidated message
    int source, i;
    char *shm_results;
    for (i = 1; i < nodes; i++){
        shm_results = readResults(MPI_COMM_WORLD, &source, &bytes);
        printSamples(shm_results, bytes);
    }
}

int calculateNumberOfNodes()
{
    // Gather all SHM rank (local rank) from all the processes
    int *shm_ranks = (int *)malloc(sizeof(int) * world_size);

    // Note: elements are ordered by process rank in MPI_COMM_WORLD communicator
    MPI_Gather (&shm_rank, 1, MPI_INT, shm_ranks, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int i = 0;
    int nodes = 0;
    if (world_rank == 0) {
        for (i = 0; i < world_size; i++)
            if (shm_ranks[i] == 0) nodes++;
    }

    return nodes;
}

int setup(int argc, char *argv[], int howmany, struct params parameters)
{
    // seedMatrix       : matrix containing the RNG seeds to be distributed to working processes.
    // localSeedMatrix  : matrix used by workers to receive RNG seeds from master.
    unsigned short *seedMatrix;
    unsigned short localSeedMatrix[3];

    if (getenv("MSPARSM_DIAGNOSE")) diagnose = 1;

    // MPI Initialization
    MPI_Init(&argc, &argv );
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // MPI_COMM_TYPE_SHARED: This type splits the communicator into subcommunicators, each of which can create a shared memory region.
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
    MPI_Comm_size(shmcomm, &shm_size);
    MPI_Comm_rank(shmcomm, &shm_rank);

    if (world_rank == 0) { // print out program parameters
        int i;
        for(i=0; i<argc; i++)
            fprintf(stdout, "%s ",argv[i]);
        fflush(stdout);
    }

    initializeSeedMatrix(argc, argv, howmany);

    int nodes = calculateNumberOfNodes();

    if (diagnose)
        fprintf(stderr, "[%d] -> SHM Rank=%d, SHM Size=%d, WORLD Size=%d\n", world_rank, shm_rank, shm_size, world_size);

    if(diagnose && world_rank == 0)
        fprintf(stderr, "[%d] -> # of nodes=%d\n", world_rank, nodes);

    return nodes;
}

void teardown() {
    MPI_Finalize();
}

void masterWorker(int argc, char *argv[], int howmany, struct params parameters, unsigned int maxsites)
{
    int nodes = setup(argc, argv, howmany, parameters);

    // Filter out workers with rank higher than howmany, meaning there are more workers than samples to be generated.
    if(world_rank < howmany) {
        if (world_size == shm_size) { // There is only one node
            int bytes;
            singleNodeProcessing(howmany, parameters, maxsites, &bytes);
        } else {
            MPI_Bcast(&nodes, 1, MPI_INT, 0, MPI_COMM_WORLD);

            int nodeSamples = howmany / nodes;
            int remainingGlobal = howmany % nodes;
            int workerSamples = nodeSamples / (shm_size - 1);
            int remainingLocal = nodeSamples % (shm_size - 1);

            if (world_rank != 0 && shm_rank != 0) {
                int bytes = 0;
                char *results = generateSamples(workerSamples, parameters, maxsites, &bytes);

                if (world_rank == shm_rank)
                    printSamples(results, bytes);
                else  // Send results to shm_rank = 0
                    sendResultsToMaster(results, bytes, shmcomm);
            } else {
                if (world_rank != 0 && shm_rank == 0) {
                    secondaryNodeProcessing(remainingLocal, parameters, maxsites);
                } else
                    principalMasterProcessing(remainingGlobal +  remainingLocal, nodes, parameters, maxsites);
            }
        }
    }

    teardown();
}

char *readResults(MPI_Comm comm, int *source, int *bytes)
{
    MPI_Status status;

    MPI_Probe(MPI_ANY_SOURCE, RESULTS_TAG, comm, &status);
    MPI_Get_count(&status, MPI_CHAR, bytes);
    *source = status.MPI_SOURCE;

    char *results = (char *) malloc(*bytes * sizeof(char));

    MPI_Recv(results, *bytes, MPI_CHAR, *source, RESULTS_TAG, comm, MPI_STATUS_IGNORE);

    if (diagnose)
        fprintf(stderr, "[%d] -> Read [%d] bytes from worker %d.\n", world_rank, *bytes, *source);

    return results;
}

char *generateSamples(int samples, struct params parameters, unsigned maxsites, int *bytes)
{
    char *results;
    char *sample;
    int length;

    results = generateSample(parameters, maxsites, &length);

    *bytes = length;

    int i;
    for (i = 1; i < samples; ++i) { // starts in 1 because a sample was already generated before the for-loop
        sample = generateSample(parameters, maxsites, &length);

        results = realloc(results, *bytes + length + 2);

        memcpy(results + *bytes, sample, length);

        *bytes += length;
        free(sample);
    }

    if (diagnose)
        fprintf(stderr, "[%d] -> Generated [%d] samples.\n", world_rank, samples);

    return results;
}

/*
 * Logic to generate a sample.
 *
 * @param samples samples to be generated
 *
 * @return the sample generated by the worker
 */
char* generateSample(struct params parameters, unsigned maxsites, int *bytes)
{
    int segsites;
    size_t offset, positionStrLength, gametesStrLenght;
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
    *bytes = offset;


    if(segsites > 0)
    {
        char *positionsStr = doPrintWorkerResultPositions(segsites, parameters.output_precision, gensamResults.positions);
        positionStrLength = strlen(positionsStr);

        char *gametesStr = doPrintWorkerResultGametes(segsites, parameters.cp.nsam, gametes);
        gametesStrLenght = strlen(gametesStr);

        results = realloc(results, offset + positionStrLength + gametesStrLenght + 2);

        memcpy(results+offset, positionsStr, positionStrLength+1);

        offset += positionStrLength;
        *bytes += positionStrLength;

        memcpy(results+offset, gametesStr, gametesStrLenght+1);

        *bytes += gametesStrLenght;

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
        if (!pars.mp.treeflag)
            asprintf(&treeOutput, "\n");

        length += strlen(treeOutput);

        if( (pars.mp.segsitesin > 0 ) && ( pars.mp.theta > 0.0 ))
            length += 17; // "prob: " + estimation of probss digits + CR/LF
    }
    results = malloc(sizeof(char)*length);

    sprintf(results, "\n//");

    if( (segsites > 0 ) || ( pars.mp.theta > 0.0 ) ) {
        sprintf(results, "%s%s", results, treeOutput);

        if( (pars.mp.segsitesin > 0 ) && ( pars.mp.theta > 0.0 ))
            sprintf(results, "%sprob: %g\n", results, probss);

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
        sprintf(gameteStr, "%s\n ", gametes[i]);
        memcpy(results+offset, gameteStr, gameteStrLength+2);
        offset += gameteStrLength;
    }

    return results;
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
}

/*
 * doInitializeRng - Initializes the Random Number Generator
 * @argc number of arguments passed to the program
 * @argv array holding the arguments passed to the program
 *
 * Reads the RGN seeds from command arguments and use them for initialize
 * the RGN that will generate the seeds that worker process will use
 * for initialize their own RGN.
 *
 * This function must be called by the master process located at the
 * main node only.
 */
int doInitializeRng(int argc, char *argv[])
{
    int arg = 0;
    int result = 0;

    while(arg < argc){
        switch(argv[arg++][1]){
        case 's':
            if(argv[arg-1][2] == 'e')
                commandlineseed(argv+arg);
            break;
        default:
            continue;
        }
    }

    return result;
}

void initializeSeedMatrix(int argc, char *argv[], int howmany) {
    int i;
    int size = world_size;

    // If there are (not likely) more processes than samples, then the process pool
    // is cut up to the number of samples. */
    if(world_size > howmany)
        size = howmany + 1; // the extra 1 is due to the master

    int dimension = 3 * size;
    unsigned short *seedMatrix = (unsigned short *) malloc(sizeof(unsigned short) * dimension);

    if (world_rank == 0) {
        doInitializeRng(argc, argv);

        for(i=0; i<dimension;i++)
            seedMatrix[i] = (unsigned short) (ran1()*100000);

    }

    unsigned short localSeedMatrix[3];
    MPI_Scatter(seedMatrix, 3, MPI_UNSIGNED_SHORT, localSeedMatrix, 3, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);

    if (world_rank != 0)
        seed48(localSeedMatrix);
}