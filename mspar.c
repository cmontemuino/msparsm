#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ms.h"
#include "mspar.h"
#include <mpi.h> /* OpenMPI library */

const int SAMPLES_NUMBER_TAG = 200;
const int RESULTS_TAG = 300;
const int GO_TO_WORK_TAG = 400;

// Following variables are with global scope in order to facilitate its sharing among routines.
// They are going to be updated in the masterWorkerSetup routine only, which is called only one, therefore there is no
// risk of race conditions or whatever other concurrency related problem.
int world_rank, shm_rank;
int world_size, shm_size;
int shm_mode = 0;  // indicates whether MPI-3 SHM is going to be applied

// **************************************  //
// MASTER
// **************************************  //

int
masterWorkerSetup(int argc, char *argv[], int howmany, struct params parameters, int maxsites)
{
    // goToWork         : used by workers to realize if there is more work to do.
    // seedMatrix       : matrix containing the RNG seeds to be distributed to working processes.
    // localSeedMatrix  : matrix used by workers to receive RNG seeds from master.
    unsigned short *seedMatrix;
    unsigned short localSeedMatrix[3];

    MPI_Comm shmcomm;

    // MPI Initialization
    MPI_Init(&argc, &argv );
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // MPI_COMM_TYPE_SHARED: This type splits the communicator into subcommunicators, each of which can create a shared memory region.
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
    MPI_Comm_size( shmcomm, &shm_size );
    MPI_Comm_rank( shmcomm, &shm_rank );

    if (shm_size != world_rank) // there are MPI process in ore than 1 computing node
    {
        shm_mode = 1;
    }

    if(world_rank == 0)
    {
        int i;
        // Only the master process prints out the application's parameters
        for(i=0; i<argc; i++)
        {
            fprintf(stdout, "%s ",argv[i]);
        }

        // If there are (not likely) more processes than samples, then the process pull
        // is cut up to the number of samples. */
        if(world_size > howmany)
        {
            world_size = howmany + 1; // the extra 1 is due to the master
        }

        int nseeds;
        doInitializeRng(argc, argv, &nseeds, parameters);
        int dimension = nseeds * world_size;
        seedMatrix = (unsigned short *) malloc(sizeof(unsigned short) * dimension);
        for(i=0; i<dimension;i++)
        {
            seedMatrix[i] = (unsigned short) (ran1()*100000);
        }
    }

    // Filter out workers with rank higher than howmany, meaning there are more workers than samples to be generated.
    if(world_rank < howmany)
    {
        MPI_Scatter(seedMatrix, 3, MPI_UNSIGNED_SHORT, localSeedMatrix, 3, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
        if(world_rank == 0) // or shm_rank = 0
        {
            // 1. if there are more nodes
            // 2. then we should split the work among nodes
            // 2.1 by doing world_size/shm_size we know how many masters are there around
            // 2.2 each shm_rank = 0 does the "masterProcessingLogic", but with reduced howmany and world_size

            // Master Processing
            masterProcessingLogic(howmany, 0);

            // if we're doing shm
            // then listen for results from other nodes
        } else
        {
            // Worker Processing
            parallelSeed(localSeedMatrix);
        }
    }

    return world_rank;
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
masterProcessingLogic(int howmany, int lastAssignedProcess) // we should get the communicator here
{
    int *processActivity = (int*) malloc(world_size * sizeof(int));
    processActivity[0] = 1; // Master does not generate replicas
    int i;

    for(i=1; i<world_size; i++)
    {
        processActivity[i] = 0;
    }

    int pendingJobs = howmany; // number of jobs already assigned but pending to be finalized

    char *results;
    char *sample;
    size_t offset, length;
    results = malloc(sizeof(char)*1000);

    while(howmany > 0)
    {
        int idleProcess = findIdleProcess(processActivity, lastAssignedProcess);
        if(idleProcess >= 0)
        {
            assignWork(processActivity, idleProcess, 1); // we should pass the communicator
            lastAssignedProcess = idleProcess;
            howmany--;
        }
        else
        {
            sample = readResultsFromWorkers(1, processActivity); // we should pass the communicator
            offset = strlen(results);
            length = strlen(sample);
            results = realloc(results, offset + length + 1);
            memcpy(results+offset, sample, length);
            free(sample);

            pendingJobs--;
        }
    }
    while(pendingJobs > 0)
    {
        char *sample = readResultsFromWorkers(0, processActivity); // we should pass the communicator
        offset = strlen(results);
        length = strlen(sample);
        results = realloc(results, offset + length + 1);
        memcpy(results+offset, sample, length);
        free(sample);

        pendingJobs--;
    }

    fprintf(stdout, "%s", results);
    free(results); // be good citizen

    // 1. at this point all of the work was done
    // 2. if there are more nodes (shm_mode?)
    // 3. if I'm not the super master, then send results to WORLD rank 0
    // 4. otherwise print the accumulated results? (or not accumulate them if not needed?)
    // don't forget to 'free' the results
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
    MPI_Status status;
    int size;
    int source;

    MPI_Probe(MPI_ANY_SOURCE, RESULTS_TAG, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_CHAR, &size);
    source = status.MPI_SOURCE;
    char * results = (char *) malloc(size*sizeof(char));

    MPI_Recv(results, size, MPI_CHAR, source, RESULTS_TAG, MPI_COMM_WORLD, &status);
    source = status.MPI_SOURCE;
    MPI_Send(&goToWork, 1, MPI_INT, source, GO_TO_WORK_TAG, MPI_COMM_WORLD);

    workersActivity[source]=0;
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
int findIdleProcess(int *processActivity, int lastAssignedProcess) {
    /*
     * Implementation note: lastAssignedProcess is used to implement a fairness policy in which every available process
     * can be assigned with some work.
     */
    int result = -1;
    int i= lastAssignedProcess+1; // master process does not generate replicas
    while(i < world_size && processActivity[i] == 1){
        i++;
    };

    if(i >= world_size){
        i=1; // master process does not generate replicas
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
    MPI_Send(&samples, 1, MPI_INT, worker, SAMPLES_NUMBER_TAG, MPI_COMM_WORLD);
    //TODO check usage of MPI_Sendv??
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

    samples = receiveWorkRequest();
    results = generateSamples(samples, parameters, maxsites);

    sendResultsToMasterProcess(results);

    free(results); // be good citizen
    return isThereMoreWork();
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
int receiveWorkRequest(){
    int samples;
    MPI_Status status;

    MPI_Recv(&samples, 1, MPI_INT, 0, SAMPLES_NUMBER_TAG, MPI_COMM_WORLD, &status);
    return samples;
}

int isThereMoreWork() {
    int goToWork;
    MPI_Status status;

    MPI_Recv(&goToWork, 1, MPI_INT, 0, GO_TO_WORK_TAG, MPI_COMM_WORLD, &status);

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
    double *positions;
    struct gensam_result gensamResults;

    if( parameters.mp.segsitesin ==  0 )
        gametes = cmatrix(parameters.cp.nsam,maxsites+1);
    else
        gametes = cmatrix(parameters.cp.nsam, parameters.mp.segsitesin+1 );

    gensamResults = gensam(gametes, &probss, &tmrca, &ttot, parameters, &segsites);
    positions = gensamResults.positions;
    results = doPrintWorkerResultHeader(segsites, probss, parameters, gensamResults.tree);
    offset = strlen(results);

    if(segsites > 0)
    {
        char *positionsStr = doPrintWorkerResultPositions(segsites, parameters.output_precision, positions);
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
        free(gensamResults.positions);
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
void sendResultsToMasterProcess(char* results)
{
    MPI_Send(results, strlen(results)+1, MPI_CHAR, 0, RESULTS_TAG, MPI_COMM_WORLD);
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
char *
append(char *lhs, const char *rhs)
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