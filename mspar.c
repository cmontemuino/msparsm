#define _GNU_SOURCE

const int SEEDS_COUNT = 3;
const int SEED_TAG = 100;
const int SAMPLES_NUMBER_TAG = 200;
const int RESULTS_TAG = 300;
const int GO_TO_WORK_TAG = 400;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "ms.h"
#include "mspar.h"
#include <mpi.h> /* OpenMPI library */

// **************************************  //
// MASTER
// **************************************  //

int
masterWorkerSetup(int argc, char *argv[], int howmany, struct params parameters, int maxsites)
{
    // global_rank           : rank of the current process in the global MPI ecosystem.
    // poolSize         : number of processes in the MPI ecosystem.
    // goToWork         : used by workers to realize if there is more work to do.
    // seedMatrix       : matrix containing the RNG seeds to be distributed to working processes.
    // localSeedMatrix  : matrix used by workers to receive RNG seeds from master.
    int global_rank;
    int poolSize;
    unsigned short *seedMatrix;
    unsigned short localSeedMatrix[3];

    MPI_Comm shmcomm;

    // MPI Initialization
    MPI_Init(&argc, &argv );
    MPI_Comm_size(MPI_COMM_WORLD, &poolSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);

    if(global_rank == 0)
    {
        int i;
        // Only the master process prints out the application's parameters
        for(i=0; i<argc; i++)
        {
          fprintf(stdout, "%s ",argv[i]);
        }

        // If there are (not likely) more processes than samples, then the process pull
        // is cut up to the number of samples. */
        if(poolSize > howmany)
        {
            poolSize = howmany + 1; // the extra 1 is due to the master
        }

        int nseeds;
        doInitializeRng(argc, argv, &nseeds, parameters);
        int dimension = nseeds * poolSize;
        seedMatrix = (unsigned short *) malloc(sizeof(unsigned short) * dimension);
        for(i=0; i<dimension;i++)
        {
          seedMatrix[i] = (unsigned short) (ran1()*100000);
        }
    }

    // Filter out workers with rank higher than howmany, meaning there are more workers than samples to be generated.
    if(global_rank <= howmany)
    {
        MPI_Scatter(seedMatrix, 3, MPI_UNSIGNED_SHORT, localSeedMatrix, 3, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
        if(global_rank == 0)
        {
            // Master Processing
            masterProcessingLogic(howmany, 0, poolSize, parameters, maxsites);
        } else
        {
            // Worker Processing
            parallelSeed(localSeedMatrix);
        }
    }

    return 0;
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
 * @param poolSize number of processes in the MPI ecosystem
 * @param parameters parameters used to generate replicas. Used when worker process is the master itself
 * @param maxsites maximum number of sites. Used when worker process is the master itself
 */
void
masterProcessingLogic(int howmany, int lastAssignedProcess, int poolSize, struct params parameters, int maxsites)
{
    char *generateSamples(int samples, struct params parameters, int maxsites);

    int *processActivity = (int*) malloc(poolSize * sizeof(int));
    int i;
    for(i=0; i<poolSize; i++)   processActivity[i] = 0;

    int pendingJobs = howmany; // number of jobs already assigned but pending to be finalized
    int assignToMySelf = 0;

    while(howmany > 0)
    {
        int idleProcess = findIdleProcess(processActivity, poolSize, lastAssignedProcess);
        if(idleProcess >= 0)
        {
          if(idleProcess == 0) {
              assignToMySelf = 1;
          } else {
            assignWork(processActivity, idleProcess, 1);
          }
          lastAssignedProcess = idleProcess;
          howmany--;
        }
        else
        {
            if(assignToMySelf == 1) {
                char *results;

                results = generateSamples(1, parameters, maxsites);
                fprintf(stdout, "%s", results);
                assignToMySelf = 0;
                free(results); // be good citizen
            } else {
                readResultsFromWorkers(1, processActivity);
            }
            processActivity[lastAssignedProcess] = 0;
            pendingJobs--;
        }
    }
    while(pendingJobs > 0)
    {
        readResultsFromWorkers(0, processActivity);
        pendingJobs--;
    }
}

/*
 *
 * Esta función realiza dos tareas: por un lado hace que el Master escuche los resultados enviados por los workers y por
 * otro lado, se envía al worker que se ha recibido la información un mensaje sobre si debe seguir esperando por
 * trabajos o si ha de finalizar su contribución al sistema.
 *
 * @param goToWork indica si el worker queda en espera de más trabajo (1) o si ya puede finalizar su ejecución (0)
 * @param workersActivity el vector con el estado de actividad de los workers
 * @return
 *
 */
void readResultsFromWorkers(int goToWork, int* workersActivity)
{
    MPI_Status status;
    int size;
    int source;

    MPI_Probe(MPI_ANY_SOURCE, RESULTS_TAG, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_CHAR, &size);
    source = status.MPI_SOURCE;
    char * results = (char *) malloc(size*sizeof(char));

    MPI_Recv(results, size, MPI_CHAR, source, RESULTS_TAG, MPI_COMM_WORLD, &status);
    MPI_Send(&goToWork, 1, MPI_INT, source, GO_TO_WORK_TAG, MPI_COMM_WORLD);

    workersActivity[source]=0;
    fprintf(stdout, "%s", results);
    free(results); // be good citizen
}

/*
 * Finds an idle process from a list of potential worker processes (that could potentially
 * include the master process as well).
 *
 * @param workersActivity status of all processes that can generate some work (0=idle; 1=busy)
 * @param poolSize number of worker processes
 * @lastAssignedProcess last process assigned with some work
 *
 * @return idle process index or -1 if all processes are busy.
 */
int findIdleProcess(int *processActivity, int poolSize, int lastAssignedProcess) {
  /*
   * Implementation note: lastAssignedProcess is used to implement a fairness policy in which every availble process
   * can be assigned with some work.
   */
  int result = -1;
  int i= lastAssignedProcess + 1;
  while(i < poolSize && processActivity[i] == 1){
    i++;
  };

  if(i >= poolSize){
    i=0;
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
workerProcess(int myRank, struct params parameters, int maxsites)
{
    char *generateSamples(int samples, struct params parameters, int maxsites);

    int samples;
    char *results;

    samples = receiveWorkRequest();
    results = generateSamples(samples, parameters, maxsites);

    sendResultsToMasterProcess(results);

    free(results); // be good citizen
    return isThereMoreWork();
}

char *generateSamples(int samples, struct params parameters, int maxsites)
{
    char *results;
    char *singleResult;

    results = generateSample(parameters, maxsites);

    for (int i = 1; i < samples; ++i) {
        singleResult = generateSample(parameters, maxsites);
        results = append(results, singleResult);
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
    struct gensam_result gensam(char **gametes, double *probss, double *ptmrca, double *pttot, struct params pars, int* segsites);
    char *doPrintWorkerResultHeader(int segsites, double probss, struct params paramters, char *results);
    char *doPrintWorkerResultPositions(int segsites, int output_precision, double *posit, char *results);

    char *doPrintWorkerResultGametes(int segsites, int nsam, char **gametes, char *results);
    int segsites;
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
    if(segsites > 0)
    {
        results = doPrintWorkerResultPositions(segsites, parameters.output_precision, positions, results);
        results = doPrintWorkerResultGametes(segsites, parameters.cp.nsam, gametes, results);
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
    char *tempString;

    char *results = malloc(4);
    sprintf(results, "\n//");

    if( (segsites > 0 ) || ( pars.mp.theta > 0.0 ) ) {
        if( pars.mp.treeflag ) {
            results = append(results, treeOutput);
        } else {
            results = append(results, "\n");
        }
        if( (pars.mp.segsitesin > 0 ) && ( pars.mp.theta > 0.0 )) {
            asprintf(&tempString, "prob: %g\n", probss);
            results = append(results, tempString);
        }
        asprintf(&tempString, "segsites: %d\n",segsites);
        results = append(results, tempString);
    }

    return results;
}

/*
 * Prints the segregation site positions:
 *      positions: 0.xxxxx 0.xxxxx .... etc.
 */
char *doPrintWorkerResultPositions(int segsites, int output_precision, double *positions, char *results){
    int i;
    char tempString[3 + output_precision]; //number+decimal point+space

    results = append(results, "positions: ");

    for(i=0; i<segsites; i++){
        sprintf(tempString, "%6.*lf ", output_precision, positions[i]);
        results = append(results, tempString);
    }
    return append(results, tempString);
}

/*
 * Prints the gametes
 */
char *doPrintWorkerResultGametes(int segsites, int nsam, char **gametes, char *results){
    int i;
    char tempString[segsites+1]; //segsites + LF/CR

    results = append(results, "\n");

    for(i=0;i<nsam; i++) {
        sprintf(tempString, "%s\n", gametes[i]);
        results = append(results, tempString);
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
	size_t len1 = strlen(lhs);
	size_t len2 = strlen(rhs) + 1; //+1 because of the terminating null

	lhs = realloc(lhs, len1 + len2);

	return strncat(lhs, rhs, len2);
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
