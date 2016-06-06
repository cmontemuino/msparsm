#include <mpi.h>

int masterWorkerSetup(int argc, char *argv[], int howmany, struct params parameters, int unsigned maxsites, int *excludeFrom, int *node_rank);
void masterWorkerTeardown();
void masterProcessingLogic(int howmany, int lastIdleWorker, struct params parameters, unsigned int maxsites);
int workerProcess(struct params parameters, unsigned int maxsites);
void doInitializeRng(int argc, char *argv[], int *seeds, struct params parameters);
void sendResultsToMasterProcess(char* results, int master);
int receiveWorkRequest(int *master);
void assignWork(int* workersActivity, int assignee, int samples);
char* readResultsFromWorkers(int goToWork, int* workersActivity);
int findIdleProcess(int *processActivity, int lastAssignedProcess, int node_offset);
char* generateSample(struct params parameters, unsigned int maxsites);
char *generateSamples(int, struct params, unsigned);
struct gensam_result gensam(char **gametes, double *probss, double *ptmrca, double *pttot, struct params pars, int* segsites);
int isThereMoreWork(int master);
unsigned short* parallelSeed(unsigned short *seedv);
char *append(char *lhs, const char *rhs);
char *doPrintWorkerResultHeader(int segsites, double probss, struct params paramters, char *treeOutput);
char *doPrintWorkerResultPositions(int segsites, int output_precision, double *posit);
char *doPrintWorkerResultGametes(int segsites, int nsam, char **gametes);

int readAckFromLocalWorker(int remaining, int *workersActivity, struct params parameters, unsigned int maxsites, int *goToWork);
int numberOfNodes(void*);
char *readResults(MPI_Comm, int*);

/* From ms.c*/
char ** cmatrix(int nsam, int len);
double ran1();

/*
void ordran(int n, double pbuf[]);
void ranvec(int n, double pbuf[]);
void order(int n, double pbuf[]);

void biggerlist(int nsam,  char **list );
int poisso(double u);
void locate(int n,double beg, double len,double *ptr);
void mnmial(int n, int nclass, double p[], int rv[]);
void usage();
int tdesn(struct node *ptree, int tip, int node );
int pick2(int n, int *i, int *j);
int xover(int nsam,int ic, int is);
int links(int c);
*/
