#include <mpi.h>

// LOCAL_FILE_NAME_SIZE: max size for local output file name, assuming the following format: out_xxxxxx.txt
// where 'xxxxxx' is to host the worker number (i.e.: max 1000000 workers).
#define LOCAL_FILE_NAME_SIZE 15

void masterWorker(int argc, char *argv[], int howmany, struct params parameters, int unsigned maxsites);
void teardown();
int setup(int argc, char *argv[], int howmany, struct params parameters);
int doInitializeRng(int argc, char *argv[]);
char *generateSample(struct params parameters, unsigned int maxsites, int *bytes);
void generateSamples(int, struct params, unsigned, int *bytes);
struct gensam_result gensam(char **gametes, double *probss, double *ptmrca, double *pttot, struct params pars, int* segsites);
char *append(char *lhs, const char *rhs);
char *doPrintWorkerResultHeader(int segsites, double probss, struct params paramters, char *treeOutput);
char *doPrintWorkerResultPositions(int segsites, int output_precision, double *posit);
char *doPrintWorkerResultGametes(int segsites, int nsam, char **gametes);
void initializeSeedMatrix(int argc, char *argv[], int howmany);
void printSamples(char *results, int bytes);

int calculateNumberOfNodes();

/* From ms.c*/
char ** cmatrix(int nsam, int len);
double ran1();
int commandlineseed(char **);
