#include <mpi.h>

void masterWorkerSetup(int argc, char *argv[], int howmany, struct params parameters, int unsigned maxsites);
void masterWorkerTeardown();
void doInitializeRng(int argc, char *argv[], int *seeds, struct params parameters);
char* generateSample(struct params parameters, unsigned int maxsites, int *bytes);
char *generateSamples(int, struct params, unsigned, int *bytes);
struct gensam_result gensam(char **gametes, double *probss, double *ptmrca, double *pttot, struct params pars, int* segsites);
char *append(char *lhs, const char *rhs);
char *doPrintWorkerResultHeader(int segsites, double probss, struct params paramters, char *treeOutput, int *bytes);
char *doPrintWorkerResultPositions(int segsites, int output_precision, double *posit);
char *doPrintWorkerResultGametes(int segsites, int nsam, char **gametes);
char *readResults(MPI_Comm comm, int* source, int *bytes);
unsigned short *initializeSeedMatrix(int argc, char *argv[], int howmany, struct params parameters);
void singleNodeProcessing(int howmany, struct params parameters, unsigned int maxsites, int *bytes);
void printSamples(char *results, int bytes);
void secondaryNodeProcessing(int remaining, struct params parameters, unsigned int maxsites);
void sendResultsToMaster(char *results, int bytes, MPI_Comm comm);

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
