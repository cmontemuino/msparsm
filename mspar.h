int masterWorkerSetup(int argc, char *argv[], int howmany, struct params parameters, int maxsites);
void masterWorkerTeardown();
void masterProcessingLogic(int howmany, int lastIdleWorker, int poolSize, struct params parameters, int maxsites);
int workerProcess(int myRank, struct params parameters, int maxsites);
char* workerProcessingLogic(int myRank, int samples, struct params parameters, unsigned maxsites);
void doInitializeRng(int argc, char *argv[], int *seeds, struct params parameters);
void sendResultsToMasterProcess(char* results);
int receiveWorkRequest();
void doInitGlobalDataStructures(int argc, char *argv[], int *howmany);
void assignWork(int* workersActivity, int assignee, int samples);
void readResultsFromWorkers(int goToWork, int* workersActivity);
int findIdleProcess(int *processActivity, int poolSize, int lastAssignedProcess);
char* generateSample(struct params parameters, unsigned maxsites);
int isThereMoreWork();
unsigned short* parallelSeed(unsigned short *seedv);
char *append(char *lhs, const char *rhs);

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
