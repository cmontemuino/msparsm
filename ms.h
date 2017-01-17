struct devent {
	double time;                    // time in the past for the demographic event change
	int popi;                       // subpopulation i (only when 'detype' in [n,g,s,m,j])
	int popj;                       // subpopulation j (only when 'detype' in [a,j])
	double paramv;                  // growth rate value
	double **mat ;                  // values for all elements of the migration matrix (only when 'detype' is 'a')
	char detype ;                   // past demographic event type (N,G,M,n,g,s,m,a,j)
	struct devent *nextde;
} ;

struct c_params {
	int npop;                       // number of subpopulations (island models)
	int nsam;                       // sample size
	int *config;                    // number of chromosomes sampled from each subpopulation
	double **mig_mat;               // migration matrix for subpopulations
	double r;                       // scaled recombination rate
	int nsites;                     // number of base pairs in the locus
	double f;                       // denote ratio g/r (g = gene conversion probability)
	double track_len;               // mean conversion track length (related to 'f')
	double *size;                   // size of subpopulations relative to 'No' (default to 1.0)
	double *alphag;                 // growth rate of subpopulations
	struct devent *deventlist ;
} ;
struct m_params {
	double theta;
	int segsitesin;
	int treeflag;
	int timeflag;
	int mfreq;
} ;
struct params {
	struct c_params cp;
	struct m_params mp;
	int commandlineseedflag ;
	int output_precision;
};

// Represents a node in the history tree of a gamete
struct node{
	int abv;        // number of node ancestral to this node
	int ndes;
	float time;     // time (in units of 4N generations) of the node
};

// Result structure returned by the gensam function
struct gensam_result {
	// positions of the segregating sites (on a scale of 0.0 - 1.0)
	double 	*positions;
	// tree output
	char	*tree;
};

// Represents the history for a given segment of a gamete
struct segl {
    int beg;                // starting point of segment
    struct node *ptree;     // points to the first node of the tree representing the history of the segment
    int next;               // index number of the next segment
};

struct seg{
    int beg;
    int end;
    int desc;
};

struct chromo{
    int nseg;
    int pop;
    struct seg  *pseg;
};

/*KRT -- prototypes added*/
void ordran(int n, double pbuf[]);
void ranvec(int n, double pbuf[]);
void order(int n, double pbuf[]);

void biggerlist(int nsam,  char **list, unsigned maxsites );
int poisso(double u);
void locate(int n,double beg, double len,double *ptr);
void mnmial(int n, int nclass, double p[], int rv[]);
void usage();
int tdesn(struct node *ptree, int tip, int node );
int pick2(int n, int *i, int *j);
void pick2_chrom(int pop,int config[], int *pc1, int *pc2, struct chromo *chrom);
int re(int nsam, int *nsegs, struct segl **seglst, struct chromo **chrom, int *maxchr, int *nchrom, int **nnodes, long *nlinks, double *cleft, double pc);
int xover(int nsam,int ic, int is, int *nsegs, struct segl **seglst, struct chromo **chrom, int *maxchr, int *nchrom, int **nnodes, long *nlinks, double *cleft, double pc);
int ca(int nsam, int nsites, int c1, int c2, double time, int nsegs, struct segl *seglst, struct chromo *chrom, int *nchrom, int *nnodes, long *nlinks, double *cleft, double pc, double *wallclock);
int cinr( int nsam, int nsites, double time, int *nsegs, struct segl **seglst, struct chromo **chrom, int *maxchr, int *nchrom, int **nnodes, long *nlinks, double *cleft, double pc, double lnpc);
int cleftr( int nsam, int *nsegs, struct segl **seglst, struct chromo **chrom, int *maxchr, int *nchrom, int **nnodes, long *nlinks, double *cleft, double pc, double lnpc);
int links(int c, struct chromo *chrom);
int isseg(int start, int c, int *psg, struct chromo *chrom);
double ran1();
