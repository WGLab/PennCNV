
typedef struct {
	int N;			/* number of states;  Q={1,2,...,N} */
	int M; 			/* number of observation symbols; V={1,2,...,M}*/
	double **A;		/* A[1..N][1..N]. a[i][j] is the transition prob
			   	of going from state i at time t to state j
			   	at time t+1 */
	double **B;		/* B[1..N][1..M]. b[j][k] is the probability of
			   	of observing symbol k in state j */
	double *pi;		/* pi[1..N] pi[i] is the initial state distribution. */
	double *B1_mean;	/* B1_mean[1..N] mean of a continuous Gaussian distribution for state 1 through N*/
	double *B1_sd;		/*B1_sd standard deviation of B1 values, which is the same for all states*/
	double B1_uf;		/*B1_uniform_fraction: the contribution of uniform distribution to the finite mixture model */
	double *B2_mean;	/* B2_mean[1..4] is the average of B_allele_freq*/
	double *B2_sd;		/* B2_sd[1..4] is the standard deviation of four B_allele_freq, B2_sd[5] is specially for state1, where B is modelled as a wide normal distribution */
	double B2_uf;		/* B2_uniform_fraction: the fraction of uniform distribution in the finite mixture model */
	
	int NP_flag;		/*flag of 1 and 0 to indicate whether Non-Polymorhpic marker information is contained with HMM file*/
	double *B3_mean;	/* B3_mean[1..N] mean of non-polymorphic probe for state 1 through N*/
	double *B3_sd;		/* B3_sd[1..4] is the standard deviation of B3 values*/
	double B3_uf;		/* B3_uniform_fraction: */
	int dist;		/* new parameter to facilitate CNV calling from resequencing data (2009 April) */
} CHMM;


/************************************
*subroutines for processing continuous HMM models
************************************/
CHMM ReadCHMM (char *filename);
void PrintCHMM(FILE *fp, CHMM *phmm);
void CopyCHMM(CHMM *phmm1, CHMM *phmm2);
void FreeCHMM(CHMM *phmm);
void testVit_CHMM (CHMM hmm, int T, double *O1, double *O2, double *pfb, int *snpdist, double *plogproba);
void tumorVit_CHMM (CHMM hmm, int T, double *O1, double *O2, double *pfb, int *snpdist, double *plogproba, double stroma);
void ViterbiLogNP_CHMM(CHMM *phmm, int T, double *O1, double *O2, double *pfb, int *snpdist, double **delta, int **psi, int *q, double *pprob);
void ViterbiLogNPStroma_CHMM(CHMM *phmm, int T, double *O1, double *O2, double *pfb, int *snpdist, double **delta, int **psi, int *q, double *pprob, double stroma);

void estHMMFromFile_CHMM (CHMM hmm, int T, FILE *fp, int *niter, double *logprobinit, double *logprobfinal);
void BaumWelchNP_CHMM(CHMM *phmm, int T, double *O1, double *O2, double *pfb, int *snpdist, double **alpha, double **beta, double **gamma, int *pniter, double *plogprobinit, double *plogprobfinal);
void ComputeXi_CHMM (CHMM* phmm, int T, double *O1, double *O2, double **biot, int *snpdist, double **alpha, double **beta, double ***xi);
void ComputeGamma_CHMM (CHMM *phmm, int T, double **alpha, double **beta, double **gamma);
void BackwardWithScale_CHMM (CHMM *phmm, int T, double *O1, double *O2, double **biot, int *snpdist, double **beta, double *scale, double *pprob);
void ForwardWithScale_CHMM (CHMM *phmm, int T, double *O1, double *O2, double **biot, int *snpdist, double **alpha, double *scale, double *pprob);

void GetStateProb_CHMM(CHMM *phmm, int T, double *O1, double *O2, double *pfb, int *snpdist, double *pprob, int state);
void convertHMMTransition (CHMM *phmm, double **A1, int dist);
void adjustBSD (CHMM *phmm, double sdo);

double b1iot (int state, double *mean, double *sd, double uf, double o);
double b2iot (int state, double *mean, double *sd, double uf, double pfb, double b);
double b5iot (int state, double *mean, double *sd, double uf, double o, double stroma);
double b6iot (int state, double *mean, double *sd, double uf, double pfb, double b, double stroma);

void testVitTrio_CHMM (CHMM *phmm, int T, double *Of1, double *Of2, double *Om1, double *Om2, double *Oo1, double *Oo2, double *pfb, int *snpdist, double *plogproba, double *trio_lrr_sd, int osex, int chrx_flag);
double trioTransition (double cimcell, double chit, double **A1, double e1, double e2, int dist, int pdn, int pifat, int pimot, int pioff, int dn, int ifat, int imot, int ioff);
void calculateCHIT (double a, double M[5][5][5][5][5][5]);
void calculateCHIT_male_x (double a, double M[5][5][5][5][5][5]);
void calculateCHIT_female_x (double a, double M[5][5][5][5][5][5]);

double *** AllocXi(int T, int N);
void FreeXi(double *** xi, int T, int N);

void callCNVFromFile_SEQ (CHMM hmm, int T, FILE *fp, int *mlstate, double gamma_k, double gamma_theta);
void ViterbiLogNP_SEQ (CHMM *phmm, int T, double k, double theta, double *O1, int *O2, double *pfb, int *length, double **delta, int **psi, int *q);
double b4iot (int state, double *mean, double *sd, double uf, double pfb, int o1, int o2);
double b3iot (int state, double *mean, double *sd, double uf, double o1, int length);
void adjustHMMExpected (CHMM *phmm, double coverage);
