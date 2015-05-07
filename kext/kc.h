
#include "EXTERN.h"   /* std perl include */
#include "perl.h"     /* std perl include */
#include "XSUB.h"     /* XSUB include */



/********************
Perl-C interface
********************/

FILE *fh_stdout ();
FILE *fh_stderr ();
FILE *fh_stdin ();
void kcerror(char *error_text);
void kcwarn(char *error_text);
SV* C2perl_vector (double *vector, int n);
SV* C2perl_vectori (int *vector, int n);
double* perl2C_vector (SV* arg);
int* perl2C_vectori (SV* arg);
SV* C2perl_matrix (double** mat, int n1, int n2);
SV* C2perl_matrixi (int** mat, int n1, int n2);
double **perl2C_matrix (SV* arg);
int** perl2C_matrixi (SV* arg);
void print_vectori (int *array, int n);
void print_vector (double *array, int n);
void print_matrixi (int** A, int m, int n);
void print_matrix (double** A, int m, int n);

int is_scalar_ref (SV* arg);
void* pack1D(SV* arg, char packtype);
void  unpack1D(SV* arg, void * var, char packtype, int n);
AV*   coerce1D ( SV* arg, int n );
void* get_mortalspace( int n, char packtype );

/***********************
Matrix memory allocation
***********************/
double **matrix_new (int row, int col);
int **matrixi_new (int row, int col);
char **matrixc_new (int row, int col);
int** matrix_toMatrixi (double** A, int m, int n);
double** matrixi_toMatrix (int** A, int m, int n);
double *vector_new (int n);
int *vectori_new (int n);
void matrix_free (double **m);
void matrixi_free (int **m);
void vector_free (double *m);
void vectori_free (int *m);
double matrix_at (double **A, int m, int n);
int matrixi_at (int **A, int m, int n);
double vector_at (double *A, int n);
int vectori_at (int *A, int n);
double* matrix_row (double** A, int m);
int* matrixi_row (int** A, int m);
void matrix_change (double **A, int m, int n, double a);
void matrixi_change (int **A, int m, int n, int a);
void vector_change (double *A, int n, double a);
void vectori_change (int *A, int n, int a);
int* vector_toVectori (double* A, int n);
double* vectori_toVector (int* A, int n);


void nrerror(char *error_text);
int *ivector(long nl, long nh);
unsigned char *cvector(long nl, long nh);
double *dvector(long nl, long nh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
double **submatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch, long newrl, long newcl);
double **convert_matrix(double *a, long nrl, long nrh, long ncl, long nch);
void free_ivector(int *v, long nl, long nh);
void free_cvector(unsigned char *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_subdmatrix(double **b, long nrl, long nrh, long ncl, long nch);


/**********************
Basic matrix operations
*********************/
double** matrix_copy (double **A, int m, int n);
int** matrixi_copy (int **A, int m, int n);
double** matrix_transpose (double **A, int m, int n);
int** matrixi_transpose (int **A, int m, int n);
double** matrix_newSeq (int m, int n);
int** matrixi_newSeq (int m, int n);
double** matrix_newRandom (int m, int n, int seed);
double** matrix_newDiagonal (double* cell, int length, int m, int n);
double** matrix_newIdentity (int m);
double** matrix_newSubrow (double** A, int m, int n, int r1, int r2);
double** matrix_newSubcol (double** A, int m, int n, int c1, int c2);
double** matrix_newSubmat (double** A, int m, int n, int r1, int r2, int c1, int c2);
double** matrix_delRow (double** A, int m, int n, int row);
double** matrix_delCol (double** A, int m, int n, int col);
double** matrix_cat (double** A, double** B, int m1, int m2, int n);
double** matrix_paste (double** A, double** B, int m, int n1, int n2);
void matrix_zero (double**A, int m, int n);
void matrix_multiply (double** A, int m, int n, double scale);
double** matrix_product (double** A, double** B, int m, int n, int o);
double** matrix_product3 (double** A, double** B, double** C, int m, int n, int o, int p);
void matrix_add (double** A, int m, int n, double a);
void matrix_sum (double** A, double** B, int m, int n);
void matrix_subtract (double** A, double** B, int m, int n);
void matrix_stdRow (double** A, int m, int n);
void matrix_stdCol (double** A, int m, int n);
double* matrix_rowMean (double** A, int m, int n);
double* matrix_colMean (double** A, int m, int n);
double* matrix_rowVar (double** A, int m, int n);
double* matrix_colVar (double** A, int m, int n);
void matrix_sortByVector(double** A, int m, int n, double* v);

void vector_zero (double* A, int n);
double* vector_copy (double *A, int n);
int* vectori_copy (int *A, int n);
double vector_innerProduct (double *A, double *B, int n);
double** vector_outerProduct (double *A, double *B, int n);
void vector_normalize (double* A, int n);
double** vector_toMatrix (double* A, int n);
double* matrix_toVector (double** A, int n);

/**************************
Advanced matrix operations
**************************/
double** matrix_corMatrix (double** A, int m, int n);
double** matrix_covMatrix (double** A, int m, int n);
void matrix_solve(double **a, int n, double **b, int m);
void matrix_ludecomp(double **a, int n, int *indx, double *d);
void matrix_lusubst(double **a, int n, int *indx, double* b);
void matrix_inverse (double **a, int n);
double matrix_det (double **a, int n);
void matrix_ludecomp2 (const int n, double **a, int *p);
void matrix_lusubst2(const int n, double **a, int *p, double *b);
void matrix_inverse2 (double **a, int n);
void matrix_tridiag(const int n, double *a, double *d, double *c, double *b);
double pythag(double a, double b);
void matrix_svd(double **a, int m, int n, double *w, double **v);
void matrix_svd2 (double *W, double *Z, int nRow, int nCol);
void matrix_eigen(double **a, int n, double *d, double **v, int *nrot);
double** matrix_pcaCor (double** A, int m, int n, double* eval);
double** matrix_pcaCov (double** A, int m, int n, double* eval);
double** matrix_pcaScore (double** A, int m, int n, double** evec, int k);
double* matrix_pcaComVar (double** A, int m, int n, double** evec);
double** matrix_reg (double** X, int m, int n, double** Y);
double matrix_regMSE (double** X, int m, int n, double** Y, double** B);
double* matrix_regModelStat (double** X, int m, int n, double** Y, double** B);
double** matrix_regCoefStat (double** X, int m, int n, double** Y, double** B);
double** matrix_regObsStat (double** X, int m, int n, double** Y, double** B);
void test1 ();



/***************************
statistical distributions
**************************/
double cdf_binomial (int n, int k, double p);
double cdf_chi2 (double df, double x);
double cdf_poisson (int k, double x);
double cdf_poisson_inc (int k, double x);
double cdf_normal (double x, double mu, double sigma);
double cdf_stdnormal (double x);
double cdf_f (double df1, double df2, double x);
double cdf_t (double df, double x);
double cdfinv_t (double df, double p);

double pdf_stdnormal (double x);
double pdf_normal (double x, double mu, double sigma);
double pdf_binomial (int n, int k, double p);
double pdf_beta (double a, double b, double x);
double pdf_poisson (int k, double x);
double pdf_geometric (int k, double p);
double pdf_hypergeometric (int a, int b, int c, int d);
double bico (int n, int k);
double lnbico (int n, int k);
double invnormal (double x);
void reg_linear (double *x, double *y, int ndata, double *a, double *b, double *F, double *P);



/*********************
summary statistics
**********************/
double mean(double *data, int n);
double mean2(double *data, int n);
void avevar(double *data, int n, double *ave, double *var);
void avevar2(double *data, int n, double *ave, double *var);
void averms(double *data, int n, double *ave, double *rms);
void moment(double *data, int n, double *ave, double *adev, double *sdev, double *var, double *skew, double *curt);
double cc(double* x, double* y, int n);
double cov(double* x, double* y, int n);

/********************
statistical tests
********************/
double bitest (int n, int k, double proportion);
void ttest_onesample (double *data, long n, double expected, double *t, double *p);
void ttest_ev(double *data1, long n1, double *data2, long n2, double *t, double *p);
void ttest_uev(double *data1, long n1, double *data2, long n2, double *t, double *p);
void ftest(double *data1, long n1, double *data2, long n2, double *f, double *p);
void chi2test(double *bins1, double *bins2, int nbins, int knstrn, double *df, double *chi2, double *p);
void chi2test_onesample(double *bins, double *ebins, int nbins, int knstrn, double *df, double *chi2, double *p);

void chi2test_trend_2by3table (double *bins, double *chi2, double *p);
void chi2test_trend_3by2table(double *bins, double *chi2, double *p);

void chi2test_2by2table(double *bins, double *chi2, double *p);
void chi2test_3by2table(double *bins, double *chi2, double *p);
void chi2test_2by3table(double *bins, double *chi2, double *p);

void kstest(double *data1, long n1, double *data2, long n2, double *d, double *prob);
void kstest_onesample(double *data, long n, double (*func)(double), double *d, double *prob);
double probks(double alam);

double fisher_exact_1sided (int a, int b, int c, int d);
double fisher_exact_2sided (int a, int b, int c, int d);

void quick_sort(int elements, double *arr);


/*******************
mathematical functions
********************/

double gammln(double x);
double lnfactorial (double x);
int factorial (int x);
double gammp(double a, double x);
double gammq(double a, double x);
void gser(double *gamser, double a, double x, double *gln);
void gcf(double *gammcf, double a, double x, double *gln);
double beta(double z, double w);
double betai(double a, double b, double x);
double betacf(double a, double b, double x);
double inverff (double x);
double erf(double x);

/****************
random number
****************/
double ran ();
double ran3(int *idum);
double random_exp(double lamda, int *idum);
double random_stdnormal (int *idum);


/******************
3rd party code
******************/
double SNPHWE(int obs_hets, int obs_hom1, int obs_hom2);
