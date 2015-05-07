%include typemaps.i

%module khmm
%{

#include "khmm.h"
#include "kc.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

%}

%{
typedef double * arr1d;
typedef float * arr1f;
typedef int * arr1i;
%}


%typemap(in) arr1d {
  $1 = (double *) pack1D($input,'d');
}
%typemap(in) arr1f {
  $1 = (float *) pack1D($input,'f');
}
%typemap(in) arr1i {
  $1 = (int *) pack1D($input,'i');
}

%typemap(argout) arr1d {
  unpack1D((SV*)$input, (void *)$1, 'd', 0);
}
%typemap(argout) arr1f {
  unpack1D((SV*)$input, (void *)$1, 'f', 0);
}
%typemap(argout) arr1i {
  unpack1D((SV*)$input, (void *)$1, 'i', 0);
}

%typemap(in) double * {
	double dvalue;
	SV* tempsv;
	if (!SvROK($input)) {
		croak("expected a reference\n");
	}
	tempsv = SvRV($input);
	if (SvTRUE(tempsv)) {
		if ((!SvNOK(tempsv)) && (!SvIOK(tempsv))) {
			croak("expected a double reference\n");
		}
	}
	dvalue = SvNV(tempsv);
	$1 = &dvalue;
}

%typemap(argout) double * {
	SV *tempsv;
	tempsv = SvRV($input);
	sv_setnv(tempsv, *$1);
}
%typemap(in) float * {
	float dvalue;
	SV* tempsv;
	if (!SvROK($input)) {
		croak("expected a reference\n");
	}
	tempsv = SvRV($input);
	if (SvTRUE(tempsv)) {
		if ((!SvNOK(tempsv)) && (!SvIOK(tempsv))) {
			croak("expected a double reference\n");
		}
	}
	dvalue = (float) SvNV(tempsv);
	$1 = &dvalue;
}

%typemap(argout) float * {
	SV *tempsv;
	tempsv = SvRV($input);
	sv_setnv(tempsv, *$1);
}
%typemap(in) int * {
	int ivalue;
	SV *tempsv;
	if (!SvROK($input)) {
		croak("expected a reference\n");
	}
	tempsv = SvRV($input);
	if (SvTRUE(tempsv)) {
		if (!SvIOK(tempsv)) {
			croak("expected a int reference\n");
		}
	}
	ivalue = SvIV(tempsv);
	$1 = &ivalue;
}

%typemap(argout) int * {
	SV *tempsv;
	tempsv = SvRV($input);
	sv_setiv(tempsv, *$1);
}


CHMM ReadCHMM (char *filename);
void estHMMFromFile_CHMM (CHMM hmm, int T, FILE *fp, int *niter, double *logprobinit, double *logprobfinal);
void testVit_CHMM (CHMM hmm, int T, arr1d O1, arr1d O2, arr1d pfb, arr1i snpdist, double *plogproba);
void PrintCHMM(FILE *fp, CHMM *pchmm);
void GetStateProb_CHMM(CHMM *phmm, int T, arr1d O1, arr1d O2, arr1d pfb, arr1i snpdist, double *pprob, int state);
void adjustBSD (CHMM *phmm, double sdo);
void FreeCHMM(CHMM *phmm);
void testVitTrio_CHMM (CHMM *phmm, int T, arr1d Of1, arr1d Of2, arr1d Om1, arr1d Om2, arr1d Oo1, arr1d Oo2, arr1d pfb, arr1i snpdist, double *plogproba, arr1d trio_lrr_sd, int osex, int chrx_flag);
void reg_linear (arr1d x, arr1d y, int ndata, double *a, double *b, double *F, double *P);


/*the following subroutines are built-in functions in C language*/
FILE *fopen (const char *filename, const char *mode);
int fclose (FILE *fp);

FILE *fh_stdout ();
double fisher_exact_2sided (int a, int b, int c, int d);
double fisher_exact_1sided (int a, int b, int c, int d);
double bitest (int n, int k, double proportion);

void callCNVFromFile_SEQ (CHMM hmm, int T, FILE *fp, arr1i mlstate, double gamma_k, double gamma_theta);
void adjustHMMExpected (CHMM *phmm, double coverage);