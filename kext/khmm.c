#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "khmm.h"
#include "kc.h"

#define STATE_CHANGE 100000.0			/*this is the expected changes (D value) in the transition matrix*/
#define VITHUGE  100000000000.0
#define FLOAT_MINIMUM 1.175494351e-38;		/*this is indeed machine dependent*/
#define DELTA 1

/*	This file was re-written from several subroutines from the UMDHMM package by Tapas Kanungo (Date: 15 December 1997), which has excellent framework of the implementation of Forward-Backward, Viterbi, and Baum-Welch algorithms.
	The original UMDHMM package was downloaded from http://www.kanungo.com/software/software.html. The citation for the UMDHMM program is "UMDHMM: Hidden Markov Model Toolkit," in "Extended Finite State Models of Language," A. Kornai (editor), Cambridge University Press, 1999."
	The basic framework (including variable name, subroutine name) is highly similar to the original UMDHMM package, but the actual implementation is completely different as no "discrete symbol emission" is used in PennCNV.
*/

void testVit_CHMM (CHMM hmm, int T, double *O1, double *O2, double *pfb, int *snpdist, double *plogproba)
{
	int i;
	int	*q;						/* state sequence q[1..T] */
	double **delta;
	int	**psi;

	q = ivector(1,T);
	delta = dmatrix(1, T, 1, hmm.N);
	psi = imatrix(1, T, 1, hmm.N);

	ViterbiLogNP_CHMM(&hmm, T, O1, O2, pfb, snpdist, delta, psi, q, plogproba);
	
	for (i=1; i<=T; i++) pfb[i]=q[i];			/*assign the most likely state value to return*/
	
	free_ivector(q, 1, T);
	free_imatrix(psi, 1, T, 1, hmm.N);
	free_dmatrix(delta, 1, T, 1, hmm.N);
}

void tumorVit_CHMM (CHMM hmm, int T, double *O1, double *O2, double *pfb, int *snpdist, double *plogproba, double stroma)
{
	int i;
	int	*q;						/* state sequence q[1..T] */
	double **delta;
	int	**psi;

	q = ivector(1,T);
	delta = dmatrix(1, T, 1, hmm.N);
	psi = imatrix(1, T, 1, hmm.N);

	ViterbiLogNPStroma_CHMM(&hmm, T, O1, O2, pfb, snpdist, delta, psi, q, plogproba, stroma);
	
	for (i=1; i<=T; i++) pfb[i]=q[i];			/*assign the most likely state value to return*/
	
	free_ivector(q, 1, T);
	free_imatrix(psi, 1, T, 1, hmm.N);
	free_dmatrix(delta, 1, T, 1, hmm.N);
}

void estHMMFromFile_CHMM (CHMM hmm, int T, FILE *fp, int *niter, double *logprobinit, double *logprobfinal)
{
	double 	**alpha; 
	double	**beta;
	double	**gamma;
	double *O1;
	double *O2;
	double *pfb;
	int *snpdist;
	int i;

	alpha = dmatrix(1, T, 1, hmm.N);
	beta = dmatrix(1, T, 1, hmm.N);
	gamma = dmatrix(1, T, 1, hmm.N);
	
	O1 = dvector (1, T);
	O2 = dvector (1, T);
	pfb = dvector (1, T);
	snpdist = ivector (1, T);
	
	for (i=1;i<=T;i++) {
		if (fscanf(fp, "%lf\t", &(O1[i])) != 1) kcerror ("ERROR: cannot read O1 from file");
		if (fscanf(fp, "%lf\t", &(O2[i])) != 1) kcerror ("ERROR: cannot read O2 from file");
		if (fscanf(fp, "%lf\n", &(pfb[i])) != 1) kcerror ("ERROR: cannot read PFB from file");
		if (fscanf(fp, "%i\n", &(snpdist[i])) != 1) kcerror ("ERROR: cannot read SNPDIST from file");
	}

	BaumWelchNP_CHMM (&hmm, T, O1, O2, pfb, snpdist, alpha, beta, gamma, niter, logprobinit, logprobfinal);

	free_dmatrix(alpha, 1, T, 1, hmm.N);
	free_dmatrix(beta, 1, T, 1, hmm.N);
	free_dmatrix(gamma, 1, T, 1, hmm.N);
	free_dvector(O1, 1, T);
	free_dvector(O2, 1, T);
	free_dvector(pfb, 1, T);
	free_ivector(snpdist, 1, T);
}

double b1iot (int state, double *mean, double *sd, double uf, double o)
{
	double p = 0;
	p = uf;
	p += (1-uf) * pdf_normal (o, mean[state], sd[state]);

	if (p==0) p=FLOAT_MINIMUM;
	return log(p);
}

double b2iot (int state, double *mean, double *sd, double uf, double pfb, double b)
{
	double p = 0;
	double mean0  = mean[1];
	double mean25 = mean[2];
	double mean33 = mean[3];
	double mean50 = mean[4];
	double mean50_state1 = mean[5];
	double sd0  = sd[1];
	double sd25 = sd[2];
	double sd33 = sd[3];
	double sd50 = sd[4];
	double sd50_state1 = sd[5];

	p = uf;
	if (state == 1) {
		if (b==0) {
			p+= (1-uf) * cdf_normal (0, mean50_state1, sd50_state1);
		} else if (b==1) {
			p+= (1-uf) * cdf_normal (0, mean50_state1, sd50_state1);
		} else {
			p+= (1-uf) * pdf_normal (b, mean50_state1, sd50_state1);
		}
	} else if (state == 2) {
		if (b==0) {
			p+= (1-uf) * (1-pfb)/2;
		} else if (b==1) {
			p+= (1-uf) * pfb/2;
		} else {
			p+= (1-uf) * (1-pfb) * pdf_normal (b, mean0, sd0);
			p+= (1-uf) * pfb     * pdf_normal (b, 1-mean0, sd0);
		}
	} else if (state == 3) {
		if (b==0) {
			p+= (1-uf) * (1-pfb)*(1-pfb)/2;
		} else if (b==1) {
			p+= (1-uf) * pfb*pfb/2;
		} else {
			p+= (1-uf) * (1-pfb)*(1-pfb) * pdf_normal (b, mean0, sd0);
			p+= (1-uf) * 2*pfb*(1-pfb)   * pdf_normal (b, mean50, sd50);
			p+= (1-uf) * pfb*pfb         * pdf_normal (b, 1-mean0, sd0);
		}
	} else if (state == 4) {
		if (b==0) {
			p+= (1-uf) * (1-pfb)/2;
		} else if (b==1) {
			p+= (1-uf) * pfb/2;
		} else {
			p+= (1-uf) * (1-pfb) * pdf_normal (b, mean0, sd0);
			p+= (1-uf) * pfb     * pdf_normal (b, 1-mean0, sd0);
		}
	} else if (state == 5) {
		if (b==0) {
			p+= (1-uf) * (1-pfb)*(1-pfb)*(1-pfb)/2;
		} else if (b==1) {
			p+= (1-uf) * pfb*pfb*pfb/2;
		} else {
			p+= (1-uf) * (1-pfb)*(1-pfb)*(1-pfb)   * pdf_normal (b, mean0, sd0);
			p+= (1-uf) * 3*(1-pfb)*(1-pfb)*pfb     * pdf_normal (b, mean33, sd33);
			p+= (1-uf) * 3*(1-pfb)*pfb*pfb         * pdf_normal (b, 1-mean33, sd33);
			p+= (1-uf) * pfb*pfb*pfb               * pdf_normal (b, 1-mean0, sd0);
		}
	} else if (state == 6) {
		if (b==0) {
			p+= (1-uf) *(1-pfb)*(1-pfb)*(1-pfb)*(1-pfb)/2;
		} else if (b==1) {
			p+= (1-uf) * pfb*pfb*pfb*pfb/2;
		} else {
			p+= (1-uf) * (1-pfb)*(1-pfb)*(1-pfb)*(1-pfb)   * pdf_normal (b, mean0, sd0);
			p+= (1-uf) * 4*(1-pfb)*(1-pfb)*(1-pfb)*pfb     * pdf_normal (b, mean25, sd25);
			p+= (1-uf) * 6*(1-pfb)*(1-pfb)*pfb*pfb         * pdf_normal (b, mean50, sd50);
			p+= (1-uf) * 4*(1-pfb)*pfb*pfb*pfb             * pdf_normal (b, 1-mean25, sd25);
			p+= (1-uf) * pfb*pfb*pfb*pfb                   * pdf_normal (b, 1-mean0, sd0);
		}
	}
	if (p==0) p=FLOAT_MINIMUM;
	return log(p);
}

double b5iot (int state, double *mean, double *sd, double uf, double o, double stroma)
/* calculate LRR */
{
	double p = 0;
	double lrrmean, lrrsd;		/* mean and SD of the LRR values for a given state and given stroma contamination value */
	lrrmean = (1-stroma)*mean[state] + stroma*mean[3];
	lrrsd = (1-stroma)*sd[state] + stroma*sd[3];
	
	p = uf;
	p += (1-uf) * pdf_normal (o, lrrmean, lrrsd);

	if (p==0) p=FLOAT_MINIMUM;
	return log(p);
}

double b6iot (int state, double *mean, double *sd, double uf, double pfb, double b, double stroma)
/* calculate BAF */
{
	double p = 0;
	double mean0  = mean[1];		/*mean1=0 mean2=.25 mean3=.33 mean4=.5 */
	double mean25 = ( (1-stroma) + stroma*mean[1] ) / ( (1-stroma)*4 + stroma * 2);
	double mean33 = ( (1-stroma) + stroma*mean[1] ) / ( (1-stroma)*3 + stroma * 2);
	double mean50 = mean[4];
	double mean50_state1 = mean[5];
	double sd0  = sd[1];
	double sd25 = (1-stroma)*sd[2]+stroma*sd[1];
	double sd33 = (1-stroma)*sd[3]+stroma*sd[1];
	double sd50 = sd[4];
	double sd50_state1 = sd[5];

	p = uf;
	if (state == 1) {
		if (b==0) {
			p+= (1-uf) * cdf_normal (0, mean50_state1, sd50_state1);
		} else if (b==1) {
			p+= (1-uf) * cdf_normal (0, mean50_state1, sd50_state1);
		} else {
			p+= (1-uf) * pdf_normal (b, mean50_state1, sd50_state1);
		}
	} else if (state == 2) {
		if (b==0) {
			p+= (1-uf) * (1-pfb)/2;
		} else if (b==1) {
			p+= (1-uf) * pfb/2;
		} else {
			p+= (1-uf) * (1-pfb) * pdf_normal (b, mean0, sd0);
			p+= (1-uf) * pfb     * pdf_normal (b, 1-mean0, sd0);
		}
	} else if (state == 3) {
		if (b==0) {
			p+= (1-uf) * (1-pfb)*(1-pfb)/2;
		} else if (b==1) {
			p+= (1-uf) * pfb*pfb/2;
		} else {
			p+= (1-uf) * (1-pfb)*(1-pfb) * pdf_normal (b, mean0, sd0);
			p+= (1-uf) * 2*pfb*(1-pfb)   * pdf_normal (b, mean50, sd50);
			p+= (1-uf) * pfb*pfb         * pdf_normal (b, 1-mean0, sd0);
		}
	} else if (state == 4) {		//the hypothesis here is that stroma do not have LOH
		if (b==0) {
			p+= (1-uf) * (1-pfb)/2;
		} else if (b==1) {
			p+= (1-uf) * pfb/2;
		} else {
			//p+= (1-uf) * (1-pfb) * pdf_normal (b, mean0, sd0);
			//p+= (1-uf) * pfb     * pdf_normal (b, 1-mean0, sd0);
			/* there are 3 possibilities, AA, AB and BB, in the stroma*/
			p+= (1-uf)* (1-pfb)*(1-pfb) * pdf_normal (b, mean0, sd0);
			p+= (1-uf) * pfb*(1-pfb)* pdf_normal (b, (1-stroma)*mean0+stroma*mean50, (1-stroma)*sd0+stroma*sd50);		//losing B allele
			p+= (1-uf) * pfb*(1-pfb)* pdf_normal (b, 1-(1-stroma)*mean0-stroma*mean50, (1-stroma)*sd0+stroma*sd50);		//losing A allele
			p+= (1-uf)* pfb*pfb * pdf_normal(b, 1-mean0, sd0);
		}
	} else if (state == 5) {		//the hypothesis here is that stroma has two copies
		if (b==0) {
			p+= (1-uf) * (1-pfb)*(1-pfb)*(1-pfb)/2;
		} else if (b==1) {
			p+= (1-uf) * pfb*pfb*pfb/2;
		} else {
			//p+= (1-uf) * (1-pfb)*(1-pfb)*(1-pfb)   * pdf_normal (b, mean0, sd0);
			//p+= (1-uf) * 3*(1-pfb)*(1-pfb)*pfb     * pdf_normal (b, mean33, sd33);
			//p+= (1-uf) * 3*(1-pfb)*pfb*pfb         * pdf_normal (b, 1-mean33, sd33);
			//p+= (1-uf) * pfb*pfb*pfb               * pdf_normal (b, 1-mean0, sd0);
			
			p+= (1-uf)*(1-pfb)*(1-pfb) * pdf_normal (b, mean0, sd0);		//AAA genotype (but we believe that it is derived from AA genotype, so it is (1-pfb)^2
			p+= (1-uf)*(1-pfb)*pfb    * pdf_normal (b, mean33, sd33);		//ABB genotype derived from AB genotype
			p+= (1-uf)*(1-pfb)*pfb    * pdf_normal (b, 1-mean33, sd33);		//ABB genotype derived from AB genotype
			p+= (1-uf)*pfb*pfb        * pdf_normal (b, 1-mean0, sd0);		//BBB genotype derived from AB genotype
			
		}
	} else if (state == 6) {
		if (b==0) {
			p+= (1-uf) *(1-pfb)*(1-pfb)*(1-pfb)*(1-pfb)/2;
		} else if (b==1) {
			p+= (1-uf) * pfb*pfb*pfb*pfb/2;
		} else {
			//p+= (1-uf) * (1-pfb)*(1-pfb)*(1-pfb)*(1-pfb)   * pdf_normal (b, mean0, sd0);
			//p+= (1-uf) * 4*(1-pfb)*(1-pfb)*(1-pfb)*pfb     * pdf_normal (b, mean25, sd25);
			//p+= (1-uf) * 6*(1-pfb)*(1-pfb)*pfb*pfb         * pdf_normal (b, mean50, sd50);
			//p+= (1-uf) * 4*(1-pfb)*pfb*pfb*pfb             * pdf_normal (b, 1-mean25, sd25);
			//p+= (1-uf) * pfb*pfb*pfb*pfb                   * pdf_normal (b, 1-mean0, sd0);
			p+= (1-uf) * (1-pfb)*(1-pfb)  * pdf_normal (b, mean0, sd0);		//AAAA genotype derived from AA genotype
			p+= (1-uf)*(1-pfb)*pfb   *0.5 * pdf_normal (b, mean25, sd25);		//ABBB genotype derived from AB genotype
			p+= (1-uf)*(1-pfb)*pfb        * pdf_normal (b, 1-mean50, sd50);		//AABB genotype derived from AB genotype
			p+= (1-uf)*(1-pfb)*pfb   *0.5 * pdf_normal (b, 1-mean25, sd25);		//AAAB genotype derived from AB genotype
			p+= (1-uf)*pfb*pfb        * pdf_normal (b, 1-mean0, sd0);		//BBBB genotype derived from BB genotype
		}
	}
	if (p==0) p=FLOAT_MINIMUM;
	return log(p);
}

void GetStateProb_CHMM(CHMM *phmm, int T, double *O1, double *O2, double *pfb, int *snpdist, double *pprob, int state)
{
	double b1iot (int state, double *mean, double *sd, double uf, double o);
	double b2iot (int state, double *mean, double *sd, double uf, double pfb, double b);
	int t;
	double logprob = 0;
	double **A1;
	
	A1 = dmatrix (1, phmm->N, 1, phmm->N);
	snpdist[0] = 5000;				/*arbitrarily set the transition probability to 5000*/
	for (t=1;t<=T;t++) {
		if (O2[t] > 1) {
			if (!phmm->NP_flag) nrerror ("FATAL ERROR: CN probe detected but HMM model does not contain parameters for them\n");
			logprob += b1iot (state, phmm->B1_mean, phmm->B1_sd, phmm->B1_uf, O1[t]);
		} else {
			logprob += b1iot (state, phmm->B1_mean, phmm->B1_sd, phmm->B1_uf, O1[t]);
			logprob += b2iot (state, phmm->B2_mean, phmm->B2_sd, phmm->B2_uf, pfb[t], O2[t]);
		}
		/*printf ("state=%i t=%d logprob=%f b1iot=%f b2iot=%f snpdist=%i transition_raw=%f pfb=%f O2=%f\n", state, t, logprob, b1iot (state, phmm->B1_mean, phmm->B1_sd, phmm->B1_uf, O1[t]), b2iot (state, phmm->B2_mean, phmm->B2_sd, phmm->B2_uf, pfb[t], O2[t]), snpdist[t-1], A1[state][state], pfb[t], O2[t]);*/
	}
	*pprob = logprob;
	free_dmatrix (A1, 1, phmm->N, 1, phmm->N);
}	

void ViterbiLogNP_CHMM(CHMM *phmm, int T, double *O1, double *O2, double *pfb, int *snpdist, double **delta, int **psi, int *q, double *pprob)
{
        int     i, j;   /* state indices */
        int     t;      /* time index */
 
 	int snp_count = 0;
 	int cn_count = 0;

        int     maxvalind;
        double  maxval, val;
	double  **biot;
	double **A1;

	A1 = dmatrix (1, phmm->N, 1, phmm->N);					/*initialize A1 matrix*/
	for (i = 1; i <= phmm->N; i++) {
		for (j = 1; j <= phmm->N; j++) {
			A1[i][j] = phmm->A[i][j];
		}
	}

	/* 0. Preprocessing */

	for (i = 1; i <= phmm->N; i++) {
		if (phmm->pi[i] == 0) phmm->pi[i] = 1e-9;			/*eliminate problems with zero probability*/
		phmm->pi[i] = log(phmm->pi[i]);
	}

	biot = dmatrix(1, phmm->N, 1, T);
	for (i = 1; i <= phmm->N; i++) {
		for (t = 1; t <= T; t++) {
			if (O2[t] > 1) {					/*normally BAF>=0 and BAF<=1; we use BAF>1 to indicate a Non-Polymorphic marker*/
				if (!phmm->NP_flag) nrerror ("FATAL ERROR: CN probe detected but HMM model does not contain parameters for them\n");
				biot[i][t] = b1iot (i, phmm->B3_mean, phmm->B3_sd, phmm->B3_uf, O1[t]);
				cn_count++;
			} else {						/*a regular SNP marker; we use both LRR and BAF information to calculate logProb for the marker*/
				biot[i][t] = b1iot (i, phmm->B1_mean, phmm->B1_sd, phmm->B1_uf, O1[t]);
				biot[i][t] += b2iot (i, phmm->B2_mean, phmm->B2_sd, phmm->B2_uf, pfb[t], O2[t]);
				snp_count++;
			}
		}
	}
	/*fprintf(stderr, "NOTICE: Encounterd %i snp probes and %i cn probes\n", snp_count, cn_count);*/

        /* 1. Initialization  */
 
        for (i = 1; i <= phmm->N; i++) {
                delta[1][i] = phmm->pi[i] + biot[i][1];
                psi[1][i] = 0;
        }
 
        /* 2. Recursion */
 
        for (t = 2; t <= T; t++) {

		if (phmm->dist != 1) convertHMMTransition (phmm, A1, snpdist[t-1]);			/*t-1 is used because the current val is calculated from previous values*/

                for (j = 1; j <= phmm->N; j++) {
                        maxval = -VITHUGE;
                        maxvalind = 1;
                        for (i = 1; i <= phmm->N; i++) {
                                val = delta[t-1][i] + log (A1[i][j]);
                                if (val > maxval) {
                                        maxval = val;
                                        maxvalind = i;
                                }
                        }
 
                        delta[t][j] = maxval + biot[j][t]; 
                        psi[t][j] = maxvalind;
                }
        }
 
        /* 3. Termination */
 
        *pprob = -VITHUGE;
        q[T] = 1;
        for (i = 1; i <= phmm->N; i++) {
                if (delta[T][i] > *pprob) {
                        *pprob = delta[T][i];
                        q[T] = i;
                }
        }
 
	/* 4. Path (state sequence) backtracking */

	for (t = T - 1; t >= 1; t--)
		q[t] = psi[t+1][q[t+1]];
	
	for (i = 1; i <= phmm->N; i++) {					/*recover the HMM model as original*/
		phmm->pi[i] = exp(phmm->pi[i]);
	}
	free_dmatrix(biot, 1, phmm->N, 1, T);
	free_dmatrix(A1, 1, phmm->N, 1, phmm->N);
}

void ViterbiLogNPStroma_CHMM(CHMM *phmm, int T, double *O1, double *O2, double *pfb, int *snpdist, double **delta, int **psi, int *q, double *pprob, double stroma)
{
        int     i, j;   /* state indices */
        int     t;      /* time index */
 
 	int snp_count = 0;
 	int cn_count = 0;

        int     maxvalind;
        double  maxval, val;
	double  **biot;
	double **A1;

	A1 = dmatrix (1, phmm->N, 1, phmm->N);					/*initialize A1 matrix*/
	for (i = 1; i <= phmm->N; i++) {
		for (j = 1; j <= phmm->N; j++) {
			A1[i][j] = phmm->A[i][j];
		}
	}

	/* 0. Preprocessing */

	for (i = 1; i <= phmm->N; i++) {
		if (phmm->pi[i] == 0) phmm->pi[i] = 1e-9;			/*eliminate problems with zero probability*/
		phmm->pi[i] = log(phmm->pi[i]);
	}

	biot = dmatrix(1, phmm->N, 1, T);
	for (i = 1; i <= phmm->N; i++) {
		for (t = 1; t <= T; t++) {
			if (O2[t] > 1) {					/*normally BAF>=0 and BAF<=1; we use BAF>1 to indicate a Non-Polymorphic marker*/
				if (!phmm->NP_flag) nrerror ("FATAL ERROR: CN probe detected but HMM model does not contain parameters for them\n");
				biot[i][t] = b5iot (i, phmm->B3_mean, phmm->B3_sd, phmm->B3_uf, O1[t], stroma);
				cn_count++;
			} else {						/*a regular SNP marker; we use both LRR and BAF information to calculate logProb for the marker*/
				biot[i][t] = b5iot (i, phmm->B1_mean, phmm->B1_sd, phmm->B1_uf, O1[t], stroma);
				biot[i][t] += b6iot (i, phmm->B2_mean, phmm->B2_sd, phmm->B2_uf, pfb[t], O2[t], stroma);
				snp_count++;
			}
		}
	}
	/*fprintf(stderr, "NOTICE: Encounterd %i snp probes and %i cn probes\n", snp_count, cn_count);*/

        /* 1. Initialization  */
 
        for (i = 1; i <= phmm->N; i++) {
                delta[1][i] = phmm->pi[i] + biot[i][1];
                psi[1][i] = 0;
        }
 
        /* 2. Recursion */
 
        for (t = 2; t <= T; t++) {

		if (phmm->dist != 1) convertHMMTransition (phmm, A1, snpdist[t-1]);			/*t-1 is used because the current val is calculated from previous values*/

                for (j = 1; j <= phmm->N; j++) {
                        maxval = -VITHUGE;
                        maxvalind = 1;
                        for (i = 1; i <= phmm->N; i++) {
                                val = delta[t-1][i] + log (A1[i][j]);
                                if (val > maxval) {
                                        maxval = val;
                                        maxvalind = i;
                                }
                        }
 
                        delta[t][j] = maxval + biot[j][t]; 
                        psi[t][j] = maxvalind;
                }
        }
 
        /* 3. Termination */
 
        *pprob = -VITHUGE;
        q[T] = 1;
        for (i = 1; i <= phmm->N; i++) {
                if (delta[T][i] > *pprob) {
                        *pprob = delta[T][i];
                        q[T] = i;
                }
        }
 
	/* 4. Path (state sequence) backtracking */

	for (t = T - 1; t >= 1; t--)
		q[t] = psi[t+1][q[t+1]];
	
	for (i = 1; i <= phmm->N; i++) {					/*recover the HMM model as original*/
		phmm->pi[i] = exp(phmm->pi[i]);
	}
	free_dmatrix(biot, 1, phmm->N, 1, T);
	free_dmatrix(A1, 1, phmm->N, 1, phmm->N);
}

void adjustBSD (CHMM *phmm, double sdo)
/*adjust the CHMM model so that the standard deviation of B1 and B2 match the observed data (by an empirical method)*/
{
	int i;
	double ratio = sdo / phmm->B1_sd[3];
	for (i = 1; i <= phmm->N; i++) {
		phmm->B1_sd[i] = phmm->B1_sd[i] * ratio;
	}
	for (i = 1; i <= 5; i++) {
		phmm->B2_sd[i] = phmm->B2_sd[i] * ratio;
	}
	if (phmm->NP_flag) {
		for (i = 1; i <= 5; i++) {
			phmm->B3_sd[i] = phmm->B3_sd[i] * ratio;
		}
	}
}

void adjustHMMExpected (CHMM *phmm, double coverage)
{
	int i;
	for (i = 1; i <= phmm->N; i++) {
		phmm->B1_mean[i] = phmm->B1_mean[i] * coverage;
		phmm->B1_sd[i] = phmm->B1_sd[i] * coverage;
	}
	if (phmm->NP_flag) {
		for (i = 1; i <= phmm->N; i++) {
			phmm->B3_mean[i] = phmm->B3_mean[i] * coverage;
			phmm->B3_sd[i] = phmm->B3_sd[i] * coverage;
		}
	}
}

void convertHMMTransition (CHMM *phmm, double **A1, int dist)
/*this subroutine convert HMM transition probabilities using the P=Pref*(1-exp(-d/D)) formula for off-diagonal cells in the matrix*/
{
	int i, j;
	double D = STATE_CHANGE;
	double offdiagonal_sum = 0;
	for (i = 1; i <= phmm->N; i++) {
		offdiagonal_sum = 0;
		for (j = 1; j <= phmm->N; j++) {
			if (i != j) {
				if (i == 4) {
					A1[i][j] = phmm->A[i][j] * (1-exp(-dist/D/1000)) / (1-exp(-5000/D/1000));
				} else {
					A1[i][j] = phmm->A[i][j] * (1-exp(-dist/D)) / (1-exp(-5000/D));
				}
				if (A1[i][j] > 1) {
					printf ("WARNING: Off-diagonal cell A[%i][%i] (%f to %f by %i) in transition matrix is over boundary of 1 (HMM model is not optimized). Assign 0.999 as the value instead.\n", i, j, phmm->A[i][j], A1[i][j], dist);
					A1[i][j] = 0.999;			/*maximum possible off-diagonal value (since state3 frequency is 0.999)*/
				}
				offdiagonal_sum += A1[i][j];
			}
		}
		if (offdiagonal_sum >= 1) {
			for (j = 1; j <= phmm->N; j++) {
				A1[i][j] /= (offdiagonal_sum/0.999);
			}
			offdiagonal_sum = 0.999;
		}
		A1[i][i] = 1 - offdiagonal_sum;
	}
}

void ForwardWithScale_CHMM (CHMM *phmm, int T, double *O1, double *O2, double **biot, int *snpdist, double **alpha, double *scale, double *pprob)
/*  pprob is the LOG probability */
{
	int	i, j; 	/* state indices */
	int	t;	/* time index */

	double sum;	/* partial sum */
	double **A1;
	A1 = dmatrix (1, phmm->N, 1, phmm->N);

	/* 1. Initialization */
	scale[1] = 0.0;	
	for (i = 1; i <= phmm->N; i++) {
		alpha[1][i] = phmm->pi[i]*biot[i][1];
		scale[1] += alpha[1][i];
	}
	for (i = 1; i <= phmm->N; i++) 
		alpha[1][i] /= scale[1]; 

	/* 2. Induction */
	for (t = 1; t <= T - 1; t++) {
		scale[t+1] = 0.0;
		
		if (phmm->dist != 1) convertHMMTransition (phmm, A1, snpdist[t]);

		for (j = 1; j <= phmm->N; j++) {
			sum = 0.0;
			for (i = 1; i <= phmm->N; i++) {
				if (A1[i][j]<0) {printf ("FATAL ERROR at i=%i j=%i t=%i\n", i, j, t); PrintCHMM (stdout, phmm);}
				sum += alpha[t][i] * A1[i][j];
			}
			alpha[t+1][j] = sum*biot[j][t+1];
			scale[t+1] += alpha[t+1][j];
		}
		for (j = 1; j <= phmm->N; j++) alpha[t+1][j] /= scale[t+1]; 
	}

	/* 3. Termination */
	*pprob = 0.0;
	free_dmatrix (A1, 1, phmm->N, 1, phmm->N);

	for (t = 1; t <= T; t++) *pprob += log(scale[t]);
	
}

void BackwardWithScale_CHMM (CHMM *phmm, int T, double *O1, double *O2, double **biot, int *snpdist, double **beta, double *scale, double *pprob)
{
	int     i, j;   /* state indices */
	int     t;      /* time index */
	double sum;
	double **A1;

	A1 = dmatrix (1, phmm->N, 1, phmm->N);
        /* 1. Initialization */
 
        for (i = 1; i <= phmm->N; i++)
                beta[T][i] = 1.0/scale[T]; 
 
        /* 2. Induction */
 
	for (t = T - 1; t >= 1; t--) {

		if (phmm->dist != 1) convertHMMTransition (phmm, A1, snpdist[t]);

		for (i = 1; i <= phmm->N; i++) {
			sum = 0.0;
                        for (j = 1; j <= phmm->N; j++) {
                        	sum += A1[i][j] * biot[j][t+1] * beta[t+1][j];
                        }
                        beta[t][i] = sum/scale[t];
		}
	}
	free_dmatrix (A1, 1, phmm->N, 1, phmm->N);
}

void BaumWelchNP_CHMM(CHMM *phmm, int T, double *O1, double *O2, double *pfb, int *snpdist, double **alpha, double **beta, double **gamma, int *pniter, double *plogprobinit, double *plogprobfinal)
{
	int	i, j;
	int	t, l = 0;

	double	logprobf, logprobb;
	double	numeratorA, denominatorA;

	double ***xi, *scale;
	double delta, deltaprev, logprobprev;
	double **A1;
	
	double param_change = 1e-9;			/*eliminate problems with zero probability in the cell*/

	/*the following paragraph prepares the biot matrix, which stores the emission probability for each i (model) at each time (t)*/
	double **biot;
	double *biotemp;
	double biosumtemp;
	double b1iot (int state, double *mean, double *sd, double uf, double o);
	double b2iot (int state, double *mean, double *sd, double uf, double pfb, double b);
	
	double inflation, inflation_s4, adjust=0.0;

	biot = dmatrix(1, phmm->N, 1, T);
	A1 = dmatrix(1, phmm->N, 1, phmm->N);
	biotemp = dvector(1, phmm->N);
	
	for (t = 1; t <= T; t++) {
		biosumtemp = 0;
		for (i = 1; i <= phmm->N; i++) {
			if (O2[t] > 1) {
				if (!phmm->NP_flag) nrerror ("FATAL ERROR: detected CN probe but HMM model does not contain info for CN probe\n");
				biotemp[i] = exp (b1iot (i, phmm->B3_mean, phmm->B3_sd, phmm->B3_uf, O1[t]));
			} else {
				biotemp[i] = exp (b1iot (i, phmm->B1_mean, phmm->B1_sd, phmm->B1_uf, O1[t]));
			}
			biosumtemp += biotemp[i];
		}
		for (i = 1; i <= phmm->N; i++) {
			biot[i][t] = biotemp[i]/biosumtemp;
		}
		biosumtemp = 0;
		for (i = 1; i <= phmm->N; i++) {
			if (O2[t] > 1) {
				biotemp[i] = 1;
			} else {
				biotemp[i] = exp (b2iot (i, phmm->B2_mean, phmm->B2_sd, phmm->B2_uf, pfb[t], O2[t]));
			}
			biosumtemp += biotemp[i];
		}
		for (i = 1; i <= phmm->N; i++) {
			biot[i][t] *= (biotemp[i] / biosumtemp);
		}
	}

	deltaprev = 10e-70;

	xi = AllocXi(T, phmm->N);
	scale = dvector(1, T);

	ForwardWithScale_CHMM (phmm, T, O1, O2, biot, snpdist, alpha, scale, &logprobf);
	*plogprobinit = logprobf;	/* log P(O |intial model) */
	BackwardWithScale_CHMM (phmm, T, O1, O2, biot, snpdist, beta, scale, &logprobb);
	ComputeGamma_CHMM (phmm, T, alpha, beta, gamma);
	ComputeXi_CHMM (phmm, T, O1, O2, biot, snpdist, alpha, beta, xi);
	logprobprev = logprobf;

	do  {	
		for (i = 1; i <= phmm->N; i++) { 
			denominatorA = 0.0;
			for (t = 1; t <= T - 1; t++) 
				denominatorA += gamma[t][i];

			for (j = 1; j <= phmm->N; j++) {
				numeratorA = 0.0;
				for (t = 1; t <= T - 1; t++) 
					numeratorA += xi[t][i][j];
				phmm->A[i][j] = param_change + (1-param_change)*numeratorA/denominatorA;	/*making sure that smallest value in each cell is param_change*/
			}
		}
		
		/* eliminate out-of-boundary error from transition matrix (i.e. making sure that when dist=100_000_000 all cells are still within 0-1 region) */
		
		inflation = floor (1e9 * (1-exp(-5000/STATE_CHANGE))) / 1e9;
		inflation_s4 = floor (1e9 * (1-exp(-5000/STATE_CHANGE/1000))) / 1e9;
		for (i = 1; i <= phmm->N; i++) {
			for (j = 1; j <= phmm->N; j++) {
				if (i !=j && ((phmm->A[i][j] >= inflation && i != 4) || (phmm->A[i][j] >= inflation_s4 && i == 4))) {
						if (i != 4) adjust = phmm->A[i][j] - inflation;
						if (i == 4) adjust = phmm->A[i][j] - inflation_s4;
						phmm->A[i][j] -= adjust;
						phmm->A[i][i] += adjust;
				}
			}
		}
		

		printf ("\n<---------------HMM: Expectation at iteration %d-------------------\n", l);
		PrintCHMM (stdout, phmm);
		printf ("-------------------------------------------------------------------->\n\n");

		ForwardWithScale_CHMM (phmm, T, O1, O2, biot, snpdist, alpha, scale, &logprobf);
		BackwardWithScale_CHMM (phmm, T, O1, O2, biot, snpdist, beta, scale, &logprobb);
		ComputeGamma_CHMM (phmm, T, alpha, beta, gamma);
		ComputeXi_CHMM (phmm, T, O1, O2, biot, snpdist, alpha, beta, xi);

		/* compute difference between log probability of two iterations, logprobf but not logprobb is used here */
		delta = logprobf - logprobprev; 
		logprobprev = logprobf;
		l++;
		printf("NOTICE: Finished Baum-Welch iteration=%i delta=%f current_logprobf=%f current_logprobb=%f\n", l, delta, logprobf, logprobb);
	}
	while (delta > DELTA);		/* if log probability does not change much, exit */ 
 
	*pniter = l;			/* number of iterations */
	*plogprobfinal = logprobf;	/* log P(O|estimated model) */
	FreeXi(xi, T, phmm->N);
	free_dvector(scale, 1, T);
	free_dvector(biotemp, 1, phmm->N);
	free_dmatrix(A1, 1, phmm->N, 1, phmm->N);
}

void ComputeGamma_CHMM (CHMM *phmm, int T, double **alpha, double **beta, double **gamma)
{

	int 	i, j;
	int	t;
	double	denominator;

	for (t = 1; t <= T; t++) {
		denominator = 0.0;
		for (j = 1; j <= phmm->N; j++) {
			gamma[t][j] = alpha[t][j]*beta[t][j];
			denominator += gamma[t][j];
		}
		for (i = 1; i <= phmm->N; i++) 
			gamma[t][i] = gamma[t][i]/denominator;
	}
}

void ComputeXi_CHMM (CHMM* phmm, int T, double *O1, double *O2, double **biot, int *snpdist, double **alpha, double **beta, double ***xi)
{
	int i, j;
	int t;
	double sum;
	double **A1;
	A1 = dmatrix (1, phmm->N, 1, phmm->N);

	for (t = 1; t <= T - 1; t++) {
		sum = 0.0;	

		if (phmm->dist != 1) convertHMMTransition (phmm, A1, snpdist[t]);

		for (i = 1; i <= phmm->N; i++) 
			for (j = 1; j <= phmm->N; j++) {
				xi[t][i][j] = alpha[t][i]*beta[t+1][j]
					*A1[i][j]
					*biot[j][t+1];
				sum += xi[t][i][j];
			}

		for (i = 1; i <= phmm->N; i++) 
			for (j = 1; j <= phmm->N; j++) 
				xi[t][i][j]  /= sum;
	}
	free_dmatrix (A1, 1, phmm->N, 1, phmm->N);
}

double *** AllocXi(int T, int N)
{
	int t;
	double ***xi;

	xi = (double ***) malloc(T*sizeof(double **));

	xi --;

	for (t = 1; t <= T; t++)
		xi[t] = dmatrix(1, N, 1, N);
	return xi;
}

void FreeXi(double *** xi, int T, int N)
{
	int t;



	for (t = 1; t <= T; t++)
		free_dmatrix(xi[t], 1, N, 1, N);

	xi ++;
	free(xi);

}

CHMM ReadCHMM (char *filename)
{
	FILE *fp;
	CHMM hmm;
	CHMM *phmm;
	int i, j, k;
	
	phmm = &hmm;
	fp = fopen (filename, "r");
	if (!fp) kcerror ("Error: cannot read from HMM file");
	
	if (fscanf(fp, "M=%d\n", &(phmm->M)) == EOF) kcerror ("khmm::ReadCHMM: cannot read M annotation from HMM file");
	if (fscanf(fp, "N=%d\n", &(phmm->N)) == EOF) kcerror ("khmm::ReadCHMM: cannot read N annotation from HMM file");

	if (fscanf(fp, "A:\n") == EOF) kcerror ("khmm::ReadCHMM: cannot read A annotation from HMM file");
	phmm->A = (double **) dmatrix(1, phmm->N, 1, phmm->N);
	for (i = 1; i <= phmm->N; i++) { 
		for (j = 1; j <= phmm->N; j++) {
			if (fscanf(fp, "%lf", &(phmm->A[i][j])) == EOF) kcerror ("khmm::ReadCHMM: cannot read A matrix from HMM file");
		}
		if (fscanf(fp,"\n") == EOF) kcerror ("khmm::ReadCHMM: cannot read return character from HMM file");
	}

	if (fscanf(fp, "B:\n") == EOF) kcerror ("khmm::ReadCHMM: cannot read B annotation from HMM file");
	phmm->B = (double **) dmatrix(1, phmm->N, 1, phmm->M);
	for (j = 1; j <= phmm->N; j++) { 
		for (k = 1; k <= phmm->M; k++) {
			if (fscanf(fp, "%lf", &(phmm->B[j][k])) == EOF) kcerror ("khmm::ReadCHMM: cannot read B matrix from HMM file");
		}
		if (fscanf(fp,"\n") == EOF) kcerror ("khmm::ReadCHMM: cannot read return character from HMM file");
	}

	if (fscanf(fp, "pi:\n") == EOF) kcerror ("khmm::ReadCHMM: cannot read PI annotation from HMM file");
	phmm->pi = (double *) dvector(1, phmm->N);
	for (i = 1; i <= phmm->N; i++) {
		if (fscanf(fp, "%lf", &(phmm->pi[i])) == EOF) kcerror ("khmm::ReadCHMM: cannot read PI vector from HMM file");
		if (phmm->pi[i] < 1e-6) phmm->pi[i] = 1e-6;
	}
	if (fscanf(fp,"\n") == EOF) kcerror ("khmm::ReadCHMM: cannot read return character from HMM file");
	
	if (fscanf(fp, "B1_mean:\n") == EOF) kcerror ("khmm::ReadCHMM: cannot read B1_mean annotation from HMM file");
	phmm->B1_mean = (double *) dvector(1, phmm->N);
	for (i = 1; i <= phmm->N; i++)
		if (fscanf(fp, "%lf", &(phmm->B1_mean[i])) == EOF) kcerror ("khmm::ReadCHMM: cannot read B1_mean vector from HMM file");
	if (fscanf(fp,"\n") == EOF) kcerror ("khmm::ReadCHMM: cannot read return character from HMM file");

	if (fscanf(fp, "B1_sd:\n") == EOF) kcerror ("khmm::ReadCHMM: cannot read B1_sd annotation from HMM file");
	phmm->B1_sd = (double *) dvector(1, phmm->N);
	for (i = 1; i <= phmm->N; i++)
		if (fscanf(fp, "%lf", &(phmm->B1_sd[i])) == EOF) kcerror ("khmm::ReadCHMM: cannot read B1_sd from HMM file");
	if (fscanf(fp,"\n") == EOF) kcerror ("khmm::ReadCHMM: cannot read return character from HMM file");
	
	if (fscanf(fp, "B1_uf:\n") == EOF) kcerror ("khmm::ReadCHMM: cannot read B1_uf annotation from HMM file");
	if (fscanf(fp, "%lf", &(phmm->B1_uf)) == EOF) kcerror ("khmm::ReadCHMM: cannot read B1_uf from HMM file");
	if (fscanf(fp,"\n") == EOF) kcerror ("khmm::ReadCHMM: cannot read return character from HMM file");

	if (fscanf(fp, "B2_mean:\n") == EOF) kcerror ("khmm::ReadCHMM: cannot read B2_mean annotation from HMM file");
	phmm->B2_mean = (double *) dvector(1, 5);
	for (i = 1; i <= 5; i++)
		if (fscanf(fp, "%lf", &(phmm->B2_mean[i])) == EOF) kcerror ("khmm::ReadCHMM: cannot read B2_mean from HMM file");
	if (fscanf(fp,"\n") == EOF) kcerror ("khmm::ReadCHMM: cannot read return character from HMM file");

	if (fscanf(fp, "B2_sd:\n") == EOF) kcerror ("khmm::ReadCHMM: cannot read B2_sd annotation from HMM file");
	phmm->B2_sd = (double *) dvector(1, 5);
	for (i = 1; i <= 5; i++)
		if (fscanf(fp, "%lf", &(phmm->B2_sd[i])) == EOF) kcerror ("khmm::ReadCHMM: cannot read B2_sd from HMM file");
	if (fscanf(fp,"\n") == EOF) kcerror ("khmm::ReadCHMM: cannot read return character from HMM file");

	if (fscanf(fp, "B2_uf:\n") == EOF) kcerror ("khmm::ReadCHMM: cannot read B2_uf annotation from HMM file");
	if (fscanf(fp, "%lf", &(phmm->B2_uf)) == EOF) kcerror ("khmm::ReadCHMM: cannot read B2_uf from HMM file");
	if (fscanf(fp,"\n") == EOF) kcerror ("khmm::ReadCHMM: cannot read return character from HMM file");
	
	if (fscanf(fp, "B3_mean:\n") != EOF) {
		phmm->NP_flag = 1;
		phmm->B3_mean = (double *) dvector (1, phmm->N);
		for (i = 1; i <= phmm->N; i++)
			if (fscanf(fp, "%lf", &(phmm->B3_mean[i])) == EOF) kcerror ("khmm::ReadCHMM: cannot read B3_mean from HMM file");
		if (fscanf(fp,"\n") == EOF) kcerror ("khmm::ReadCHMM: cannot read return character from HMM file");
		if (fscanf(fp, "B3_sd:\n") == EOF) kcerror ("khmm::ReadCHMM: cannot read B3_sd annotation from HMM file");
		phmm->B3_sd = (double *) dvector (1, phmm->N);
		for (i = 1; i <= phmm->N; i++)
			if (fscanf(fp, "%lf", &(phmm->B3_sd[i])) == EOF) kcerror ("khmm::ReadCHMM: cannot read B3_sd from HMM file");
		if (fscanf(fp,"\n") == EOF) kcerror ("khmm::ReadCHMM: cannot read return character from HMM file");
		if (fscanf(fp, "B3_uf:\n") == EOF) kcerror ("khmm::ReadCHMM: cannot read B3_uf annotation from HMM file");
		if (fscanf(fp, "%lf", &(phmm->B3_uf)) == EOF) kcerror ("khmm::ReadCHMM: cannot read B3_uf from HMM file");
		if (fscanf(fp,"\n") == EOF) kcerror ("khmm::ReadCHMM: cannot read return character from HMM file");
	} else {
		phmm->NP_flag = 0;
	}
	
	if (fscanf(fp, "DIST:\n") != EOF) {
		if (fscanf(fp, "%d", &(phmm->dist)) == EOF) kcerror ("khmm:ReadCHMM: cannot read DIST from HMM file");
	} else {
		phmm->dist = STATE_CHANGE;
	}
	
	fclose (fp);
	return hmm;
}

void FreeCHMM(CHMM *phmm)
{
	free_dmatrix(phmm->A, 1, phmm->N, 1, phmm->N);
	free_dmatrix(phmm->B, 1, phmm->N, 1, phmm->M);
	free_dvector(phmm->pi, 1, phmm->N);
	free_dvector(phmm->B1_mean, 1, phmm->N);
	free_dvector(phmm->B1_sd, 1, phmm->N);
	free_dvector(phmm->B2_mean, 1, 5);
	free_dvector(phmm->B2_sd, 1, 5);
	
	if (phmm->NP_flag) {
		free_dvector(phmm->B3_mean, 1, phmm->N);
		free_dvector(phmm->B3_sd, 1, phmm->N);
	}
}

void CopyCHMM(CHMM *phmm1, CHMM *phmm2)
{
        int i, j, k;
 
        phmm2->M = phmm1->M;
        phmm2->N = phmm1->N;
 
        phmm2->A = (double **) dmatrix(1, phmm2->N, 1, phmm2->N);
 
        for (i = 1; i <= phmm2->N; i++)
                for (j = 1; j <= phmm2->N; j++)
                        phmm2->A[i][j] = phmm1->A[i][j];
 
        phmm2->B = (double **) dmatrix(1, phmm2->N, 1, phmm2->M);
        for (j = 1; j <= phmm2->N; j++)
                for (k = 1; k <= phmm2->M; k++)
                        phmm2->B[j][k] = phmm1->B[j][k];
 
        phmm2->pi = (double *) dvector(1, phmm2->N);
        for (i = 1; i <= phmm2->N; i++)
                phmm2->pi[i] = phmm1->pi[i]; 

	phmm2->B1_mean = (double *) dvector(1, phmm2->N);
	for (i=1;i<=phmm2->N;i++) phmm2->B1_mean[i] = phmm1->B1_mean[i];
	phmm2->B1_sd = (double *) dvector(1, phmm2->N);
	for (i=1;i<=phmm2->N;i++) phmm2->B1_sd[i] = phmm1->B1_sd[i];
	phmm2->B1_uf = phmm1->B1_uf;
	phmm2->B2_mean = (double *) dvector(1, 5);
	for (i=1;i<=phmm2->N;i++) phmm2->B2_mean[i] = phmm1->B2_mean[i];
	phmm2->B2_sd = (double *) dvector(1, 5);
	for (i=1;i<=phmm2->N;i++) phmm2->B2_sd[i] = phmm1->B2_sd[i];
	phmm2->B2_uf = phmm1->B2_uf;
	
	if (phmm1->NP_flag) {
		phmm2->NP_flag = phmm1->NP_flag;
		phmm2->B3_mean = (double *) dvector(1, phmm2->N);
		phmm2->B3_sd = (double *) dvector(1, phmm2->N);
		for (i=1; i<=phmm1->N; i++) phmm2->B3_mean[i] = phmm1->B3_mean[i];
		for (i=1; i<= phmm1->N; i++) phmm2->B3_sd[i] = phmm1->B3_sd[i];
		phmm2->B3_uf = phmm1->B3_uf;
	} else {
		phmm2->NP_flag = 0;
	}
}

void PrintCHMM(FILE *fp, CHMM *phmm)
{
        int i, j, k;

	fprintf(fp, "M=%d\n", phmm->M); 
	fprintf(fp, "N=%d\n", phmm->N); 
 
	fprintf(fp, "A:\n");
        for (i = 1; i <= phmm->N; i++) {
                for (j = 1; j <= phmm->N; j++) {
                        fprintf(fp, "%.9f ", phmm->A[i][j] );
		}
		fprintf(fp, "\n");
	}
 
	fprintf(fp, "B:\n");
        for (j = 1; j <= phmm->N; j++) {
                for (k = 1; k <= phmm->M; k++){
                        fprintf(fp, "%f ", phmm->B[j][k]);
		}
		fprintf(fp, "\n");
	}
 
	fprintf(fp, "pi:\n");
        for (i = 1; i <= phmm->N; i++) {
		fprintf(fp, "%f ", phmm->pi[i]);
	}
	fprintf(fp, "\n");
	
	fprintf(fp, "B1_mean:\n");
	for (i = 1; i <= phmm->N; i++) {
		fprintf(fp, "%f ", phmm->B1_mean[i]);
	}
	fprintf(fp, "\n");
	
	fprintf(fp, "B1_sd:\n");
	for (i = 1; i <= phmm->N; i++) {
		fprintf(fp, "%f ", phmm->B1_sd[i]);
	}
	fprintf(fp, "\n");
	
	fprintf(fp, "B1_uf:\n");
	fprintf(fp, "%f", phmm->B1_uf);
	fprintf(fp, "\n");
	
	fprintf(fp, "B2_mean:\n");
	for (i = 1; i<= 5; i++) {
		fprintf(fp, "%f ", phmm->B2_mean[i]);
	}
	fprintf(fp, "\n");
	
	fprintf(fp, "B2_sd:\n");
	for (i = 1; i<= 5; i++) {
		fprintf(fp, "%f ", phmm->B2_sd[i]);
	}
	fprintf(fp, "\n");
	
	fprintf(fp, "B2_uf:\n");
	fprintf(fp, "%f", phmm->B2_uf);
	fprintf(fp, "\n");
	
	if (phmm->NP_flag) {
		fprintf(fp, "B3_mean:\n");
		for (i = 1; i <= phmm->N; i++)
			fprintf(fp, "%f ", phmm->B3_mean[i]);
		fprintf(fp, "\n");
		fprintf(fp, "B3_sd:\n");
		for (i = 1; i <= phmm->N; i++)
			fprintf(fp, "%f ", phmm->B3_sd[i]);
		fprintf(fp, "\n");
		fprintf(fp, "B3_uf:\n");
		fprintf(fp, "%f", phmm->B3_uf);
		fprintf(fp, "\n");
	}
}

double b4iot (int state, double *mean, double *sd, double uf, double pfb, int o1, int o2)
{
	double p = 0;
	
	p = uf/o1;
	if (p>uf) p=uf;
	if (pfb<0.01) pfb = 0.01;				/*set up the minimum PFB value as 0.01*/
		
	
	if (state == 1) {
		if (o1) p += pdf_binomial (o1, o2, 0.25);
	} else if (state == 2) {
		p+= (1-uf) * (1-pfb) * pdf_binomial (o1, o2, 0.04);
		p+= (1-uf) * pfb * pdf_binomial (o1, o2, 0.96);
	} else if (state == 3) {
		p += (1-uf) * (1-pfb) * (1-pfb) * pdf_binomial (o1, o2, 0.04);
		p += (1-uf) * 2 * pfb * (1-pfb) * pdf_binomial (o1, o2, 0.5);
		p += (1-uf) * pfb * pfb * pdf_binomial (o1, o2, 0.96);
	} else if (state == 4) {
		p += (1-uf) * (1-pfb) * (1-pfb)        * pdf_binomial (o1, o2, 0.04);
		p += (1-uf) * pfb * pfb                * pdf_binomial (o1, o2, 0.96);
	} else if (state == 5) {
		p+= (1-uf) * (1-pfb)*(1-pfb)*(1-pfb)   * pdf_binomial (o1, o2, 0.04);
		p+= (1-uf) * 3*(1-pfb)*(1-pfb)*pfb     * pdf_binomial (o1, o2, 0.3333);
		p+= (1-uf) * 3*(1-pfb)*pfb*pfb         * pdf_binomial (o1, o2, 0.6667);
		p+= (1-uf) * pfb*pfb*pfb               * pdf_binomial (o1, o2, 0.96);
	} else if (state == 6) {
		p+= (1-uf) * (1-pfb)*(1-pfb)*(1-pfb)*(1-pfb)   * pdf_binomial (o1, o2, 0.04);
		p+= (1-uf) * 4*(1-pfb)*(1-pfb)*(1-pfb)*pfb     * pdf_binomial (o1, o2, 0.25);
		p+= (1-uf) * 6*(1-pfb)*(1-pfb)*pfb*pfb         * pdf_binomial (o1, o2, 0.5);
		p+= (1-uf) * 4*(1-pfb)*pfb*pfb*pfb             * pdf_binomial (o1, o2, 0.75);
		p+= (1-uf) * pfb*pfb*pfb*pfb                   * pdf_binomial (o1, o2, 0.96);
	}
	if (p==0) p=FLOAT_MINIMUM;
	if (p>1.0) {
		fprintf(stderr, "state=%i o1=%i o2=%i p=%f uf=%f mean3=%f sd3=%f total is %f uffraction=%f\n", state, o1, o2, p, uf, mean[3], sd[3], pdf_binomial (o1, o2, 0.75), uf / sqrt (mean[3]) / sd[3] / 4);
		nrerror("ERROR");
	}
	return log(p);
}

double b3iot (int state, double *mean, double *sd, double uf, double o1, int length)
{
	double p = 0;
	double k, theta;
	
	k = mean[state];
	theta = sd[state];
	
	
	if (length > 36) {
		k = k*(double)length/36;			/*36 is the current default read length (will change in the future)*/
		theta = theta/(double)length*36;
	}
	
	p = uf/o1;
	if (p>uf) p=uf;
	p += exp (log(1-uf) + (k-1)*log(o1) - o1/theta - k*log(theta) - gammln(k));

	if (p==0) p=FLOAT_MINIMUM;
	return log(p);
}

void callCNVFromFile_SEQ (CHMM hmm, int T, FILE *fp, int *mlstate, double gamma_k, double gamma_theta)
{
	double *O1;
	int *O2;
	double *pfb;
	int *length;
	int i;
	int	*q;						/* state sequence q[1..T] */
	double **delta;
	int	**psi;
	char string[1000];
	double bpcovcount=0, bpcount=0, bplncovcount=0;
	double s, k, theta;

	O1 = dvector (1, T);
	O2 = ivector (1, T);
	pfb = dvector (1, T);
	length = ivector (1, T);
	q = ivector(1,T);
	delta = dmatrix(1, T, 1, hmm.N);
	psi = imatrix(1, T, 1, hmm.N);
	
	if (fgets(string, 1000, fp) == NULL) kcerror ("ERROR: cannot read first line from file");
	for (i=1;i<=T;i++) {
		if (fscanf(fp, "%s\t", string) != 1) kcerror ("ERROR: cannot read string from file");
		if (fscanf(fp, "%lf\t", &(O1[i])) != 1) kcerror ("ERROR: cannot read O1 from file");
		if (fscanf(fp, "%i\t", &(O2[i])) != 1) kcerror ("ERROR: cannot read O2 from file");
		if (fscanf(fp, "%lf\t", &(pfb[i])) != 1) kcerror ("ERROR: cannot read PFB from file");
		if (fscanf(fp, "%i\n", &(length[i])) != 1) kcerror ("ERROR: cannot read LENGTH from file");
		if (i<5) 
			fprintf(stderr, "found marker=%s o1=%f o2=%i pfb=%f length=%i\n", string, O1[i], O2[i], pfb[i], length[i]);
		

		bpcount++;
		bpcovcount += O1[i];
		if (O1[i]) {
			bplncovcount += log(O1[i]);
		} else {
			bplncovcount += log(1e-9);
		}
	}
	s = log (bpcovcount/ bpcount) - bplncovcount/ bpcount;
	k = (3-s+sqrt((s-3)*(s-3)+24*s)) / (12*s);
	theta = bpcovcount/ bpcount / k;
	
	
	if (gamma_k > 0) {
		k = gamma_k;
		theta = gamma_theta;
	} else {
		fprintf(stderr, "NOTICE: Inferring inaccurate gamma distribution from data: bpcount=%f bpcovcount=%f bplncovcount=%f realcov=%f s=%f k=%f theta=%f\n", bpcount, bpcovcount, bplncovcount, bpcovcount/ bpcount, s, k, theta);
	}
	
	ViterbiLogNP_SEQ (&hmm, T, k, theta, O1, O2, pfb, length, delta, psi, q);
	
	for (i=1; i<=T; i++) {
		mlstate[i]=q[i];			/*assign the most likely state value to return*/
	}

	free_ivector(q, 1, T);
	free_imatrix(psi, 1, T, 1, hmm.N);
	free_dmatrix(delta, 1, T, 1, hmm.N);
	free_dvector(O1, 1, T);
	free_ivector(O2, 1, T);
	free_dvector(pfb, 1, T);
	free_ivector(length, 1, T);
}


void ViterbiLogNP_SEQ (CHMM *phmm, int T, double k, double theta, double *O1, int *O2, double *pfb, int *length, double **delta, int **psi, int *q)
{
        int     i, j;   /* state indices */
        int     t;      /* time index */
 
 	int snp_count = 0;
 	int cn_count = 0;

        int     maxvalind;
        double  maxval, val;
	double  **biot;
	double **A1;
	double *pprob;
	
	A1 = dmatrix (1, phmm->N, 1, phmm->N);					/*initialize A1 matrix*/
	for (i = 1; i <= phmm->N; i++) {
		for (j = 1; j <= phmm->N; j++) {
			A1[i][j] = phmm->A[i][j];
		}
	}
	
	phmm->B1_mean[1] = k; phmm->B1_sd[1] = 1/k;
	phmm->B1_mean[2] = k; phmm->B1_sd[2] = theta/2.0;
	phmm->B1_mean[3] = k; phmm->B1_sd[3] = theta;
	phmm->B1_mean[4] = k; phmm->B1_sd[4] = theta;
	phmm->B1_mean[5] = k; phmm->B1_sd[5] = theta*1.5;
	phmm->B1_mean[6] = k; phmm->B1_sd[6] = theta*2.0;
	phmm->B3_mean[1] = k; phmm->B3_sd[1] = 1/k;
	phmm->B3_mean[2] = k; phmm->B3_sd[2] = theta/2.0;
	phmm->B3_mean[3] = k; phmm->B3_sd[3] = theta;
	phmm->B3_mean[4] = k; phmm->B3_sd[4] = theta;
	phmm->B3_mean[5] = k; phmm->B3_sd[5] = theta*1.5;
	phmm->B3_mean[6] = k; phmm->B3_sd[6] = theta*2.0;
	
		
	/* 0. Preprocessing */

	for (i = 1; i <= phmm->N; i++) {
		if (phmm->pi[i] == 0) phmm->pi[i] = 1e-9;			/*eliminate problems with zero probability*/
		phmm->pi[i] = log(phmm->pi[i]);
	}

	biot = dmatrix(1, phmm->N, 1, T);
	for (t = 1; t <= T; t++) {
		for (i = 1; i <= phmm->N; i++) {
			if (pfb[t] > 1) {
				if (t<10) {
					fprintf(stderr, "NOTICE: i=%i t=%i l=%i b3iot=%f o1=%f o2=%i pfb=%f mena=%f sd=%f uf=%f\n", i, t, length[t], b3iot (i, phmm->B1_mean, phmm->B1_sd, phmm->B1_uf, O1[t], length[t]), O1[t], O2[t], pfb[t], phmm->B3_mean[i], phmm->B3_sd[i], phmm->B3_uf);
				}
				if (!phmm->NP_flag) nrerror ("FATAL ERROR: CN probe detected but HMM model does not contain parameters for them\n");
				biot[i][t] = b3iot (i, phmm->B3_mean, phmm->B3_sd, phmm->B3_uf, O1[t], length[t]);
				cn_count++;
			} else {						/*a regular SNP marker; we use both LRR and BAF information to calculate logProb for the marker*/
				if (t<10) {
					fprintf(stderr, "NOTICE: i=%i t=%i l=%i b3iot=%f b4iot=%f o1=%f o2=%i pfb=%f mena=%f sd=%f uf=%f\n", i, t, length[t], b3iot (i, phmm->B1_mean, phmm->B1_sd, phmm->B1_uf, O1[t], length[t]), b4iot (i, phmm->B3_mean, phmm->B3_sd, phmm->B3_uf, pfb[t], (int) O1[t], O2[t]), O1[t], O2[t], pfb[t], phmm->B3_mean[i], phmm->B3_sd[i], phmm->B3_uf);
				}
				biot[i][t] = b3iot (i, phmm->B1_mean, phmm->B1_sd, phmm->B1_uf, O1[t], length[t]);
				biot[i][t] += b4iot (i, phmm->B1_mean, phmm->B1_sd, phmm->B1_uf, pfb[t], (int) O1[t], O2[t]);
				snp_count++;
				
			}
		}
	}
	
        /* 1. Initialization  */
 
        for (i = 1; i <= phmm->N; i++) {
                delta[1][i] = phmm->pi[i] + biot[i][1];
                psi[1][i] = 0;
                fprintf(stderr, "i=%i t=%i delta[t][i]=%f psi=%i pi=%f biot=%f\n", i, 1, delta[1][i],psi[1][i], phmm->pi[i], biot[i][1]);
        }
  	
        /* 2. Recursion */
 
        for (t = 2; t <= T; t++) {


                for (j = 1; j <= phmm->N; j++) {
                        maxval = -VITHUGE;
                        maxvalind = 1;
                        for (i = 1; i <= phmm->N; i++) {
                                val = delta[t-1][i] + log (A1[i][j]);
                                if (val > maxval) {
                                        maxval = val;
                                        maxvalind = i;
                                }
                        }
 
                        delta[t][j] = maxval + biot[j][t]; 
                        psi[t][j] = maxvalind;
                        if (t<10) fprintf(stderr, "i=%i t=%i delta[t][i]=%f psi=%i\n", j, t,  delta[t][j],psi[t][j]);
                }
        }
         
        /* 3. Termination */
 
        *pprob = -VITHUGE;
        q[T] = 1;
        for (i = 1; i <= phmm->N; i++) {
                if (delta[T][i] > *pprob) {
                        *pprob = delta[T][i];
                        q[T] = i;
                        fprintf(stderr, "state i=%i is better than prob=%f\n", i, *pprob);

                }
        }


	/* 4. Path (state sequence) backtracking */

	for (t = T - 1; t >= 1; t--)
		q[t] = psi[t+1][q[t+1]];
	
	for (i = 1; i <= phmm->N; i++) {					/*recover the HMM model as original*/
		phmm->pi[i] = exp(phmm->pi[i]);
	}
	
	free_dmatrix(biot, 1, phmm->N, 1, T);
	free_dmatrix(A1, 1, phmm->N, 1, phmm->N);
	
}
