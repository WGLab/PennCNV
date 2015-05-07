#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "kc.h"

/************************   DESCRIPTION OF THE kc.c FILE   *******************

$Revision: 298 $
$LastChangedDate: 2009-09-01 07:17:15 -0700 (Tue, 01 Sep 2009) $

The kc.c file contains C subroutines that collectively form the core of kc.pm. 
The kc.pm is a light-weight Perl module that I developed to provide fast 
computation for commonly used mathematical and statistical functions. A Perl 
programmer can easily use the kc.pm module to run mathematically intensive tasks 
implemented in C, which may be 10-50 times faster than the counterparts 
written in native Perl.

The source code is divided into several major sections:
1. macros and definitions
2. Perl-C interface
3. matrix memory allocations
4. basic matrix operations
5. advanced matrix operations
6. statistical distributions
7. summary statistics
8. statistical test
9. mathematical functions
10. random number generator
11. third-party code with modification

The code is under constant development. The majority of the functions are 
obtained from googling with certain modifications. Only some of the functions 
have been tested rigorously, while other functions need more testing and 
debugging. For questions, bug reports or various other concerns, email 
kai@openbioinformatics.org.

*****************************************************************************/



/********************   EXAMPLES ON USING THE kc.pm MODULE *******************

Several notes on using the kc.pm is described below:

The kc.pm module is a simple perl module that is not object-oriented, that does 
not export any function or variable. As a result, when using the functions in 
kc.pm, one must always prefix the function name with "kc::".

For the majority of functions (except matrix functions), appropriate care has 
been taken such that a Perl reference to scalar variable directly corresponds to 
a pointer in the C function, and that a Perl array directly corresponds to a 
pointer to double array in the C function. For matrix operations, one must 
manually convert a Perl array or array-of-array to a C vector or a C matrix 
(vector-to-vector), or vice versa. Eight functions, including four that start 
with C2perl and four that start with perl2C, are designed for this manual 
conversion of data structures.

Some examples are given below to show how kc.pm can be used from a Perl program:

#EXAMPLE START
require "kc.pm";
my @data = qw/1.2 2.5 3.4 7.8/;
my ($mean, $variance);
kc::avevar (\@data, 4, \$mean, \$variance);
print "input data contains four samples with mean=$mean; variance=$variance\n";
#EXAMPLE END

The above Perl code calculates the mean and variance of a Perl array, which is 
transmitted to the kc.pm module through a Perl reference to array, and returns 
the results are returned to the $mean and $variance variable via a reference to 
these two scalars. The implicit conversion of data structure is performed in a 
black box not visible to users.

#EXAMPLE START
require "kc.pm";
my $output = kc::pdf_hypergeometric (10, 20, 30, 40);
#EXAMPLE END

The above code simply calculate the hypergeometric distribution (probability 
mass function, despite the PDF name!) and returns the results to $output Perl 
scalar variable, as if the kc::pdf_hypergeometrix were a native Perl function.

#EXAMPLE START
require "kc.pm";
my @sex = qw/0 1 1 1 0 1/;
my @height = qw/5.8 5.9 4.2 4.7 5.5 6.0/;
my @weight = qw/130 120 140 155 134 148/;
my $x = [[(1)x6], \@sex, \@height];
my $y = [\@weight];
my $xc = kc::perl2C_matrix ($x, 3, 6);		#generate C matrix
my $yc = kc::perl2C_matrix ($y, 1, 6);
my $xt = kc::matrix_transpose ($x, 3, 6);	#tranpose C matrix
my $yt = kc::matrix_tranpose ($y, 1, 6);
my $b = kc::matrix_reg ($xt, 3, 6, $yt);	#regression
kc::print_matrix ($b, 3, 1);			#print coefficient
my $c = kc::matrix_regCoefStat ($xt, 3, 6, $yt, $b);
kc::print_matrix ($c, 3, 5);			#print statistic
my $reg_coef = kc::C2perl_matrix ($b, 3, 1);	#return to Perl
#EXAMPLE END

The above code first build a X matrix and a Y matrix, then perform standard 
multiple linear regression, print out the regression coefficient and their test 
statistics, finally returns the coefficient as a Perl reference to 
array-of-array. In this example, the @sex, @height, @weight, $x, $y and 
$reg_coef are all Perl scalar or array variables, but the $xc, $yc, $xt, $yt, $b 
and $c are all references to C data structures.

Unlike most other functions in the kc.pm module, the matrix operations in kc.pm 
requires manual conversion of Perl and C data structures, since the automatic 
conversion does not work well for two-dimensional arrays, and since usually a 
procedure of matrix operations (like multiple regression) requires multiple 
steps, and that I think there is much less overhead if all temporary variables 
are kept as C data structure, rather than Perl data structure.

#EXAMPLE START
require "kc.pm";
my @matrix = [[1,2,3],[1,2,3],[2,4,4]];
my $mat = kc::perl2C_matrix (\@matrix, 3, 3);
eval {
	kc::matrix_ludecomp ($mat, 3, 3);
}
if ($@ and $@=~m/singular/i) {
	print STDERR "WARNING: singular matrix encounterd\n";
	print STDERR "ERROR MESSAGE from kc.pm: $@\n";
}
#EXAMPLE END

The above code demonstrate how to catch the error message when executing the C 
subroutines and bypass the normal program exit in C, and handle the error code 
to Perl for further analysis and remedy. This "error-catching" mechanism is what 
Perl is good for, and is not something that you could easily do in C.

Altogether, these above examples should give a hint on what kc.pm is and what it 
does. The code is under continuous update and development and will incorporate 
more functionalities in the future.

******************************************************************************/


/******************************************************************************
the following section contains definitions for macros (this section will be gradually eliminated in future versions)
******************************************************************************/

#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
#define PI 3.141592653579893
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
#define TINY 1.0e-20;


/*******************************************************************************
the following section contains helper functions to handle Perl-C interface

the four functions starting with "C2perl" and four functions starting with 
"perl2C" can be used to convert variables between perl and C. However, little 
error-checking functionality has been implemented in these subroutines yet, as 
it is quite non-trivial for me. Therefore, users must exercise caution when 
using these functions. For example, make sure that all Perl arrays contains 
"numbers" before sending the data to C double array.

To understand how these functions work or to modify the functionality of these 
functions, it is recommended to read the perlguts manual thoroughly. Besides the 
perlguts documentation embeded in the Perl language, this linke 
(http://www.perl.org/tpc/1998/Perl_Language_and_Modules/Perl%20Illustrated/) 
provides very good illustration on the relationships between different variables 
and functions and it is one of the best references that one can find.

*******************************************************************************/

FILE *fh_stdout ()
{
	return stdout;
}

FILE *fh_stderr ()
{
	return stderr;
}

FILE *fh_stdin ()
{
	return stdin;
}

void kcerror(char *error_text)
{	
	fprintf (stderr, "kc.pm module run-time error (use eval{} in Perl to catch the error) ...\n");
	
	/*This is the XSUB-writer's interface to Perl's die function.*/
	croak (error_text);
}

void kcwarn(char *error_text)
{	
	fprintf (stderr,"kc.pm module warning (use eval{} in Perl to catch the warning message) ...\n");
	
	/*This is the XSUB-writer's interface to Perl's warn function.*/
	warn (error_text);
}

SV* C2perl_vector (double *vector, int n)
/* convert a C double vector of size n into a perl array, return the array reference */
{
	AV* newvec;
	SV* newref;
	int i;
	
	newvec = newAV ();
	for (i=0; i<n; i++) {
		av_push (newvec, newSVnv (vector[i]));
	}
	
	newref = newRV ((SV*) newvec);
	return newref;
}

SV* C2perl_vectori (int *vector, int n)
/* convert a C int vector of size n into a perl array, return the array reference */
{
	AV* newvec;
	SV* newref;
	int i;
	
	newvec = newAV ();
	for (i=0; i<n; i++) {
		/* AN IV CAST SHOULD BE USED HERE: typedef IV is a simple signed integer type that is guaranteed to be large enough to hold a pointer (as well as an integer) */
		av_push (newvec, newSViv ( (IV) vector[i]));
	}
	
	newref = newRV ((SV*) newvec);
	return newref;
}

double* perl2C_vector (SV* arg)
/* convert a Perl array (size is implictely defined in the array object) into a double C array, return the pointer to the C array */
{
	double *newvec;
	AV* array;
	SV** work;
	I32 i, n;
	
	if (!SvROK (arg))
		kcerror ("Error: arg is not a reference varable");

	if (SvTYPE (arg) != SVt_PVGV && SvTYPE (SvRV (arg)) != SVt_PVAV)
		kcerror ("Error: arg is not a reference varable to arrays");

	if (SvTYPE(arg)==SVt_PVGV) {
		array = GvAVn((GV*) arg);		/* glob */
	} else{
		array = (AV *) SvRV(arg);		/* reference */
	}
	
	n = av_len (array);
	newvec = vector_new (n+1);
	
	for (i=0; i<=n; i++) {
		work = av_fetch (array, i, 0);
		newvec[i] = SvNV (*work);		/* automatically coerce the actual scalar type into an NV */
	}
	return newvec;
}

int* perl2C_vectori (SV* arg)
/* convert a Perl array (size is implictely defined in the array object) into an int C array, return the pointer to the C array */
{
	int *newvec;
	AV* array;
	SV** work;
	I32 i, n;
	
	if (!SvROK (arg))
		kcerror ("Error: arg is not a reference varable");

	if (SvTYPE (arg) != SVt_PVGV && SvTYPE (SvRV (arg)) != SVt_PVAV)
		kcerror ("Error: arg is not a reference varable to arrays");

	if (SvTYPE(arg)==SVt_PVGV) {
		array = GvAVn((GV*) arg);		/* glob */
	} else{
		array = (AV *) SvRV(arg);		/* reference */
	}
	
	n = av_len (array);
	newvec = vectori_new (n+1);
	
	for (i=0; i<=n; i++) {
		work = av_fetch (array, i, 0);
		newvec[i] = SvIV (*work);		/* automatically coerce the actual scalar type into an IV */
	}
	return newvec;
}

SV* C2perl_matrix (double** mat, int n1, int n2)
/* convert a double C matrix (pointer-to-pointer) into a perl array-of-array */
{
	AV* newmat;
	AV* newvec;
	SV* newref;
	int i, j;
	
	newmat = newAV ();
	for (i=0; i<n1; i++) {
		newvec = newAV ();
		for (j=0; j<n2; j++) {
			av_push (newvec, newSVnv (mat[i][j]));
		}
		av_push (newmat, newRV ((SV*) newvec));
	}
	newref = newRV ((SV*) newmat);
	return newref;
}

SV* C2perl_matrixi (int** mat, int n1, int n2)
/* convert a int C matrix (pointer-to-point) into a perl array-of-array */
{
	AV* newmat;
	AV* newvec;
	SV* newref;
	int i, j;
	
	newmat = newAV ();
	for (i=0; i<n1; i++) {
		newvec = newAV ();
		for (j=0; j<n2; j++) {
			av_push (newvec, newSViv ( (IV) mat[i][j]));
		}
		av_push (newmat, newRV ((SV*) newvec));
	}
	newref = newRV ((SV*) newmat);
	return newref;
}

double** perl2C_matrix (SV* arg)
/*convert Perl reference (array-of-array) to double C two-dimensional arrays (pointer-to-pointer).
The new matrix can be subsequently used in Perl (by calling additional C subroutines)
once the matrix is no longer useful, one should use matrix_free to empty it
*/
{	
	double **newmat=0;
	SV** work2;
	AV* array;
	AV* array2;
	I32 i, j, n, m=-1;
	
	if (!SvROK (arg))
		kcerror ("kc::perl2C_matrix: arg is not a reference varable");

	if (SvTYPE (arg) != SVt_PVGV && SvTYPE (SvRV (arg)) != SVt_PVAV)
		kcerror ("kc::perl2C_matrix: arg is not a reference varable to arrays");

	if (SvTYPE(arg)==SVt_PVGV) {
		array = GvAVn((GV*) arg);		/* glob */
	} else{
		array = (AV *) SvRV(arg);		/* reference */
	}
	
	n = av_len (array);				/* n=the highest index value in array (length-1) */
	
	for (i=0; i<=n; i++) {
		work2 = av_fetch (array, i, 0);
		if (work2 == NULL || !SvROK (*work2))
			kcerror ("Error in perl2C_matrix: second dimension of arg is not a reference to arrays");
		array2 = (AV *) SvRV(*work2);		/* array of 2nd dimension */
		
		if (m >= 0 && m != av_len (array2)) {
			kcerror ("Error: array-to-array need to have the same elements in each sub-array");
		} else if (m < 0) {
			m = av_len (array2);
			newmat = matrix_new (n+1, m+1);	/*construct the return matrix*/
		}
		
		for (j=0; j<=m; j++) {
			work2 = av_fetch (array2, j, 0);
			if (SvROK (*work2))
				kcerror ("Error: the array-to-array is more than two dimensions");
			newmat[i][j] = SvNV (*work2);	/*work2 might be a string, but SvNV automatically coerce it into numerical value*/
		}
	}
	
	return newmat;
}

int** perl2C_matrixi (SV* arg)
/*convert Perl reference (array-of-array) to int C two-dimensional arrays (pointer-to-pointer).
The new matrix can be subsequently used in Perl (by calling additional C subroutines)
once the matrix is no longer useful, one should use matrixi_free to empty it
*/
{
	int **newmat=0;
	SV** work2;
	AV* array;
	AV* array2;
	I32 i, j, n, m=-1;
	
	if (!SvROK (arg))
		kcerror ("kc::perl2C_matrixi: arg is not a reference varable");

	if (SvTYPE (arg) != SVt_PVGV && SvTYPE (SvRV (arg)) != SVt_PVAV)
		kcerror ("kc::perl2C_matrixi: arg is not a reference varable to arrays");

	if (SvTYPE(arg)==SVt_PVGV) {
		array = GvAVn((GV*) arg);		/* glob */
	} else{
		array = (AV *) SvRV(arg);		/* reference */
	}
	
	n = av_len (array);				/* the highest index value in array (length-1) */
	
	for (i=0; i<=n; i++) {
		work2 = av_fetch (array, i, 0);
		if (work2 == NULL || !SvROK (*work2))
			kcerror ("kc::perl2C_matrixi: second dimension of arg is not a reference to arrays");
		array2 = (AV *) SvRV(*work2);		/* array of 2nd dimension */
		
		if (m >= 0 && m != av_len (array2)) {
			/*fprintf (stderr, "number of elements in previous (%i) row is %i but in current (%i) row is %i\n", (int) i-1, (int) m, (int) i, (int) av_len (array2));*/
			kcerror ("kc::perl2C_matrixi: Perl array-of-array need to have the same elements in each sub-array");
		} else if (m < 0) {
			m = av_len (array2);
			newmat = matrixi_new (n+1, m+1);	/*construct the return matrix*/
		}
		
		for (j=0; j<=m; j++) {
			work2 = av_fetch (array2, j, 0);
			if (SvROK (*work2))
				kcerror ("kc::perl2C_matrixi: the array-of-array is more than two dimensions");
			newmat[i][j] = SvIV (*work2);	/*work2 might be a string, but SvNV automatically coerce it into integer value*/
		}
	}
	
	return newmat;
}


void print_vectori (int *array, int n)
{
	int i;
	printf ("[ %i", array[0]);
	for (i=1; i<n; i++) {
		printf ("\t%i", array[i]);
	}
	printf (" ]\n");
}

void print_vector (double *array, int n)
{
	int i;
	printf ("[ %f", array[0]);
	for (i=1; i<n; i++) {
		printf ("\t%f", array[i]);
	}
	printf (" ]\n");
}

void print_matrixi (int** A, int m, int n)
{
	int i, j;
	for (i=0; i<m; i++) {
		printf ("[ %i", A[i][0]);
		for (j=1; j<n; j++) {
			printf ("\t%i", A[i][j]);
		}
		printf (" ]\n");
	}
}

void print_matrix (double** A, int m, int n)
{
	int i, j;
	for (i=0; i<m; i++) {
		printf ("[ %f", A[i][0]);
		for (j=1; j<n; j++) {
			printf ("\t%f", A[i][j]);
		}
		printf (" ]\n");
	}
}


/* the following five subroutines were copied and modified from the arrays.c 
file in the Perl Data Language. in older versions (pre-2009August) of kc.c, 
these subroutines were in separate files in the include/ directory, however, 
they have created multiple issues with respect to compilation in different 
operating systems. Therefore, I decided that I include them directly into the 
kc.c program with appropriate modification to the source code

The beauty of the packing/unpacking function is that they have used the mortal 
space provided by Perl. Therefore, memory allocation and free-up is performed 
automatically by Perl, and users do not have to worry about these tedious work 

*/


int is_scalar_ref (SV* arg) {
/* Utility to determine if ref to scalar */
	SV* foo;
	if (!SvROK(arg)) return 0;
	foo = SvRV(arg);
	if (SvPOK(foo)) {
		return 1;
	} else {
		return 0;
	}
}

void* pack1D ( SV* arg, char packtype ) {
/* ####################################################################################

   pack1D - argument is perl scalar variable and one char pack type. 
   If it is a reference to a 1D array pack it and return pointer.
   If it is a glob pack the 1D array of the same name.
   If it is a scalar pack as 1 element array.  
   If it is a reference to a scalar then assume scalar is prepacked binary data

   [1D-ness is checked - routine croaks if any of the array elements
   themselves are references.] 


*/
	int iscalar;
	float scalar;
	double dscalar;
	short sscalar;
	unsigned char uscalar;
	AV* array;
	I32 i,n;
	SV* work;
	SV** work2;
	double nval;
	STRLEN len;

	if (is_scalar_ref(arg))                 /* Scalar ref */
		return (void*) SvPV(SvRV(arg), len);

	if (packtype!='f' && packtype!='i' && packtype!='d' && packtype!='s' && packtype != 'u')
		croak("ERROR: valid type conversion for pack1D are f, i, d, s and u only");

	work = sv_2mortal(newSVpv("", 0));
   
	/* Is arg a scalar? Return scalar*/
	if (!SvROK(arg) && SvTYPE(arg)!=SVt_PVGV) {
		if (packtype=='f') {
			scalar = (float) SvNV(arg);             /* Get the scalar value */
			sv_setpvn(work, (char *) &scalar, sizeof(float)); /* Pack it in */
		}
		if (packtype=='i') {
			iscalar = (int) SvNV(arg);             /* Get the scalar value */
			sv_setpvn(work, (char *) &iscalar, sizeof(int)); /* Pack it in */
		}
		if (packtype=='d') {
			dscalar = (double) SvNV(arg);		/*Get the scalar value */
			sv_setpvn(work, (char *) &dscalar, sizeof(double)); /* Pack it in */
		}
		if (packtype=='s') {
			sscalar = (short) SvNV(arg);		/*Get the scalar value */
			sv_setpvn(work, (char *) &sscalar, sizeof(short)); /* Pack it in */
		}
		if (packtype=='u') {
			uscalar = (unsigned char) SvNV(arg);	/*Get the scalar value */
			sv_setpvn(work, (char *) &uscalar, sizeof(char)); /* Pack it in */
		}
		return (void *) SvPV(work, PL_na);        /* Return the pointer */
	}
   
	/* Is it a glob or reference to an array? */
	if (SvTYPE(arg)==SVt_PVGV || (SvROK(arg) && SvTYPE(SvRV(arg))==SVt_PVAV)) {
		if (SvTYPE(arg)==SVt_PVGV) {
			array = (AV *) GvAVn((GV*) arg);   /* glob */
		} else{
			array = (AV *) SvRV(arg);   /* reference */
		}

		n = av_len(array);
   
		if (packtype=='f')
			SvGROW( work, sizeof(float)*(n+1) );  /* Pregrow for efficiency */
		if (packtype=='i')
			SvGROW( work, sizeof(int)*(n+1) );   
		if (packtype=='d')
			SvGROW( work, sizeof(double)*(n+1) );
		if (packtype=='s')
			SvGROW( work, sizeof(short)*(n+1) );   
		if (packtype=='u')
			SvGROW( work, sizeof(char)*(n+1) );
      

		/* Pack array into string */
		for(i=0; i<=n; i++) {
			work2 = av_fetch( array, i, 0 ); /* Fetch */
			if (work2==NULL) 
				nval = 0.0;   /* Undefined */
			else {
				if (SvROK(*work2)) 
					croak("Routine can only handle scalar values or refs to 1D arrays of scalars");	/*  Checking 1D-ness: Croak if reference [i.e. not 1D] */
				nval = SvNV(*work2);               
			}   
			
			if (packtype=='f') {
				scalar = (float) nval;
				sv_catpvn( work, (char *) &scalar, sizeof(float));
			}
			if (packtype=='i') {
				iscalar = (int) nval;
				sv_catpvn( work, (char *) &iscalar, sizeof(int));
			}
			if (packtype=='d') {
				dscalar = (double) nval;
				sv_catpvn( work, (char *) &dscalar, sizeof(double));
			}
			if (packtype=='s') {
				sscalar = (short) nval;
				sv_catpvn( work, (char *) &sscalar, sizeof(short));
			}
			if (packtype=='u') {
				uscalar = (unsigned char) nval;
				sv_catpvn( work, (char *) &uscalar, sizeof(char));
			}
		}

		/* Return a pointer to the byte array */
		return (void *) SvPV(work, PL_na);
	}

	croak("Routine can only handle scalar values or refs to 1D arrays of scalars");
}

void unpack1D ( SV* arg, void * var, char packtype, int n ) {
/* ##################################################################################

   unpack1D - take packed string (C array) and write back into perl 1D array.
   If 1st argument is a reference, unpack into this array.
   If 1st argument is a glob, unpack into the 1D array of the same name.

   Can only be used in a typemap if the size of the array is known
   in advance or is the size of a preexisting perl array (n=0). If it
   is determined by another variable you may have to put in in some
   direct CODE: lines in the XSUB file.

*/
/* n is the size of array var[] (n=1 for 1 element, etc.) If n=0 take
      var[] as having the same dimension as array referenced by arg */
   
	int* ivar=0;
	float* fvar=0;
	double* dvar=0;
	short* svar=0;
	unsigned char* uvar=0;
	AV* array;
	I32 i,m;
	
	/* Note in ref to scalar case data is already changed */
	
	if (is_scalar_ref(arg)) /* Do nothing */
		return;
	
	if (packtype!='f' && packtype!='i' && packtype!= 'd' && packtype!='u' && packtype!='s')
		croak("ERROR: valid type conversion specified to unpack1D are f, i, d, u and s only");
   
	m=n;  
	array = coerce1D( arg, m );   /* Get array ref and coerce */
	
	if (m==0) 
		m = av_len( array )+1;  
	
	if (packtype=='i')        /* Cast void array var[] to appropriate type */
		ivar = (int *) var;
	if (packtype=='f') 
		fvar = (float *) var;
	if (packtype=='d') 
		dvar = (double *) var;
	if (packtype=='u') 
		uvar = (unsigned char *) var;
	if (packtype=='s') 
		svar = (short *) var;
 
	/* Unpack into the array */
	for(i=0; i<m; i++) {
		if (packtype=='i') 
			av_store( array, i, newSViv( (IV)ivar[i] ) );
		if (packtype=='f') 
			av_store( array, i, newSVnv( (double)fvar[i] ) );
		if (packtype=='d') 
			av_store( array, i, newSVnv( (double)dvar[i] ) );
		if (packtype=='u') 
			av_store( array, i, newSViv( (IV)uvar[i] ) );
		if (packtype=='s') 
			av_store( array, i, newSViv( (IV)svar[i] ) );
	}
	return;
}




AV* coerce1D ( SV* arg, int n ) {
/* #################################################################################

   coerce1D - utility function. Make sure arg is a reference to a 1D array 
   of size at least n, creating/extending as necessary. Fill with zeroes.
   Return reference to array. If n=0 just returns reference to array,
   creating as necessary.
*/
/* n is the size of array var[] (n=1 for 1 element, etc.) */
   
	AV* array;
	I32 i,m;
   
   /* In ref to scalar case we can do nothing - we can only hope the
      caller made the scalar the right size in the first place  */

	if (is_scalar_ref(arg)) /* Do nothing */
		return (AV*)NULL;
   
   /* Check what has been passed and create array reference whether it
      exists or not */

	if (SvTYPE(arg)==SVt_PVGV) {
		array = GvAVn((GV*)arg);                             /* glob */
	} else if (SvROK(arg) && SvTYPE(SvRV(arg))==SVt_PVAV) {
		array = (AV *) SvRV(arg);                           /* reference */
	} else{
		array = newAV();                                    /* Create */
		sv_setsv(arg, newRV((SV*) array));                            
	}
   
	m = av_len(array);
	
	for (i=m+1; i<n; i++) {
		av_store( array, i, newSViv( (IV) 0 ) );
	}
	return array;
}

void* get_mortalspace( int n, char packtype ) {
/* ################################################################################

   get_mortalspace - utility to get temporary memory space. Uses
   a mortal *SV for this so it is automatically freed when the current
   context is terminated. Useful in typemap's for OUTPUT only arrays.

*/
/* n is the number of elements of space required, packtype is 'f' or 'i' */
   
	SV* work;
	
	if (packtype!='f' && packtype!='i' && packtype!='d' && packtype!='u' && packtype!='s')
		croak("ERROR: valid type conversion specified to get_mortalspace are f, i, d, u and s only");
	
	work = sv_2mortal(newSVpv("", 0));
	
	if (packtype=='f')
		SvGROW( work, sizeof(float)*n );  /* Pregrow for efficiency */
	if (packtype=='i')
		SvGROW( work, sizeof(int)*n );  
	if (packtype=='d')
		SvGROW( work, sizeof(double)*n);
	if (packtype=='u')
		SvGROW( work, sizeof(char)*n);
	if (packtype=='s')
		SvGROW( work, sizeof(short)*n);
	
	return (void *) SvPV(work, PL_na);
}


/***********************************************************************************

the following section contains helper functions for memory allocation of arrays 
and matrices There are two major subsections: (1) subroutines written by myself 
(2) subroutines from Numerical Recipe util.c for historical reasons

The NR_END in NR functions is defined as 1, which means that we always waste one 
spot when allocating memory

Note that matrix dimension is represented by int, rather than long, since in 
most applications in 32-bit or 64-bit computers, int is more than sufficient to 
hold array indexes.

***********************************************************************************/

void MATRIX_MEMORY ()
{
	printf ("This section contains matrix/vector memory allocation functions\n");
}

double **matrix_new (int row, int col)
/* allocate a double matrix with subscript range m[0..row-1][0..col-1] 
Note that all the data in the matrix are arranged consecutively in memory.
The idea is borrowed from Numerical Recipe, where one vector and one matrix is both allocated in memory.
*/
{
	int i;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((row)*sizeof(double*)));
	if (!m) kcerror("kc::matrix_new: allocation failure 1");

	/* allocate rows and set pointers to them */
	m[0]=(double *) malloc((size_t)((row*col)*sizeof(double)));
	if (!m[0]) kcerror("kc::matrix_new: allocation failure 2");

	for(i=1;i<row;i++) m[i]=m[i-1]+col;
	return m;
}

int **matrixi_new (int row, int col)
/* allocate a int matrix with subscript range m[0..row-1][0..col-1] 
*/
{
	int i;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((row)*sizeof(int*)));
	if (!m) kcerror("kc::matrix_new: allocation failure 1");

	/* allocate rows and set pointers to them */
	m[0]=(int *) malloc((size_t)((row*col)*sizeof(int)));
	if (!m[0]) kcerror("kc::matrix_new: allocation failure 2");

	for(i=1;i<row;i++) m[i]=m[i-1]+col;
	return m;
}

char **matrixc_new (int row, int col)
/* allocate a int matrix with subscript range m[0..row-1][0..col-1] 
*/
{
	int i;
	char **m;

	/* allocate pointers to rows */
	m=(char **) malloc((size_t)((row)*sizeof(char*)));
	if (!m) kcerror("kc::matrixc_new: allocation failure 1");

	/* allocate rows and set pointers to them */
	m[0]=(char *) malloc((size_t)((row*col)*sizeof(char)));
	if (!m[0]) kcerror("kc::matrixc_new: allocation failure 2");

	for(i=1;i<row;i++) m[i]=m[i-1]+col;
	return m;
}

int** matrix_toMatrixi (double** A, int m, int n)
/* convert a double matrix to int matrix
 */
{
	int **B;
	int i, j;
	
	B = matrixi_new (m, n);
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++) {
			B[i][j] = (double) A[i][j];
		}
	}
	return B;
}

double** matrixi_toMatrix (int** A, int m, int n)
/* convert an int matrix to a double matrix
 */
{
	double** B;
	int i, j;
	
	B = matrix_new (m, n);
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++) {
			B[i][j] = (int) floor (A[i][j]+0.5);
		}
	}
	return B;
}

double *vector_new (int n)
/* allocate a double vector with subscript range v[0..n-1] */
{
	double *v;

	v=(double *)malloc((size_t) ((n)*sizeof(double)));
	if (!v) kcerror("kc::vector_new: allocation failure");
	return v;
}

int *vectori_new (int n)
/* allocate a int vector with subscript range v[0..n-1] */
{
	int *v;

	v=(int *)malloc((size_t) ((n)*sizeof(int)));
	if (!v) kcerror("kc::vector_new: allocation failure");
	return v;
}

void matrix_free (double **m)
/* free a double matrix allocated by matrix_new () */
{
	free((void *) m[0]);
	free((void *) m);
}

void matrixi_free (int **m)
/* free a double matrix allocated by matrix_new () */
{
	free((void *) m[0]);
	free((void *) (m));
}

void vector_free (double *v)
/* free a double vector allocated with vector_new() */
{
	free((void *) v);
}

void vectori_free (int *v)
/* free a double vector allocated with vectori_new() */
{
	free((void *) v);
}

double matrix_at (double **A, int m, int n)
/* return the value at row m and column n of a matrix */
{
	return A[m][n];
}

int matrixi_at (int **A, int m, int n)
/* return the value at row m and column n of a matrix */
{
	return A[m][n];
}

double vector_at (double *A, int n)
{
	return A[n];
}

int vectori_at (int *A, int n)
{
	return A[n];
}

double* matrix_row (double** A, int m)
{
	return A[m];
}

int* matrixi_row (int** A, int m)
{
	return A[m];
}

void matrix_change (double **A, int m, int n, double a)
/* change the element of matrix A at row m and column n to a */
{
	A[m][n] = a;
}

void matrixi_change (int **A, int m, int n, int a)
/* change the element of matrix A at row m and column n to a */
{
	A[m][n] = a;
}

void vector_change (double *A, int n, double a)
{
	A[n] = a;
}

void vectori_change (int *A, int n, int a)
{
	A[n] = a;
}

int* vector_toVectori (double* A, int n)
/* convert a double vector to int vector */
{
	int* B;
	int i;
	
	B = vectori_new (n);
	for (i=0; i<n; i++) {
		B[i] = (int) floor (A[i]+0.5);
	}
	return B;
}

double* vectori_toVector (int* A, int n)
/* convert an int vector to double vector */
{
	double *B;
	int i;
	
	B = vector_new (n);
	for (i=0; i<n; i++) {
		B[i] = (double) A[i];
	}
	return B;
}

#define NR_END 1
#define FREE_ARG char*

void nrerror(char *error_text)
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	void kcerror(char *error_text);
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) kcerror ("kc::ivector: allocation failure in ivector()");
	return v-nl+NR_END;
}

unsigned char *cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	void kcerror(char *error_text);
	unsigned char *v;

	v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) kcerror("kc::cvector: allocation failure in cvector()");
	return v-nl+NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	void kcerror(char *error_text);
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) kcerror("kc::dvector: allocation failure in dvector()");
	return v-nl+NR_END;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	void kcerror(char *error_text);
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) kcerror("kc::dmatrix: allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) kcerror("kc::dmatrix: allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	void kcerror(char *error_text);
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) kcerror("kc::imatrix: allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) kcerror("kc::imatrix: allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **submatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch, long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] 
this is modified from the NR version due to change of float to double
*/
{
	void kcerror(char *error_text);
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	double **m;

	/* allocate array of pointers to rows */
	m=(double **) malloc((size_t) ((nrow+NR_END)*sizeof(double*)));
	if (!m) kcerror("kc::submatrix: allocation failure in subdmatrix()");
	m += NR_END;
	m -= newrl;

	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **convert_matrix(double *a, long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	void kcerror(char *error_text);
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t) ((nrow+NR_END)*sizeof(double*)));
	if (!m) kcerror("kc::convert_matrix: allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;

	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(unsigned char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_subdmatrix(double **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

#undef NR_END
#undef FREE_ARG



/***********************************************************************************

This section contains basic matrix manipulation functions, such as copy, 
transpose and generation of various new matrices.

***********************************************************************************/

void MATRIX_BASIC ()
{
	printf ("This section contains basic matrix operations\n");
}

double** matrix_copy (double **A, int m, int n)
/* copy A and generate a new matrix B */
{
	double** B;
	int i, j;
	B = matrix_new (m, n);
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++) {
			B[i][j] = A[i][j];
		}
	}
	return B;
}

int** matrixi_copy (int **A, int m, int n)
/* copy A and generate a new matrix B */
{
	int** B;
	int i, j;
	B = matrixi_new (m, n);
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++) {
			B[i][j] = A[i][j];
		}
	}
	return B;
}

double** matrix_transpose (double **A, int m, int n)
/* transpose A and generate a new matrix B */
{
	double** B;
	int i, j;
	B = matrix_new (n, m);
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++) {
			B[j][i] = A[i][j];
		}
	}
	return B;
}

int** matrixi_transpose (int **A, int m, int n)
/* transpose A and generate a new matrix B */
{
	int** B;
	int i, j;
	B = matrixi_new (n, m);
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++) {
			B[j][i] = A[i][j];
		}
	}
	return B;
}

double** matrix_newSeq (int m, int n)
/* create a new double matrix filled with sequence */
{
	double** A = matrix_new (m, n);
	long i, j, count = 1;
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++) {
			A[i][j] = count++;
		}
	}
	return A;
}

int** matrixi_newSeq (int m, int n)
/* create a new int matrix filled with sequence */
{
	int** A = matrixi_new (m, n);
	long i, j, count = 1;
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++) {
			A[i][j] = count++;
		}
	}
	return A;
}

double** matrix_newRandom (int m, int n, int seed)
/* create a mxn matrix filled with random numbers from 0 to 1 */
{
	double** A;
	int i, j;
	
	A = matrix_new (m, n);
	srand (seed);
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++) {
			A[i][j] = ran ();
		}
	}
	return A;
}

double** matrix_newDiagonal (double* cell, int length, int m, int n)
/* create a new diagonal matrix from a given array of values, not necessarily a square matrix */
{
	double** A;
	int i, j, k=0;
	
	A = matrix_new (m, n);
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++) {
			if (i==j) {
				if (k<length) {
					A[i][j] = cell[k];
					k++;
				} else {
					A[i][j] = 0;
				}
			} else {
				A[i][j] = 0;
			}
		}
	}
	return A;
}

double** matrix_newIdentity (int m)
/* create a new identity matrix (mxm), where diagonal entries equal to 1 and all other entries equal to zero */
{
	double** A;
	int i, j;
	
	A = matrix_new (m, m);
	for (i=0; i<m; i++) {
		for (j=0; j<m; j++) {
			if (i==j) {
				A[i][j] = 1;
			} else {
				A[i][j] = 0;
			}
		}
	}
	return A;
}

double** matrix_newSubrow (double** A, int m, int n, int r1, int r2)
/* create a new matrix composing sub-rows of previous matrix A (mxn matrix)
 * r1 and r2 are 1-based counts of subrows
 */
{
	int nrow, i, j;
	double** B;
	
	if (r2<r1) kcerror ("kc::matrix_newsubrow: start row (r1) should be less than or equal to end row (r2)");
	if (r1>m) kcerror ("kc::matrix_newsubrow: start row (r1) should be less than or equal to total row (m)");
	if (r2>m) kcerror ("kc::matrix_newsubrow: end row (r2) should be less than or equal to total row (m)");
	nrow = r2-r1+1;
	B = matrix_new (nrow, n);
	
	for (i=0; i<nrow; i++) {
		for (j=0; j<n; j++) {
			B[i][j] = A[r1+i-1][j];
		}
	}
	return B;
}

double** matrix_newSubcol (double** A, int m, int n, int c1, int c2)
/* create a new matrix composing sub-columns of previous matrix A (mxn matrix) 
 * c1 and c2 are 1-based counts of subcolumns
 */
{
	int ncol, i, j;
	double** B;
	
	if (c2<c1) kcerror ("kc::matrix_newsubcol: start col (c1) should be less than or equal to end col (c2)");
	if (c1>n) kcerror ("kc::matrix_newsubcol: start col (c1) should be less than or equal to total col (n)");
	if (c2>n) kcerror ("kc::matrix_newsubcol: end col (c2) should be less than or equal to total col (n)");
	ncol = c2-c1+1;
	B = matrix_new (m, ncol);
	
	for (i=0; i<m; i++) {
		for (j=0; j<ncol; j++) {
			B[i][j] = A[i][c1+j-1];
		}
	}
	return B;
}

double** matrix_newSubmat (double** A, int m, int n, int r1, int r2, int c1, int c2)
/* create a new matrix composing sub-matrix of previous matrix A (mxn matrix) 
 * r1,r2,c1,c2 are 1-based counts of subrows and subcolumns
 */
{
	int nrow, ncol, i, j;
	double** B;
	
	if (r2<r1) kcerror ("kc::matrix_newsubmat: start row (r1) should be less than or equal to end row (r2)");
	if (r1>m) kcerror ("kc::matrix_newsubmat: start row (r1) should be less than or equal to total row (m)");
	if (r2>m) kcerror ("kc::matrix_newsubmat: end row (r2) should be less than or equal to total row (m)");
	if (c2<c1) kcerror ("kc::matrix_newsubmat: start col (c1) should be less than or equal to end col (c2)");
	if (c1>n) kcerror ("kc::matrix_newsubmat: start col (c1) should be less than or equal to total col (n)");
	if (c2>n) kcerror ("kc::matrix_newsubmat: end col (c2) should be less than or equal to total col (n)");
	nrow = r2-r1+1;
	ncol = c2-c1+1;
	B = matrix_new (nrow, ncol);
	
	for (i=0; i<nrow; i++) {
		for (j=0; j<ncol; j++) {
			B[i][j] = A[r1+i-1][c1+j-1];
		}
	}
	return B;
}

double** matrix_delRow (double** A, int m, int n, int row)
/* delete a row from the matrix and return the new matrix */
{
	int i, j;
	double** B;;
	
	if (row>m) kcerror ("kc::matrix_delRow: the deleted row should be less than the total number of rows of matrix");
	B = matrix_new (m-1, n);
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++) {
			if (i+1>row) {
				B[i-1][j] = A[i][j];
			} else if (i+1<row) {
				B[i][j] = A[i][j];
			}
		}
	}
	return B;
}

double** matrix_delCol (double** A, int m, int n, int col)
/* delete a column from the matrix and return the new matrix */
{
	int i, j;
	double** B;;
	
	if (col>n) kcerror ("kc::matrix_delCol: the deleted column should be less than the total number of columns of matrix");
	B = matrix_new (m, n-1);
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++) {
			if (j+1>col) {
				B[i][j-1] = A[i][j];
			} else if (j+1<col) {
				B[i][j] = A[i][j];
			}
		}
	}
	return B;
}	

double** matrix_cat (double** A, double** B, int m1, int m2, int n)
/* concatenate two matrix together, into a new matrix with more rows */
{
	int i, j;
	double** C;
	
	C = matrix_new (m1+m2, n);
	for (i=0; i<m1; i++) {
		for (j=0; j<n; j++) {
			C[i][j] = A[i][j];
		}
	}
	for (i=m1; i<m1+m2; i++) {
		for (j=0; j<n; j++) {
			C[i][j] = B[i-m1][j];
		}
	}
	return C;
}

double** matrix_paste (double** A, double** B, int m, int n1, int n2)
/* paste two matrix together into a new matrix with more columns */
{
	int i, j;
	double** C;
	
	C = matrix_new (m, n1+n2);
	for (i=0; i<m; i++) {
		for (j=0; j<n1; j++) {
			C[i][j] = A[i][j];
		}
	}
	for (i=0; i<m; i++) {
		for (j=n1; j<n1+n2; j++) {
			C[i][j] = B[i][j-n1];
		}
	}
	return C;
}


void matrix_zero (double**A, int m, int n)
/* make all cells in a matrix as zero*/
{
	int i, j;
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++) {
			A[i][j] = 0;
		}
	}
}


void matrix_multiply (double** A, int m, int n, double scale)
/* multiply a matrix by a constant in place */
{
	int i,j;
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++) {
			A[i][j] *= scale;
		}
	}
}

double** matrix_product (double** A, double** B, int m, int n, int o)
/* caculate the product of two matrices, one with mxn dimension, one with nxk dimension*/
{
	int i, j, k;
	double** C;
	
	C = matrix_new (m, o);
	matrix_zero (C, m, o);
	for (i=0; i<m; i++) {
		for (j=0; j<o; j++) {
			for (k=0; k<n; k++) {
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
	return C;
}

double** matrix_product3 (double** A, double** B, double** C, int m, int n, int o, int p)
/* caculate the product of three matrices, one with mxn dimension, one with nxn dimension, one with nxk dimension*/
{
	double **D, **E;
	
	D = matrix_product (A, B, m, n, o);
	E = matrix_product (D, C, m, o, p);
	matrix_free (D);
	return E;
}

void matrix_add (double** A, int m, int n, double a)
/* add all element of mxn matrix by a */
{
	int i, j;
	
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++) {
			A[i][j] += a;
		}
	}
}

void matrix_sum (double** A, double** B, int m, int n)
/* calculate the sum of two matrices, both with mxn dimensions, and put the results in the first matrix */
{
	int i, j;
	
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++) {
			A[i][j] += B[i][j];
		}
	}
}

void matrix_subtract (double** A, double** B, int m, int n)
/* calculate the difference of two matrices, both with mxn dimensions, and put the results in the first matrix */
{
	int i, j;
	
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++) {
			A[i][j] -= B[i][j];
		}
	}
}

void matrix_stdRow (double** A, int m, int n)
/* standardize each row of the matrix by (x-mean)/SD formula 
 * the matrix is changed in place
 */
{
	int i, j;
	double mean, sd;
	
	for (i=0; i<m; i++) {
		avevar(A[i], n, &mean, &sd);
		sd = sqrt(sd);
		for (j=0; j<n; j++) {
			if (sd) {
				A[i][j] = (A[i][j]-mean)/sd;
			} else {
				A[i][j] = 0;			/* when all element in the row are the same */
			}
		}
	}
}

void matrix_stdCol (double** A, int m, int n)
/* standardize each column of the matrix by (x-mean)/SD formula 
 * the matrix is changed in place
 */
{
	void matrix_stdRow (double** A, int m, int n);
	int i, j;
	double** B;
	
	B = matrix_transpose (A, m, n);
	matrix_stdRow (B, n, m);
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++) {
			A[i][j] = B[j][i];
		}
	}
	matrix_free (B);
}

double* matrix_rowMean (double** A, int m, int n)
/* calculate the mean of each row and return as a vector */
{
	double* S;
	int i;
	
	S = vector_new (m);
	for (i=0; i<m; i++) {
		S[i] = mean (matrix_row (A, i), n);
	}
	return S;
}

double* matrix_colMean (double** A, int m, int n)
/* calculate the mean of each column and return as a vector */
{
	double *S, *T;
	int i, j;
	
	S = vector_new (n);
	T = vector_new (m);
	for (j=0; j<n; j++) {
		for (i=0; i<m; i++) {
			T[i] = A[i][j];
		}
		S[j] = mean (T, m);
	}
	
	vector_free (T);
	return S;
}

double* matrix_rowVar (double** A, int m, int n)
/* calculate the mean of each row and return as a vector */
{
	double* S;
	double s;
	int i;
	
	S = vector_new (m);
	for (i=0; i<m; i++) {
		avevar (matrix_row(A, i), n, &s, &S[i]);
	}
	return S;
}

double* matrix_colVar (double** A, int m, int n)
/* calculate the mean of each column and return as a vector */
{
	double *S, *T;
	double s;
	int i, j;
	
	S = vector_new (n);
	T = vector_new (m);
	for (j=0; j<n; j++) {
		for (i=0; i<m; i++) {
			T[i] = A[i][j];
		}
		avevar (T, m, &s, &S[j]);
	}
	vector_free (T);
	return S;
}

void matrix_sortByVector(double**A, int m, int n, double* v)
/* sort the column of a matrix, based on a vector of values 
 * the values are sorted from largest one to smallest one */
{
	int k,j,i;
	double p;

	for (i=0;i<n-1;i++) {
		p=v[k=i];
		for (j=i+1;j<n;j++)
			if (v[j] >= p) p=v[k=j];
		if (k != i) {
			v[k]=v[i];
			v[i]=p;
			for (j=0;j<m;j++) {
				p=A[j][i];
				A[j][i]=A[j][k];
				A[j][k]=p;
			}
		}
	}
}

void vector_zero (double* A, int n)
/* make all cells in a matrix as zero*/
{
	int i;
	for (i=0; i<n; i++) {
		A[i] = 0;
	}
}

double* vector_copy (double *A, int n)
/* create a copy of a given vector */
{
	double *B;
	int i;
	
	B = vector_new (n);
	for (i=0; i<n; i++) {
		B[i] = A[i];
	}
	return B;
}

int* vectori_copy (int *A, int n)
/* create a copy of a given vector */
{
	int *B;
	int i;
	
	B = vectori_new (n);
	for (i=0; i<n; i++) {
		B[i] = A[i];
	}
	return B;
}

double vector_innerProduct (double *A, double *B, int n)
/* calculate the inner product of two vectors of dimension n */
{
	int i;
	double p=0;
	
	for (i=0; i<n; i++) {
		p += A[i]*B[i];
	}
	return p;
}

double** vector_outerProduct (double *A, double *B, int n)
/* calculate the outer product of two vectors of dimension n */
{
	int i, j;
	double** C;
	
	C = matrix_new (n, n);
	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			C[i][j] = A[i]*B[j];
		}
	}
	return C;
}

void vector_normalize (double* A, int n)
/* normalize a vector to unit Euclidean length */
{
	double len = 0.0;
	int i;
	
	for (i=0; i<n; i++) {
		len += A[i]*A[i];
	}
	len = sqrt(len);
	if (len) {
		for (i=0; i<n; i++) {
			A[i] /= len;
		}
	} else {
		kcerror ("kc::vector_normalize: length of input vector is zero");
	}
}

double** vector_toMatrix (double* A, int n)
/* convert a vector with n elements to 1xn matrix */
{
	double** B;
	int i;
	
	B = matrix_new (1, n);
	for (i=0; i<n; i++) {
		B[0][i] = A[i];
	}
	return B;
}

double* matrix_toVector (double** A, int n)
/* convert a 1xn matrix to a vector with n element */
{
	double* B;
	int i;
	
	B = vector_new (n);
	for (i=0; i<n; i++) {
		B[i] = A[0][i];
	}
	return B;
}


/**********************************************************************************

the following section contains functions for commonly used matrix operations

Notes on certain terminology:

1. some function involves interchanging rows (partial pivoting) or rows and 
columns (full pivoting)

2. to use simple language to summarize the linear equation problem: A x = b, 
where A is a mxn (m rows, n columns) matrix, and b is a 1xn matrix (column 
vector). Suppose that the A matrix is not row or column degenerative, or not 
singular (and in fact not remotely close to singular despite round-off error), 
then (1) if m<n, there are either no solution or many solutions, and the 
solution space can be defined by singular value decomposition (2) if m>n, this 
reduce to a simple multiple linear regression problem (least square solution), 
(3) if m=n, then there are multiple ways to give unique answer to the equation, 
as shown in this section.

3. The matrix is called as row by column. For example, a matrix with 5 rows and 
3 columns were referred to as 5x3 matrix. This is the convention in linear 
algebra. Unfortunately, this directly contradict the natural writting of matrix 
in computer language: when inputting matrix, we first write all numbers in one 
row, then putting several rows together to form a matrix. This means that the 
first "dimension" of a matrix is the column, while the secon "dimension" is the 
row. Typically, this creates a problem for 3-dimensional or multi-dimensional 
data (for example, in Perl Data Language), but this issue is not considered 
here; instead we just follow the convention of rowxcolumn in refering to 
matrices.

More specific functions are briefly described below:

1. some standard matrix manipulation functions, such as transposing matrix or 
generating new matrices of various properties.

2. For solving linear equations, there are several methods

(1) Gauss-Jordan elimination with full pivoting, which is a robust solution, is 
3 times slower than LU decomposition, but it generates inverse matrix. In 
principle, since it generates inverse matrix, given a new right-side vector, one 
can immediately give a new solution by the product of inverse matrix and the new 
vector; howver, this is susceptible to round-off error. The implementation is 
adapted from NR.

(2) Gaussian elimination, which is 3 times faster than Gauss-Jordan (according 
to NR). It essentially reduces A to a upper triangular matrix via partial 
pivoting (row exchanges), then do a backward substitution. This method is 
implemented but it has not been thoroughly tested!

(3) LU decomposition, which is fast and can be applied to varying right-side 
matrices. It essentially reduce A to L x U, and L and U are both triangular 
matrices, so that we can do forward substitution on L and backward substitution 
on U to give linear solution. The ludecomp program implements 
partial pivoting (row permutation), and the lusubst implements forward and 
backward substitution.

3. Eigen systems

Given a linear transformation A, a non-zero vector x is defined to be an 
eigenvector of the transformation if it satisfies the eigenvalue equation Ax = 
lx for some scalar l (lamda). In this situation, the scalar l is called an 
eigenvalue of A corresponding to the eigenvector x. The eigenvalue equation 
means that under the transformation A eigenvectors experience only changes in 
magnitude and sign — the direction of Ax is the same as that of x

If x is an eigenvector of the linear transformation A with eigenvalue l, then 
any scalar multiple x is also an eigenvector of A with the same eigenvalue. 
Similarly if more than one eigenvector share the same eigenvalue l, any linear 
combination of these eigenvectors will itself be an eigenvector with eigenvalue 
l. Together with the zero vector, the eigenvectors of A with the same 
eigenvalue form a linear subspace of the vector space called an eigenspace.

4. Principle component analysis

Principle Component Analysis is just an eigen decomposition (or a singular value 
decomposition) of a covariance or a correlation matrix. It is a commonly used 
dimension reduction technique.

**********************************************************************************/

void MATRIX_ADVANCED () {
	printf ("This section contains advanced matrix operations\n");
}

double** matrix_corMatrix (double** A, int m, int n)
/* generate a correlation matrix for mxn matrix A,
 * correlations are calculated for each pairs of columns, returning nxn matrix B
 * this is an expensive (yet elegant) way of calculating correlation matrix!
 */
{
	double **B, **C, **D;
	
	B = matrix_copy (A, m, n);
	matrix_stdCol (B, m, n);
	C = matrix_transpose (B, m, n);
	D = matrix_product (C, B, n, m, n);
	matrix_multiply (D, n, n, 1.0/(double) (m-1));
	matrix_free (B);
	matrix_free (C);
	return D;
}

double** matrix_covMatrix (double** A, int m, int n)
/* generate a covariance matrix for mxn matrix A
 */
{
	double **B, **C, **D;
	double *col_mean;
	int i, j;
	
	B = matrix_copy (A, m, n);
	col_mean = matrix_colMean (B, m, n);
	
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++) {
			B[i][j] -= col_mean[j];
		}
	}
	C = matrix_transpose (B, m, n);
	D = matrix_product (C, B, n, m, n);
	matrix_multiply (D, n, n, 1.0/(double) (m-1));
	matrix_free (B);
	matrix_free (C);
	vector_free (col_mean);
	return D;
}


void matrix_solve (double **a, int n, double **b, int m)
/* Linear equation solution by Gauss-Jordan elimination with FULL pivoting
 * the function can (1) solve simultaneous linear equations, and (2) generates inverse matrix
 * A is a nxn matrix, B is a nxm matrix
 * after the solution
 * A becomes the inverse matrix, while B becomes the solution matrix
*/
{
	int *indxc,*indxr,*ipiv;
	int i,j,k,l,ll;
	int icol=0, irow=0;
	double big,dum,pivinv,temp;

	indxc=vectori_new (n);
	indxr=vectori_new (n);
	ipiv=vectori_new (n);
	for (j=0;j<n;j++) ipiv[j]=0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if (ipiv[j] != 1)
				for (k=0;k<n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1) kcerror("kc::solve: singular matrix-1");
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=0;l<m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) kcerror("kc::solve: singular matrix-2");
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=0;l<n;l++) a[icol][l] *= pivinv;
		for (l=0;l<m;l++) b[icol][l] *= pivinv;
		for (ll=0;ll<n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n-1;l>=0;l--) {
		if (indxr[l] != indxc[l])
			for (k=0;k<n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	vectori_free(ipiv);
	vectori_free(indxr);
	vectori_free(indxc);
}

void matrix_ludecomp(double **a, int n, int *indx, double *d)
/* perform LU decomposition on nxn matrix A
indx is an output vector that records the row permutation effected by the partial
pivoting; d is output as ±1 depending on whether the number of row interchanges was even
or odd, respectively.
*/
{
	int i,j,k;
	int imax=0;
	double big,dum,sum,temp;
	double *vv;

	vv=vector_new (n);
	*d=1.0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) kcerror("kc::matrix_ludecomp: singular matrix");
		vv[i]=1.0/big;
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i][j];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<n;i++) {
			sum=a[i][j];
			for (k=0;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<n;i++) a[i][j] *= dum;
		}
	}
	vector_free (vv);
}

void matrix_lusubst(double **a, int n, int *indx, double* b)
/* solve linear equations by forward and backward substitution.
A and indx are the output of ludcmp
A is indeed a hybrid of L and U matrix, where L has diagonal of 1 (not shown in A)
*/
{
	int i,ii=-1,ip,j;
	double sum;
	
	for (i=0;i<n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii>=0)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	
	for (i=n-1;i>=0;i--) {
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}

void matrix_inverse (double **a, int n)
/* use LU decomposition to find the inverse of a matrix */
{
	void matrix_ludecomp(double **a, int n, int *indx, double *d);
	void matrix_lusubst(double **a, int n, int *indx, double* b);
	int *indx, i, j;
	double d;
	double* b;
	double** newmat;
	
	indx = vectori_new (n);
	b = vector_new (n);
	newmat = matrix_new (n, n);
	
	matrix_ludecomp (a, n, indx, &d);		/*first perform LU decomposition*/
	for (j=0; j<n; j++) {
		for (i=0; i<n; i++) b[i] = 0.0;
		b[j] = 1.0;
		matrix_lusubst (a, n, indx, b);
		for (i=0; i<n; i++) newmat[i][j] = b[i];
	}
	
	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			a[i][j] = newmat[i][j];
		}
	}
	vectori_free (indx);
	vector_free (b);
	matrix_free (newmat);
}

double matrix_det (double **a, int n)
/* calculate the determinant of a matrix by LU decomposition */
{
	void matrix_ludecomp(double **a, int n, int *indx, double *d);
	double** matrix_copy (double **A, int m, int n);
	int* p;
	double d;
	int j;
	double** b;		/* copy of the a matrix */
	
	p = vectori_new (n);
	b = matrix_copy (a, n, n);
	matrix_ludecomp (b, n, p, &d);
	for(j=0; j<n; j++)
		d *= b[j][j];
	matrix_free (b);
	vectori_free (p);
	return d;
}

void matrix_ludecomp2 (const int n, double **a, int *p)
/* Modified From Perl Data Language:

 * matrix_ludecomp2 and matrix_lusubst2 are plain LU decomposition and substitution routines, 
 * The version here is Gaussian elimination 
 * with scaled row pivoting, and it is based on the algorithm given in Numerical Analysis,
 *       D. Kincaid and W. Cheney,
 *       Brooks/Cole Publishing Company, 1991.
 *
 * The parameters are used:
 *   n            the dimension of the matrix.
 *   a            the matrix; it contains both L and U at termination.
 *   p            permutation index.  
 *   b            the constant vector, and at termination it will contain
 *                the solution.
 *
 */
{
	int      i, j, k;             /* counters           */
	double   z;                   /* temporary real     */
	double  *s;                   /* pivot elements     */
	int      not_finished;        /* loop control var.  */
	int      i_swap;              /* swap var.          */
	double   temp;                /* another temp. real */
	
	s=vector_new(n);
	for(i=0; i<n; i++) {
		p[i]=i;
		s[i]=0.0;
		for(j=0; j<n; j++) {
			z=fabs(a[i][j]);
			if (s[i]<z)
			s[i]=z;
		} /* for j */
	} /* for i */
	
	for(k=0; k<(n-1); k++) {
		j=k-1;           /* select j>=k so ... */
		not_finished=1;
		while (not_finished) {
			j++;
			temp=fabs(a[p[j]][k]/s[p[j]]);
			for(i=k; i<n; i++)  
				if (temp>=(fabs(a[p[i]][k])/s[p[i]]))
					not_finished=0;          /* end loop */
		} /* while */
		i_swap=p[k];
		p[k]=p[j];
		p[j]=i_swap;
		temp=1.0/a[p[k]][k];
		for(i=(k+1); i<n; i++) {
			z=a[p[i]][k]*temp;
			a[p[i]][k]=z;
			for(j=(k+1); j<n; j++) 
			a[p[i]][j]-=z*a[p[k]][j];
		} /* for i */
	} /* for k */
	
	vector_free(s);
}

void matrix_lusubst2(const int n, double **a, int *p, double *b)
{
	int        i, j, k;           /* counters               */
	double     sum;               /* temporary sum variable */
	double    *x;                 /* solution               */
	
	x=vector_new(n);
	
	for(k=0; k<(n-1); k++)        /* forward subst */
		for(i=(k+1); i<n; i++)
			b[p[i]]-=a[p[i]][k]*b[p[k]];

	for(i=(n-1); i>=0; i--) {     /* back subst */
		sum=b[p[i]];
		for(j=(i+1); j<n; j++)
			sum-=a[p[i]][j]*x[j];
		x[i]=sum/a[p[i]][i];
	}
	for(i=0; i<n; i++)             /* copy solution */
		b[i]=x[i];
	vector_free (x);
}

void matrix_inverse2 (double **a, int n)
/* use LU decomposition to find the inverse of a matrix */
{
	void matrix_ludecomp2 (const int n, double **a, int *p);
	void matrix_lusubst2(const int n, double **a, int *p, double *b);
	int *indx, i, j;
	double* b;
	double** newmat;
	
	indx = vectori_new (n);
	b = vector_new (n);
	newmat = matrix_new (n, n);
	
	matrix_ludecomp2 (n, a, indx);		/*first perform LU decomposition*/
	for (j=0; j<n; j++) {
		for (i=0; i<n; i++) b[i] = 0.0;
		b[j] = 1.0;
		matrix_lusubst2 (n, a, indx, b);
		for (i=0; i<n; i++) newmat[i][j] = b[i];
	}
	

	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			a[i][j] = newmat[i][j];
		}
	}
	vectori_free (indx);
	vector_free (b);
	matrix_free (newmat);
}

void matrix_tridiag(const int n, double *a, double *d, double *c, double *b)
/* Tridiag solves a tridiagonal linear system of equations. A row in the 
 * matrix is (a_i, d_i, c_i). The vector b is the 
 * constant vector, and it will contain the solution at termination.
 * Function Adapted From Perl Data Language
 */
{
	
	int        i;            /* counter */
	double    *x;            /* solution */
	
	x = vector_new (n);
	
	for(i=1; i<n; i++) {
		d[i]-=(a[i-1]/d[i-1])*c[i-1];
		b[i]-=(a[i-1]/d[i-1])*b[i-1];
	}
	x[n-1]=b[n-1]/d[n-1];
	for(i=(n-2); i>=0; i--) 
		x[i]=(b[i]-c[i]*b[i+1])/d[i];
	
	for(i=0; i<n; i++)
		b[i]=x[i];
	vector_free (x);
}

double pythag(double a, double b)
/* used by the matrix_svd function */
{
    double at = fabs(a), bt = fabs(b), ct, result;

    if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
    else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
    else result = 0.0;
    return(result);
}

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
void matrix_svd(double **a, int m, int n, double *w, double **v)
/* Singular Value Decomposition
 * Code adapted from Numerial Recipe
 Given a matrix a[0..m-1][0..n-1], this routine computes its singular value decomposition, A =
U·W·V^T. The matrix U (mxn) replaces a on output. The diagonal matrix of singular values W is output
as a vector w[0..n-1]. The matrix V (not the transpose V^T ) is output as v[0..n-1][0..n-1].

 * see the matrix_svd2 function for another implementation from Perl Data Language

*/
{
    double pythag(double a, double b);
    int flag, i, its, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    double *rv1;
  
    if (m < n) kcerror ("kc::matrix_svd: #rows should be greater than or equal to #cols (please create more rows and fill in zero)");
  
    rv1 = vector_new (n);

	/* Householder reduction to bidiagonal form */
    for (i = 0; i < n; i++) 
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m) 
        {
            for (k = i; k < m; k++) 
                scale += fabs(a[k][i]);
            if (scale) 
            {
                for (k = i; k < m; k++) 
                {
                    a[k][i] /= scale;
                    s += a[k][i] * a[k][i];
                }
                f = a[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][i] = f - g;
                if (i != n - 1) 
                {
                    for (j = l; j < n; j++) 
                    {
                        for (s = 0.0, k = i; k < m; k++) 
                            s += a[k][i] * a[k][j];
                        f = s / h;
                        for (k = i; k < m; k++) 
                            a[k][j] += f * a[k][i];
                    }
                }
                for (k = i; k < m; k++) 
                    a[k][i] *= scale;
            }
        }
        w[i] = scale * g;
    
        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != n - 1) 
        {
            for (k = l; k < n; k++) 
                scale += fabs(a[i][k]);
            if (scale) 
            {
                for (k = l; k < n; k++) 
                {
                    a[i][k] /= scale;
                    s += a[i][k] * a[i][k];
                }
                f = a[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][l] = f - g;
                for (k = l; k < n; k++) 
                    rv1[k] = a[i][k] / h;
                if (i != m - 1) 
                {
                    for (j = l; j < m; j++) 
                    {
                        for (s = 0.0, k = l; k < n; k++) 
                            s += a[j][k] * a[i][k];
                        for (k = l; k < n; k++) 
                            a[j][k] += s * rv1[k];
                    }
                }
                for (k = l; k < n; k++) 
                    a[i][k] *= scale;
            }
        }
        if (anorm < fabs(w[i]) + fabs(rv1[i])) anorm = fabs(w[i]) + fabs(rv1[i]);
        /*anorm = FMAX(anorm, (fabs(w[i]) + fabs(rv1[i])));*/
    }
  
    /* accumulate the right-hand transformation */
    for (i = n - 1; i >= 0; i--) 
    {
        if (i < n - 1) 
        {
            if (g) 
            {
                for (j = l; j < n; j++)
                    v[j][i] = (a[i][j] / a[i][l]) / g;
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < n; k++) 
                        s += a[i][k] * v[k][j];
                    for (k = l; k < n; k++) 
                        v[k][j] += s * v[k][i];
                }
            }
            for (j = l; j < n; j++) 
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }
  
    /* accumulate the left-hand transformation */
    for (i = n - 1; i >= 0; i--) 
    {
        l = i + 1;
        g = (double)w[i];
        if (i < n - 1) 
            for (j = l; j < n; j++) 
                a[i][j] = 0.0;
        if (g) 
        {
            g = 1.0 / g;
            if (i != n - 1) 
            {
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < m; k++) 
                        s += a[k][i] * a[k][j];
                    f = (s / a[i][i]) * g;
                    for (k = i; k < m; k++) 
                        a[k][j] += f * a[k][i];
                }
            }
            for (j = i; j < m; j++) 
                a[j][i] *= g;
        }
        else 
        {
            for (j = i; j < m; j++) 
                a[j][i] = 0.0;
        }
        ++a[i][i];
    }

    /* diagonalize the bidiagonal form */
    for (k = n - 1; k >= 0; k--) 
    {                             /* loop over singular values */
        for (its = 0; its < 100; its++) 
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--) 
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm) 
                {
                    flag = 0;
                    break;
                }
                if (fabs(w[nm]) + anorm == anorm) 
                    break;
            }
            if (flag) 
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++) 
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm) 
                    {
                        g = w[i];
                        h = pythag(f, g);
                        w[i] = h; 
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++) 
                        {
                            y = a[j][nm];
                            z = a[j][i];
                            a[j][nm] = y * c + z * s;
                            a[j][i] = z * c - y * s;
                        }
                    }
                }
            }
            z = w[k];
            if (l == k) 
            {                  /* convergence */
                if (z < 0.0) 
                {              /* make singular value nonnegative */
                    w[k] = -z;
                    for (j = 0; j < n; j++) 
                        v[j][k] = -v[j][k];
                }
                break;
            }
            if (its >= 100) {
                vector_free (rv1);
                kcerror ("kc::matrix_svd: No convergence after 100 iterations");
            }
    
            /* shift from bottom 2 x 2 minor */
            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = pythag(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
          
            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++) 
            {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = pythag(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < n; jj++) 
                {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = x * c + z * s;
                    v[jj][i] = z * c - x * s;
                }
                z = pythag(f, h);
                w[j] = z;
                if (z) 
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++) 
                {
                    y = a[jj][j];
                    z = a[jj][i];
                    a[jj][j] = y * c + z * s;
                    a[jj][i] = z * c - y * s;
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }
    vector_free (rv1);
}
#undef SIGN

#define TOLERANCE 1.0e-22
void matrix_svd2 (double *W, double *Z, int nRow, int nCol)
/* This is copied from Perl Data Language, the validity of this function has not been verified yet!!!

This SVD routine is based on pgs 30-48 of "Compact Numerical Methods
   for Computers" by J.C. Nash (1990), used to compute the pseudoinverse.
   Modifications include:
        Translation from Pascal to ANSI C.
        Array indexing from 0 rather than 1.
        Float replaced by double everywhere.
        Support for the Matrix structure.
        I changed the array indexing so that the matricies (float [][])
           could be replaced be a single list (double *) for more
           efficient communication with Mathematica.
*/
{
  int i, j, k, EstColRank, RotCount, SweepCount, slimit;
  double eps, e2, tol, vt, p, x0, y0, q, r, c0, s0, d1, d2;  
  eps = TOLERANCE;
  slimit = nCol/4;
  if (slimit < 6.0)
    slimit = 6;
  SweepCount = 0;
  e2 = 10.0*nRow*eps*eps;
  tol = eps*.1;
  EstColRank = nCol;
  for (i=0; i<nCol; i++) {
    for (j=0; j<nCol; j++) {
      W[nCol*(nRow+i)+j] = 0.0;
    }
    W[nCol*(nRow+i)+i] = 1.0;  /* rjrw 7/7/99: moved this line out of j loop */
  }
  RotCount = EstColRank*(EstColRank-1)/2;
  while (RotCount != 0 && SweepCount <= slimit)
    {
      RotCount = EstColRank*(EstColRank-1)/2;
      SweepCount++;
      for (j=0; j<EstColRank-1; j++)
        {
          for (k=j+1; k<EstColRank; k++)
            {
              p = q = r = 0.0;
              for (i=0; i<nRow; i++)
                {
                  x0 = W[nCol*i+j]; y0 = W[nCol*i+k];
                  p += x0*y0; q += x0*x0; r += y0*y0;
                }
              Z[j] = q; Z[k] = r;
              if (q >= r)
                {
                  if (q<=e2*Z[0] || fabs(p)<=tol*q) RotCount--;
                  else
                    {
                      p /= q; r = 1 - r/q; vt = sqrt(4*p*p+r*r);
                      c0 = sqrt(fabs(.5*(1+r/vt))); s0 = p/(vt*c0);
                      for (i=0; i<nRow+nCol; i++)
                        {
                          d1 = W[nCol*i+j]; d2 = W[nCol*i+k];
                          W[nCol*i+j] = d1*c0+d2*s0; W[nCol*i+k] = -d1*s0+d2*c0;
                        }
                    }
                }
              else
                {
                  p /= r; q = q/r-1; vt = sqrt(4*p*p+q*q);
                  s0 = sqrt(fabs(.5*(1-q/vt)));
                  if (p<0) s0 = -s0;
                  c0 = p/(vt*s0);
                  for (i=0; i<nRow+nCol; i++)
                    {
                      d1 = W[nCol*i+j]; d2 = W[nCol*i+k];
                      W[nCol*i+j] = d1*c0+d2*s0; W[nCol*i+k] = -d1*s0+d2*c0;
                    }
                }
            }
        }
      while (EstColRank>=3 && Z[(EstColRank-1)]<=Z[0]*tol+tol*tol)
        EstColRank--;
    }
}
#undef TOLERANCE

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);
void matrix_eigen (double **a, int n, double *d, double **v, int *nrot)
/* eigen value and eigen vector calculation by Jacobi method 
 */
{
	int j,iq,ip,i;
	double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

	b=vector_new(n);
	z=vector_new(n);
	for (ip=0;ip<n;ip++) {
		for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=0;ip<n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	*nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=0;ip<=n-2;ip++) {
			for (iq=ip+1;iq<n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
			vector_free(z);
			vector_free(b);
			return;		/*quadratic convergence to machine underflow*/
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=0;ip<=n-2;ip++) {
			for (iq=ip+1;iq<n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && (fabs(d[ip])+g) == fabs(d[ip])
					&& (fabs(d[iq])+g) == fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((fabs(h)+g) == fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=0;j<=ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=0;j<n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=0;ip<n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	vector_free (b);
	vector_free (z);
	kcerror("kc::matrix_eigen: No convergence after 50 iterations");
}
#undef ROTATE

double** matrix_pcaCor (double** A, int m, int n, double *eval)
/* perform Principle Component Analysis on mxn matrix (for example, representing m samples, n variables)
 * using correltaion matrix method
 * returns a nxn eigen vector matrix, where each column contains a vector
 * the eigen values are stored in the eval array
 * the procedure is briefly described below:
 *         1. normalize the column (variable) by (X-mean)/SD
 *         2. calculate correlation matrix for columns
 *         3. perform eigen decomposition on the correlation matrix
 *         4. sort the eigen values and eigen vector, change signs
 */
{
	double **B, *C, **D;
	double *col_mean;
	int nrot, i, j;
	
	B = matrix_corMatrix (A, m, n);
	C = vector_new (n);
	D = matrix_new (n, n);

	matrix_eigen (B, n, C, D, &nrot);
	matrix_sortByVector (D, n, n, C);
	col_mean = matrix_colMean (D, n, n);
	for (i=0; i<n; i++) {
		if (col_mean[i] < 0) {
			for (j=0; j<n; j++) {
				D[j][i] = -D[j][i];
			}
		}
		eval[i] = C[i];
	}
	
	matrix_free (B);
	vector_free (C);
	vector_free (col_mean);
	return D;
}

double** matrix_pcaCov (double** A, int m, int n, double* eval)
/* perform Principle Component Analysis on mxn matrix (for example, representing m samples, n variables)
 * by the covariance matrix method
 * returns a nxn eigen vector matrix, where each column contains an eigen vector
 * the eigen values are stored in the eval array
 * the procedure is briefly described below:
 *         1. normalize the column by (X-mean)
 *         2. calculate covariance matrix for columns
 *         3. perform eigen decomposition on the covariance matrix
 *         4. sort the eigen vectors and values
 */
{
	double **B, *C, **D;
	double *col_mean;
	int nrot, i;
	
	B = matrix_covMatrix (A, m, n);
	C = vector_new (n);
	D = matrix_new (n, n);
	
	matrix_eigen (B, n, C, D, &nrot);
	matrix_sortByVector (D, n, n, C);
	col_mean = matrix_colMean (D, n, n);
	for (i=0; i<n; i++) {
		eval[i] = C[i];
	}
	
	/* TO BE IMPLEMENTED */
	/* I have not had time to write the code here (July 2008) */
	
	matrix_free (B);
	vector_free (C);
	vector_free (col_mean);
	return D;
}

double** matrix_pcaScore (double** A, int m, int n, double** evec, int k)
/* calculate the scores (or re-mapping of the original data) by the first k components of PCA
 * the eigen vector is assumed to be known and have been computed
 */
{
	double **evec2, **score;
	
	evec2 = matrix_newSubcol (evec, n, n, 1, k);
	score = matrix_product (A, evec2, m, n, k);
	matrix_free (evec2);
	return score;
}

double* matrix_pcaComVar (double** A, int m, int n, double** evec)
/* calculate the variance that each principle component can explain 
 * numerically should be identical to eigen value
 */
{
	double **B, **C;
	
	B = matrix_copy (A, m, n);
	matrix_stdCol (B, m, n);
	C = matrix_product (B, evec, m, n, n);
	return matrix_colVar (C, m, n);
}

double** matrix_reg (double** X, int m, int n, double** Y)
/* perform multiple least-square linear regression using matrix algebra
 * Note: for simple linear regression, the reg_linear function can be directly applied on Perl arrays
 * Note: this function here operates on C matrices
 * X is a predictor matrix, with m samples and n-1 variables (the first column is always 1)
 * Y is the response matrix, with m samples and 1 response column
 * return the coefficient matrix (nx1 dimension)
 */
{
	double **Xt, **XtX, **XtXiXt, **B;
	
	if (m<=n) kcerror ("kc::matrix_reg: the rows must more than columns for regression");
	
	Xt = matrix_transpose (X, m, n);
	XtX = matrix_product (Xt, X, n, m, n);
	matrix_inverse (XtX, n);
	XtXiXt = matrix_product (XtX, Xt, n, n, m);
	B = matrix_product (XtXiXt, Y, n, m, 1);
	
	matrix_free (Xt);
	matrix_free (XtX);
	matrix_free (XtXiXt);
	return B;
}

double matrix_regMSE (double** X, int m, int n, double** Y, double** B)
/* calculate mean squared error post regression
 * it turns out that the residue can be easily calculated as Y-XB
 */
{
	double **P;
	double mse=0;
	int i;
	
	P = matrix_product (X, B, m, n, 1);
	for (i=0; i<m; i++) {
		mse += (Y[i][0]-P[i][0])*(Y[i][0]-P[i][0]);
	}
	mse /= (double) (m-n);
	matrix_free (P);
	return mse;
}

double* matrix_regModelStat (double** X, int m, int n, double** Y, double** B)
/* evaluate the model fit to the data
 * return a vector that contains pre-defined list of model statistics
 * check the source code to know what each field is, I'll write a more detailed description in the future
 */
{
	double mse, msm, mst, sse=0, ssm=0, sst=0, y=0, f, p, r2, r2_adj;
	double **P, *stat;
	int i;
	
	for (i=0; i<m; i++) {
		y += Y[i][0];
	}
	y /= m;
	
	P = matrix_product (X, B, m, n, 1);
	for (i=0; i<m; i++) {
		sse += (Y[i][0]-P[i][0])*(Y[i][0]-P[i][0]);
		ssm += (y-P[i][0])*(y-P[i][0]);
		sst += (y-Y[i][0])*(y-Y[i][0]);
	}
	mse = sse / (double)(m-n);
	msm = ssm / (double)(n-1);
	mst = sst / (double)(m-1);
	
	f = msm/mse;
	p = 1-cdf_f (n-1, m-n, f);
	r2 = ssm/sst;
	r2_adj = 1-(1-r2)*((double)(m-1))/((double)(m-n));
	stat = vector_new (10);
	stat[0] = p;
	stat[1] = f;
	stat[2] = ssm;
	stat[3] = sse;
	stat[4] = sst;
	stat[5] = n-1;
	stat[6] = m-n;
	stat[7] = m-1;
	stat[8] = r2;
	stat[9] = r2_adj;
	return stat;
}

double** matrix_regCoefStat (double** X, int m, int n, double** Y, double** B)
/* inference of individual coefficient in regression, including the SE, t and P-value for each coefficient
 * the P-value tests the signicance of a variable given that the other variables are already in the model
 * the procedure is below:
 *         1. first create the MSEx(X'X)^-1 matrix, where diagonal line is the standard error for coefficient
 *         2. next calculate the t=b/SE, then df=m-n, then P-value
 *         3. next calculate confidence interval
 * the resulting matrix contains n rows, with columns being t, p, lower 95% CI, upper 95% CI, respectively
 */
{
	double mse, tthres;
	double **Xt, **XtX, **I;
	int i;

	I = matrix_new (n, 5);

	mse = matrix_regMSE (X, m, n, Y, B);	
	tthres = cdfinv_t ((double)(m-n), 0.975);			/*the critical t-value for the CDF=0.975 for confidence interval calculation*/	
	Xt = matrix_transpose (X, m, n);
	XtX = matrix_product (Xt, X, n, m, n);
	matrix_inverse (XtX, n);
	
	matrix_multiply (XtX, n, n, mse);

	for (i=0; i<n; i++) {
		I[i][0] = sqrt(XtX[i][i]);				/*standard error of coefficient*/
		I[i][1] = B[i][0]/I[i][0];				/*t value of coefficient*/
		I[i][2] = 2*cdf_t((double)(m-n), fabs(I[i][1]));	/*p value of coefficient*/
		if (I[i][2]>1) I[i][2] = 2-I[i][2];
		
		I[i][3] = B[i][0] - tthres*I[i][0];
		I[i][4] = B[i][0] + tthres*I[i][0];
	}
	matrix_free (Xt);
	matrix_free (XtX);
	return I;
}

double** matrix_regObsStat (double** X, int m, int n, double** Y, double** B)
/* inference of individual observation in regression
 * the results are returned in a mx5 matrix
 * the first column is expected value of response variable
 * the 2-4 columns are the SE and confidence interval for the subpopulation mean corresponding to the set of explanatory variables X
 * the 5-7 columns are the SE and confidence interval that expresses the uncertainty in our prediction of a single given observation
 * the 8-10 columns are the residue, the standardized residue and the studentized residue, respectively
 * the 11 column is the P-value testing for outliers using studentized residue (which follows m-1-n df t-distribution)
 * the 12 column is the leverage of the data point (diagonal element of HAT matrix)
 */
{
	double **P, **stat, **Xt, **XtXi, **H, *vec;
	double tthres, mse, msed, s2_mean, s2_obs;
	int i, j, k;
	
	P = matrix_product (X, B, m, n, 1);
	stat = matrix_new (m, 12);
	for (i=0; i<m; i++) {
		stat[i][0] = P[i][0];
	}
	
	mse = matrix_regMSE (X, m, n, Y, B);
	tthres = cdfinv_t ((double)(m-n), 0.975);			/*the critical t-value for the CDF=0.975 for confidence interval calculation*/
	Xt = matrix_transpose (X, m, n);
	XtXi = matrix_product (Xt, X, n, m, n);
	matrix_inverse (XtXi, n);
	
	H = matrix_product3 (X, XtXi, Xt, m, n, n, m);

	vec = vector_new (n);
	for (i=0; i<m; i++) {
		vector_zero (vec, n);
		for (j=0; j<n; j++) {
			for (k=0; k<n; k++) {
				vec[j] += X[i][k]*XtXi[k][j];
			}
		}
		s2_mean = 0;
		for (j=0; j<n; j++) {
			s2_mean += vec[j]*X[i][j];
		}
		s2_mean *= mse;
		
		stat[i][1] = sqrt(s2_mean);
		stat[i][2] = stat[i][0] - sqrt(s2_mean)*tthres;
		stat[i][3] = stat[i][0] + sqrt(s2_mean)*tthres;
		
		s2_obs = s2_mean + mse;
		stat[i][4] = sqrt(s2_obs);
		stat[i][5] = stat[i][0] - sqrt(s2_obs)*tthres;
		stat[i][6] = stat[i][0] + sqrt(s2_obs)*tthres;
		
		stat[i][7] = Y[i][0] - stat[i][0];
		stat[i][8] = stat[i][7] / sqrt(mse*(1-H[i][i]));
		
		msed = ((m-n)*mse - stat[i][7]*stat[i][7]/(1-H[i][i])) / (m-n-1);	/* MSE of deleting this observation*/
		stat[i][9] = stat[i][7] / sqrt(msed*(1-H[i][i]));
		
		stat[i][10] = 2*cdf_t (m-n-1, stat[i][9]);
		if (stat[i][10]>1) stat[i][10]=2-stat[i][10];
		
		stat[i][11] = H[i][i];
	}
	vector_free (vec);
	matrix_free (P);
	matrix_free (Xt);
	matrix_free (XtXi);
	matrix_free (H);
	return stat;
}
	
	

void test1 ()
{
	printf ("INFO: size of int=%lu size of long=%lu size of float=%lu size of double=%lu\n", (long unsigned int) sizeof(int), (long unsigned int) sizeof(long), (long unsigned int) sizeof(float), (long unsigned int) sizeof(double));
	kcerror ("kc::test1: testing error handling functionality");
}




/***********************************************************************************
the following section contains statistical distributions, including PDF/PMF and CDF
***********************************************************************************/

void STATISTICAL_DISTRIBUTION () {
	printf ("This section contains functions for statistical distributions\n");
}

double cdf_binomial (int n, int k, double p)
/* returns the probability of k-1 or less successes in n trials when the
        probability of a success on a single trial is p.  Here n > k > 0 must be non-negative
        integers, and p > 0.
 * see bitest() for another implementation specifically for statistical testing
*/
{
	double betai(double a, double b, double x);
	return 1-betai (k, n-k+1, p);
}

double cdf_chi2 (double df, double x)
/*cdf_chi2(n,x) returns the cumulative chi-squared distribution with n degrees of freedom for n > 0; chi2(df,x) = gammap(df/2,x/2).  Returns 0 unless x>0.
*/
{
	double gammp(double a, double x);
	double output;
	output = gammp(df/2.0, x/2.0);
	return output;
}

double cdf_poisson (int k, double x)
/* cumulative distribution function for the Poisson distribution
Usually the "x" is referred to as an integer lamda (for example, http://en.wikipedia.org/wiki/Poisson_distribution)
However, here we use "x" instead, to indicate that the mean does not have to be an integer!
the value is the probability that the number of Poisson random events occurring will be between 0 and k -1 inclusive, when expected mean number is x
*/
{
	double gammp(double a, double x);
	int factorial (int i);
	double sum = 0;
	int i;
	if (k <= 12) {
		for (i=0; i<k; i++) {
			sum += exp(-x)*pow(x, i)/factorial(i);
		}
		return sum;
	} else {
		return 1-gammp (k, x);
	}
}

double cdf_poisson_inc (int k, double x)
/* this function differs from cdf_poisson by the inclusion of the "observed" event in the CDF calculation
the value is the probability that the number of Poisson random events occurring will be between 0 and k inclusive, when expected mean number is x
*/
{
	double cdf_poisson (int k, double x);
	double pdf_poisson (int k, double x);
	return cdf_poisson (k, x) + pdf_poisson (k, x);
}




double cdf_normal (double x, double mu, double sigma)
/* Returns the cumulative probability density function for a normal distribution with mean as mu and standard deviation as sigma
cumulative normal distribution
                  x    2
         1      /   -t  / 2
     ---------- |  e         dt
     sqrt(2 pi) /
                 -inf
*/
{
	return (1+erf((x-mu)/(sigma*sqrt(2))))/2;
}

double cdf_stdnormal (double x)
/*Returns the cumulative density function for a standard normal distribution (see http://en.wikipedia.org/wiki/Normal_distribution)*/
{
	return (1+erf(x/sqrt(2)))/2;
}


double cdf_f (double df1, double df2, double x)
/*F(n1,n2,x) returns the cumulative F distribution with n1 numerator and n2 denominator degrees
        of freedom, for n1, n2 > 0; returns 0 unless x > 0.
*/
{
	double betai(double a, double b, double x);
	double output;
	if (x <= 0) return 0;
	output = 1 - betai (df2/2.0, df1/2.0, df2/(df2+df1*x));
	return output;
}

double cdf_t (double df, double x)
/* returns the cumulative Student's t distribution with n > 0 degrees of freedom.
 * the new formula is based on http://en.wikipedia.org/wiki/Student's_t-distribution
 */
{
	double betai(double a, double b, double x);
	double output;
	//output = 1 - betai (0.5*df, 0.5, df/(df+x*x)) / 2;
	output = betai(0.5*df, 0.5*df, (x+sqrt(x*x+df))/(2.0*sqrt(x*x+df)));
	return output;
}

#define EPSILON 1e-10
double cdfinv_t (double df, double p)
/* returns the inverse cumulative student's t distribution
 * in general, there is no tractable solution: http://en.wikipedia.org/wiki/Quantile_function
 * however, the approximate solution can be found by iterative methods.
 */
{
	double tthres, alpha, q;
	double h=1000000, l=-1000000, ph, pl, tp;				/*high and low threshold value*/
	if (df == 1) {
		tthres = tan(PI*(p-0.5));
	} else if (df == 2) {
		tthres = (2*p-1)/sqrt(2*p*(1-p));
	} else if (df == 4) {
		alpha = 4*p*(1-p);
		q = 4.0 / sqrt(alpha) * cos(acos(sqrt(alpha))/3.0);
		tthres = sqrt(q-4);
		if (p<0.5) tthres = -tthres;
	} else {
		ph = cdf_t(df, h);
		pl = cdf_t(df, l);
		if (p<pl) {
			tthres = l;
		} else if (p>ph) {
			tthres = h;
		} else {
			tthres = (h+l)/2.0;
			tp = cdf_t(df, tthres);
			while(fabs(tp-p)>EPSILON) {
				if (cdf_t(df, (h+l)/2.0) < p) {
					l=(h+l)/2.0;
				} else {
					h=(h+l)/2.0;
				}
				tp = cdf_t(df, (h+l)/2.0);
			}
			tthres = (l+h)/2.0;
		}
	}
	return tthres;
}
#undef EPSILON

double pdf_stdnormal (double x)
/*Returns the probability density function for a standard normal distribution (see http://en.wikipedia.org/wiki/Normal_distribution)*/
{
	return exp(-x*x/2)/sqrt(2*PI);
}

double pdf_normal (double x, double mu, double sigma)
/*Returns the probability density function for a normal distribution with mean as mu and standard deviation as sigma*/
{
	return exp(-(x-mu)*(x-mu)/(2*sigma*sigma)) / (sigma * sqrt(2*PI));
}

double pdf_binomial (int n, int k, double p)
/* calculate the PDF for binomial distribution
*/
{
	double bico (int n, int k);
	if (k > n) return 0;
	if (p == 0) {
		return k==0 ? 1 : 0;
	} else if (p == 1) {
		return k==n ? 1 : 0;
	} else {
		return exp(lnbico(n, k) + k * log(p) + (n-k) * log(1-p));
	}
}

double pdf_beta (double a, double b, double x)
/* calculate the PDF of beta distribution (http://en.wikipedia.org/wiki/Beta_distribution)
Note the difference between "beta distribution" and "beta function"! the beta function is used in the calculation below.
a, b > 0; returns 0 unless 0 <= x <= 1.
*/
{
	double beta(double z, double w);
	if (x<0 || x>1) return 0;			/*the beta distribution has support interval of [0,1]*/
	return pow(x, a-1) * pow(1-x, b-1) / beta (a, b);
}

double pdf_poisson (int k, double x)
/* calculate the probability mass function of poisson distribution
*/
{
	int factorial (int x);
	double lnfactorial (double x);
	if (k <= 12) {
		return exp(-x) * pow(x, k) / factorial(k);
	} else {
		return exp(-x) * pow(x, k) / exp(lnfactorial(k));
	}
}

double pdf_geometric (int k, double p)
/*the probability distribution of the number X of Bernoulli trials needed to get one success, supported on the set { 1, 2, 3, ...}
*/
{
	return pow (1-p, k-1) * p;
}

double pdf_hypergeometric (int a, int b, int c, int d)
/* hypergeometric distribution is a discrete probability distribution that describes the number of successes in a sequence of n draws from a finite population without replacement.
NOTE THAT IT IS DRAWN WITHOUT REPLACEMENT!
the four numbers are four cells in a 2x2 contingency table
*/
{
	double bico (int n, int k);
	return bico (a+b, a) * bico(c+d, c) / bico (a+b+c+d, a+c);
}

double bico (int n, int k)
/*calculate the binomial coefficient (n, k)
*/
{
	double gammln(double x);
        return floor(0.5+exp(gammln(n+1)-gammln(k+1)-gammln(n-k+1)));
}

double lnbico (int n, int k)
{
	double gammln(double x);
        return gammln(n+1)-gammln(k+1)-gammln(n-k+1);
}

double invnormal (double x)
/*inverse standard normal cumulative distribution function, or quantile function, can be expressed in terms of the inverse error function*/
{
	void kcerror(char *error_text);
	double inverff (double x);
	if (x>1 || x<0) kcerror ("kc::invnormal: x (probability value) should be between 0 and 1");
	return sqrt(2) * inverff (2*x-1);
}


void reg_linear (double *x, double *y, int ndata, double *a, double *b, double *F, double *P)
/* perform linear regression of two variables with F and P value reported
 * for multiple linear regression, a temporary solution is to use a series of matrix operation subroutines, including matrix_reg
 * in the future, I may add a subroutine specifically for multiple linear regression
*/
{
	int i;
	double t,sxoss,syoss,sx=0.0,sy=0.0, st2=0.0;
	double mss=0.0,rss=0.0,tss=0.0;			/*model sum of squares, residual sum of squares, total sum of squares*/
	*b=0.0;
	
	for (i=0;i<ndata;i++) {
		sx += x[i];
		sy += y[i];
	}
	sxoss=sx/(double) ndata;
	syoss=sy/(double) ndata;
	for (i=0;i<ndata;i++) {
		t=x[i]-sxoss;
		st2 += t*t;
		*b += t*y[i];
	}
	*b /= st2;					/*Solve for a, b*/
	*a=(sy-sx*(*b))/(double) ndata;

	for (i=0;i<ndata;i++) {
		tss += (y[i]-syoss)*(y[i]-syoss);
		rss += (y[i]-*a-*b*x[i])*(y[i]-*a-*b*x[i]);
	}
	mss = tss-rss;
	*F = mss / (rss/(double) (ndata-2));
	*P = 1-cdf_f(1.0, (double) (ndata-2), *F);
}



/********************************************************************************
the following section contains summary statistics
*********************************************************************************/

void SUMMARY_STATISTICS ()
{
	printf ("This section contains summary statistics\n");
}

double mean(double *data, int n)
/*Given array data[0..n-1], returns its mean as ave and its variance as var.
see mean2 for running statistics (continuously updating statistics) calculation for mean
*/
{
	int j;
	double ave = 0.0;
	for (j=0;j<n;j++) ave += data[j];
	ave /= n;
	return ave;
}

double mean2(double *data, int n)
/*Given array data[0..n-1], returns its mean as ave and its variance as var.
unlike the mean subroutine, this subroutine continuously update the mean by calculating running statistics
*/
{
	double ave;
	int i;
	ave = data[0];
	for (i=1; i<n; i++) {
		ave += (data[i] - ave) / (i + 1);
	}
	return ave;
}

void avevar(double *data, int n, double *ave, double *var)
/* Given array data[0..n-1], returns its mean as ave and its variance as var.
this subroutine uses idea from Numerical Recipe.
it uses a two-pass formula to correct roundoff errors for variance calculation
see avevar2 subroutine for running statistics (continuously updating statistics) for mean and avariance calculation
*/
{
	int j;
	double s,ep;
	for (*ave=0.0,j=0;j<n;j++) *ave += data[j];
	*ave /= n;
	*var=ep=0.0;
	for (j=0;j<n;j++) {
		s=data[j]-(*ave);
		ep += s;
		*var += s*s;
	}
	*var=(*var-ep*ep/n)/(n-1); 			//Corrected two-pass formula
}

void avevar2(double *data, int n, double *ave, double *var)
/* stably updating mean and variance
D.H.D. West, Updating mean and variance estimates: an improved method, Comm ACM 22:9, 532 (1979)  
it assumes that you get a weight and a data value (W_i and X_i) that you use to update the estimates XBAR and S2:

    SUMW = W_1
    M = X_1
    T = 0
    For i=2,3,...,n 
    {
       Q = X_i - M
       TEMP = SUM + W_i    // typo: He meant SUMW (I think so)
       R = Q*W_i/TEMP
       M = M + R
       T = T + R*SUMW*Q
       SUMW = TEMP
    }
    XBAR = M
    S2 = T*n/((n-1)*SUMW)
*/
{
	double sumw, m, t, q, temp, r;
	int i;
	
	sumw = 1;			//all weights are treated as 1
	m = data[0];
	t = 0;
	for (i=1; i<n; i++) {
		q = data[i] - m;
		temp = sumw + 1;	//temp is the sum of all weights for previous records
		r = q / temp;		//r is the contribution of new entry to the mean
		m += r;			//m is the current mean
		t += r * sumw * q;
		sumw = temp;
	}
	*ave = m;
	*var = t*n/((n-1)*sumw);
}

void averms(double *data, int n, double *ave, double *rms)
/* Given array data[0..n-1], returns its mean as ave and its Root Mean Square devaition from the mean
 * RMS = sqrt(sum( (x-mean(x))^2 )/N)
 * RMS is also known as the root-mean-square deviation, or the square root of the variance
 * The RMS measure differ from Standard Deviation by the denominator (the former has n, but the latter has n-1)
*/
{
	int j;
	double s,ep;
	for (*ave=0.0,j=0;j<n;j++) *ave += data[j];
	*ave /= n;
	*rms=ep=0.0;
	for (j=0;j<n;j++) {
		s=data[j]-(*ave);
		ep += s;
		*rms += s*s;
	}
	*rms=(*rms-ep*ep/n)/(double) n; 			/* Corrected two-pass formula */
	*rms = sqrt (*rms);
}

void moment(double *data, int n, double *ave, double *adev, double *sdev, double *var, double *skew, double *curt)
/*Given an array of data[0..n-1], this routine returns its mean ave, average deviation adev,
standard deviation sdev, variance var, skewness skew, and kurtosis curt.
Note: it seems that the definition of skewness and kurtosis differ between softwares. The values reported by Numerical Recipe differ from those reported by Microsoft Excel and by STATA (which does not agree with Excel either). Excel seems to agree with many online calculators
I decided to take the more commonly referred definition, where
skewness: b1=sum(xi-mean)^3/(N-1)*s^3 (for Normal Distribution, b1=0)
kurtosis: b2=sum(xi-mean)^4/(N-1)*s^4 (for Normal Distribution, b2=3)
excess kurtosis: b3=b2-3
*/
{
	void kcerror(char *error_text);
        int j;
        double ep=0.0,s,p;

        if (n <= 1) kcerror("kc::moment: n must be at least 2 in moment");
        s=0.0;
        for (j=0;j<n;j++) s += data[j];
        *ave=s/n;
        *adev=(*var)=(*skew)=(*curt)=0.0;
        for (j=0;j<n;j++) {
                *adev += fabs(s=data[j]-(*ave));
                ep += s;
                *var += (p=s*s);
                *skew += (p *= s);
                *curt += (p *= s);
        }
        *adev /= n;
        *var=(*var-ep*ep/n)/(n-1);			//corrected formula. normally ep=0
        *sdev=sqrt(*var);
        if (*var) {
                *skew /= ((n-1)*(*var)*(*sdev));
                *curt=(*curt)/((n-1)*(*var)*(*var));
        } else
        	kcerror("kc::moment: No skew/kurtosis when variance = 0 (in moment)");
}

double cc(double* x, double* y, int n)
{
	int j;
	double xt, yt, syy=0.0,sxy=0.0,sxx=0.0,ay=0.0,ax=0.0;

	for (j=0;j<n;j++) {
		ax += x[j];
		ay += y[j];
	}
	ax /= n;
	ay /= n;
	for (j=0;j<n;j++) {
		xt=x[j]-ax;
		yt=y[j]-ay;
		sxx += xt*xt;
		syy += yt*yt;
		sxy += xt*yt;
	}
	if (! sxx*syy) kcerror ("kc::cc: correlation cannot be calculated due to lack of variation in data");
	return sxy/sqrt(sxx*syy);
}

double cov(double* x, double* y, int n)
/* calculate the covariance of two variables */
{
	int i;
	double ax=0, ay=0, sxy=0;
	
	for (i=0; i<n; i++) {
		ax += x[i];
		ay += y[i];
	}
	ax /= n;
	ay /= n;
	for (i=0; i<n; i++) {
		sxy += (x[i]-ax)*(y[i]-ay);
	}
	return sxy/(double)(n-1);
}

/****************************************************************************
the following section contains subroutines for performing statistcal tests
****************************************************************************/

void STATISTICAL_TEST () {
	printf ("This section contains functions for performing statistical tests\n");
}

double bitest (int n, int k, double proportion)
/* this subroutine calculates one-sided binomial test to see whether the proportion is different from specified value
*/
/* similar to fisher's exact test, the word two-sided is generally difficult to define; in this case, the return value is only one-sided as a result
* for two-sided P-value, generally one can multiply by 2 and then if the result is higher than 1, assign 1 to the P-value instead.
*/
{
	double p=0;
	int i;
	double obs;
	
	if (n<=0) kcerror ("Error in kc::bitest: n should be an positive integer");
	if (n<k) kcerror ("Error in kc::bitest: n should be more than or equal to k");
	
	obs = (double) k / (double) n;
	for (i=0; i<=n; i++) {
		if (obs > proportion) {
			if (i>=k)
				p += bico(n, i)*pow(proportion, i)*pow(1-proportion, n-i);
		} else {
			if (i<=k)
				p += bico(n, i)*pow(proportion, i)*pow(1-proportion, n-i);
		}
	}
	if (p>1) p=1;
	return p;
}

void ttest_onesample (double *data, long n, double expected, double *t, double *p)
/* this subroutine perform one sample t-test to test whether mean is different from expected value
*/
{
	void avevar(double *data, int n, double *ave, double *var);
	double cdf_t (double df, double x);
	double ave, var;
	
	avevar (data, n, &ave, &var);
	*t = (ave - expected) / (sqrt(var) / sqrt(n));
	*p = 2 * (1 - cdf_t (n-1, fabs(*t)));
}
	

void ttest_ev(double *data1, long n1, double *data2, long n2, double *t, double *p)
/* this subroutine calculate 2-sided t-test assuming equal variance of two populations. the pooled variance estimate are used
*/
{
	void avevar(double *data, int n, double *ave, double *var);
	double cdf_t (double df, double x);
	double var1,var2,pvar,df,ave1,ave2;

	avevar(data1,n1,&ave1,&var1);
	avevar(data2,n2,&ave2,&var2);
	df=n1+n2-2;
	pvar=(((n1-1)*var1+(n2-1)*var2)/df) * (1.0/n1 + 1.0/n2);	//pooled variance
	*t=(ave1-ave2)/sqrt(pvar);
	*p = 2 * (1 - cdf_t (df, fabs(*t)));
}

void ttest_uev(double *data1, long n1, double *data2, long n2, double *t, double *p)
/* this subroutine calculates the 2-sided t-test assuming unequal variance of two populations
*/
{
	void avevar(double *data, int n, double *ave, double *var);
	double cdf_t (double df, double x);
	double var1,var2,df,ave1,ave2;
	
	avevar(data1,n1,&ave1,&var1);
	avevar(data2,n2,&ave2,&var2);
	df=(var1/n1+var2/n2)*(var1/n1+var2/n2)/((var1/n1)*(var1/n1)/(n1-1)+(var2/n2)*(var2/n2)/(n2-1));
	*t=(ave1-ave2)/sqrt(var1/n1+var2/n2);
	*p = 2 * (1 - cdf_t (df, fabs(*t)));
}

void ftest(double *data1, long n1, double *data2, long n2, double *f, double *p)
/* this subroutine test whether two "normally distributed" populations have equal variance
PLEASE NOTICE THAT THE CALCULATED f VALUE IS ALWAYS LARGER THAN 1 ! SO THE *f VALUE MAY NOT ACTUALLY BE THE STANDARD DEVIATION OF FIRST DATA DIVIDED BY SECOND DATA
*/
{
	void avevar(double *data, int n, double *ave, double *var);
	double cdf_f (double df1, double df2, double x);
	double var1,var2,ave1,ave2,df1,df2;
	
	avevar(data1,n1,&ave1,&var1);
	avevar(data2,n2,&ave2,&var2);
	if (var1 > var2) {
		*f=var1/var2;
		df1=n1-1;
		df2=n2-1;
	} else {
		*f=var2/var1;
		df1=n2-1;
		df2=n1-1;
	}
	*p = 2.0 * cdf_f (df1, df2, *f);
	if (*p > 1.0) *p=2.0-*p;
}

void chi2test(double *bins1, double *bins2, int nbins, int knstrn, double *df, double *chi2, double *p)
/* Given the arrays bins1[0..nbins-1] and bins2[0..nbins-1], containing two sets of binned data, perform chi-squared test.
bins1 and bins2 should be considered as observed data and expected data, respectively
the knstrn is dependent on whether there is a priori constraint on the equal number of total elements in two arrays. If yes, then knstrn=1; otherwise knstrn=0.
*/ 
{
	double cdf_chi2 (double df, double x);
	int j;

	*df=nbins-knstrn;
	*chi2=0.0;
	for (j=0;j<nbins;j++)
		if (bins1[j] == 0.0 && bins2[j] == 0.0)
			--(*df);
		else {
			*chi2 += (bins1[j]-bins2[j])*(bins1[j]-bins2[j]) / bins2[j];
		}
	*p = 1 - cdf_chi2 (*df, *chi2);
}

void chi2test_onesample(double *bins, double *ebins, int nbins, int knstrn, double *df, double *chi2, double *p)
/* Given the array bins[1..nbins] containing the observed numbers of events, and an array
ebins[1..nbins] containing the expected numbers of events, perform chi-squared test.
*/
{
	double cdf_chi2 (double df, double x);
	void kcerror(char *error_text);
	int j;

	*df=nbins-knstrn;
	*chi2=0.0;
	for (j=0;j<nbins;j++) {
		if (ebins[j] <= 0.0) kcerror("kc::chi2test_onesample: invalid expected number encountered in chi2test_onesample: cell is less than or equal to zero");
		*chi2 += (bins[j]-ebins[j])*(bins[j]-ebins[j]) / ebins[j];
	}
	*p = 1 - cdf_chi2 (*df, *chi2);
}

void chi2test_trend_2by3table (double *bins, double *chi2, double *p)
/* calculate Cochran-Armitage trend association test for 2x3 contingency table
the bins[] contains six elements, representing the 3 elements in the first row, followed by 3 elements in the second row
*/
{
	double row1, row2, col1, col2, col3, total;
	double e0, e1;	/*expected value in the first and second cell in contingency table*/
	double cdf_chi2 (double df, double x);
	
	row1 = bins[0]+bins[1]+bins[2];
	row2 = bins[3]+bins[4]+bins[5];
	col1 = bins[0]+bins[3];
	col2 = bins[1]+bins[4];
	col3 = bins[2]+bins[5];
	total = row1+row2;
	
	if (! (row1*row2*(total * (col2 + 4*col1) - (col2 + 2*col1)*(col2 + 2*col1)))) {
		*chi2 = -1;
		*p = -1;
		return;
	}
	
	e0 = row1*col1/total;
	e1 = row1*col2/total;
	
	*chi2 = (bins[1] - e1) + 2 * (bins[0] - e0);
	*chi2 *= *chi2;
	*chi2 *= total * total * (total-1) / row1 / row2 / (total * (col2 + 4*col1) - (col2 + 2*col1)*(col2 + 2*col1));
	*p = 1 - cdf_chi2 (1, *chi2);
}

void chi2test_trend_3by2table(double *bins, double *chi2, double *p)
{
	double *newbins = dvector (0, 5);
	newbins[0] = bins[0];
	newbins[1] = bins[2];
	newbins[2] = bins[4];
	newbins[3] = bins[1];
	newbins[4] = bins[3];
	newbins[5] = bins[5];
	chi2test_trend_2by3table (newbins, chi2, p);
	free_dvector (newbins, 0, 5);
}

void chi2test_2by2table(double *bins, double *chi2, double *p)
/* calculate chi2test for a 2x2 contingency table
*/
{
	double row1, row2, col1, col2, total;
	double cdf_chi2 (double df, double x);
	
	row1 = bins[0]+bins[1];
	row2 = bins[2]+bins[3];
	col1 = bins[0]+bins[2];
	col2 = bins[1]+bins[3];
	total = row1+row2;
	
	if (! (row1*row2*col1*col2)) {
		*chi2 = -1;
		*p = -1;
		return;
	}

/*	This method (calculating expected value for each cell before calculating chi2 value) is too complicated so I use a shortcut formula instead
	void chi2test_onesample(double bins[], double ebins[], int nbins, int knstrn, double *df, double *chi2, double *p);
	double *ebins; 	
 	ebins = dvector (0, 3);
	ebins[0] = row1*col1/total;
	ebins[1] = row1*col2/total;
	ebins[2] = row2*col1/total;
	ebins[3] = row2*col2/total;
	
	chi2test_onesample (bins, ebins, 4, 3, &df, chi2, p);
	free_dvector (ebins, 0, 3);
*/

	*chi2 = (bins[0]*bins[3]-bins[1]*bins[2]);
	*chi2 *= *chi2;
	*chi2 *= total / row1 / row2 / col1 / col2;
	*p = 1 - cdf_chi2 (1, *chi2);
}

void chi2test_3by2table(double *bins, double *chi2, double *p)
{
	void chi2test_2by3table(double *bins, double *chi2, double *p);
	double *newbins = dvector (0, 5);
	newbins[0] = bins[0];
	newbins[1] = bins[2];
	newbins[2] = bins[4];
	newbins[3] = bins[1];
	newbins[4] = bins[3];
	newbins[5] = bins[5];
	chi2test_2by3table (newbins, chi2, p);
	free_dvector (newbins, 0, 5);
}

void chi2test_2by3table(double *bins, double *chi2, double *p)
/* calculate chi2test for a 2x3 contingency table (for example, the table used in genotypeic association test, or commonly referred to as "2df test")
this table has 2 rows and 3 columns, the input in bins[] array are a11, a12, a13, a21, a22, a23, respectively
*/
{
	void chi2test_onesample(double *bins, double *ebins, int nbins, int knstrn, double *df, double *chi2, double *p);
	double *ebins;
	double row1, row2, col1, col2, col3, total, df;
	
	row1 = bins[0]+bins[1]+bins[2];
	row2 = bins[3]+bins[4]+bins[5];
	col1 = bins[0]+bins[3];
	col2 = bins[1]+bins[4];
	col3 = bins[2]+bins[5];
	total = row1+row2;

	if (! (row1*row2*col1*col2*col3)) {
		*chi2 = -1;
		*p = -1;
		return;
	}
	
	ebins = dvector (0, 5);
	ebins[0] = row1*col1/total;
	ebins[1] = row1*col2/total;
	ebins[2] = row1*col3/total;
	ebins[3] = row2*col1/total;
	ebins[4] = row2*col2/total;
	ebins[5] = row2*col3/total;
	
	chi2test_onesample (bins, ebins, 6, 4, &df, chi2, p);
	free_dvector (ebins, 0, 5);
}

void kstest(double *data1, long n1, double *data2, long n2, double *d, double *prob)
/* perform the Kolmogorov-Smirnov test to examine whether two data sets are drawn from the same distribution
the code is adapted from numerical recipe, but it does not handle ties well. I may change it in the future
*/
{
	double probks(double alam);
	void quick_sort(int elements, double *arr);
	long j1=0,j2=0;
	double d1,d2,dt,en1,en2,en,fn1=0.0,fn2=0.0;

	quick_sort(n1,data1);
	quick_sort(n2,data2);
	en1=n1;
	en2=n2;
	*d=0.0;
	while (j1 < n1 && j2 < n2) {
		if ((d1=data1[j1]) <= (d2=data2[j2])) fn1=j1++/en1;
		if (d2 <= d1) fn2=j2++/en2;
		if ((dt=fabs(fn2-fn1)) > *d) *d=dt;
	}
	en=sqrt(en1*en2/(en1+en2));
	*prob=probks((en+0.12+0.11/en)*(*d));
}

void kstest_onesample(double *data, long n, double (*func)(double), double *d, double *prob)
/* perform the Kolmogorov-Smirnov test to examine whether two data sets are drawn from the same distribution
*/
{
	double probks(double alam);
	void quick_sort(int elements, double *arr);
	long j;
	double dt,en,ff,fn,fo=0.0;

	quick_sort(n,data);
	en=n;
	*d=0.0;
	for (j=0;j<n;j++) {
		fn=j/en;
		ff=(*func)(data[j]);
		dt=fabs(fo-ff);
		if (dt<fabs(fn-ff)) dt=fabs(fn-ff);
		/*dt=FMAX(fabs(fo-ff),fabs(fn-ff));*/
		if (dt > *d) *d=dt;
		fo=fn;
	}
	en=sqrt(en);
	*prob=probks((en+0.12+0.11/en)*(*d));
}

#define EPS1 0.001
#define EPS2 1.0e-8
double probks(double alam)
{
	int j;
	double a2,fac=2.0,sum=0.0,term,termbf=0.0;

	a2 = -2.0*alam*alam;
	for (j=1;j<=1000;j++) {
		term=fac*exp(a2*j*j);
		sum += term;
		if (fabs(term) <= EPS1*termbf || fabs(term) <= EPS2*sum) return sum;
		fac = -fac;
		termbf=fabs(term);
	}
	return 1.0;
}
#undef EPS1
#undef EPS2

double fisher_exact_1sided (int a, int b, int c, int d)
/* calculate 1-sided Fisher's exact test
I used lnfactorial, rather than relying on the more familiar hypergeometric distribution, since lnfactorial is more direct to compute
*/
{
	double lnfactorial (double x);
	int lownum=a;
	int temp, i, newa, newb, newc, newd;
	double prob = 0, firstprob, currentprob;
	
	if (a<b && a<c && a<d) {
		lownum=a;
	} else if (b<c && b<d) {
		temp=b; b=a; a=temp;
		temp=d; d=c; c=temp;
	} else if (c<d) {
		temp=c; c=a; a=temp;
		temp=d; d=b; b=temp;
	} else {
		temp=d; d=a; a=temp;
		temp=c; c=b; b=temp;
	}

	firstprob = lnfactorial(a+b)+lnfactorial(a+c)+lnfactorial(b+d)+lnfactorial(c+d)-lnfactorial(a+b+c+d);

	for (i=0; i<a; i++) {
		newa = i;
		newb = a+b-i;
		newc = a+c-i;
		newd = -a+d+i;
		currentprob = exp(firstprob - lnfactorial(newa) - lnfactorial(newb) - lnfactorial(newc) - lnfactorial(newd));
		prob += currentprob;
	}

	if (a*d < b*c) {			//the a cell is under-represented; calculate the sum of prob where the cell range from 0 to a
		prob += exp(firstprob - lnfactorial(a) - lnfactorial(b) - lnfactorial(c) - lnfactorial(d));
		if (prob > 1) prob = 1; if (prob < 0) prob = 0;
		return prob;
	} else {
		if (prob > 1) prob = 1; if (prob < 0) prob = 0;
		return 1-prob;
	}
}

double fisher_exact_2sided (int a, int b, int c, int d)
/* calculate 2-sided Fisher's exact test. I first re-arrange the contingency table and make sure that the first cell is the smallest cell
next determine whether the first cell is under-represented or over-represented and then decide the direction of iteration
I used lnfactorial, rather than relying on the more familiar hypergeometric distribution, since lnfactorial is more direct to compute
*/
{
	double lnfactorial (double x);
	int lownum=a;
	int temp, i, newa, newb, newc, newd, maxa;
	double prob = 0, firstprob, currentprob, tableprob;
	
	if (a<b && a<c && a<d) {
		lownum=a;
	} else if (b<c && b<d) {
		temp=b; b=a; a=temp;
		temp=d; d=c; c=temp;
	} else if (c<d) {
		temp=c; c=a; a=temp;
		temp=d; d=b; b=temp;
	} else {
		temp=d; d=a; a=temp;
		temp=c; c=b; b=temp;
	}

	firstprob = lnfactorial(a+b)+lnfactorial(a+c)+lnfactorial(b+d)+lnfactorial(c+d)-lnfactorial(a+b+c+d);
	tableprob = exp(firstprob - lnfactorial(a) - lnfactorial(b) - lnfactorial(c) - lnfactorial(d));			//prob of the original table

	maxa = a+b;
	if (a+c<maxa) maxa=a+c;
	if (a*d < b*c) {			//the a cell is under-represented; calculate the sum of prob where the cell range from a to large numbers until prob is less than tableprob
		for (i=a+1; i<=maxa; i++) {
			newa = i;
			newb = a+b-i;
			newc = a+c-i;
			newd = -a+d+i;
			currentprob = exp(firstprob - lnfactorial(newa) - lnfactorial(newb) - lnfactorial(newc) - lnfactorial(newd));
			if (currentprob < tableprob) break;
			prob += currentprob;
		}
	} else {				//the a cell is over-represented; calculate the sum of prob where the cell range from a to smaller number until prob is less than tableprob
		for (i=a-1; i>=0; i--) {
			newa = i;
			newb = a+b-i;
			newc = a+c-i;
			newd = -a+d+i;
			currentprob = exp(firstprob - lnfactorial(newa) - lnfactorial(newb) - lnfactorial(newc) - lnfactorial(newd));
			if (currentprob < tableprob) break;
			prob += currentprob;
		}
	}
	if (prob > 1) prob = 1; if (prob < 0) prob = 0;
	return 1-prob;
}


void quick_sort(int elements, double *arr) 
/* use quick sort algorithm to sort an array of double numbers
code adapted from http://alienryderflex.com/quicksort/
*/
{
#define  MAX_LEVELS  300
	double piv;
	int  beg[MAX_LEVELS], end[MAX_LEVELS], i=0, L, R, swap ;
	
	beg[0]=0; end[0]=elements;
	while (i>=0) {
		L=beg[i]; R=end[i]-1;
		if (L<R) {
			piv=arr[L];
			while (L<R) {
				while (arr[R]>=piv && L<R) R--; if (L<R) arr[L++]=arr[R];
				while (arr[L]<=piv && L<R) L++; if (L<R) arr[R--]=arr[L]; 
			}
			arr[L]=piv; beg[i+1]=L+1; end[i+1]=end[i]; end[i++]=L;
			if (end[i]-beg[i]>end[i-1]-beg[i-1]) {
				swap=beg[i]; beg[i]=beg[i-1]; beg[i-1]=swap;
				swap=end[i]; end[i]=end[i-1]; end[i-1]=swap; 
			}
		} else {
			i--; 
		}
	}
#undef MAX_LEVELS
}


/*******************************************************************************
the following section contains mathematical functions
********************************************************************************/

void MATH_FUNCTION () {
	printf ("This section contains mathematical functions\n");
}

double gammln(double x)
/*calculate the natural logarithm of the Gamma function (http://en.wikipedia.org/wiki/Gamma_function)

The gamma function is defined by
                +inf
                /   - t  (z - 1)
     gamma(z) = |  e    t       dt
                /
                0

If z is a positive integer, then

    Gamma(z)=(z-1)!,

this subroutine was modified from jdhedden's formula that I found at http://www.nr.com/forum/showthread.php?t=606
*/
{
    double tmp, ser;
    tmp = x + 4.5 - (x - 0.5) * log(x + 4.5);
    ser = 1.000000000190015 + (76.18009172947146 / x) - (86.50532032941677 / (x + 1.0)) + (24.01409824083091 / (x + 2.0)) - (1.231739572450155 / (x + 3.0)) + (0.1208650973866179e-2 / (x + 4.0)) - (0.5395239384953e-5 / (x + 5.0));
    return (log(2.5066282746310005 * ser) - tmp);
}

double lnfactorial (double x)
/*calculate the natural logarithm of the factorial of x, this is equivalent to gammln(x+1)
*/
{
	return gammln (x+1.0);
}

int factorial (int x)
{
	int i;
	int total=1;
	if (x==0) return 1;
	if (x>12) return exp(lnfactorial((double) x));
	for (i=1; i<=x; i++) {
		total *= i;
	}
	return total;
}

double gammp(double a, double x)
/* calculate incomplete gamma function (http://en.wikipedia.org/wiki/Incomplete_gamma_function)

                             x
                      1     /   - t  (a - 1)
gammainc(a, x) = ---------- |  e    t       dt
                  gamma(a)  /
                             0

this subroutine is modified from Numerical Recipe
there are two ways to calculate incomplete gamma function: by series representation (gsea) and by continued fraction representation (gcf)
when x<a+1, gsea converges faster; otherwise gcf converges faster
	
*/
{
        void gcf(double *gammcf, double a, double x, double *gln);
        void gser(double *gamser, double a, double x, double *gln);
        void kcerror(char *error_text);
        double gamser,gammcf,gln;

        if (x < 0.0 || a <= 0.0) kcerror("kc::gammp: invalid x or a value");
        if (x < (a+1.0)) {
                gser(&gamser,a,x,&gln);
                return gamser;
        } else {
                gcf(&gammcf,a,x,&gln);
                return 1.0-gammcf;
        }
}

double gammq(double a, double x)
{
	double gammp(double a, double x);
	return 1.0 - gammp(a, x);
}

void gser(double *gamser, double a, double x, double *gln)
/* the incomplete gamma function P(a, x) evaluated by its series representation as gamser.
*/
{
        double gammln(double x);
        void kcerror(char *error_text);
        int n;
        double sum,del,ap;

        *gln=gammln(a);
        if (x <= 0.0) {
                if (x < 0.0) kcerror("kc::gser: x must be greater than or equal to zero");
                *gamser=0.0;
                return;
        } else {
                ap=a;
                del=sum=1.0/a;
                for (n=1;n<=ITMAX;n++) {
                        ++ap;
                        del *= x/ap;
                        sum += del;
                        if (fabs(del) < fabs(sum)*EPS) {
                                *gamser=sum*exp(-x+a*log(x)-(*gln));
                                return;
                        }
                }
                kcerror("kc::gser: a too large or ITMAX too small");
                return;
        }
}

void gcf(double *gammcf, double a, double x, double *gln)
/* calculate the incomplete gamma function Q(a, x) evaluated by its continued fraction representation as gammcf.

*/
{
        double gammln(double xx);
        int i;
        double an,b,c,d,del,h;

        *gln=gammln(a);
        b=x+1.0-a;
        c=1.0/FPMIN;
        d=1.0/b;
        h=d;
        for (i=1;i<=ITMAX;i++) {
                an = -i*(i-a);
                b += 2.0;
                d=an*d+b;
                if (fabs(d) < FPMIN) d=FPMIN;
                c=b+an/c;
                if (fabs(c) < FPMIN) c=FPMIN;
                d=1.0/d;
                del=d*c;
                h *= del;
                if (fabs(del-1.0) < EPS) break;
        }
        if (i > ITMAX) kcerror("kc::gcf: a too large or ITMAX too small");
        *gammcf=exp(-x+a*log(x)-(*gln))*h;
}


double beta(double z, double w)
/* calculate the beta function (http://en.wikipedia.org/wiki/Beta_function)
The beta function is defined by
                   1
                  /   (z - 1)       (w - 1)
     beta(z, w) = |  t       (1 - t)       dt
                  /
                   0

z and w are non-negative shape parameters
*/
{
        double gammln(double x);
        return exp(gammln(z)+gammln(w)-gammln(z+w));
}

double betai(double a, double b, double x)
/* calculate the incomplete beta function 
                             x
                      1     /   (a - 1)       (b - 1)
betai(a, b, x) = ---------- |  t       (1 - t)       dt
                 beta(a, b) /
                             0
*/
{
        double betacf(double a, double b, double x);
        double gammln(double xx);
        void kcerror(char *error_text);
        double bt;

        if (x < 0.0 || x > 1.0) kcerror("kc::betai: Invalid x value (x value must be bewteen 0 and 1)");
        if (x == 0.0 || x == 1.0) bt=0.0;
        else
                bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
        if (x < (a+1.0)/(a+b+2.0))
                return bt*betacf(a,b,x)/a;
        else
                return 1.0-bt*betacf(b,a,1.0-x)/b;
}

double betacf(double a, double b, double x)
/* continued fraction for incomplete beta function by modified Lentz’s method
*/
{
        int m,m2;
        double aa,c,d,del,h,qab,qam,qap;

        qab=a+b;
        qap=a+1.0;
        qam=a-1.0;
        c=1.0;
        d=1.0-qab*x/qap;
        if (fabs(d) < FPMIN) d=FPMIN;
        d=1.0/d;
        h=d;
        for (m=1;m<=ITMAX;m++) {
                m2=2*m;
                aa=m*(b-m)*x/((qam+m2)*(a+m2));
                d=1.0+aa*d;
                if (fabs(d) < FPMIN) d=FPMIN;
                c=1.0+aa/c;
                if (fabs(c) < FPMIN) c=FPMIN;
                d=1.0/d;
                h *= d*c;
                aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
                d=1.0+aa*d;
                if (fabs(d) < FPMIN) d=FPMIN;
                c=1.0+aa/c;
                if (fabs(c) < FPMIN) c=FPMIN;
                d=1.0/d;
                del=d*c;
                h *= del;
                if (fabs(del-1.0) < EPS) break;
        }
        if (m > ITMAX) kcerror("kc::betacf: a or b too big (no convergene after ITMAX iteration)");
        return h;
}



double inverff (double x)
/*inverse error function (see http://en.wikipedia.org/wiki/Error_function)
the formula is adapted from http://www.theorie.physik.uni-muenchen.de/~serge/erf-approx.pdf
*/
{
	double y;
	double a = 8887.0/63473.0;
	y = -2.0/PI/a - log (1-x*x)/2.0 + sqrt ((2.0/PI/a + log(1-x*x)/2.0)*(2.0/PI/a + log(1-x*x)/2.0) - log(1-x*x)/a);
	y = sqrt (y);
	if (x<0) y = -y;
	return y;
}

double erf(double x)
/* calculate the error function (http://en.wikipedia.org/wiki/Error_function)
The error function is defined by
                        x     2
                  2    /    -t
     erf(x) = -------- |  e    dt
              sqrt(pi) /
                        0
*/
{
        double gammp(double a, double x);
        return (x < 0.0) ? (-gammp(0.5,x*x)) : gammp(0.5,x*x);
}




/*******************************************************************************
random number generators

Several functions for generating random numbers are given.

Besides uniform random numbers, we can generate random numbers from any given 
distribution by inverse transform sampling, which generates random number from 
any probability distribution given its cumulative distribution function (cdf).

*******************************************************************************/

void RANDOM_NUMBER ()
{
	printf ("Random number generators\n");
}

double ran ()
{
	return (double) rand() / (RAND_MAX+1.0);
}


#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

double ran3(int *idum)
/* random number generator that generate [0,1] uniform distribution
According to Knuth, any large MBIG, and any smaller (but still large) MSEED can be substituted
for the above values.
*/
{
	static int inext,inextp;
	static int ma[56];
	static int iff=0;
	int mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC



double random_exp(double lamda, int *idum)
/* random number generator for exponential distribution, with mean of lamda
*/
{
	double ran3(int *idum);
	double dum;

	do
		dum=ran3(idum);
	while (dum == 0.0);
	return -log(dum) / lamda;
}

double random_stdnormal (int *idum)
/* random number generator for standard normal distribution, with mean of zero and variance of 1
*/
{
	double ran3(int *idum);
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*ran3(idum)-1.0;
			v2=2.0*ran3(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}






/*********************************************************************************
the following section contains third-party subroutines with slight modifications
*********************************************************************************/


/*
// This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
// Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of 
// Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000  
//
// Written by Jan Wigginton
*/
/*
// This code has been modified by Kai to conform to ISO C90 standard to prevent warning messages during compilation by modern compilers
*/
double SNPHWE(int obs_hets, int obs_hom1, int obs_hom2)
{
	int obs_homc, obs_homr, rare_copies, genotypes, i, mid, curr_hets, curr_homr, curr_homc;
	double *het_probs;
	double sum, p_hwe;
        void kcerror(char *error_text);

	if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0) {
		fprintf (stderr, "NOTICE: the three genotype counts are %d  %d %d\n", obs_hets, obs_hom1, obs_hom2);
		kcerror ("kc::SNPHWE: Current genotype configuration includes a negative count");
	}
	
	obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
	obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;
	
	rare_copies = 2 * obs_homr + obs_hets;
	genotypes   = obs_hets + obs_homc + obs_homr;
	
	het_probs = (double *) malloc((size_t) (rare_copies + 1) * sizeof(double));
	if (het_probs == NULL) {
		kcerror ("kcerror::SNPHWE: Unable to allocate array for heterozygote probabilities" );
	}
   

	for (i = 0; i <= rare_copies; i++)
		het_probs[i] = 0.0;

	/* start at midpoint */
	mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);
	
	/* check to ensure that midpoint and rare alleles have same parity */
	if ((rare_copies & 1) ^ (mid & 1))
		mid++;

	curr_hets = mid;
	curr_homr = (rare_copies - mid) / 2;
	curr_homc = genotypes - curr_hets - curr_homr;

	het_probs[mid] = 1.0;
	sum = het_probs[mid];
	for (curr_hets = mid; curr_hets > 1; curr_hets -= 2) {
		het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0)
	                       / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
		sum += het_probs[curr_hets - 2];
	
		/* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
		curr_homr++;
		curr_homc++;
	}

	curr_hets = mid;
	curr_homr = (rare_copies - mid) / 2;
	curr_homc = genotypes - curr_hets - curr_homr;
	for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2) {
		het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc
		                    /((curr_hets + 2.0) * (curr_hets + 1.0));
		sum += het_probs[curr_hets + 2];
		
		/* add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote */
		curr_homr--;
		curr_homc--;
	}

	for (i = 0; i <= rare_copies; i++)
		het_probs[i] /= sum;

	p_hwe = 0.0;
	/*  p-value calculation for p_hwe  */
	for (i = 0; i <= rare_copies; i++) {
		if (het_probs[i] > het_probs[obs_hets])
			continue;
		p_hwe += het_probs[i];
	}
	
	p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;
	
	free(het_probs);
	
	return p_hwe;
}
