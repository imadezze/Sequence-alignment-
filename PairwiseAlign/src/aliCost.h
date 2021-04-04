#ifndef _ALICOST_H_
#define _ALICOST_H_


struct cost {
	/* for SW (linear cost), use indelOpen and ignore indelExtend;
	   for AE both indelOpen and indelExtend are needed */
	double indelOpen;
	double indelExtend;
	/* subst: pointer to function returning the substitution score of
	   its args. 
	   Of course the function will be different for DNA and proteins.
	   Precondition: args are both uppercase letters */
	double (*subst)(char,char);
};


/* return a cost structure for DNA alignment.
   This cost structure is designed for an affine gap cost model (AE),
   but implementations of the linear model (SW) should also use
   it and simply ignore indelExtend.
   Currently the subst function uses match==+5 and mismatch==-4, this
   can be easily customized by editing the static dnaSubst function.
   Requires indelOpen <= indelExtend < 0.
*/
struct cost *costDna(double indelOpen, double indelExtend);


/* return a cost structure for protein alignment.
   Again the struct is designed for AE but should also be used for SW.
   The returned structure currently implements a BLOSUM62 substitution 
   cost.
   Requires indelOpen <= indelExtend < 0.
*/
struct cost *costProt(double indelOpen, double indelExtend);

#endif
