#include <stdio.h>
#include <stdlib.h>

#include "aliCost.h"
#include "aliGetSeq.h"
#include "aliCalcAE.h"
#include "aliOutAE.h"

int main(void)
{
	char *s1 ;
	while((s1 = getSeq(1)) == NULL) {
		// nothing to do
	}
	char *s2 ;
	while((s2 = getSeq(1)) == NULL) {
	}
	printf("Sequences read:\ns1\t%s\ns2\t%s\n\n", s1, s2) ;

	/* BLOSUM62 prot subst cost with affine cost for short indels */
	struct cost *cost = costProt(-10,-0.5);
	struct matrix *mat = aliInitMat(s1,s2);
	double bestScore = aliFillMat(mat,cost,s1,s2);
	/* for debugging you can uncomment:
	   aliPrintMat(mat); */
	aliPrintBestAlis(mat,bestScore,cost,s1,s2);

	aliFreeMat(mat);
	free(cost);
	free(s1);
	free(s2);
	return(0);
}
