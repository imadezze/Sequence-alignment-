#include <stdio.h>
#include <stdlib.h>

#include "aliCost.h"
#include "aliGetSeq.h"
#include "aliCalcSW.h"
#include "aliOutSW.h"

int main(void)
{
	char *s1 ;
	while((s1 = getSeq(0)) == NULL) {
		// nothing to do
	}
	char *s2 ;
	while((s2 = getSeq(0)) == NULL) {
	}
	printf("Sequences read:\ns1\t%s\ns2\t%s\n\n", s1, s2) ;

	/* affine cost for long gaps eg spliced RNA on genome */
	struct cost *cost = costDna(-100,-0.05);
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
