#ifndef _ALIOUTSW_H_
#define _ALIOUTSW_H_

#include "aliCalcSW.h"
#include "aliCost.h"

/* Find highest scoring local alignment(s) in mat, and print to stdout
   the corresponding best alignments.
   mat must have been filled with scores and prevs.
   bestScore is the highest score in mat (provided for convenience).
   cost is provided so mismatches with negative scores can be lowercased.
*/
void aliPrintBestAlis(struct matrix *mat, double bestScore, struct cost *cost, char *s1, char *s2) ;

#endif
