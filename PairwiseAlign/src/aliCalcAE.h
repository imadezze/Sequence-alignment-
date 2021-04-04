//
// Created by imad on 17/03/2021.
//

#ifndef PAIRWISEALIGN_ALICALCAE_H
#define PAIRWISEALIGN_ALICALCAE_H



#include <stdint.h>
#include "aliCost.h"


struct cell {
    double scoreV;
    double scoreD;
    double scoreH;
    uint16_t prevs;
};

struct matrix {
    unsigned int w;
    unsigned int h;
    struct cell *cells; /* pointer to array of w*h cells
			       cells[w*i+j] contains cell (i,j) */
};


/* allocate and initialize (first row and col) a matrix for
   AE alignment of strings s1 and s2 */
struct matrix *aliInitMat(char *s1, char *s2);


/* Fill the mat matrix, using SW with a linear indel model
   using cost->indelOpen, or AE with an affine indel model using
   cost->indelOpen and cost->indelExtend.
   Return the best score found.
   Preconditions:
   - mat is correctly allocated and initialized (by aliInitMat)
   - cost->subst is defined for each pair of letters in s1 and s2
*/
double aliFillMat(struct matrix *mat, struct cost *cost, char *s1, char *s2) ;


/* free all allocated memory in mat */
void aliFreeMat(struct matrix *mat);

/* print contents of matrix, for debugging */
void aliPrintMat(struct matrix *mat);

/* fills the three scores */
void fillPrevs(struct matrix *mat, struct cost *cost, char *s1, char *s2, int i, int j);

#endif //PAIRWISEALIGN_ALICALCAE_H