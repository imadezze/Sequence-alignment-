#ifndef _ALICALCSW_H_
#define _ALICALCSW_H_

#include <stdint.h>
#include "aliCost.h"

/* defines the matrix datatype and methods for Smith-Waterman 
   sequence alignment.
   In order to implement Altschul-Erickson the cell and matrix 
   types will need to be modified, but all function prototypes
   should remain the same. */


/* in struct cell, the three last bits of prevs are used to 
   say whether best paths come from Top, Left, or Diag (in that order).
   So eg if (prevs&4) then a best path comes from the
   cell above, and if (prevs&1) a best path comes from diag.
   The three bits are non-exclusive (can have multiple best paths).
   If prevs==0 then score must be 0.
*/
struct cell {
	double score;
	uint8_t prevs;
};

struct matrix {
	unsigned int w;
	unsigned int h;
	struct cell *cells; /* pointer to array of w*h cells
			       cells[w*i+j] contains cell (i,j) */
};


/* allocate and initialize (first row and col) a matrix for SW 
   or AE alignment of strings s1 and s2 */
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

#endif
