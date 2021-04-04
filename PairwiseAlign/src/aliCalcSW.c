//
// Created by ensimag on 03/03/21.
//

#include "aliCalcSW.h"
#include <string.h>
#include "mem.h"
#include <math.h>
#include <stdio.h>
#include <stdbool.h>


/* allocate and initialize (first row and col) a matrix for SW
   or AE alignment of strings s1 and s2 */
struct matrix *aliInitMat(char *s1, char *s2){
    unsigned int w = strlen(s2)+1;
    unsigned int h = strlen(s1)+1;
    char errmsg[20] = "allocation failed";
    struct cell *cells = mallocOrDie(sizeof(struct cell)*w*h, errmsg);
    struct matrix *matrice = mallocOrDie(sizeof(struct matrix), errmsg);
    matrice->w = w;
    matrice->h = h;
    matrice->cells = cells;
    for (int i = 0; i < h; i++){
        cells[w*i].score = 0;
        cells[w*i].prevs = 0;
    }
    for (int j = 0; j < w; j++){
        cells[j].score = 0;
        cells[j].prevs = 0;
    }
    return matrice;
}

/* Fill the mat matrix, using SW with a linear indel model
   using cost->indelOpen, or AE with an affine indel model using
   cost->indelOpen and cost->indelExtend.
   Return the best score found.
   Preconditions:
   - mat is correctly allocated and initialized (by aliInitMat)
   - cost->subst is defined for each pair of letters in s1 and s2
*/
double aliFillMat(struct matrix *mat, struct cost *cost, char *s1, char *s2){
    unsigned int w = mat->w;
    unsigned int h = mat->h;
    struct cell *cells = mat->cells;
    double scoremax = 0;
    bool wasZero = false;
    for (int i = 1; i < h; i++){
       for (int j = 1; j < w; j++) {
           double max = fmax(fmax(cells[(i-1)*w+j-1].score + cost->subst(s2[j-1], s1[i-1]), 0),
                             fmax(cells[(i-1)*w+j].score + cost->indelOpen,cells[i*w+j-1].score + cost->indelOpen));
           cells[w*i+j].score = max;
           if (max == cells[(i-1)*w+j-1].score + cost->subst(s2[j-1], s1[i-1])){
                cells[w*i+j].prevs |= 1;
                wasZero = true;
           }
           if (max == cells[(i-1)*w+j].score + cost->indelOpen){
               cells[w*i+j].prevs |= 4;
               wasZero = true;
           }
           if (max == cells[i*w+j-1].score + cost->indelOpen){
               cells[w*i+j].prevs |= 2;
               wasZero = true;
           }
           if ((!wasZero) && (max == 0)){
               cells[w*i+j].prevs = 0;
           }
           scoremax = fmax(max, scoremax);
       }
   }
    return scoremax;
}
/* print contents of matrix, for debugging */
void aliPrintMat(struct matrix *mat){
    unsigned int w = mat->w;
    unsigned int h = mat->h;
    struct cell *cells = mat->cells;
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            printf("%0.f   ", cells[w*i+j].score);
        }
        printf("\n");
    }
}

/* free all allocated memory in mat */
void aliFreeMat(struct matrix *mat){
    free(mat->cells);
    free(mat);
}






