//
// Created by imad on 17/03/2021.
//

#include "aliCalcAE.h"
#include <string.h>
#include "mem.h"
#include <math.h>
#include <stdio.h>

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
        cells[w*i].scoreV = 0;
        cells[w*i].scoreD = 0;
        cells[w*i].scoreH = 0;
        cells[w*i].prevs = 0;
    }
    for (int j = 0; j < w; j++){
        cells[j].scoreV = 0;
        cells[j].scoreD = 0;
        cells[j].scoreH = 0;
        cells[j].prevs = 0;
    }
    return matrice;
}

double aliFillMat(struct matrix *mat, struct cost *cost, char *s1, char *s2){
    unsigned int w = mat->w;
    unsigned int h = mat->h;
    struct cell *cells = mat->cells;
    double scoremax = 0;
    for (int i = 1; i < h; i++){
        for (int j = 1; j < w; j++) {
            double maxD = fmax(fmax(cells[(i-1)*w+j-1].scoreD + cost->subst(s2[j-1], s1[i-1]), 0),
                              fmax(cells[(i-1)*w+j-1].scoreH + cost->subst(s2[j-1], s1[i-1]),cells[(i-1)*w+j-1].scoreV + cost->subst(s2[j-1], s1[i-1])));

            double maxV = fmax(fmax(cells[(i-1)*w+j].scoreD + cost->indelOpen, cells[(i-1)*w+j].scoreV + cost->indelExtend),
                               fmax(cells[(i-1)*w+j].scoreH + cost->indelOpen,0));

            double maxH = fmax(fmax(cells[i*w+j-1].scoreD + cost->indelOpen, cells[i*w+j-1].scoreV + cost->indelOpen),
                               fmax(cells[i*w+j-1].scoreH + cost->indelExtend,0));

            cells[w*i+j].scoreD = maxD;
            cells[w*i+j].scoreH = maxH;
            cells[w*i+j].scoreV = maxV;

            fillPrevs(mat, cost, s1, s2, i, j);

            scoremax = fmax(fmax(scoremax, maxD), fmax(maxH, maxV));
        }
    }
    return scoremax;
}

/* fills the three scores */
void fillPrevs(struct matrix *mat, struct cost *cost, char *s1, char *s2, int i, int j) {
    unsigned int w = mat->w;
    struct cell *cells = mat->cells;
    double maxD = cells[w*i+j].scoreD;
    double maxV = cells[w*i+j].scoreV;
    double maxH = cells[w*i+j].scoreH;

    /* prevD */
    if(maxD == cells[(i-1)*w+j-1].scoreD + cost->subst(s2[j-1], s1[i-1])) {
        cells[w*i+j].prevs |= 1;
    }
    if(maxD == cells[(i-1)*w+j-1].scoreH + cost->subst(s2[j-1], s1[i-1])) {
        cells[w*i+j].prevs |= 2;
    }
    if(maxD == cells[(i-1)*w+j-1].scoreV + cost->subst(s2[j-1], s1[i-1])) {
        cells[w*i+j].prevs |= 4;
    }
    if(maxD == 0) {
        cells[w*i+j].prevs &= ~(0b111);
    }

    /* prevH */
    if(maxH == cells[i*w+j-1].scoreD + cost->indelOpen) {
        cells[w*i+j].prevs |= 8;
    }
    if(maxH ==cells[i*w+j-1].scoreH + cost->indelExtend) {
        cells[w*i+j].prevs |= 16;
    }
    if(maxH == cells[i*w+j-1].scoreV + cost->indelOpen) {
        cells[w*i+j].prevs |= 32;
    }
    if(maxH == 0) {
        cells[w*i+j].prevs &= ~(0b111000);
    }

    /* prevV */
    if(maxV == cells[(i-1)*w+j].scoreD + cost->indelOpen) {
        cells[w*i+j].prevs |= 64;
    }
    if(maxV == cells[(i-1)*w+j].scoreH + cost->indelOpen) {
        cells[w*i+j].prevs |= 128;
    }
    if(maxV == cells[(i-1)*w+j].scoreV + cost->indelExtend) {
        cells[w*i+j].prevs |= 256;
    }
    if(maxV == 0) {
        cells[w*i+j].prevs &= ~(0b111000000);
    }

}

/* print contents of matrix, for debugging */
void aliPrintMat(struct matrix *mat) {
    unsigned int w = mat->w;
    unsigned int h = mat->h;
    struct cell *cells = mat->cells;
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            printf("%0.f / ", cells[w*i+j].scoreD);
            printf("%0.f / ", cells[w*i+j].scoreH);
            printf("%0.f   ", cells[w*i+j].scoreV);
        }
        printf("\n");
    }
}

/* free all allocated memory in mat */
void aliFreeMat(struct matrix *mat){
    free(mat->cells);
    free(mat);
}