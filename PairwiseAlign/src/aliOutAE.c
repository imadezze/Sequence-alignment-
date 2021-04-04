//
// Created by imad on 17/03/2021.
//

#include "aliOutAE.h"
#include "aliCost.h"
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>


void aliPrintBestAlis(struct matrix *mat, double bestScore, struct cost *cost, char *s1, char *s2){
    unsigned int w = mat->w;
    unsigned int h = mat->h;
    struct cell *cells = mat->cells;
    unsigned int positionbestscorelocal = 0;
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            double maxScore = fmax(fmax(cells[w*i+j].scoreD, cells[w*i+j].scoreH), cells[w*i+j].scoreV);
            if(maxScore == bestScore){
                positionbestscorelocal = w*i+j;
                aliPrintBestAlisprime(mat,bestScore, cost, s1, s2, positionbestscorelocal);
            }
        }
    }
}

unsigned int findBestScore(struct matrix *mat, double bestScore){
    unsigned int w = mat->w;
    unsigned int h = mat->h;
    struct cell *cells = mat->cells;
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            double maxScore = fmax(fmax(cells[w*i+j].scoreD, cells[w*i+j].scoreH), cells[w*i+j].scoreV);
            if(maxScore == bestScore){
                return w*i+j;
            }
        }
    }
}

/* Find highest scoring local alignment(s) in mat, and print to stdout
        the corresponding best alignments.
mat must have been filled with scores and prevs.
bestScore is the highest score in mat (provided for convenience).
cost is provided so mismatches with negative scores can be lowercased.
*/

void aliPrintBestAlisOld(struct matrix *mat, double bestScore, struct cost *cost, char *s1, char *s2) {
    unsigned int w = mat->w;
    unsigned int h = mat->h;
    struct cell *cells = mat->cells;

    /* allocation de deux strings == tableaux de char */
    char out1[h + w];
    char out2[h + w];

    /* pointeurs vers les caractères courants de out1 et out2, on
    commence a la fin des chaines et on va remonter en remplissant */
    char *out1_cur = out1 + h + w -1;
    char *out2_cur = out2 + h + w -1;

    // une string finit par \0
    *out1_cur = '\0' ;
    *out2_cur = '\0' ;

    unsigned int posBestScore = findBestScore(mat, bestScore);
    unsigned int i = posBestScore / w, j = posBestScore % w;
    uint8_t ultimatePrev = 1;
    while((cells[w*i+j].scoreD != 0) || (cells[w*i+j].scoreH != 0) || (cells[w*i+j].scoreV != 0)) {
        /* we came from the diag */
        if(ultimatePrev == 1) {
            if ((cells[w * i + j].scoreD - cells[w * (i - 1) + j - 1].scoreD) <= 0) {
                    *(--out1_cur) = tolower(s1[i - 1]);
                    *(--out2_cur) = tolower(s2[j - 1]);
            } else {
                    *(--out1_cur) = toupper(s1[i - 1]);
                    *(--out2_cur) = toupper(s2[j - 1]);
            }

            if ((cells[w * i + j].prevs & 1) == 1) {
                ultimatePrev = 1;
            } else if ((cells[w * i + j].prevs & 2) == 2) {

                ultimatePrev = 2;
            } else if ((cells[w * i + j].prevs & 4) == 4) {

                ultimatePrev = 4;
            }
            i--;
            j--;
        }

        /* we came from the side */
        else if(ultimatePrev == 2) {
            *(--out1_cur) = '-';
            *(--out2_cur) = tolower(s2[j - 1]);
            if ((cells[w * i + j].prevs & 8) == 8) {
                ultimatePrev = 1;
            } else if ((cells[w * i + j].prevs & 16) == 16) {
                ultimatePrev = 2;
            } else if ((cells[w * i + j].prevs & 32) == 32) {
                ultimatePrev = 4;
            }
            j--;
        }

        /* we came from below */
        else if(ultimatePrev == 4) {
            *(--out2_cur) = '-';
            *(--out1_cur) = tolower(s1[i - 1]);
            if ((cells[w * i + j].prevs & 64) == 64) {
                ultimatePrev = 1;
            } else if ((cells[w * i + j].prevs & 128) == 128) {
                ultimatePrev = 2;
            } else if ((cells[w * i + j].prevs & 256) == 256) {
                ultimatePrev = 4;
            }
            i--;
        }

    }

    printf("Best score is %f, the best-scoring alignements are : \n",bestScore);
    printf("\n");
    printf("s1 alignment starts at coord %u, s2 starts at coord %u \n", i+1, j+1);

    printf("s1          %s\ns2          %s\n\n", out1_cur, out2_cur);
}


void recursivCall(struct matrix *mat, struct cost *cost, char *s1, char *s2, unsigned int i, unsigned int j,
                  char *out1, char *out2, unsigned int start, uint8_t ultimatePrev) {

    unsigned int w = mat->w;
    unsigned int h = mat->h;
    struct cell *cells = mat->cells;
    char *out1_cur = out1;
    char *out2_cur = out2;
    bool recursiv = false;

    if ((i != 0 && j != 0) && (start == 1)) {
        if ((cells[w * i + j].scoreD - cells[w * (i - 1) + j - 1].scoreD) <= 0) {
            *(--out1_cur) = tolower(s1[i - 1]);
            *(--out2_cur) = tolower(s2[j - 1]);
        } else {
            *(--out1_cur) = toupper(s1[i - 1]);
            *(--out2_cur) = toupper(s2[j - 1]);
        }
        ultimatePrev = 0;
        if (cells[w * i + j].prevs == 0) {
            ultimatePrev = 0;
        }
        if ((cells[w * i + j].prevs & 1) == 1) {
            ultimatePrev |= 1;
        }
        if ((cells[w * i + j].prevs & 2) == 2) {
            ultimatePrev |= 2;
        }
        if ((cells[w * i + j].prevs & 4) == 4) {
            ultimatePrev |= 4;
        }
        i--;
        j--;
    } else if ((i != 0 && j != 0) && (start == 2)) {
        *(--out1_cur) = '-';
        *(--out2_cur) = tolower(s2[j - 1]);
        ultimatePrev = 0;
        if (cells[w * i + j].prevs == 0) {
            ultimatePrev = 0;
        }
        if ((cells[w * i + j].prevs & 8) == 8) {
            ultimatePrev |= 1;
        }
        if ((cells[w * i + j].prevs & 16) == 16) {
            ultimatePrev |= 2;
        }
        if ((cells[w * i + j].prevs & 32) == 32) {
            ultimatePrev |= 4;
        }
        j--;
    } else if ((i != 0 && j != 0) && (start == 3)) {
        *(--out2_cur) = '-';
        *(--out1_cur) = tolower(s1[i - 1]);
        ultimatePrev = 0;
        if (cells[w * i + j].prevs == 0) {
            ultimatePrev = 0;
        }
        if ((cells[w * i + j].prevs & 64) == 64) {
            ultimatePrev |= 1;
        }
        if ((cells[w * i + j].prevs & 128) == 128) {
            ultimatePrev |= 2;
        }
        if ((cells[w * i + j].prevs & 256) == 256) {
            ultimatePrev |= 4;
        }
        i--;
    }

    // The elexir to BUGS
    if(cells[w*i + j].prevs == 0) ultimatePrev = 0;

    while (ultimatePrev != 0) {
        if ((cells[w * i + j].scoreD == 0) && (cells[w * i + j].scoreH == 0) && (cells[w * i + j].scoreV == 0)) {
            printf("s1 alignment starts at coord %u, s2 starts at coord %u \n", i + 1, j + 1);
            printf("out1 %s\nout2 %s\n\n", out1_cur, out2_cur);
        }
        if ((ultimatePrev & 7) == 7) {
//            if(i == 0 || j == 0) printf("\nprev of cell is : %u\n", cells[w*i + j].prevs);
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 1, ultimatePrev);
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 2, ultimatePrev);
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 3, ultimatePrev);
            recursiv = true;
            break;
        } else if ((ultimatePrev & 6) == 6) {
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 2, ultimatePrev);
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 3, ultimatePrev);
            recursiv = true;
            break;
        } else if ((ultimatePrev & 5) == 5) {
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 1, ultimatePrev);
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 3, ultimatePrev);
            recursiv = true;
            break;
        } else if ((ultimatePrev & 3) == 3) {
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 1, ultimatePrev);
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 2, ultimatePrev);
            recursiv = true;
            break;
        }

            /* we came from the diag */
        else if ((ultimatePrev & 1) == 1) {
            ultimatePrev = 0;
            if ((cells[w * i + j].scoreD - cells[w * (i - 1) + j - 1].scoreD) <= 0) {
                *(--out1_cur) = tolower(s1[i - 1]);
                *(--out2_cur) = tolower(s2[j - 1]);
            } else {
                *(--out1_cur) = toupper(s1[i - 1]);
                *(--out2_cur) = toupper(s2[j - 1]);
            }

            if (cells[w * i + j].prevs == 0) {
                ultimatePrev = 0;
            }
            if ((cells[w * i + j].prevs & 1) == 1) {
                ultimatePrev |= 1;
            }
            if ((cells[w * i + j].prevs & 2) == 2) {
                ultimatePrev |= 2;
            }
            if ((cells[w * i + j].prevs & 4) == 4) {
                ultimatePrev |= 4;
            }
            i--;
            j--;
        }

            /* we came from the side */
        else if ((ultimatePrev & 2) == 2) {
            ultimatePrev = 0;
            *(--out1_cur) = '-';
            *(--out2_cur) = tolower(s2[j - 1]);
            if (cells[w * i + j].prevs == 0) {
                ultimatePrev = 0;
            }
            if ((cells[w * i + j].prevs & 8) == 8) {
                ultimatePrev |= 1;
            }
            if ((cells[w * i + j].prevs & 16) == 16) {
                ultimatePrev |= 2;
            }
            if ((cells[w * i + j].prevs & 32) == 32) {
                ultimatePrev |= 4;
            }
            j--;
        }

            /* we came from below */
        else if ((ultimatePrev & 4) == 4) {
            ultimatePrev = 0;
            *(--out2_cur) = '-';
            *(--out1_cur) = tolower(s1[i - 1]);
            if (cells[w * i + j].prevs == 0) {
                ultimatePrev = 0;
            }
            if ((cells[w * i + j].prevs & 64) == 64) {
                ultimatePrev |= 1;
            }
            if ((cells[w * i + j].prevs & 128) == 128) {
                ultimatePrev |= 2;
            }
            if ((cells[w * i + j].prevs & 256) == 256) {
                ultimatePrev |= 4;
            }
            i--;
        }

    }
}


void aliPrintBestAlisprime(struct matrix *mat, double bestScore, struct cost *cost, char *s1, char *s2, unsigned int posBestScore) {
    unsigned int w = mat->w;
    unsigned int h = mat->h;
    struct cell *cells = mat->cells;

    /* allocation de deux strings == tableaux de char */
    char out1[h + w];
    char out2[h + w];

    /* pointeurs vers les caractères courants de out1 et out2, on
    commence a la fin des chaines et on va remonter en remplissant */
    char *out1_cur = out1 + h + w -1;
    char *out2_cur = out2 + h + w -1;

    // une string finit par \0
    *out1_cur = '\0' ;
    *out2_cur = '\0' ;

//    unsigned int posBestScore = findBestScore(mat, bestScore);
    printf("Best score is %f, the best-scoring alignements are : \n", bestScore);
    printf("\n");

    unsigned int i = posBestScore / w, j = posBestScore % w;
    uint8_t ultimatePrev = 1;
    bool recursiv = false;
    while(ultimatePrev != 0) {
        if((cells[w*i+j].scoreD == 0) && (cells[w*i+j].scoreH == 0) && (cells[w*i+j].scoreV == 0)) {
            printf("s1 alignment starts at coord %u, s2 starts at coord %u \n", i+1, j+1);
            printf("out1 %s\nout2 %s\n\n", out1_cur, out2_cur);
        }
        if((ultimatePrev & 7) == 7) {
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 1, ultimatePrev);
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 2, ultimatePrev);
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 3, ultimatePrev);
            recursiv = true;
            break;
        }
        else if((ultimatePrev & 6) == 6) {
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 2, ultimatePrev);
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 3, ultimatePrev);
            recursiv = true;
            break;
        }
        else if((ultimatePrev & 5) == 5) {
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 1, ultimatePrev);
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 3, ultimatePrev);
            recursiv = true;
            break;
        }
        else if((ultimatePrev & 3) == 3) {
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 1, ultimatePrev);
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 2, ultimatePrev);
            recursiv = true;
            break;
        }

        /* we came from the diag */
        else if((ultimatePrev & 1) == 1) {
            ultimatePrev = 0;
            if ((cells[w * i + j].scoreD - cells[w * (i - 1) + j - 1].scoreD) <= 0) {
                *(--out1_cur) = tolower(s1[i - 1]);
                *(--out2_cur) = tolower(s2[j - 1]);
            } else {
                *(--out1_cur) = toupper(s1[i - 1]);
                *(--out2_cur) = toupper(s2[j - 1]);
            }

            if (cells[w * i + j].prevs == 0) {
                ultimatePrev = 0;
            } if ((cells[w * i + j].prevs & 1) == 1) {
                ultimatePrev |= 1;
            } if ((cells[w * i + j].prevs & 2) == 2) {
                ultimatePrev |= 2;
            } if ((cells[w * i + j].prevs & 4) == 4) {
                ultimatePrev |= 4;
            }
            i--;
            j--;
        }

            /* we came from the side */
        else if((ultimatePrev & 2) == 2) {
            ultimatePrev = 0;
            *(--out1_cur) = '-';
            *(--out2_cur) = tolower(s2[j - 1]);
            if (cells[w * i + j].prevs == 0) {
                ultimatePrev = 0;
            } if ((cells[w * i + j].prevs & 8) == 8) {
                ultimatePrev |= 1;
            } if ((cells[w * i + j].prevs & 16) == 16) {
                ultimatePrev |= 2;
            } if ((cells[w * i + j].prevs & 32) == 32) {
                ultimatePrev |= 4;
            }
            j--;
        }

            /* we came from below */
        else if((ultimatePrev & 4) == 4) {
            ultimatePrev = 0;
            *(--out2_cur) = '-';
            *(--out1_cur) = tolower(s1[i - 1]);
            if (cells[w * i + j].prevs == 0) {
                ultimatePrev = 0;
            } if ((cells[w * i + j].prevs & 64) == 64) {
                ultimatePrev |= 1;
            } if ((cells[w * i + j].prevs & 128) == 128) {
                ultimatePrev |= 2;
            } if ((cells[w * i + j].prevs & 256) == 256) {
                ultimatePrev |= 4;
            }
            i--;
        }

    }

    if(!recursiv) {
        printf("s1 alignment starts at coord %u, s2 starts at coord %u \n", i + 1, j + 1);
        printf("s1          %s\ns2          %s\n\n", out1_cur, out2_cur);
    }
}
