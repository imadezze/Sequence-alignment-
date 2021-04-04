
#include "aliCalcSW.h"
#include "aliCost.h"
#include "aliOutSW.h"
#include <ctype.h>
#include <stdio.h>
#include <stdbool.h>

unsigned int findBestScore(struct matrix *mat, double bestScore){
    unsigned int w = mat->w;
    unsigned int h = mat->h;
    struct cell *cells = mat->cells;
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            if(cells[w*i+j].score == bestScore){
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
    while(cells[w*i+j].score != 0) {
        if ((cells[w * i + j].prevs & 1) == 1) {
            if ((cells[w * i + j].score - cells[w * (i - 1) + j - 1].score) <= 0) {
                *(--out1_cur) = tolower(s1[i - 1]);
                *(--out2_cur) = tolower(s2[j - 1]);
            } else {
                *(--out1_cur) = toupper(s1[i - 1]);
                *(--out2_cur) = toupper(s2[j - 1]);
            }
            i--;
            j--;
        } else if ((cells[w * i + j].prevs & 2) == 2) {
            *(--out1_cur) = '-';
            *(--out2_cur) = tolower(s2[j - 1]);
            j--;
        } else if ((cells[w * i + j].prevs & 4) == 4) {
            *(--out2_cur) = '-';
            *(--out1_cur) = tolower(s1[i - 1]);
            i--;
        }

    }
    printf("%f \n",bestScore);

    printf("out1 %s\nout2 %s\n\n", out1_cur, out2_cur);
}


void recursivCall(struct matrix *mat, struct cost *cost, char *s1, char *s2, unsigned int i, unsigned int j,
                    char *out1, char *out2, unsigned int start) {
    unsigned int w = mat->w;
    unsigned int h = mat->h;
    struct cell *cells = mat->cells;
    char *out1_cur = out1;
    char *out2_cur = out2;
    bool recursiv = false;

    if(start == 1) {
        if ((cells[w * i + j].score - cells[w * (i - 1) + j - 1].score) <= 0) {
            *(--out1_cur) = tolower(s1[i - 1]);
            *(--out2_cur) = tolower(s2[j - 1]);
        } else {
            *(--out1_cur) = toupper(s1[i - 1]);
            *(--out2_cur) = toupper(s2[j - 1]);
        }
        i--;
        j--;
    }
    else if(start == 2) {
        *(--out1_cur) = '-';
        *(--out2_cur) = tolower(s2[j - 1]);
        j--;
    }
    else if(start == 3) {
        *(--out2_cur) = '-';
        *(--out1_cur) = tolower(s1[i - 1]);
        i--;
    }
    while(cells[w*i+j].prevs != 0) {
        if(cells[w*i+j].score == 0) {
            printf("s1 alignment starts at coord %u, s2 starts at coord %u \n", i+1, j+1);
            printf("out1 %s\nout2 %s\n\n", out1_cur, out2_cur);
        }
        if((cells[w * i + j].prevs & 7) == 7) {
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 1);
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 2);
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 3);
            recursiv = true;
            break;
        }
        else if((cells[w * i + j].prevs & 6) == 6) {
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 2);
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 3);
            recursiv = true;
            break;
        }
        else if((cells[w * i + j].prevs & 5) == 5) {
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 1);
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 3);
            recursiv = true;
            break;
        }
        else if((cells[w * i + j].prevs & 3) == 3) {
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 1);
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 2);
            recursiv = true;
            break;
        }
        else if ((cells[w * i + j].prevs & 1) == 1) {
            if ((cells[w * i + j].score - cells[w * (i - 1) + j - 1].score) <= 0) {
                *(--out1_cur) = tolower(s1[i - 1]);
                *(--out2_cur) = tolower(s2[j - 1]);
            } else {
                *(--out1_cur) = toupper(s1[i - 1]);
                *(--out2_cur) = toupper(s2[j - 1]);
            }
            i--;
            j--;
        } else if ((cells[w * i + j].prevs & 2) == 2) {
            *(--out1_cur) = '-';
            *(--out2_cur) = tolower(s2[j - 1]);
            j--;
        } else if ((cells[w * i + j].prevs & 4) == 4) {
            *(--out2_cur) = '-';
            *(--out1_cur) = tolower(s1[i - 1]);
            i--;
        }

    }
    if(!recursiv) {
        printf("s1 alignment starts at coord %u, s2 starts at coord %u \n", i+1, j+1);
        printf("out1 %s\nout2 %s\n\n", out1_cur, out2_cur);
    }

}

void aliPrintBestAlis(struct matrix *mat, double bestScore, struct cost *cost, char *s1, char *s2) {
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
    bool recursiv = false;

    printf("%f \n",bestScore);

    while(cells[w*i+j].prevs != 0) {
        if(cells[w*i+j].score == 0) {
            printf("s1 alignment starts at coord %u, s2 starts at coord %u \n", i+1, j+1);
            printf("out1 %s\nout2 %s\n\n", out1_cur, out2_cur);
        }
        if((cells[w * i + j].prevs & 7) == 7) {
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 1);
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 2);
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 3);
            recursiv = true;
            break;
        }
        else if((cells[w * i + j].prevs & 6) == 6) {
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 2);
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 3);
            recursiv = true;
            break;
        }
        else if((cells[w * i + j].prevs & 5) == 5) {
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 1);
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 3);
            recursiv = true;
            break;
        }
        else if((cells[w * i + j].prevs & 3) == 3) {
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 1);
            recursivCall(mat, cost, s1, s2, i, j, out1_cur, out2_cur, 2);
            recursiv = true;
            break;
        }
        else if ((cells[w * i + j].prevs & 1) == 1) {
            if ((cells[w * i + j].score - cells[w * (i - 1) + j - 1].score) <= 0) {
                *(--out1_cur) = tolower(s1[i - 1]);
                *(--out2_cur) = tolower(s2[j - 1]);
            } else {
                *(--out1_cur) = toupper(s1[i - 1]);
                *(--out2_cur) = toupper(s2[j - 1]);
            }
            i--;
            j--;
        } else if ((cells[w * i + j].prevs & 2) == 2) {
            *(--out1_cur) = '-';
            *(--out2_cur) = tolower(s2[j - 1]);
            j--;
        } else if ((cells[w * i + j].prevs & 4) == 4) {
            *(--out2_cur) = '-';
            *(--out1_cur) = tolower(s1[i - 1]);
            i--;
        }

    }

    if(!recursiv) {
        printf("s1 alignment starts at coord %u, s2 starts at coord %u \n", i+1, j+1);
        printf("out1 %s\nout2 %s\n\n", out1_cur, out2_cur);
    }

}