#include "stdlib.h"
#include "colamd.h"
#define A_NNZ 11
#define A_NROW 5
#define A_NCOL 4
#define ALEN 150

int main (void)
{
    int A [ALEN] = {
        0, 1, 4,             /* row indices of nonzeros in column 0 */
        2, 4,                    /* row indices of nonzeros in column 1 */
        0, 1, 2, 3,      /* row indices of nonzeros in column 2 */
        1, 3} ;                 /* row indices of nonzeros in column 3 */

    int p [ ] = {
        0,                      /* column 0 is in A [0..2] */
        3,                      /* column 1 is in A [3..4] */
        5,                      /* column 2 is in A [5..8] */
        9,                      /* column 3 is in A [9..10] */
        A_NNZ                   /* number of nonzeros in A */
    };

    int stats [COLAMD_STATS];
    int ok;

    /* compute the column ordering */
    ok = colamd (A_NROW, A_NCOL, ALEN, A, p, (double *) NULL, stats);
    if (!ok)
    {
        printf ("colamd error!\n") ;
        exit (1) ;
    }

    printf ("colamd column ordering:\n") ;
    printf ("1st column: %d\n", p [0]) ;
    printf ("2nd column: %d\n", p [1]) ;
    printf ("3rd column: %d\n", p [2]) ;
    printf ("4th column: %d\n", p [3]) ;
}
