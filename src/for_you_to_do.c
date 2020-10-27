#include "../include/for_you_to_do.h"
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

int get_block_size()
{
    //return the block size you'd like to use
    /*add your code here */
    return 66;
}

/**
 * 
 * this function computes LU factorization
 * for a square matrix
 * 
 * syntax 
 *  
 *  input : 
 *      A     n by n , square matrix
 *      ipiv  1 by n , vector
 *      n            , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf(double *A, int *ipiv, int n)
{
    int i, j, max_idx, k;
    double max_v;
    double *tmpv = (double *)malloc(sizeof(double) * n);

    for (i = 0; i < n; i++)
    {

        // Pivot
        max_v = fabs(A[i * n + i]);
        max_idx = i;
        for (j = i + 1; j < n; j++)
        {
            k = j * n + i;
            if (fabs(A[k]) > max_v)
            {
                max_v = fabs(A[k]);
                max_idx = j;
            }
        }
        if (max_v < 1e-3)
            return -1; // sigular
        max_v = A[max_idx * n + i];

        // swap
        if (max_idx != i)
        {
            int tmp = ipiv[i];
            ipiv[i] = ipiv[max_idx];
            ipiv[max_idx] = tmp;
            memcpy(tmpv, &A[i * n], n * sizeof(double));
            memcpy(&A[i * n], &A[max_idx * n], n * sizeof(double));
            memcpy(&A[max_idx * n], tmpv, n * sizeof(double));
        }

        // pivot = A[ipiv[i], i] = max_v
        for (j = i + 1; j < n; j++)
        {
            A[j * n + i] /= max_v;
        }

        // update matrix
        for (j = i + 1; j < n; j++)
        {
            for (k = i + 1; k < n; k++)
            {
                A[j * n + k] -= A[j * n + i] * A[i * n + k];
            }
        }
    }
    free(tmpv);
    return 0;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix . according to lecture slides, this
 * function computes forward AND backward subtitution in the
 * same function.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      none
 * 
 **/
void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    /* add your code here */
    double *y = (double *)malloc(n * sizeof(double));
    int i, j;
    double sum;
    if (UPLO == 'L')
    {
        y[0] = B[ipiv[0]];
        for (i = 1; i < n; i++)
        {
            sum = 0.0;
            for (j = 0; j < i; j++)
            {
                sum += y[j] * A[i * n + j];
            }
            y[i] = B[ipiv[i]] - sum;
        }
    }
    else if (UPLO == 'U')
    {
        y[n - 1] = B[n - 1] / A[(n - 1) * n + n - 1];
        for (i = n - 2; i >= 0; i--)
        {
            sum = 0;
            for (j = i + 1; j < n; j++)
            {
                sum += y[j] * A[i * n + j];
            }
            y[i] = (B[i] - sum) / A[i * n + i];
        }
    }
    memcpy(B, y, sizeof(double) * n);
    free(y);

    return;
}

/**
 * 
 * Same function as what you used in lab1, cache_part4.c : optimal( ... ).
 * 
 **/
void mydgemm(double *A, double *B, double *C, int n, int i, int j, int k, int b)
{
    /* add your code here */
    /* please just copy from your lab1 function optimal( ... ) */
    int i1, j1, k1, l, ni, nj, nk;
    ni = i + b;
    if (ni > n)
        ni = n;
    nj = j + b;
    if (nj > n)
        nj = n;
    nk = k + b;
    if (nk > n)
        nk = n;

    for (i1 = i; i1 < ni; i1 += 3)
    {
        for (j1 = j; j1 < nj; j1 += 3)
        {
            int t = i1 * n + j1;
            int tt = t + n;
            int ttt = tt + n;
            register double c00 = C[t];
            register double c01 = C[t + 1];
            register double c02 = C[t + 2];
            register double c10 = C[tt];
            register double c11 = C[tt + 1];
            register double c12 = C[tt + 2];
            register double c20 = C[ttt];
            register double c21 = C[ttt + 1];
            register double c22 = C[ttt + 2];

            for (k1 = k; k1 < nk; k1 += 3)
            {
                for (l = 0; l < 3; l++)
                {
                    int ta = i1 * n + k1 + l;
                    int tta = ta + n;
                    int ttta = tta + n;
                    int tb = k1 * n + j1 + l * n;
                    register double a00 = A[ta];
                    register double a10 = A[tta];
                    register double a20 = A[ttta];
                    register double b00 = B[tb];
                    register double b01 = B[tb + 1];
                    register double b02 = B[tb + 2];

                    c00 -= a00 * b00;
                    c01 -= a00 * b01;
                    c02 -= a00 * b02;
                    c10 -= a10 * b00;
                    c11 -= a10 * b01;
                    c12 -= a10 * b02;
                    c20 -= a20 * b00;
                    c21 -= a20 * b01;
                    c22 -= a20 * b02;
                }
            }
            C[t] = c00;
            C[t + 1] = c01;
            C[t + 2] = c02;
            C[tt] = c10;
            C[tt + 1] = c11;
            C[tt + 2] = c12;
            C[ttt] = c20;
            C[ttt + 1] = c21;
            C[ttt + 2] = c22;
        }
    }
    return;
}

/**
 * 
 * this function computes LU decomposition
 * for a square matrix using block gepp introduced in course
 * lecture .
 * 
 * just implement the block algorithm you learned in class.
 * 
 * syntax 
 *  
 *  input :
 *      
 * 
 *      A     n by n     , square matrix
 * 
 *    
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *     
 *      b                , block size   
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf_block(double *A, int *ipiv, int n, int b)
{
    int ib, end, i, j, k, max_idx;
    double max_v, sum;
    double *tmpv = (double *)malloc(sizeof(double) * n);

    for (ib = 0; ib < n; ib += b)
    {
        end = ib + b;
        if(end > n) end = n;

        // A(ib:n , ib:end) = P’ * L’ * U’
        for (i = ib; i < end; i++)
        {
            // Pivot
            max_v = fabs(A[i * n + i]);
            max_idx = i;
            for (j = i + 1; j < n; j++)
            {
                k = j * n + i;
                if (fabs(A[k]) > max_v)
                {
                    max_v = fabs(A[k]);
                    max_idx = j;
                }
            }
            if (max_v < 1e-3)
                return -1; // sigular
            max_v = A[max_idx * n + i];

            // swap
            if (max_idx != i)
            {
                int tmp = ipiv[i];
                ipiv[i] = ipiv[max_idx];
                ipiv[max_idx] = tmp;
                memcpy(tmpv, &A[i * n], n * sizeof(double));
                memcpy(&A[i * n], &A[max_idx * n], n * sizeof(double));
                memcpy(&A[max_idx * n], tmpv, n * sizeof(double));
            }

            // LU factorization (Blue in slides)
            for (j = i + 1; j < n; j++)
            {
                A[j * n + i] /= max_v;
                for (k = i + 1; k < end; k++)
                {
                    A[j * n + k] -= A[j * n + i] * A[i * n + k];
                }
            }
        }

        // update A(ib:end , end+1:n) ( pink in slides)
        for (i = ib; i < end; i++)
        {
            for (j = end; j < n; j++)
            {
                sum = 0;
                for (k = ib; k < i; k++)
                {
                    sum += A[i * n + k] * A[k * n + j];
                }
                A[i * n + j] -= sum;
            }
        }

        for (i = end; i < n; i += b)
        {
            for (j = end; j < n; j += b)
            {
                mydgemm(A, A, A, n, i, j, ib, b);
            }
        }
    }
    free(tmpv);
    return 0;
}