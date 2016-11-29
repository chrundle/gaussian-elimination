#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include "../linearalgebra.h"

/* -------------------------- LUpivot -------------------------- */
/*  Given a matrix and its dimension, this function computes 
    a permutation matrix P, an upper triangular matrix U and a 
    lower triangular matrix L with ones on the main diagonal 
    such that PA = LU. 

    Input variables:
        a : pointer to array whose elements are the columns
             of matrix A.
        l : pointer to array in which the columns of L are to 
             be stored.
        p : pointer to the array {0, 1, 2, ..., m - 1} which 
             will be updated to represent the permutation 
             matrix P.
        m : number of rows and number of columns in matrix A.

    Features: This algorithm has time complexity ~(2/3)m^3
    flops and requires O(1) additional memory.                   */

void LUpivot (double ** a, double ** l, int * p, int m) {
    int i, j, k, index;
    
    for(i = 0; i < m - 1; i++) {
        index = subinfnorm_index(a[i], i, m - i);
        if(index > i) {
            /* interchange rows a[:][i] and a[:][index] */
            row_swap(a, m, index, i);

            /* interchange rows l[:][i] and l[:][index] */
            row_swap(l, m, index, i);

            /* interchange columns p[i] and p[index] */
            k = p[i];
            p[i] = p[index];
            p[index] = k;
        }

        for(j = i + 1; j < m; j++) {
            l[i][j] = a[i][j]/a[i][i];
            rowsubrow(a, l[i][j], m - i, i, j, i);   
        }
    }
}

int main () {
    int i, j, m, sym;
    double x;

    /* let user set the dimension of matrix A */
    printf("Enter the dimension m (where A is a m by m matrix): ");
    scanf("%i", &m);
    printf("Enter either 0 to test a nonsymmetric matrix\n"
           "          or 1 to test a symmetric matrix: ");
    scanf("%i", &sym);

    /* allocate memory for A, L, and P  */
    double ** a = new double*[m];
    double ** l = new double*[m];
    int * p = new int[m];
    for(i = 0; i < m; i++) {
        a[i] = new double[m];
        l[i] = new double[m];
    }

    /* initialize the values in matrix A */
    for(i = 0; i < m; i++) {
        for(j = 0; j < m; j++) {
            if(j < i) {
                if(sym) {
                    a[i][j] = i - j + 1;
                }
                else {
                    a[i][j] = i * j + 1;
                }
            }
            else {
                a[i][j] = j - i + 1; // this choice of values was arbitrary
            }
        }
    }

    /* initialize the values in matrix P */
    for(i = 0; i < m; i++) {
        p[i] = i;
    }

    /* print the matrix A before calling LUpivot */
    printf("A = \n");
    for(i = 0; i < m; i++) {
        for(j = 0; j < m; j++) {

            printf("%9.6g ", a[j][i]);
        }
        printf("\n");
    }
    printf("\n");

    /* execute householder recudtion to hessenberg form */
    LUpivot(a, l, p, m);

    /* print the matrix U (stored in A) after calling LUpivot */
    printf("R = \n");
    for(i = 0; i < m; i++) {
        for(j = 0; j < m; j++) {
            printf("%9.6g ", a[j][i]);
        }
        printf("\n");
    }
    printf("\n");

    /* print the matrix L after calling LUpivot */
    printf("L = \n");
    for(i = 0; i < m; i++) {
        for(j = 0; j < i; j++) {
            printf("%9.6g ", l[j][i]);
        }
        printf("%9.6g ", 1.0);
        j++;
        for(; j < m; j++) {
            printf("%9.6g ", l[j][i]);
        }
        printf("\n");
    }
    printf("\n");

    /* print the array P after calling LUpivot */
    printf("P = \n");
    for(i = 0; i < m; i++) {
        for(j = 0; j < m; j++) {
            if(j == p[i]) {
                printf("%9.6g ", 1.0); 
            }
            else {
                printf("%9.6g ", 0.0); 
            }
        }
        printf("\n");
    }

    /* free memory */
    for(i = 0; i < m; i++) {
        delete[] a[i];
        delete[] l[i];
    }
    delete[] a;
    delete[] l;
    delete[] p;
    return 0;
}
