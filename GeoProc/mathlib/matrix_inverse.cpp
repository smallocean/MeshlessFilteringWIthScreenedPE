//
// matrix_inverse.cpp
//

/************************************************************
*
*   Copyritht (C) 2005    PAMI Lab, SJTU
*   All rights reserved
*
*   filename: matrix_inverse.h
*   abstract: implement a kind of matrix class
*
*   coder: DeepBlue
*   date: 2005/11/12
*
************************************************************/

//#include "stdafx.h"
#include <cmath>
#include "matrix_inverse.h"

NAMESPACE_MATHLIB_START

bool matrix_inverse(mat_f8& src)
{
    assert( src.rows() == src.cols() );
    if ( src.rows() != src.cols() )
    {
        printf("matrix is not a square matrix!\n");
        return false;
    }

    if ( compute_rank(mat_f8(src)) < src.rows() )
    {
        printf("failed inverse matrix!\n");
        return false;
    }

    int i, j;
    int n = src.rows();

    vec_int indx(n);
    double  d = 0;
    mat_f8  y(n, n);
    vec_f8  col(n);

    if ( ! ludcmp(src, indx, d))
    {
        printf("failed inverse matrix!\n");
        return false;
    }

    for ( j = 0; j < n; j++ )
    {
        for ( i = 0; i < n; i++ )
            col[i] = 0;
        col[j] = 1.0;
        lubksb(src, indx, col);
        for ( i = 0; i < n; i++ )
            y(i, j) = col[i];
    }

    for ( i = 0; i < n; i++ )
    {
        for ( j = 0; j < n; j++ )
            src(i, j) = y(i, j);
    }

    return true;
}

bool ludcmp(mat_f8& src, vec_int& indx, double& d)
{
    assert( src.rows() == src.cols() );
    if ( src.rows() != src.cols() )
    {
        printf("matrix is not a square matrix!\n");
        return false;
    }

    int i, imax, j, k;
    int n = src.rows();
    double  big, dum, sum, temp;
    vec_f8  vv(n);
    d = 1.0;

    for ( i = 0; i < n; i++ )
    {
        big = 0.0;
        for ( j = 0; j < n; j++ )
        {
            if ( (temp = fabs(src(i, j))) > big )
                big = temp;
        }
        if ( fabs(big) < ZERO )
        //if ( big == 0.0f )
        {
            printf("Singular matrix in routine ludcmp\n");
            return false;
        }
        vv[i] = 1.0 / big;
    }

    for ( j = 0; j < n; j++ )
    {
        for (i=0; i < j; i++ )
        {
            sum = src(i, j);
            for ( k = 0; k < i; k++ )
                sum -= src(i, k) * src(k, j);
            src(i, j) = sum;
        }

        big = 0.0;
        for ( i = j; i < n; i++ )
        {
            sum = src(i, j);
            for ( k = 0; k < j; k++ )
                sum -= (src(i, k) * src(k, j));
            src(i, j) = sum;
            if ( (dum = vv[i] * fabs(sum) ) >= big )
            {
                big=dum;
                imax=i;
            }
        }

        if ( j != imax )
        {
            for ( k = 0; k < n; k++ )
            {
                dum = src(imax, k);
                src(imax, k) = src(j, k);
                src(j, k) = dum;
            }
            d = (0 - d);
            vv[imax] = vv[j];
        }

        indx[j] = imax;
        if ( fabs(src(j, j)) < ZERO )
            src(j, j) = ZERO;
        if (j != (n - 1) )
        {
            dum = 1.0 / (src(j, j));
            for ( i = j + 1; i < n; i++ )
                src(i, j) *= dum;
        }
    }

    return true;
}

void lubksb(mat_f8& src, vec_int& indx, vec_f8& b)
{
    int i, ii = 0, ip, j;
    double  sum;

    int n = src.rows();
    for ( i = 0; i < n; i++ )
    {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if ( ii != 0 )
        {
            for ( j = ii - 1; j < i; j++ )
                sum -= (src(i, j) * b[j]);
        }
        else if ( fabs(sum) > ZERO )
        //else if ( sum != 0.0f )
        {
            ii = i + 1;
        }

        b[i] = sum;
    }

    for ( i = n - 1; i >= 0; i-- )
    {
        sum = b[i];
        for ( j = i + 1; j < n; j++ )
            sum -= (src(i, j) * b[j]);
        b[i] = sum / src(i, i);
    }
}

int compute_rank(mat_f8& src)
{
    int rmax = 0;
    int j = 0;
    int rank = 0;

    int rows = src.rows();
    int cols = src.cols();

    for ( int i = 0; i < rows; i++ )
    {
        j = i;
        rmax = find_pivot_in_column(src, i, j);

        while ( fabs(src(rmax, j)) < ZERO )
        //while ( fabs(src(rmax, j)) == 0.0f )
        {
            if ( j == (cols - 1) )
                break;
            j++;
            rmax = find_pivot_in_column(src, i, j);
        }

        printf("%e\n", src(rmax, j));
        if ( (fabs(src(rmax, j)) < ZERO) && (j == (cols - 1)) )
        {
            rank = i;
            return rank;
        }
        else if ( rmax != i )
        {
            exchange(src, i, j, rmax);
        }

        elimination(src, i, j);
        rank++;
    }

    return rank;
}

// 将第i行与列主元素行交换
void exchange(mat_f8& src, const int row_i, int j, int row_j)
{
    int     p = 0;
    double  t = 0;

    int rows = src.rows();
    int cols = src.cols();
    for ( p = j; p < cols; p++ ) 
    {
        t = src(row_i, p);
        src(row_i, p) = src(row_j, p);
        src(row_j, p) = t;
    }
}

// 找出列主元素过程
int find_pivot_in_column(mat_f8& src, int i, int j)
{
     // 记录列主元素的行号，rmax
    int rmax = i;
    int rows= src.rows();
    for ( int k = i + 1; k < rows; k++ ) 
    {
        if ( fabs(src(k, j)) > fabs(src(i, j)) )
            rmax = k;
    }

    return rmax;
}

/*消元过程*/
void elimination(mat_f8& src, int i, int j)
{
    int s, k;
    int rows = src.rows();
    int cols = src.cols();

    double t = 0;
    for ( s = i + 1; s < rows; s++ )
    {
        t = src(s, j);
        for ( k = i; k < cols; k++ )
            src(s, k) = src(s, k) - t / src(i, j) * src(i, k);
    }
}

NAMESPACE_MATHLIB_END
