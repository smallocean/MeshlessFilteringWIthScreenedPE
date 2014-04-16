//
// math_util.cpp
//

/************************************************************
*
*   Copyritht (C) 2005    PAMI Lab, SJTU
*   All rights reserved
*
*   filename: math_util.h
*   abstract: basic marco define
*
*   coder: zhuww
*   date: 2005/08/18
*
************************************************************/

//#include "stdafx.h"

#include "math_util.h"
#include <cmath>

NAMESPACE_MATHLIB_START

// [m, n] --> [1, n]
mat_f8 sum_col(const mat_f8& src)
{
    assert( 0 != src.size() );

    if ( 0 == src.size() )
        return mat_f8();

    mat_f8 rlt(1, src.cols());

    int i, j;
    for ( j = 0; j < src.cols(); j++ )
    {
        rlt(0, j) = 0;
        for ( i = 0; i < src.rows(); i++ )
        {
            rlt(0, j) += src(i, j);
        }
    }

    return rlt;
}

// [m, n] --> [m, 1]
mat_f8 sum_row(const mat_f8& src)
{
    assert( 0 != src.size() );

    if ( 0 == src.size() )
        return mat_f8();

    mat_f8 rlt(src.rows(), 1);

    int i, j;
    for ( i = 0; i < src.rows(); i++ )
    {
        rlt(i, 0) = 0;
        for ( j = 0; j < src.cols(); j++ )
        {
            rlt(i, 0) += src(i, j);
        }
    }

    return rlt;
}

void cov(mat_f8& rlt, const mat_f8& src)
{
    //rlt.resize(src.cols(), src.cols());

    mat_f8 mat1(src);
    vec_f8 avr(src.cols());

    int rows = src.rows();
    int cols = src.cols();

    int i, j;
    for ( j = 0; j < cols; j++ )
    {
        avr[j] = 0;
        for ( i = 0; i < rows; i++ )
        {
            avr[j] += src(i, j);
        }
        avr[j] /= rows;

        for ( i = 0; i < rows; i++ )
        {
            mat1(i, j) -= avr[j];
        }
    }

    mat_f8 mat2 = mat1;
    mat2.transpose();

    rlt = mat2 * mat1;

    rows = rlt.rows();
    cols = rlt.cols();
    int N = src.rows() - 1;
    for ( i = 0; i < rows; i++ )
    {
        for ( j = 0; j < cols; j++ )
        {
            rlt(i, j) /= N;
        }
    }
}

mat_f8 cov(const mat_f8& src)
{
    mat_f8 result;//(src.cols(), src.cols());

    mat_f8 mat1(src);
    vec_f8 avr(src.cols());

    int rows = src.rows();
    int cols = src.cols();

    int i, j;
    for ( j = 0; j < cols; j++ )
    {
        avr[j] = 0;
        for ( i = 0; i < rows; i++ )
        {
            avr[j] += src(i, j);
        }
        avr[j] /= rows;

        for ( i = 0; i < rows; i++ )
        {
            mat1(i, j) -= avr[j];
        }
    }

    mat_f8 mat2 = mat1;
    mat2.transpose();

    result = mat2 * mat1;

    rows = result.rows();
    cols = result.cols();
    int N = src.rows() - 1;
    for ( i = 0; i < rows; i++ )
    {
        for ( j = 0; j < cols; j++ )
        {
            result(i, j) /= N;
        }
    }

    return result;
}

// 对列标准化，列模为1
void normalize_col(mat_f8& mat, const double sigma)
{
    mat_f8 temp(mat);
    temp *= mat;

    int rows = mat.rows();
    int cols = mat.cols();

    double sum;
    int i, j;
    for ( j = 0; j < cols; j++ )
    {
        sum = 0.0f;
        for ( i = 0; i < rows; i++ )
        {
            sum += temp(i, j);
        }
        sum = sqrt(sum);
        sum /= sigma;

        for ( i = 0; i < rows; i++ )
        {
            mat(i, j) /= sum;
        }
    }
}

// 对行标准化，行模为1
void normalize_row(mat_f8& mat, const double sigma)
{
    mat_f8 temp(mat);
    temp *= mat;

    int rows = mat.rows();
    int cols = mat.cols();

    double sum;
    int i, j;
    for ( i = 0; i < rows; i++ )
    {
        sum = 0.0f;
        for ( j = 0; j < cols; j++ )
        {
            sum += temp(i, j);
        }
        sum = sqrt(sum);
        sum /= sigma;

        for ( j = 0; j < cols; j++ )
        {
            mat(i, j) /= sum;
        }
    }
}

// 映射到0-1之间，列
void normalize2_col(mat_f8& mat, const double sigma)
{
    int i, j;
    double min, max;
    for ( j = 0; j < mat.cols(); j++ )
    {
        min = mat(0, j);
        max = mat(0, j);
        for ( i = 1; i < mat.rows(); i++ )
        {
            if ( min > mat(i, j) )
                min = mat(i, j);
            if ( max < mat(i, j) )
                max = mat(i, j);
        }
        double d = max - min;
        for ( i = 0; i < mat.rows(); i++ )
        {
            mat(i, j) = (mat(i, j) - min) / d * sigma;
        }
    }
}

// 映射到0-1之间，行
void normalize2_row(mat_f8& mat, const double sigma)
{
    int i, j;
    double min, max;
    for ( i = 0; i < mat.rows(); i++ )
    {
        min = mat(i, 0);
        max = mat(i, 0);
        for ( j = 1; j < mat.cols(); j++ )
        {
            if ( min > mat(i, j) )
                min = mat(i, j);
            if ( max < mat(i, j) )
                max = mat(i, j);
        }
        double d = max - min;
        for ( j = 0; j < mat.cols(); j++ )
        {
            mat(i, j) = (mat(i, j) - min) / d * sigma;
        }
    }
}

mat_f8 get_sub_cols(const mat_f8& src, const int begin, const int num)
{
    assert( (begin >= 0) && ((begin + num) <= src.cols()) );

    mat_f8 rlt;
    if ( (begin < 0) || ((begin + num) > src.cols()) || (0 == num) )
        return rlt;

    int row = src.rows();
    int col = num;
    rlt.resize(row, col);

    int i, j;
    for ( i = 0; i < row; i++ )
    {
        for ( j = 0; j < col; j++ )
        {
            rlt(i, j) = src(i, j + begin);
        }
    }

    return rlt;
}

mat_f8 get_sub_rows(const mat_f8& src, const int begin, const int num)
{
    assert( (begin >= 0) && ((begin + num) <= src.rows()) );

    mat_f8 rlt;
    if ( (begin < 0) || ((begin + num) > src.rows()) || (0 == num) )
        return rlt;

    int row = num;
    int col = src.cols();
    rlt.resize(row, col);

    int i, j;
    for ( i = 0; i < row; i++ )
    {
        for ( j = 0; j < col; j++ )
        {
            rlt(i, j) = src(i + begin, j);
        }
    }

    return rlt;
}

int max_in_vec(const vec_f8& src)
{
    double  max = src[0];
    int     size = src.size();
    int     pos = 0;
    for (int i = 1; i < size; i++)
    {
        if (src[i] > max)
        {
            max = src[i];
            pos = i;
        }
    }

    return pos;
}

int mix_in_vec(const vec_f8& src)
{
    double  min = src[0];
    int     size = src.size();
    int     pos = 0;
    for (int i = 1; i < size; i++)
    {
        if (src[i] < min)
        {
            min = src[i];
            pos = i;
        }
    }

    return pos;
}

void sum_row2sum(mat_f8& src, double sum)
{
    int i, j;
    int row = src.rows();
    int col = src.cols();
    double temp;
    for ( i = 0; i < row; i++ )
    {
        temp = 0.0;
        for ( j = 0; j < col; j++ )
        {
            temp += src(i, j);
        }
        for ( j  = 0; j < col; j++ )
        {
            src(i, j) /= (temp * sum);
        }
    }
}

void sum_col2sum(mat_f8& src, double sum)
{
}

NAMESPACE_MATHLIB_END
