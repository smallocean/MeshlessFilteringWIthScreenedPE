//
// mathlib.h
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

#pragma once

#ifndef _MATH_UTIL_H_
#define _MATH_UTIL_H_

#include "math_type.h"
#include <cmath>

NAMESPACE_MATHLIB_START

template <typename T_>
// ����Ԫ����exp
matrix<T_> mat_exp(const matrix<T_>& rhs)
{
    matrix<T_> temp(rhs.rows(), rhs.cols());

    int i, j;
    for ( i = 0; i < rhs.rows(); i++ )
    {
        for ( j = 0; j < rhs.cols(); j++ )
        {
            temp(i, j) = exp(rhs(i, j));
        }
    }

    return temp;
}

// ����Ԫ����pow
template <typename T_>
matrix<T_> mat_pow(const matrix<T_>& rhs, const int power)
{
    matrix<T_> temp(rhs.rows(), rhs.cols());

    int i, j;
    for ( i = 0; i < rhs.rows(); i++ )
    {
        for ( j = 0; j < rhs.cols(); j++ )
        {
            temp(i,j) = pow(rhs(i, j), power);
        }
    }

    return temp;
}

// ����Ԫ����cos
template <class T_>
matrix<T_> mat_cos(const matrix<T_> rhs)
{
    matrix<T_> temp(rhs.rows(), rhs.cols());

    int i, j;
    for ( i = 0; i < rhs.rows(); i++ )
    {
        for ( j = 0; j < rhs.cols(); j++ )
        {
            temp(i, j) = cos(rhs(i, j));
        }
    }

    return temp;
}

// ����Ԫ����sin
template <class T_>
matrix<T_> mat_sin(const matrix<T_> rhs)
{
    matrix<T_> temp(rhs.rows(), rhs.cols());

    int i, j;
    for ( i = 0; i < rhs.rows(); i++ )
    {
        for ( j = 0; j < rhs.cols(); j++ )
        {
            temp(i, j) = sin(rhs(i, j));
        }
    }

    return temp;
}

// ������Э�������
void cov(mat_f8& rlt, const mat_f8& src);

mat_f8 cov(const mat_f8& src);


// ���б�׼������ģΪ1
void normalize_col(mat_f8& mat, const double sigma = 1.0f);

// ���б�׼������ģΪ1
void normalize_row(mat_f8& mat, const double sigma = 1.0f);

// ӳ�䵽0-1֮��
void normalize2_col(mat_f8& mat, const double sigma = 1.0f);

void normalize2_row(mat_f8& mat, const double sigma = 1.0f);

mat_f8 sum_col(const mat_f8& src);

mat_f8 sum_row(const mat_f8& src);

mat_f8 get_sub_cols(const mat_f8& src, const int begin, const int num);

mat_f8 get_sub_rows(const mat_f8& src, const int begin, const int num);

int max_in_vec(const vec_f8& src);

int mix_in_vec(const vec_f8& src);

void sum_row2sum(mat_f8& src, double sum);

void sum_col2sum(mat_f8& src, double sum);

NAMESPACE_MATHLIB_END

#endif // _MATH_UTIL_H_
