//
// matrix_inverse.h
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

#pragma once

#ifndef _MATRIX_INVERSE_H_
#define _MATRIX_INVERSE_H_

#include "math_type.h"

NAMESPACE_MATHLIB_START

// 矩阵求逆
bool matrix_inverse(mat_f8& src);

// 矩阵LU分解
bool ludcmp(mat_f8& src, vec_int& indx, double& d);

// 
void lubksb(mat_f8& src, vec_int& indx, vec_f8& b);

// 计算矩阵的秩
int compute_rank(mat_f8& src);

// 将第i行与列主元素行交换
void exchange(mat_f8& src, const int row_i, int j, int row_j);

// 找出列主元素过程
int find_pivot_in_column(mat_f8& src, int i, int j);

// 消元过程
void elimination(mat_f8& src, int i, int j);

NAMESPACE_MATHLIB_END

#endif    // _MATRIX_INVERSE_H_
