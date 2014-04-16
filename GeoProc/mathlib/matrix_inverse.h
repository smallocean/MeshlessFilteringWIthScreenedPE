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

// ��������
bool matrix_inverse(mat_f8& src);

// ����LU�ֽ�
bool ludcmp(mat_f8& src, vec_int& indx, double& d);

// 
void lubksb(mat_f8& src, vec_int& indx, vec_f8& b);

// ����������
int compute_rank(mat_f8& src);

// ����i��������Ԫ���н���
void exchange(mat_f8& src, const int row_i, int j, int row_j);

// �ҳ�����Ԫ�ع���
int find_pivot_in_column(mat_f8& src, int i, int j);

// ��Ԫ����
void elimination(mat_f8& src, int i, int j);

NAMESPACE_MATHLIB_END

#endif    // _MATRIX_INVERSE_H_
