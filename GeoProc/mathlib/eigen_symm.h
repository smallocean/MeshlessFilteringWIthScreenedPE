//
// eigen_symm.h
//

/************************************************************
*
*   Copyritht (C) 2005    PAMI Lab, SJTU
*   All rights reserved
*
*   filename: eigen_symm.h
*   abstract: 
*
*   coder: DeepBlue
*   date: 2005/11/12
*
************************************************************/

#pragma once

#ifndef _EIGEN_SYSMM_H_
#define _EIGEN_SYSMM_H_

#include "math_type.h"

NAMESPACE_MATHLIB_START

///////////////////////////////////////////////////////////////////////////////
// Function:
//  int eigen_symm(vec_f8& eig_val, mat_f8& eig_vec, mat_f8 src)
//
// Params:
//  vec_f8& eig_val     -- 分解得到的特征值，降序排列
//  mat_f8& eig_vec     -- 分解得到的特征向量，与特征值对应，列向量
//  mat_f8 src          -- 待分解矩阵，实对称矩阵
//
// Return:
//  int nrot            -- 旋转次数
//
///////////////////////////////////////////////////////////////////////////////
int eigen_symm(vec_f8& eig_val, mat_f8& eig_vec, mat_f8 src);
//void jacobi(mat_f8& src, vec_f8& d, mat_f8& v, int& nrot);

NAMESPACE_MATHLIB_END

#endif //    _EIGEN_SYSMM_H_
