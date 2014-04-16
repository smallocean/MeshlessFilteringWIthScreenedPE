//
// math_type.h
//

/************************************************************
*
*   Copyritht (C) 2005    PAMI Lab, SJTU
*   All rights reserved
*
*   filename: math_type.h
*   abstract: basic marco define
*
*   coder: zhuww
*   date: 2005/08/18
*
************************************************************/

#pragma once

#ifndef _MATH_TYPE_H_
#define _MATH_TYPE_H_

#include "matrix.h"
#include "vector.h"

#ifndef NAMESPACE_MATHLIB
#define NAMESPACE_MATHLIB_START namespace MATHLIB { // namespace mathlib start
#define NAMESPACE_MATHLIB_END } // namespace mathlib end
#endif // NAMESPACE_MATHLIB

NAMESPACE_MATHLIB_START

const double ZERO = 1.0e-12;
//const double TINY = 1.0e-12;

typedef unsigned long int uint;

typedef matrix<double> mat_f8;

typedef matrix<int> mat_int;

typedef matrix<bool> mat_bool;

typedef matrix<unsigned char> mat_byte;

typedef vector<mat_f8> vec_mf8;

typedef vector<double> vec_f8;

typedef vector<int> vec_int;

typedef vector<vec_f8> vec_vf8;

NAMESPACE_MATHLIB_END

#endif // _MATH_TYPE_H_
