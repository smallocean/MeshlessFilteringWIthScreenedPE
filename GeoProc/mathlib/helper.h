//
// helper.h
//

/************************************************************
*
*   Copyritht (C) 2005    PAMI Lab, SJTU
*   All rights reserved
*
*   filename: helper.h
*   abstract: function declaration
*
*   coder: zhuww
*   date: 2005/11/21
*
************************************************************/

#pragma once

#ifndef _HELPER_H_
#define _HELPER_H_

#include "math_type.h"

NAMESPACE_MATHLIB_START

template <typename T_>
void rand_init(vector<T_>& rhs)
{
    for ( int i = 0; i < rhs.size(); i++ )
    {
        rhs[i] = (double)rand() / 0x7fff * 10.0;
    }
}

template <typename T_>
void display(const vector<T_>& vec)
{
    printf("\n");
    for ( int i = 0; i < vec.size(); i++ )
        printf("%.4f\t", vec[i]);
    
    printf("\n");
}

//template <>
//void display<int>(const MATHLIB::vector<int>& vec)
//{
//    printf("\n");
//    for ( int i = 0; i < vec.size(); i++ )
//        printf("%d\t", vec[i]);
//    
//    printf("\n");
//}

template <typename T_>
void rand_init(matrix<T_>& rhs)
{
    int i, j;
    for ( i = 0; i < rhs.rows(); i++ )
    {
        for ( j = 0; j < rhs.cols(); j++ )
        {
            rhs(i, j) = (double)rand() / 0x7fff * 10.0;
        }
    }
}

template <typename T_>
void display(const matrix<T_>& rhs)
{
    printf("\n");
    if ( !rhs.empty() )
    {
        int i, j;
        for ( i = 0; i < rhs.rows(); i++ )
        {
            for ( j = 0; j < rhs.cols(); j++ )
            {
                printf("%.4f", rhs(i, j));
                if ( j != rhs.cols() - 1 )
                    printf("\t");
            }
            printf("\n");
        }
    }
    printf("\n");
}

NAMESPACE_MATHLIB_END

#endif //_HELPER_H_
