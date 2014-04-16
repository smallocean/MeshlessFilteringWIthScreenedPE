//
// eigen_symm.cpp
//

/************************************************************
*
*   Copyritht (C) 2005    PAMI Lab, SJTU
*   All rights reserved
*
*   filename: eigen_symm.cpp
*   abstract: 
*
*   coder: DeepBlue
*   date: 2005/11/12
*
************************************************************/

//#include "stdafx.h"

#include "eigen_symm.h"

NAMESPACE_MATHLIB_START

// rotate
inline void eigen_rot(mat_f8& src, const double s, const double tau, const int i,
                const int j, const int k, const int l);

// sort eigenvalue from big to small
void eigen_sort(vec_f8& eig_val, mat_f8& eig_vec);

///////////////////////////////////////////////////////////////////////////////
// Function:
//  int eigen_symm(vec_f8& eig_val, mat_f8& eig_vec, mat_f8 src)
//
// Params:
//  vec_f8& eig_val -- 分解得到的特征值，降序排列
//  mat_f7& eig_vec -- 分解得到的特征向量，与特征值对应，列向量
//  mat_f8 src      -- 待分解矩阵，实对称矩阵
//
// Return:
//  int nrot        -- 旋转次数
//
///////////////////////////////////////////////////////////////////////////////
int eigen_symm(vec_f8& eig_val, mat_f8& eig_vec, mat_f8 src)
{
    if ( src.rows() != src.cols() )
        return 0;

    int i, j, ip, iq;
    double tresh, theta, tau, t, sm, s, h, g, c;

    int n = src.rows();
    vec_f8 b(n), z(n);
    if (n != eig_val.size())
        eig_val.resize(n);
    if ((n != eig_vec.rows()) || (n != eig_vec.cols()))
        eig_vec.resize(n, n);

    for ( ip = 0; ip < n; ip++ )
    {
        for ( iq = 0; iq < n; iq++ )
            eig_vec(ip, iq)=0.0;
        eig_vec(ip, ip)=1.0;
    }

    for ( ip = 0; ip < n; ip++ )
    {
        b[ip] = eig_val[ip] = src(ip, ip);
        z[ip] = 0.0;
    }

    int nrot=0;
    for ( i = 1; i <= 50; i++ )
    {
        sm=0.0;

        for ( ip = 0; ip < n - 1; ip++ )
        {
            for ( iq = ip + 1; iq < n; iq++ )
                sm += fabs(src(ip, iq));
        }

        //if ( sm == 0.0 )
        if ( fabs(sm) < ZERO )
        {
            eigen_sort(eig_val, eig_vec);
            return nrot;
        }

        if ( i < 4 )
            tresh = 0.2 * sm / (n * n);
        else
            tresh = 0.0;

        for ( ip = 0; ip < n - 1; ip++ )
        {
            for ( iq = ip + 1; iq < n; iq++ )
            {
                g = 100.0 * fabs(src(ip, iq));

                if ( i > 4 && (fabs(eig_val[ip]) + g) == fabs(eig_val[ip])
                    && (fabs(eig_val[iq]) + g) == fabs(eig_val[iq]) )
                {
                    src(ip, iq) = 0.0;
                }
                else if ( fabs(src(ip, iq)) > tresh )
                {
                    h = eig_val[iq] - eig_val[ip];

                    if ((fabs(h) + g) == fabs(h))
                    {
                        t = (src(ip, iq)) / h;
                    }
                    else
                    {
                        theta = 0.5 * h / (src(ip, iq));
                        t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
                        if ( theta < 0.0 )
                            t = -t;
                    }

                    c = 1.0 / sqrt(1 + (t * t));
                    s = t * c;
                    tau = s / (1.0 + c);
                    h = t * src(ip, iq);
                    z[ip] -= h;
                    z[iq] += h;
                    eig_val[ip] -= h;
                    eig_val[iq] += h;
                    src(ip, iq) = 0.0;

                    for ( j = 0; j < ip; j++ )
                        eigen_rot(src, s, tau, j, ip, j, iq);

                    for ( j = ip + 1; j < iq; j++ )
                        eigen_rot(src, s, tau, ip, j, j, iq);

                    for ( j = iq + 1; j < n; j++ )
                        eigen_rot(src, s, tau, ip, j, iq, j);

                    for ( j = 0; j < n; j++ )
                        eigen_rot(eig_vec, s, tau, j, ip, j, iq);

                    ++nrot;
                }
            }
        }

        for ( ip = 0; ip < n; ip++ )
        {
            b[ip] += z[ip];
            eig_val[ip] = b[ip];
            z[ip] = 0.0;
        }
    }

    eigen_sort(eig_val, eig_vec);

    return nrot;
}

inline void eigen_rot(mat_f8& src, const double s, const double tau, const int i,
                const int j, const int k, const int l)
{
    double g, h;

    g = src(i, j);
    h = src(k, l);
    src(i, j) = g - s * (h + g * tau);
    src(k, l) = h + s * (g - h * tau);
}

void eigen_sort(vec_f8& eig_val, mat_f8& eig_vec)
{
    int i, j, k;
    double p;

    int n = eig_val.size();
    for ( i = 0; i < n - 1; i++ )
    {
        p = eig_val[i];
        k = i;
        for ( j = i; j < n; j++ )
        {
            if ( eig_val[j] >= p )
            {
                p = eig_val[j];
                k = j;
            }
        }
        if ( k != i )
        {
            eig_val[k] = eig_val[i];
            eig_val[i] = p;
            for ( j = 0; j < n; j++ )
            {
                p = eig_vec(j, i);
                eig_vec(j, i) = eig_vec(j, k);
                eig_vec(j, k) = p;
            }
        }
    }
}

NAMESPACE_MATHLIB_END
