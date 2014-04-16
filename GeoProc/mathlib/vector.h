//
// vector.h
//

/************************************************************
*
*   Copyritht (C) 2005    PAMI Lab, SJTU
*   All rights reserved
*
*   filename: vector.h
*   abstract: implement a kind of vector class
*
*   coder: DeepBlue
*   date: 2005/08/18
*
************************************************************/

#pragma once

#ifndef _VECTOR_H_
#define _VECTOR_H_

#include "matrix.h"
#include <cassert>
#include <cmath>

//#ifndef NAMESPACE_MATHLIB
//#define NAMESPACE_MATHLIB_START namespace MATHLIB { // namespace mathlib start
//#define NAMESPACE_MATHLIB_END } // namespace mathlib end
//#endif // NAMESPACE_MATHLIB

namespace MATHLIB
{
    #ifndef NULL
    #define NULL 0
    #endif

    template <typename T_> class vector;
    template <typename T_> const vector<T_> operator + (const vector<T_>& Left_, const vector<T_>& Right_);
    template <typename T_> const vector<T_> operator - (const vector<T_>& Left_, const vector<T_>& Right_);
    template <typename T_> const T_ operator * (const vector<T_>& Left_, const vector<T_>& Right_);

    template <typename T_> class matrix;
    template <typename T_> const vector<T_> operator * (const matrix<T_>& Left_, const vector<T_>& Right_);
    template <typename T_> const vector<T_> operator * (const vector<T_>& Left_, const matrix<T_>& Right_);

    // declare vector class
    template <typename T_> class matrix;
    template <typename T_>
    class vector
    {
        // constructors
    public:
        vector()
            : m_data(NULL)
            , m_len(0)
        {    };

        vector(const vector& rhs)
            : m_data(NULL)
            , m_len(0)
        { *this = rhs; };

        vector(const matrix<T_>& rhs)
            : m_data(NULL)
            , m_len(0)
        { *this = rhs; };

        explicit vector(const int Len_)
            : m_data(new T_[Len_]())
            , m_len(Len_)
        {    };

        explicit vector(const int Len_, const T_& Val_)
            : m_data(new T_[Len_]())
            , m_len(Len_)
        {
            for ( int i = 0; i < m_len; i++ )
                m_data[i] = Val_;
        };

        // destructor
    public:
        ~vector()
        { freemem(); }

        // get element for vector
    public:
        T_& operator () (const int Len_)
        {
            assert( (Len_ >= 0) && (Len_ < m_len) );

            if ( (Len_ >= 0) && (Len_ < m_len) )
                return m_data[Len_];
            else if ( Len_ < 0 )
                return m_data[0];
            else // if ( Len_ >= m_len )
                return m_data[m_len - 1];
        };

        const T_& operator () (const int Len_) const
        {
            //assert( (Len_ >= 0) && (Len_ < m_len) );

            if ( (Len_ >= 0) && (Len_ < m_len) )
                return m_data[Len_];
            else if ( Len_ < 0 )
                return m_data[0];
            else // if ( Len_ >= m_len )
                return m_data[m_len - 1];
        };

        T_& operator [] (const int Len_)
        {
            //assert( (Len_ >= 0) && (Len_ < m_len) );

            if ( (Len_ >= 0) && (Len_ < m_len) )
                return m_data[Len_];
            else if ( Len_ < 0 )
                return m_data[0];
            else // if ( Len_ >= m_len )
                return m_data[m_len - 1];
        };

        const T_& operator [] (const int Len_) const
        {
            //assert( (Len_ >= 0) && (Len_ < m_len) );

            if ( (Len_ >= 0) && (Len_ < m_len) )
                return m_data[Len_];
            else if ( Len_ < 0 )
                return m_data[0];
            else // if ( Len_ >= m_len )
                return m_data[m_len - 1];
        };

        // resize the vector
    public:
        const vector& resize(const int Len_)
        {
            assert( 0 != Len_ );

            if ( 0 == Len_ )
                return *this;
            else
                return create(Len_);
        };

        const vector& resize(const int Len_, const T_& Val_)
        {
            assert( 0 != Len_ );

            if ( 0 == Len_ )
                return *this;
            else
                return create(Len_, Val_);
        };

        // clear vector
    public:
        const vector& clear() {freemem(); return *this};

        // overload operators
    public:
        const vector& operator = (const vector& Right_)
        {
            create(Right_.m_len);

            for ( int i = 0; i < m_len; i++ )
                m_data[i] = Right_.m_data[i];

            return *this;
        };

        const vector& operator = (const matrix<T_>& Right_)
        {
            assert( (1 == Right_.rows()) || (1 == Right_.cols()) );

            if ( (1 != Right_.rows()) && (1 != Right_.cols()))
                return *this;

            int i;
            if ( 1 == Right_.rows() )
            {
                resize(Right_.cols());
                for ( i = 0; i < m_len; i++ )
                    m_data[i] = Right_(0, i);
            }
            else
            {
                resize(Right_.rows());
                for ( i = 0; i < m_len; i++ )
                    m_data[i] = Right_(i, 0);
            }

            return *this;
        };

    public:
        const vector& operator += (const vector& Right_)
        {
            assert( m_len == Right_.m_len );

            if ( m_len != Right_.m_len )
                return *this;

            for ( int i = 0; i < m_len; i++ )
                m_data[i] += Right_.m_data[i];

            return *this;
        };

        const vector& operator -= (const vector& Right_)
        {
            assert( m_len == Right_.m_len );

            if ( m_len != Right_.m_len )
                return *this;

            for ( int i = 0; i < m_len; i++ )
                m_data[i] -= Right_.m_data[i];

            return *this;
        };

        const vector& operator -= (const T_& Right_)
        {
            for ( int i = 0; i < m_len; i++ )
                m_data[i] -= Right_;

            return *this;
        };

        const vector& operator *= (const vector& Right_)
        {
            assert( m_len == Right_.m_len );

            if ( m_len != Right_.m_len )
                return *this;

            for ( int i = 0; i < m_len; i++ )
                m_data[i] *= Right_.m_data[i];

            return *this;
        };

        const vector& operator *= (const T_& Right_)
        {
            for ( int i = 0; i < m_len; i++ )
                m_data[i] *= Right_;

            return *this;
        };

        const vector& operator /= (const vector& Right_)
        {
            assert( m_len == Right_.m_len );

            if ( m_len != Right_.m_len )
                return *this;

            for ( int i = 0; i < m_len; i++ )
                m_data[i] /= Right_.m_data[i];

            return *this;
        };

        const vector& operator /= (const T_& Right_)
        {
            for ( int i = 0; i < m_len; i++ )
                m_data[i] /= Right_;

            return *this;
        };

    public:
        friend const vector<T_> operator + (const vector<T_>& Left_, const vector<T_>& Right_)
        {
            assert( Left_.m_len == Right_.m_len );

            vector<T_> rlt;
            if ( Left_.m_len != Right_.m_len )
                return rlt;

            rlt.resize(Left_.m_len);
            for ( int i = 0; i < rlt.m_len; i++ )
                rlt.m_data[i] = Left_.m_data[i] + Right_.m_data[i];

            return rlt;
        };

        friend const vector<T_> operator - (const vector<T_>& Left_, const vector<T_>& Right_)
        {
            assert( Left_.m_len == Right_.m_len );

            vector<T_> rlt;
            if ( Left_.m_len != Right_.m_len )
                return rlt;

            rlt.resize(Left_.m_len);
            for ( int i = 0; i < rlt.m_len; i++ )
                rlt[i] = Left_[i] - Right_[i];

            return rlt;
        };

        friend const T_ operator * (const vector<T_>& Left_, const vector<T_>& Right_)
        {
            assert( Left_.m_len == Right_.m_len );

            T_ rlt = T_();
            if ( Left_.m_len != Right_.m_len )
                return rlt;

            for ( int i = 0; i < Left_.m_len; i++ )
                rlt += Left_.m_data[i] * Right_.m_data[i];

            return rlt;
        };

        friend const vector<T_> operator * (const vector<T_>& Left_, const T_& Right_)
        {
            vector<T_> rlt(Left_.size());
            for ( int i = 0; i < Left_.m_len; i++ )
                rlt[i] = Left_.m_data[i] * Right_;

            return rlt;
        };

    public:
        friend class matrix<T_>;
        friend const vector<T_> operator * (const matrix<T_>& Left_, const vector<T_>& Right_)
        {
            assert( Left_.cols() == Right_.m_len );

            vector<T_> rlt;
            if ( Left_.cols() != Right_.m_len )
                return rlt;

            int i, j;
            rlt.resize(Left_.rows());
            for ( i = 0; i < Left_.rows(); i++ )
            {
                for ( j = 0; j < Left_.cols(); j++ )
                    rlt.m_data[i] += Left_(i, j) * Right_.m_data[j];
            }

            return rlt;
        };

        friend const vector<T_> operator * (const vector<T_>& Left_, const matrix<T_>& Right_)
        {
            assert( Left_.m_len == Right_.rows() );

            vector<T_> rlt;
            if ( Left_.m_len != Right_.rows() )
                return rlt;

            int i, j;
            rlt.resize(Right_.cols());
            for ( j = 0; j < Right_.cols(); j++ )
            {
                for ( i = 0; i < Right_.rows(); i++ )
                    rlt.m_data[j] += Left_.m_data[i] * Right_(i, j);
            }

            return rlt;
        };

        // some math functions
    public:
        const vector& pow(double Expo_)
        {
            assert( 0 != m_len );

            if ( 0 == m_len )
                return *this;

            for ( int i = 0; i < m_len; i++ )
                m_data[i] = pow(m_data[i], Expo_);

            return *this;
        };

        const T_ mod()  // T_ == double
        {
            assert( 0 != m_len );

            if (0 == m_len)
                return 0.0f;

            T_ rlt = 0.0f;
            for (int i = 0; i < m_len; i++)
                rlt += (m_data[i] * m_data[i]);

            rlt = sqrt(double(rlt));
        }

        // get vector's attributes
    public:
        const int   size() const { return m_len; };
        const bool  empty() const { return (NULL == m_data); };

        // set value to vector
    public:
        void set_val(const T_ Val_)
        {
            for ( int i = 0; i < m_len; i++ )
                m_data[i] = Val_;
        }

        // create and free memory for vector
    private:
        const vector& create(const int Len_)
        {
            freemem();

            m_len = Len_;
            m_data = new T_[m_len]();

            return *this;
        };

        const vector& create(const int Len_, const T_& Val_)
        {
            freemem();

            m_len = Len_;
            m_data = new T_[m_len]();
            for ( int i = 0; i < m_len; i++ )
                m_data[i] = Val_;

            return *this;
        };

        void freemem()
        {
            if ( NULL != m_data )
            {
                delete [] m_data;
                m_data = NULL;
                m_len = 0;
            }
        };

        // attributes
    private:
        T_*     m_data;
        int     m_len;
    };
}

#endif // _VECTOR_H_
