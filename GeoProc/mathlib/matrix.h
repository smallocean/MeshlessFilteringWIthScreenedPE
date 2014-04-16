//
// matrix.h
//

/************************************************************
*
*   Copyritht (C) 2005    PAMI Lab, SJTU
*   All rights reserved
*
*   filename: matrix.h
*   abstract: implement a kind of matrix class
*
*   coder: DeepBlue
*   date: 2005/08/18
*
************************************************************/

#pragma once

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <cassert>
#include <cstdlib>
#include <cstdio>
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

    template <typename T_> class matrix;
    template <typename T_> const matrix<T_> operator + (const matrix<T_>& Left_, const matrix<T_>& Right_);
    template <typename T_> const matrix<T_> operator - (const matrix<T_>& Left_, const matrix<T_>& Right_);
    template <typename T_> const matrix<T_> operator * (const matrix<T_>& Left_, const matrix<T_>& Right_);

    template <typename T_> const matrix<T_> operator * (const matrix<T_>& Left_, const T_& Right_);
    template <typename T_> const matrix<T_> operator * (const T_& Left_, const matrix<T_>& Right_);
    template <typename T_> const matrix<T_> operator / (const matrix<T_>& Left_, const T_& Right_);

    // declare matrix class
    template <typename T_> class vector;
    template <typename T_> 
    class matrix    
    {
        // constructors
    public:
        matrix()
            : m_data(NULL)
            , m_rows(0)
            , m_cols(0)
        {    };

        matrix(const matrix& rhs)
            : m_data(NULL)
            , m_rows(0)
            , m_cols(0)
        { *this = rhs; };

        matrix(const vector<T_>& rhs)
            : m_data(NULL)
            , m_rows(0)
            , m_cols(0)
        { *this = rhs; };

        explicit matrix(const int Rows_, const int Cols_)
            : m_data(new T_*[Rows_])
            , m_rows(Rows_)
            , m_cols(Cols_)
        {
            m_data[0] = new T_[m_rows * m_cols]();
            for ( int i = 1; i < m_rows; i++ )
                m_data[i] = m_data[i - 1] + m_cols;
        };

        explicit matrix(const int Rows_, const int Cols_, const T_& Val_)
            : m_data(new T_*[Rows_])
            , m_rows(Rows_)
            , m_cols(Cols_)
        {
            int len = m_rows * m_cols;
            m_data[0] = new T_[len]();

            int i, j;
            for ( i = 1; i < m_rows; i++ )
                m_data[i] = m_data[i - 1] + m_cols;

            for ( j = 0; j < len; j++ )
                m_data[0][j] = Val_;
        };

        // destructor
    public:
        ~matrix()
        { freemem(); };

        // get element from matrix
    public:
        T_& operator () (/*const*/ int Row_, /*const*/ int Col_)
        {
            //assert( (Row_ < m_rows) && (Col_ < m_cols) );

            if ( Row_ < 0 )
                Row_ = 0;
            else if ( Row_ >= m_rows )
                Row_ = m_rows - 1;

            if ( Col_ < 0 )
                Col_ = 0;
            else if ( Col_ >= m_cols )
                Col_ = m_cols - 1;

            return m_data[Row_][Col_];
        };

        const T_& operator () (/*const*/ int Row_, /*const*/ int Col_) const
        {
            //assert( (Row_ < m_rows) && (Col_ < m_cols) );

            if ( Row_ < 0 )
                Row_ = 0;
            else if ( Row_ >= m_rows )
                Row_ = m_rows - 1;

            if ( Col_ < 0 )
                Col_ = 0;
            else if ( Col_ >= m_cols )
                Col_ = m_cols - 1;

            return m_data[Row_][Col_];
        };

        T_* operator [] (/*const*/ int Row_)
        {
            assert((Row_ >= 0) && (Row_ < m_rows));
            if ( Row_ < 0 )
                Row_ = 0;
            else if ( Row_ >= m_rows )
                Row_ = m_rows - 1;

            return m_data[Row_];
        };

        const T_* operator [] (/*const*/ int Row_) const
        {
            assert((Row_ >= 0) && (Row_ < m_rows));
            if ( Row_ < 0 )
                Row_ = 0;
            else if ( Row_ >= m_rows )
                Row_ = m_rows - 1;

            return m_data[Row_];
        };

        // get row or column
    public:
        vector<T_> get_row(const int Pos_) const
        {
            vector<T_> rlt;

            assert( Pos_ < m_rows );
            if ( Pos_ >= m_rows )
                return rlt;

            rlt.resize(m_cols);
            for ( int i = 0; i < m_cols; i++ )
                rlt[i] = m_data[Pos_][i];

            return rlt;
        };

        vector<T_> get_col(const int Pos_) const
        {
            vector<T_> rlt;

            assert( Pos_ < m_cols );
            if ( Pos_ >=  m_cols)
                return rlt;

            rlt.resize(m_rows);
            for ( int i = 0; i < m_rows; i++ )
                rlt[i] = m_data[i][Pos_];

            return rlt;
        };

        // resize the matrix
    public:
        const matrix& resize(const int Rows_, const int Cols_)
        {
            //assert( 0 != Rows_ * Cols_ );

            if ( 0 == Rows_ * Cols_ )
                return this->clear();
                //return *this;
            else
                return create(Rows_, Cols_);
        };

        const matrix& resize(const int Rows_, const int Cols_, const T_& Val_)
        {
            //assert( 0 != Rows_ * Cols_ );

            if ( 0 == Rows_ * Cols_ )
                return this->clear();
                //return *this;
            else
                return create(Rows_, Cols_, Val_);
        };

        // clear matrix
    public:
        const matrix& clear() {freemem(); return *this;};

        // overload operators
    public:
        const matrix& operator = (const matrix& Right_)
        {
            if ( 0 == Right_.size() )
            {
                this->freemem();
                return *this;
            }

            if ((m_rows != Right_.m_rows) || (m_cols != Right_.m_cols))
                create(Right_.m_rows, Right_.m_cols);

            int i, j;
            for ( i = 0; i < m_rows; i++ )
            {
                for ( j = 0; j < m_cols; j++ )
                    m_data[i][j] = Right_.m_data[i][j];
            }

            return *this;
        };

        const matrix& operator = (const vector<T_>& Right_)
        {
            assert( 0 != Right_.size() );

            if ( 0 == Right_.size() )
            {
                this->freemem();
                return *this;
            }

            if ((m_rows != Right_.size()) || (m_cols != 1))
                resize(Right_.size(), 1);
            for ( int i = 0; i < m_rows; i++ )
            {
                m_data[i][0] = Right_[i];
            }

            return *this;
        };

    public:
        const matrix& operator += (const matrix& Right_)
        {
            assert( (m_rows == Right_.m_rows) && (m_cols == Right_.m_cols) );

            if ( (m_rows != Right_.m_rows) || (m_cols != Right_.m_cols) )
                return *this;

            int i, j;
            for ( i = 0; i < m_rows; i++ )
            {
                for ( j = 0; j < m_cols; j++ )
                    m_data[i][j] += Right_.m_data[i][j];
            }

            return *this;
        };

        const matrix& operator -= (const matrix& Right_)
        {
            assert( (m_rows == Right_.m_rows) && (m_cols == Right_.m_cols) );

            if ( (m_rows != Right_.m_rows) || (m_cols != Right_.m_cols) )
                return *this;

            int i, j;
            for ( i = 0; i < m_rows; i++ )
            {
                for ( j = 0; j < m_cols; j++ )
                    m_data[i][j] -= Right_.m_data[i][j];
            }

            return *this;
        };

        const matrix& operator *= (const matrix& Right_)
        {
            assert( (m_rows == Right_.m_rows) && (m_cols == Right_.m_cols) );

            if ( (m_rows != Right_.m_rows) || (m_cols != Right_.m_cols) )
                return *this;

            int i, j;
            for ( i = 0; i < m_rows; i++ )
            {
                for ( j = 0; j < m_cols; j++ )
                    m_data[i][j] *= Right_.m_data[i][j];
            }

            return *this;
        };

        const matrix& operator /= (const matrix& Right_)
        {
            assert( (m_rows == Right_.m_rows) && (m_cols == Right_.m_cols) );

            if ( (m_rows != Right_.m_rows) || (m_cols != Right_.m_cols) )
                return *this;

            int i, j;
            for ( i = 0; i < m_rows; i++ )
            {
                for ( j = 0; j < m_cols; j++ )
                    m_data[i][j] /= Right_.m_data[i][j];
            }

            return *this;
        };

    public:
        const matrix& operator += (const T_& Right_)
        {
            int i, j;
            for ( i = 0; i < m_rows; i++ )
            {
                for ( j = 0; j < m_cols; j++ )
                    m_data[i][j] += Right_;
            }
            return *this;
        }

        const matrix& operator -= (const T_& Right_)
        {
            int i, j;
            for ( i = 0; i < m_rows; i++ )
            {
                for ( j = 0; j < m_cols; j++ )
                    m_data[i][j] -= Right_;
            }
            return *this;
        }

        const matrix& operator *= (const T_& Right_)
        {
            int i, j;
            for ( i = 0; i < m_rows; i++ )
            {
                for ( j = 0; j < m_cols; j++ )
                    m_data[i][j] *= Right_;
            }
            return *this;
        }

        const matrix& operator /= (const T_& Right_)
        {
            int i, j;
            for ( i = 0; i < m_rows; i++ )
            {
                for ( j = 0; j < m_cols; j++ )
                    m_data[i][j] /= Right_;
            }
            return *this;
        }

    public:
        friend const matrix<T_> operator + (const matrix<T_>& Left_, const matrix<T_>& Right_)
        {
            assert( (Left_.m_rows == Right_.m_rows) && (Left_.m_cols == Right_.m_cols) );

            matrix<T_> rlt;
            if ( (Left_.m_rows != Right_.m_rows) || (Left_.m_cols != Right_.m_cols) )
                return rlt;

            int i, j;
            rlt.resize(Left_.m_rows, Left_.m_cols);
            for ( i = 0; i < rlt.m_rows; i++ )
            {
                for ( j = 0; j < rlt.m_cols; j++ )
                    rlt.m_data[i][j] = Left_.m_data[i][j] + Right_.m_data[i][j];
            }

            return rlt;
        };

        friend const matrix<T_> operator - (const matrix<T_>& Left_, const matrix<T_>& Right_)
        {
            assert( (Left_.m_rows == Right_.m_rows) && (Left_.m_cols == Right_.m_cols) );

            matrix<T_> rlt;
            if ( (Left_.m_rows != Right_.m_rows) || (Left_.m_cols != Right_.m_cols) )
                return rlt;

            int i, j;
            rlt.resize(Left_.m_rows, Left_.m_cols);
            for ( i = 0; i < rlt.m_rows; i++ )
            {
                for ( j = 0; j < rlt.m_cols; j++ )
                    rlt.m_data[i][j] = Left_.m_data[i][j] - Right_.m_data[i][j];
            }

            return rlt;
        };

        friend const matrix<T_> operator * (const matrix<T_>& Left_, const matrix<T_>& Right_)
        {
            assert( Left_.m_cols == Right_.m_rows );

            matrix<T_> rlt;
            if ( Left_.m_cols != Right_.m_rows )
                return rlt;

            int i, j, k;
            rlt.resize(Left_.m_rows, Right_.m_cols);
            for ( i = 0; i < rlt.m_rows; i++ )
            {
                for ( j = 0; j < rlt.m_cols; j++ )
                {
                    rlt.m_data[i][j] = T_();
                    for ( k = 0; k < Left_.m_cols; k++ )
                    {
                        rlt.m_data[i][j] += (Left_.m_data[i][k] * Right_.m_data[k][j]);
                    }
                }
            }

            return rlt;
        };

    public:
        friend const matrix<T_> operator + (const matrix<T_>& Left_, const T_& Right_)
        {
            matrix<T_> rlt(Left_);
            int i, j;
            for ( i = 0; i < rlt.m_rows; i++ )
            {
                for ( j = 0; j < rlt.m_cols; j++ )
                    rlt.m_data[i][j] += Right_;
            }
            return rlt;
        };

        friend const matrix<T_> operator + (const T_& Left_, const matrix<T_>& Right_)
        { return (Right_ + Left); };

        friend const matrix<T_> operator - (const matrix<T_>& Left_, const T_& Right_)
        {
            matrix<T_> rlt(Left_);
            int i, j;
            for ( i = 0; i < rlt.m_rows; i++ )
            {
                for ( j = 0; j < rlt.m_cols; j++ )
                    rlt.m_data[i][j] -= Right_;
            }
            return rlt;
        };

        friend const matrix<T_> operator - (const T_& Left_, const matrix<T_>& Right_)
        {
            matrix<T_> rlt(Left_);
            int i, j;
            for ( i = 0; i < rlt.m_rows; i++ )
            {
                for ( j = 0; j < rlt.m_cols; j++ )
                    rlt.m_data[i][j] = Right_ - rlt.m_data[i][j];
            }
            return rlt;
        };

        friend const matrix<T_> operator * (const matrix<T_>& Left_, const T_& Right_)
        {
            matrix<T_> rlt(Left_);
            int i, j;
            for ( i = 0; i < rlt.m_rows; i++ )
            {
                for ( j = 0; j < rlt.m_cols; j++ )
                    rlt.m_data[i][j] *= Right_;
            }
            return rlt;
        };

        friend const matrix<T_> operator * (const T_& Left_, const matrix<T_>& Right_)
        { return (Right_ * Left_); };

        friend const matrix<T_> operator / (const matrix<T_>& Left_, const T_& Right_)
        {
            matrix<T_> rlt(Left_);
            int i, j;
            for ( i = 0; i < rlt.m_rows; i++ )
            {
                for ( j = 0; j < rlt.m_cols; j++ )
                    rlt.m_data[i][j] /= Right_;
            }
            return rlt;
        };
        friend const matrix<T_> operator / (const T_& Left_, const matrix<T_>& Right_)
        {
            matrix<T_> rlt(Left_);
            int i, j;
            for ( i = 0; i < rlt.m_rows; i++ )
            {
                for ( j = 0; j < rlt.m_cols; j++ )
                    rlt.m_data[i][j] = Right_ / rlt.m_data[i][j];
            }
            return rlt;
        };



        // some math functions
    public:
        const matrix& pow(double Expo_)
        {
            assert( 0 != m_rows * m_cols );

            int len = m_rows * m_cols;
            if ( 0 == len )
                return *this;

            for ( int i = 0; i < len; i++ )
                m_data[0][i] = pow(m_data[0][i], Expo_);

            return *this;
        };

    public:
        const matrix& transpose()
        {
            // copy this matrix, and recreate this matrix
            matrix temp(*this);
            create(m_cols, m_rows);

            int i, j;
            for ( i = 0; i < m_rows; i++ )
            {
                for ( j = 0; j < m_cols; j++ )
                    m_data[i][j] = temp.m_data[j][i];
            }

            return *this;
        };

    public:
        const T_ trace()
        {
            assert( m_rows == m_cols );

            if ( m_rows != m_cols )
                return T_();

            T_    rlt = m_data[0][0];
            for ( int i = 1; i < m_rows; i++ )
                rlt += m_data[i][i];

            return rlt;
        };

        // get matrix's attributes
    public:
        const int   rows() const { return m_rows; };
        const int   cols() const { return m_cols; };
        const int   size() const { return (m_rows * m_cols); };
        const bool  empty() const { return (NULL == m_data); };

        // get matrix's minimum and maximum value
    public:
        const T_ min_val()
        {
            T_ val = m_dada[0];
            int size = m_rows * m_cols;
            for (int i = 1; i < size; i++)
            {
                if (val > m_data[i])
                    val = m_data[i];
            }

            return val;
        }

        const T_ max_val()
        {
            T_ val = m_dada[0];
            int size = m_rows * m_cols;
            for (int i = 1; i < size; i++)
            {
                if (val < m_data[i])
                    val = m_data[i];
            }

            return val;
        }

        // set value to matrix
    public:
        void set_val(const T_ Val_)
        {
            int i, j;
            for ( i = 0; i < m_rows; i++ )
            {
                for ( j = 0; j < m_cols; j++ )
                    m_data[i][j] = Val_;
            }
        }

        // create and free memory for matrix
    private:
        const matrix& create(const int Rows_, const int Cols_)
        {
            freemem();

            int len = Rows_ * Cols_;
            if ( 0 == len )
                return *this;

            int i;
            m_rows = Rows_;
            m_cols = Cols_;
            m_data = new T_*[m_rows];
            m_data[0] = new T_[len]();
            for ( i = 1; i < m_rows; i++ )
                m_data[i] = m_data[i - 1] + m_cols;

            return *this;
        };

        const matrix& create(const int Rows_, const int Cols_, const T_& Val_)
        {
            freemem();

            int len = Rows_ * Cols_;
            if ( 0 == len )
                return *this;

            int i;
            m_rows = Rows_;
            m_cols = Cols_;
            m_data = new T_*[m_rows];
            m_data[0] = new T_[len]();
            for ( i = 1; i < m_rows; i++ )
            {
                m_data[i] = m_data[i - 1] + m_cols;
            }
            for ( i = 0; i < len; i++ )
            {
                m_data[0][i] = Val_;
            }

            return *this;
        };

        void freemem()
        {
            if ( NULL != m_data )
            {
                delete [] m_data[0];
                delete [] m_data;

                m_data = NULL;
                m_rows = 0;
                m_cols = 0;
            }
        };

        // attributes
    private:
        T_**    m_data;
        int     m_rows;
        int     m_cols;
    };
}

#endif // _MATRIX_H_
