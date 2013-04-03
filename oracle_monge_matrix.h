//
//  oracle_monge_matrix.h
//  SubmatrixQueries
//
//  Created by Raphael Bost on 03/04/13.
//  Copyright (c) 2013 Raphael Bost. All rights reserved.
//

#ifndef __SubmatrixQueries__oracle_monge_matrix__
#define __SubmatrixQueries__oracle_monge_matrix__

#include <cstdlib>
#include <algorithm>


#include "matrix.h"

namespace matrix {
    
    class BadTemplateException: public std::exception{
        virtual const char* what() const throw();
    };
    class ImmutableMatrixException: public std::exception
    {
        virtual const char* what() const throw();
    };

    template <typename T>
    class OracleMongeMatrix : public Matrix<T>{
    private:
        size_t _rows, _cols;
        T *_rowSlopes;
        T *_colSlopes;
        
    public:
        OracleMongeMatrix(size_t r, size_t c): _rows(r), _cols(c)
        {
            BadTemplateException ex;
            throw ex;
        }

        OracleMongeMatrix(size_t r, size_t c, T *rowSlopes, T *colSlopes): _rows(r), _cols(c), _rowSlopes(rowSlopes), _colSlopes(colSlopes)
        {}
        
        ~OracleMongeMatrix()
        {
            delete [] _rowSlopes;
            delete [] _colSlopes;
        }
        
        size_t rows() const { return _rows; }
        size_t cols() const { return _cols; }
        
        T operator()(size_t i, size_t j) const
        {
            return (_rowSlopes[i]*j - i) + (_colSlopes[j]*i - j);
        }
        
        T& operator()(size_t i, size_t j)
        {
            ImmutableMatrixException ex;
            throw ex;
        }

    };
    
    template <> OracleMongeMatrix<double>::OracleMongeMatrix(size_t r, size_t c);
}

#endif /* defined(__SubmatrixQueries__oracle_monge_matrix__) */
