//
//  oracle_monge_matrix.cpp
//  SubmatrixQueries
//
//  Created by Raphael Bost on 03/04/13.
//  Copyright (c) 2013 Raphael Bost. All rights reserved.
//

#include "oracle_monge_matrix.h"

namespace matrix
{
    double dRand(double fMin, double fMax)
    {
        double f = (double)rand() / RAND_MAX;
        return fMin + f * (fMax - fMin);
    };
    
    const char* BadTemplateException::what() const throw()
    {
        return "Invalid template exception";
    }
    const char* ImmutableMatrixException::what() const throw()
    {
        return "Immutable matrix exception";
    }

    template <> OracleMongeMatrix<double>::OracleMongeMatrix(size_t r, size_t c): _rows(r), _cols(c)
    {
        _rowSlopes = new double[r];
        _colSlopes = new double[c];
        
        srand ( time(NULL) );
        
        for (size_t i = 0; i < r; i++) {
            _rowSlopes[i] = tan(dRand(-0.5*M_PI, 0.5*M_PI));
        }
        std::sort(_rowSlopes,_rowSlopes+r);
        
        for (size_t i = 0; i < c; i++) {
            _colSlopes[i] = tan(dRand(-0.5*M_PI, 0.5*M_PI));
        }
        std::sort(_colSlopes,_colSlopes+r);
        
    }
    

}