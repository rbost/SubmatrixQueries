//
//  tests.h
//  SubmatrixQueries
//
//  Created by Raphael Bost on 08/02/13.
//  Copyright (c) 2013 Raphael Bost. All rights reserved.
//

#ifndef __SubmatrixQueries__tests__
#define __SubmatrixQueries__tests__

#include <iostream>
#include "matrix.h"
#include "range.h"
#include "envelope_tree.h"

using namespace matrix;

/*
 * class SubmatrixQueriesTest
 *
 * This is a simple class to test the data structures we built on big matrices.
 * Naive maximum search functions are provided and their results are compared 
 * to the ones returned by our data structures.
 *
 * We also provide a function to generate inverse Monge matrices.
 */

class SubmatrixQueriesTest {
    const Matrix<double> *_testMatrix;
    
    SubmatrixQueriesDataStructure<double> *_queryDS;
    
public:
    SubmatrixQueriesTest(Matrix<double> *m);
    SubmatrixQueriesTest(size_t rows, size_t cols);
    ~SubmatrixQueriesTest();

    
    bool testColumnQuery(Range rowRange, size_t col, clock_t *naiveTime, clock_t *queryTime);
    bool testColumnQuery(clock_t *naiveTime, clock_t *queryTime);
    bool testRowQuery(Range colRange, size_t row, clock_t *naiveTime, clock_t *queryTime);
    bool testRowQuery(clock_t *naiveTime, clock_t *queryTime);
    bool testSubmatrixQuery(Range rowRange, Range colRange, clock_t *naiveTime, clock_t *queryTime);
    bool testSubmatrixQuery(clock_t *naiveTime, clock_t *queryTime);
    
    bool multipleColumnQueryTest(size_t n);
    bool multipleRowQueryTest(size_t n);
    bool multipleSubmatrixQueryTest(size_t n);
    
    static double naiveMaximumInColumn(const Matrix<double> *m, Range rowRange, size_t col);
    static double naiveMaximumInRow(const Matrix<double> *m, Range colRange, size_t row);
    static double naiveMaximumInSubmatrix(const Matrix<double> *m, Range rowRange, Range colRange);
    
    static Matrix<double>* generateInverseMongeMatrix(size_t rows, size_t cols);
    static Matrix<double>* generateInverseMongeMatrix2(size_t rows, size_t cols);
    static Matrix<double>* generateInverseMongeMatrix3(size_t rows, size_t cols);
};

#endif /* defined(__SubmatrixQueries__tests__) */
