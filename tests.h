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

#ifdef __MACH__
    typedef clock_t bench_time_t;
#else
    typedef timespec bench_time_t;
#endif
double benchTimeAsMiliSeconds(bench_time_t t);
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

    
    bool testColumnQuery(Range rowRange, size_t col, bench_time_t *naiveTime, bench_time_t *queryTime);
    bool testColumnQuery(bench_time_t *naiveTime, bench_time_t *queryTime);
    bool testCascadingColQuery(Range rowRange, size_t col, bench_time_t *queryTime, bench_time_t *cascadingTime, bench_time_t *simpleCascadingTime);
    bool testCascadingColQuery(bench_time_t *queryTime, bench_time_t *cascadingTime, bench_time_t *simpleCascadingTime);
    bool testRowQuery(Range colRange, size_t row, bench_time_t *naiveTime, bench_time_t *queryTime);
    bool testRowQuery(bench_time_t *naiveTime, bench_time_t *queryTime);
    bool testCascadingRowQuery(Range colRange, size_t row, bench_time_t *queryTime, bench_time_t *cascadingTime, bench_time_t *simpleCascadingTime);
    bool testCascadingRowQuery(bench_time_t *queryTime, bench_time_t *cascadingTime, bench_time_t *simpleCascadingTime);
    bool testSubmatrixQuery(Range rowRange, Range colRange, bench_time_t *naiveTime, bench_time_t *queryTime);
    bool testSubmatrixQuery(bench_time_t *naiveTime, bench_time_t *queryTime);
    
    void benchmarkAllRowQueries(Range colRange, size_t row, bench_time_t *naiveTime, bench_time_t *queryTime, bench_time_t *cascadingTime, bench_time_t *simpleCascadingTime);
    void benchmarkAllRowQueries(bench_time_t *naiveTime, bench_time_t *queryTime, bench_time_t *cascadingTime, bench_time_t *simpleCascadingTime);
    void benchmarkAllColQueries(Range rowRange, size_t col, bench_time_t *naiveTime, bench_time_t *queryTime, bench_time_t *cascadingTime, bench_time_t *simpleCascadingTime);
    void benchmarkAllColQueries(bench_time_t *naiveTime, bench_time_t *queryTime, bench_time_t *cascadingTime, bench_time_t *simpleCascadingTime);

    void benchmarkAllSubmatrixQueries(Range rowRange, Range colRange, bench_time_t *naiveTime, bench_time_t *explicitNodesTime, bench_time_t *implicitNodesTime);
    void benchmarkAllSubmatrixQueries(bench_time_t *naiveTime, bench_time_t *explicitNodesTime, bench_time_t *implicitNodesTime);

    static void multiBenchmarksPositionQueries(size_t maxNRows, size_t maxNCols, size_t nSamples, bench_time_t *naiveTime, bench_time_t *queryTime, bench_time_t *cascadingTime, bench_time_t *simpleCascadingTime);
    static bench_time_t** multiSizeBenchmarksPositionQueries(size_t maxNRows, size_t maxNCols, size_t nSampleSize, size_t nSamplePerSize);

    static void multiBenchmarksSubmatrixQueries(size_t nRows, size_t nCols, size_t nSamples, bench_time_t *naiveTime, bench_time_t *explicitNodesTime, bench_time_t *implicitNodesTime);
    static bench_time_t** multiSizeBenchmarksSubmatrixQueries(size_t maxNRows, size_t maxNCols, size_t nSampleSize, size_t nSamplePerSize);

    bool multipleColumnQueryTest(size_t n);
    bool multipleRowQueryTest(size_t n);
    bool multipleSubmatrixQueryTest(size_t n);
    bool multipleRowQueryTestVsCascading(size_t n);
    bool multipleColQueryTestVsCascading(size_t n);
    
    void multipleBenchmarksRowQueries(size_t n);
    void multipleBenchmarksRowQueries(size_t n, bench_time_t *naiveTime, bench_time_t *queryTime, bench_time_t *cascadingTime, bench_time_t *simpleCascadingTime);
    void multipleBenchmarksColQueries(size_t n);
    void multipleBenchmarksColQueries(size_t n, bench_time_t *naiveTime, bench_time_t *queryTime, bench_time_t *cascadingTime, bench_time_t *simpleCascadingTime);
    void multipleBenchmarksSubmatrixQueries(size_t n);
    void multipleBenchmarksSubmatrixQueries(size_t n,bench_time_t *naiveTime, bench_time_t *explicitNodesTime, bench_time_t *implicitNodesTime);

    static double naiveMaximumInColumn(const Matrix<double> *m, Range rowRange, size_t col);
    static double naiveMaximumInRow(const Matrix<double> *m, Range colRange, size_t row);
    static double naiveMaximumInSubmatrix(const Matrix<double> *m, Range rowRange, Range colRange);
    
    static Matrix<double>* generateInverseMongeMatrixStrip1(size_t rows, size_t cols);
    static Matrix<double>* generateInverseMongeMatrixStrip2(size_t rows, size_t cols);
    static Matrix<double>* generateInverseMongeMatrixSlope(size_t rows, size_t cols);
};

#endif /* defined(__SubmatrixQueries__tests__) */
