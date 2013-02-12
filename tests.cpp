//
//  tests.cpp
//  SubmatrixQueries
//
//  Created by Raphael Bost on 08/02/13.
//  Copyright (c) 2013 Raphael Bost. All rights reserved.
//

#include "tests.h"

#include <cstdlib>
#include <cassert>
#include <ctime>
#include <algorithm>
#include <climits>

#define PRINT_TEST_MATRIX false
#define BENCHMARK true

SubmatrixQueriesTest::SubmatrixQueriesTest(Matrix<double> *m)
{
    assert(m->isInverseMonge());
    // copy the matrix
    _testMatrix = new Matrix<double>(m);
    
#if BENCHMARK
    clock_t time = clock();
#endif
    _queryDS = new SubmatrixQueriesDataStructure<double>(*_testMatrix);
#if BENCHMARK
    time = clock() - time;
    cout << "Building Data Structure: " << 1000*((double)time)/CLOCKS_PER_SEC << " ms" << endl;
#endif
}

SubmatrixQueriesTest::SubmatrixQueriesTest(size_t rows, size_t cols)
{
#if BENCHMARK
    clock_t time = clock();
#endif
    _testMatrix = generateInverseMongeMatrix(rows, cols);
#if BENCHMARK
    time = clock() - time;
    cout << "Building Matrix: " << 1000*((double)time)/CLOCKS_PER_SEC << " ms" << endl;
    time = clock();
#endif
    
    _queryDS = new SubmatrixQueriesDataStructure<double>(*_testMatrix);
#if BENCHMARK
    time = clock() - time;
    cout << "Building Data Structure: " << 1000*((double)time)/CLOCKS_PER_SEC << " ms" << endl;
#endif
    
#if PRINT_TEST_MATRIX
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            cout << (*_testMatrix)(i,j) << " | \t ";
        }
        cout << endl;
    }
#endif
}

SubmatrixQueriesTest::~SubmatrixQueriesTest()
{
    delete _testMatrix;
    delete _queryDS;
}

bool SubmatrixQueriesTest::testColumnQuery(Range rowRange, size_t col)
{
    double naiveMax, queryMax;
    
    queryMax = _queryDS->rowsTree()->maxForColumnInRange(col, rowRange.min, rowRange.max);
    naiveMax = SubmatrixQueriesTest::naiveMaximumInColumn(_testMatrix, rowRange, col);
    
    return queryMax == naiveMax;
}

bool SubmatrixQueriesTest::testColumnQuery()
{
    size_t col = rand() % (_testMatrix->cols());
    
    size_t r1, r2;
    r1 = rand() % (_testMatrix->rows());
    r2 = rand() % (_testMatrix->rows());
    
    Range r = Range(min(r1,r2),max(r1,r2));
    
    return testColumnQuery(r, col);
}

bool SubmatrixQueriesTest::testRowQuery(Range colRange, size_t row)
{
    double naiveMax, queryMax;
    
    queryMax = _queryDS->columnTree()->maxForRowInRange(row, colRange.min, colRange.max);
    naiveMax = SubmatrixQueriesTest::naiveMaximumInRow(_testMatrix, colRange, row);
    
    return queryMax == naiveMax;
}

bool SubmatrixQueriesTest::testRowQuery()
{
    size_t row = rand() % (_testMatrix->rows());
    
    size_t c1,c2;
    c1 = rand() % (_testMatrix->cols());
    c2 = rand() % (_testMatrix->cols());
    
    Range c = Range(min(c1,c2),max(c1,c2));
    
    return testRowQuery(c, row);
}

bool SubmatrixQueriesTest::testSubmatrixQuery(Range rowRange, Range colRange, clock_t *naiveTime, clock_t *queryTime)
{
    double naiveMax, queryMax;
    
    clock_t clock1, clock2, clock3;
    
    clock1 = clock();
    queryMax = _queryDS->maxInRange(rowRange,colRange);
    clock2 = clock();
    naiveMax = SubmatrixQueriesTest::naiveMaximumInSubmatrix(_testMatrix, rowRange, colRange);
    clock3 = clock();
    
    if (naiveTime) {
        *naiveTime += clock3 - clock2;
    }
    if (queryTime) {
        *queryTime += clock2 - clock1;
    }
    
    if (queryMax != naiveMax) {
        cout << "query ranges : row = (" << rowRange.min <<","<<rowRange.max<<")";
        cout << " col = (" << colRange.min <<","<<colRange.max<<")"<<endl;
        cout << "queryMax: " << queryMax;
        cout << " ; naiveMax: " << naiveMax << endl;
        
    }
    return queryMax == naiveMax;    
}

bool SubmatrixQueriesTest::testSubmatrixQuery(clock_t *naiveTime, clock_t *queryTime)
{
    size_t r1, r2;
    r1 = rand() % (_testMatrix->rows());
    r2 = rand() % (_testMatrix->rows());
    
    size_t c1,c2;
    c1 = rand() % (_testMatrix->cols());
    c2 = rand() % (_testMatrix->cols());

    
    Range r = Range(min(r1,r2),max(r1,r2));
    Range c = Range(min(c1,c2),max(c1,c2));
    
    return testSubmatrixQuery(r, c, naiveTime, queryTime);
}

bool SubmatrixQueriesTest::multipleColumnQueryTest(size_t n)
{
    bool result = true;
    
    for (size_t i = 0; i < n && result; i++) {
        result = result && testColumnQuery();
    }
    return result;
}

bool SubmatrixQueriesTest::multipleRowQueryTest(size_t n)
{
    bool result = true;
    
    for (size_t i = 0; i < n && result; i++) {
        result = result && testRowQuery();
    }
    return result;
}

bool SubmatrixQueriesTest::multipleSubmatrixQueryTest(size_t n)
{
    bool result = true;
    clock_t naiveTime = 0, queryTime = 0;
    
    for (size_t i = 0; i < n && result; i++) {
        result = result && testSubmatrixQuery(&naiveTime, &queryTime);
    }
    
#if BENCHMARK
    cout << "Benchmark for " << n << " queries:" <<endl;
    cout << "Naive algorithm: " << 1000*((double)naiveTime)/CLOCKS_PER_SEC << " ms" << endl;
    cout << "Submatrix queries: " << 1000*((double)queryTime)/CLOCKS_PER_SEC << " ms" << endl;
#endif
    return result;
}


double SubmatrixQueriesTest::naiveMaximumInColumn(const Matrix<double> *m, Range rowRange, size_t col)
{
    MaxValue<double> max = MaxValue<double>();
    
    for (size_t i = rowRange.min; i <= rowRange.max; i++) {
        max.updateMax((*m)(i,col));
    }
    
    return max.value();
}

double SubmatrixQueriesTest::naiveMaximumInRow(const Matrix<double> *m, Range colRange, size_t row)
{
    MaxValue<double> max = MaxValue<double>();
    
    for (size_t j = colRange.min; j <= colRange.max; j++) {
        max.updateMax((*m)(row,j));
    }
    
    return max.value();
}

double SubmatrixQueriesTest::naiveMaximumInSubmatrix(const Matrix<double> *m, Range rowRange, Range colRange)
{
    MaxValue<double> max = MaxValue<double>();
    
    for (size_t j = colRange.min; j <= colRange.max; j++) {
        max.updateMax(naiveMaximumInColumn(m, rowRange, j));
    }
    
    return max.value();
}

#define LINE_DISTANCE 5

Matrix<double>* SubmatrixQueriesTest::generateInverseMongeMatrix(size_t rows, size_t cols)
{
    Matrix<double> *m = new Matrix<double>(rows,cols);
    
    int *rowsAbscissa, *colsAbscissa;
    
    rowsAbscissa = new int[rows];
    colsAbscissa = new int[cols];
    
    srand ( time(NULL) );
    
    // we define these values to avoid overflows that will lead to a non inverse Monge matrix 
    int max_abscissa = (int)sqrtf(INT_MAX/3) - LINE_DISTANCE;
    int rowInterval = max<int>(max_abscissa/rows,1);
    int colInterval = max<int>(max_abscissa/cols,1);
    
    assert(rowInterval > 0);
    assert(colInterval > 0);
    
    int accumulator =0;
    
    for (size_t i = 0; i < rows; i++) {
        accumulator += rand() % rowInterval;
        rowsAbscissa[i] = accumulator;
    }
//    cout << "row max : " << accumulator << endl;
    
    accumulator = 0;
    for (size_t j = 0; j < cols; j++) {
        accumulator += rand() % colInterval;
        colsAbscissa[cols-1-j] = accumulator;
    }
//    cout << "col max : " << accumulator << endl;

    
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            double abscissaDiff = rowsAbscissa[i]-colsAbscissa[j];

            (*m)(i,j) = sqrt(LINE_DISTANCE + abscissaDiff*abscissaDiff);
        }
    }
    
    delete [] rowsAbscissa;
    delete [] colsAbscissa;
    
    return m;
}