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
#include <set>
#include <cfloat>


double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}


#define PRINT_TEST_MATRIX false
#define BENCHMARK false

SubmatrixQueriesTest::SubmatrixQueriesTest(Matrix<double> *m)
{
    assert(m->isInverseMonge());
    // copy the matrix
    _testMatrix = new ComplexMatrix<double>(m);
    
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
    _testMatrix = generateInverseMongeMatrixSlope(rows, cols);
#if BENCHMARK
    time = clock() - time;
    cout << "Building Matrix: " << 1000*((double)time)/CLOCKS_PER_SEC << " ms" << endl;
    time = clock();
#endif
    assert(_testMatrix->isInverseMonge());
    
    _queryDS = new SubmatrixQueriesDataStructure<double>(*_testMatrix);
#if BENCHMARK
    time = clock() - time;
    cout << "Building Data Structure: " << 1000*((double)time)/CLOCKS_PER_SEC << " ms" << endl;
    cout << "Max/min envelope size for rows: " << _queryDS->rowsTree()->maxEnvelopeSize() << " / "<< _queryDS->rowsTree()->minEnvelopeSize() << endl;
    cout << "Max/min envelope size for cols: " << _queryDS->columnTree()->maxEnvelopeSize()  << " / "<< _queryDS->columnTree()->minEnvelopeSize() << endl;
#endif
    
#if PRINT_TEST_MATRIX
    cout << "\n";
    _testMatrix->print();
#endif
#if BENCHMARK
    cout << "\n";
#endif
}

SubmatrixQueriesTest::~SubmatrixQueriesTest()
{
    delete _testMatrix;
    delete _queryDS;
}

bool SubmatrixQueriesTest::testColumnQuery(Range rowRange, size_t col, clock_t *naiveTime, clock_t *queryTime)
{
    double naiveMax, queryMax;
    clock_t clock1, clock2, clock3;

    clock1 = clock();
    queryMax = _queryDS->rowsTree()->maxForColumnInRange(col, rowRange.min, rowRange.max);
    clock2 = clock();
    naiveMax = SubmatrixQueriesTest::naiveMaximumInColumn(_testMatrix, rowRange, col);
    
    clock3 = clock();
    
    if (naiveTime) {
        *naiveTime += clock3 - clock2;
    }
    if (queryTime) {
        *queryTime += clock2 - clock1;
    }
    
    return queryMax == naiveMax;
}

bool SubmatrixQueriesTest::testColumnQuery(clock_t *naiveTime, clock_t *queryTime)
{
    size_t col = rand() % (_testMatrix->cols());
    
    size_t r1, r2;
    r1 = rand() % (_testMatrix->rows());
    r2 = rand() % (_testMatrix->rows());
    
    Range r = Range(min(r1,r2),max(r1,r2));
    
    return testColumnQuery(r, col, naiveTime, queryTime);
}

bool SubmatrixQueriesTest::testCascadingColQuery(Range rowRange, size_t col, clock_t *queryTime, clock_t *cascadingTime, clock_t *simpleCascadingTime)
{
    double cascadingMax, simpleCascadingMax, queryMax;
    clock_t clock1, clock2, clock3, clock4;
    
    clock1 = clock();
    cascadingMax = _queryDS->rowsTree()->cascadingMaxInRange(col, rowRange);
    
    clock2 = clock();
    simpleCascadingMax = _queryDS->rowsTree()->simpleCascadingMaxInRange(col, rowRange);
    
    clock3 = clock();
    queryMax = _queryDS->rowsTree()->maxForColumnInRange(col, rowRange.min, rowRange.max);
    
    clock4 = clock();
    
    if (queryTime) {
        *queryTime += clock4 - clock3;
    }
    if (cascadingTime) {
        *cascadingTime += clock2 - clock1;
    }
    if (simpleCascadingTime) {
        *simpleCascadingTime += clock3 - clock2;
    }
    
    if (queryMax != cascadingMax) {
        cout << "Test failed: "<<endl;
        cout << "\tquery ranges : row = (" << rowRange.min <<","<<rowRange.max<<")";
        cout << " col = " << col <<endl;
        cout << "\tqueryMax: " << queryMax;
        cout << " ; cascading: " << cascadingMax << endl;
        
    }

    return queryMax == cascadingMax;
}

bool SubmatrixQueriesTest::testCascadingColQuery(clock_t *queryTime, clock_t *cascadingTime, clock_t *simpleCascadingTime)
{
    size_t col = rand() % (_testMatrix->cols());
    
    size_t r1,r2;
    r1 = rand() % (_testMatrix->rows());
    r2 = rand() % (_testMatrix->rows());
    
    Range r = Range(min(r1,r2),max(r1,r2));
    
    return testCascadingColQuery(r, col,queryTime,cascadingTime,simpleCascadingTime);
}


bool SubmatrixQueriesTest::testRowQuery(Range colRange, size_t row, clock_t *naiveTime, clock_t *queryTime)
{
    double naiveMax, queryMax;
    clock_t clock1, clock2, clock3;
    
    clock1 = clock();
    queryMax = _queryDS->columnTree()->maxForRowInRange(row, colRange.min, colRange.max);
    
    clock2 = clock();
    naiveMax = SubmatrixQueriesTest::naiveMaximumInRow(_testMatrix, colRange, row);
    
    clock3 = clock();

    if (naiveTime) {
        *naiveTime += clock3 - clock2;
    }
    if (queryTime) {
        *queryTime += clock2 - clock1;
    }

    return queryMax == naiveMax;
}

bool SubmatrixQueriesTest::testRowQuery(clock_t *naiveTime, clock_t *queryTime)
{
    size_t row = rand() % (_testMatrix->rows());
    
    size_t c1,c2;
    c1 = rand() % (_testMatrix->cols());
    c2 = rand() % (_testMatrix->cols());
    
    Range c = Range(min(c1,c2),max(c1,c2));
    
    return testRowQuery(c, row,naiveTime,queryTime);
}

bool SubmatrixQueriesTest::testCascadingRowQuery(Range colRange, size_t row, clock_t *queryTime, clock_t *cascadingTime, clock_t *simpleCascadingTime)
{
    double cascadingMax, queryMax, simpleCascadingMax;
    clock_t clock1, clock2, clock3, clock4;
    
    clock1 = clock();
    cascadingMax = _queryDS->columnTree()->cascadingMaxInRange(row, colRange);
    
    clock2 = clock();
    simpleCascadingMax = _queryDS->columnTree()->simpleCascadingMaxInRange(row, colRange);

    clock3 = clock();
    queryMax = _queryDS->columnTree()->maxForRowInRange(row, colRange.min, colRange.max);

    clock4 = clock();

    if (queryTime) {
        *queryTime += clock4 - clock3;
    }
    if (cascadingTime) {
        *cascadingTime += clock2 - clock1;
    }
    if (simpleCascadingTime) {
        *simpleCascadingTime += clock3 - clock2;
    }

    return queryMax == cascadingMax;
}

bool SubmatrixQueriesTest::testCascadingRowQuery(clock_t *queryTime, clock_t *cascadingTime, clock_t *simpleCascadingTime)
{
    size_t row = rand() % (_testMatrix->rows());
    
    size_t c1,c2;
    c1 = rand() % (_testMatrix->cols());
    c2 = rand() % (_testMatrix->cols());
    
    Range c = Range(min(c1,c2),max(c1,c2));
    
    return testCascadingRowQuery(c, row,queryTime,cascadingTime,simpleCascadingTime);
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
        cout << "Test failed: "<<endl;
        cout << "\tquery ranges : row = (" << rowRange.min <<","<<rowRange.max<<")";
        cout << " col = (" << colRange.min <<","<<colRange.max<<")"<<endl;
        cout << "\tqueryMax: " << queryMax;
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

void SubmatrixQueriesTest::benchmarkAllRowQueries(Range colRange, size_t row, clock_t *naiveTime, clock_t *queryTime, clock_t *cascadingTime, clock_t *simpleCascadingTime)
{
    clock_t clock1, clock2, clock3, clock4, clock5;
    
    clock1 = clock();
    _queryDS->columnTree()->cascadingMaxInRange(row, colRange);
    
    clock2 = clock();
    _queryDS->columnTree()->simpleCascadingMaxInRange(row, colRange);

    clock3 = clock();
     _queryDS->columnTree()->maxForRowInRange(row, colRange.min, colRange.max);
    
    clock4 = clock();
    SubmatrixQueriesTest::naiveMaximumInRow(_testMatrix, colRange, row);

    clock5 = clock();

    if (queryTime) {
        *queryTime += clock4 - clock3;
    }
    if (cascadingTime) {
        *cascadingTime += clock2 - clock1;
    }
    if (simpleCascadingTime) {
        *simpleCascadingTime += clock3 - clock2;
    }
    if (naiveTime) {
        *naiveTime += clock5 - clock4;
    }

}

void SubmatrixQueriesTest::benchmarkAllRowQueries(clock_t *naiveTime, clock_t *queryTime, clock_t *cascadingTime, clock_t *simpleCascadingTime)
{
    size_t row = rand() % (_testMatrix->rows());
    
    size_t c1,c2;
    c1 = rand() % (_testMatrix->cols());
    c2 = rand() % (_testMatrix->cols());
    
    Range c = Range(min(c1,c2),max(c1,c2));
    
    benchmarkAllRowQueries(c, row,naiveTime,queryTime,cascadingTime,simpleCascadingTime);
}

void SubmatrixQueriesTest::benchmarkAllColQueries(Range rowRange, size_t col, clock_t *naiveTime, clock_t *queryTime, clock_t *cascadingTime, clock_t *simpleCascadingTime)
{
    clock_t clock1, clock2, clock3, clock4, clock5;
    
    clock1 = clock();
    _queryDS->rowsTree()->cascadingMaxInRange(col, rowRange);
    
    clock2 = clock();
    _queryDS->rowsTree()->simpleCascadingMaxInRange(col, rowRange);
    
    clock3 = clock();
    _queryDS->rowsTree()->maxForColumnInRange(col, rowRange.min, rowRange.max);
    
    clock4 = clock();
    SubmatrixQueriesTest::naiveMaximumInColumn(_testMatrix, rowRange, col);
    
    clock5 = clock();
    
    if (queryTime) {
        *queryTime += clock4 - clock3;
    }
    if (cascadingTime) {
        *cascadingTime += clock2 - clock1;
    }
    if (simpleCascadingTime) {
        *simpleCascadingTime += clock3 - clock2;
    }
    if (naiveTime) {
        *naiveTime += clock5 - clock4;
    }

}


void SubmatrixQueriesTest::benchmarkAllColQueries(clock_t *naiveTime, clock_t *queryTime, clock_t *cascadingTime, clock_t *simpleCascadingTime)
{
    size_t col = rand() % (_testMatrix->cols());
    
    size_t r1,r2;
    r1 = rand() % (_testMatrix->rows());
    r2 = rand() % (_testMatrix->rows());
    
    Range r = Range(min(r1,r2),max(r1,r2));
    
    benchmarkAllColQueries(r, col,naiveTime,queryTime,cascadingTime,simpleCascadingTime);
    
}

bool SubmatrixQueriesTest::multipleColumnQueryTest(size_t n)
{
    bool result = true;
    clock_t naiveTime = 0, queryTime = 0;
    
    for (size_t i = 0; i < n && result; i++) {
        result = result && testColumnQuery(&naiveTime, &queryTime);
    }
#if BENCHMARK
    cout << "Benchmark for " << n << " column queries:" <<endl;
    cout << "Naive algorithm: " << 1000*((double)naiveTime)/CLOCKS_PER_SEC << " ms" << endl;
    cout << "Submatrix queries: " << 1000*((double)queryTime)/CLOCKS_PER_SEC << " ms" << endl;
#endif
    return result;
}

bool SubmatrixQueriesTest::multipleRowQueryTest(size_t n)
{
    bool result = true;
    clock_t naiveTime = 0, queryTime = 0;
    
    for (size_t i = 0; i < n && result; i++) {
        result = result && testRowQuery(&naiveTime, &queryTime);
    }

#if BENCHMARK
    cout << "Benchmark for " << n << " row queries:" <<endl;
    cout << "Naive algorithm: " << 1000*((double)naiveTime)/CLOCKS_PER_SEC << " ms" << endl;
    cout << "Submatrix queries: " << 1000*((double)queryTime)/CLOCKS_PER_SEC << " ms" << endl;
#endif
    return result;
}

bool SubmatrixQueriesTest::multipleRowQueryTestVsCascading(size_t n)
{
    bool result = true;
    clock_t queryTime = 0, cascadingTime = 0, simpleCascadingTime = 0;
    
    for (size_t i = 0; i < n && result; i++) {
        result = result && testCascadingRowQuery(&queryTime, &cascadingTime, &simpleCascadingTime);
    }
    
#if BENCHMARK
    cout << "Benchmark for " << n << " row queries:" <<endl;
    cout << "Submatrix queries: " << 1000*((double)queryTime)/CLOCKS_PER_SEC << " ms" << endl;
    cout << "Simple Cascading queries: " << 1000*((double)simpleCascadingTime)/CLOCKS_PER_SEC << " ms" << endl;
    cout << "Cascading queries: " << 1000*((double)cascadingTime)/CLOCKS_PER_SEC << " ms" << endl;
#endif
    return result;
}


bool SubmatrixQueriesTest::multipleColQueryTestVsCascading(size_t n)
{
    bool result = true;
    clock_t queryTime = 0, cascadingTime = 0, simpleCascadingTime = 0;
    
    for (size_t i = 0; i < n && result; i++) {
        result = result && testCascadingColQuery(&queryTime, &cascadingTime,  &simpleCascadingTime);
    }
    
#if BENCHMARK
    cout << "Cascading Benchmark for " << n << " row queries:" <<endl;
    cout << "Submatrix queries: " << 1000*((double)queryTime)/CLOCKS_PER_SEC << " ms" << endl;
    cout << "Simple Cascading queries: " << 1000*((double)simpleCascadingTime)/CLOCKS_PER_SEC << " ms" << endl;
    cout << "Cascading queries: " << 1000*((double)cascadingTime)/CLOCKS_PER_SEC << " ms" << endl;
#endif
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
    cout << "Cascading Benchmark for " << n << " submatrix queries:" <<endl;
    cout << "Naive algorithm: " << 1000*((double)naiveTime)/CLOCKS_PER_SEC << " ms" << endl;
    cout << "Submatrix queries: " << 1000*((double)queryTime)/CLOCKS_PER_SEC << " ms" << endl;
#endif
    return result;
}

void SubmatrixQueriesTest::multipleBenchmarksRowQueries(size_t n)
{
    clock_t naiveTime = 0, queryTime = 0, cascadingTime = 0, simpleCascadingTime = 0;
    
    multipleBenchmarksRowQueries(n,&naiveTime, &queryTime, &cascadingTime,  &simpleCascadingTime);

    cout << "Benchmark for " << n << " row queries:" <<endl;
    cout << "Naive algorithm: " << 1000*((double)naiveTime)/CLOCKS_PER_SEC << " ms" << endl;
    cout << "Canonical nodes: " << 1000*((double)queryTime)/CLOCKS_PER_SEC << " ms" << endl;
    cout << "Cascading: " << 1000*((double)cascadingTime)/CLOCKS_PER_SEC << " ms" << endl;
    cout << "Simple cascading: " << 1000*((double)simpleCascadingTime)/CLOCKS_PER_SEC << " ms" << endl;
    
}

void SubmatrixQueriesTest::multipleBenchmarksRowQueries(size_t n, clock_t *naiveTime, clock_t *queryTime, clock_t *cascadingTime, clock_t *simpleCascadingTime)
{
    for (size_t i = 0; i < n; i++) {
        benchmarkAllRowQueries(naiveTime, queryTime, cascadingTime,  simpleCascadingTime);
    }    
}

void SubmatrixQueriesTest::multipleBenchmarksColQueries(size_t n)
{
    clock_t naiveTime = 0, queryTime = 0, cascadingTime = 0, simpleCascadingTime = 0;
    
    multipleBenchmarksColQueries(n,&naiveTime, &queryTime, &cascadingTime,  &simpleCascadingTime);

    cout << "Benchmark for " << n << " row queries:" <<endl;
    cout << "Naive algorithm: " << 1000*((double)naiveTime)/CLOCKS_PER_SEC << " ms" << endl;
    cout << "Canonical nodes: " << 1000*((double)queryTime)/CLOCKS_PER_SEC << " ms" << endl;
    cout << "Cascading: " << 1000*((double)cascadingTime)/CLOCKS_PER_SEC << " ms" << endl;
    cout << "Simple cascading: " << 1000*((double)simpleCascadingTime)/CLOCKS_PER_SEC << " ms" << endl;
    
}

void SubmatrixQueriesTest::multipleBenchmarksColQueries(size_t n, clock_t *naiveTime, clock_t *queryTime, clock_t *cascadingTime, clock_t *simpleCascadingTime)
{
    for (size_t i = 0; i < n; i++) {
        benchmarkAllColQueries(naiveTime, queryTime, cascadingTime,  simpleCascadingTime);
    }
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

Matrix<double>* SubmatrixQueriesTest::generateInverseMongeMatrixStrip1(size_t rows, size_t cols)
{
    Matrix<double> *m;
    
    cout << "Initializing matrix " << rows << "x" << cols << " ... ";
    try {
        m = new ComplexMatrix<double>(rows,cols);
        cout << "Done" << endl;
    } catch (std::bad_alloc& ba) {
        cout << "\nbad_alloc caught: " << ba.what() << endl;
        cout << "Try to build an other matrix ... ";
        m = new SimpleMatrix<double>(rows,cols);
        cout << "Done" << endl;
    }
    
    cout << "Fill the Inverse Monge Matrix ... " << endl;
    
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
    
    accumulator = 0;
    for (size_t j = 0; j < cols; j++) {
        accumulator += rand() % colInterval;
        colsAbscissa[cols-1-j] = accumulator;
    }
    
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

Matrix<double>* SubmatrixQueriesTest::generateInverseMongeMatrixStrip2(size_t rows, size_t cols)
{
    Matrix<double> *m;
    
#if BENCHMARK
    cout << "Initializing matrix " << rows << "x" << cols << " ... ";
#endif
    try {
        m = new ComplexMatrix<double>(rows,cols);
#if BENCHMARK
        cout << "Done" << endl;
#endif
    } catch (std::bad_alloc& ba) {
#if BENCHMARK
        cout << "\nbad_alloc caught: " << ba.what() << endl;
        cout << "Try to build an other matrix ... ";
#endif
        m = new SimpleMatrix<double>(rows,cols);
#if BENCHMARK
        cout << "Done" << endl;
#endif
    }
    
#if BENCHMARK
    cout << "Fill the Inverse Monge Matrix ... " << endl;
#endif
    
    int *rowsAbscissa, *colsAbscissa;
    
    rowsAbscissa = new int[rows];
    colsAbscissa = new int[cols];
    
    srand ( time(NULL) );
    
    // we define these values to avoid overflows that will lead to a non inverse Monge matrix
    // our goal is to have values in rowsAbscissa and colsAbscissa that are between -max_abscissa and +max_abscissa 
    int max_abscissa;
//    max_abscissa = (int)sqrtf(INT_MAX/3) - LINE_DISTANCE;
    max_abscissa = 0.5*sqrt(INT_MAX-LINE_DISTANCE);
    
    int rowInterval = max<int>(2*max_abscissa/rows,1);
    int colInterval = max<int>(2*max_abscissa/cols,1);
    
    assert(rowInterval > 0);
    assert(colInterval > 0);
    
    int accumulator = -max_abscissa;
    
    for (size_t i = 0; i < rows; i++) {
        accumulator += rand() % rowInterval +1;
        rowsAbscissa[i] = accumulator;
    }
    
    accumulator = -max_abscissa;
    for (size_t j = 0; j < cols; j++) {
        accumulator += rand() % colInterval +1;
        colsAbscissa[cols-1-j] = accumulator;
    }
    
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            int abscissaDiff = rowsAbscissa[i]-colsAbscissa[j];
            int diffSquare = abscissaDiff*abscissaDiff;
            int distSquare = LINE_DISTANCE + diffSquare;
            double v = sqrt(distSquare);
            (*m)(i,j) = v;
        }
    }
    
    delete [] rowsAbscissa;
    delete [] colsAbscissa;
    
    return m;
}

Matrix<double>* SubmatrixQueriesTest::generateInverseMongeMatrixSlope(size_t rows, size_t cols)
{
    Matrix<double> *m;
    
#if BENCHMARK
    cout << "Initializing matrix " << rows << "x" << cols << " ... ";
#endif
    try {
        m = new ComplexMatrix<double>(rows,cols);
#if BENCHMARK
        cout << "Done" << endl;
#endif
    } catch (std::bad_alloc& ba) {
#if BENCHMARK
        cout << "\nbad_alloc caught: " << ba.what() << endl;
        cout << "Try to build an other matrix ... ";
#endif
        m = new SimpleMatrix<double>(rows,cols);
#if BENCHMARK
        cout << "Done" << endl;
#endif
    }
    
#if BENCHMARK
    cout << "Fill the Inverse Monge Matrix ... " << endl;
#endif
    
    vector<double> slope(rows);
    
    srand ( time(NULL) );

    for (size_t i = 0; i < rows; i++) {
        slope[i] = fRand(-0.5*M_PI, 0.5*M_PI);
    }
    std::sort(slope.begin(), slope.end());
    
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            double v;
            v = tan(slope[i])*j + rows -1 -i;
            
            (*m)(i,j) = v;
        }
    }
    
    slope = vector<double>(cols);
    for (size_t i = 0; i < cols; i++) {
        slope[i] = fRand(-0.5*M_PI, 0.5*M_PI);
    }
    std::sort(slope.begin(), slope.end());

    for (size_t i = 0; i < cols; i++) {
        for (size_t j = 0; j < rows; j++) {
            double v;
            v = tan(slope[i])*j + rows -1 -i;
            
            (*m)(j,i) += v;
        }
    }

    return m;
}

typedef clock_t* clock_ptr;

void SubmatrixQueriesTest::multiBenchmarksPositionQueries(size_t nRows, size_t nCols, size_t nSamples, clock_t *naiveTime, clock_t *queryTime, clock_t *cascadingTime, clock_t *simpleCascadingTime)
{
    cout << "\n";
    for (size_t i = 0; i < nSamples; i++) {
        SubmatrixQueriesTest *t = new SubmatrixQueriesTest(nRows,nCols);
        cout<< "|";
        cout.flush();
        t->multipleBenchmarksColQueries(100, naiveTime, queryTime, cascadingTime, simpleCascadingTime);
        t->multipleBenchmarksRowQueries(100, naiveTime, queryTime, cascadingTime, simpleCascadingTime);
        
        delete t;
    }
}

clock_t** SubmatrixQueriesTest::multiSizeBenchmarksPositionQueries(size_t maxNRows, size_t maxNCols, size_t nSampleSize, size_t nSamplePerSize)
{
    clock_t **results = new clock_ptr [nSampleSize];
    for (size_t i = 0; i < nSampleSize; i++) {
        results[i] = new clock_t [4];
        results[i][0] = 0;
        results[i][1] = 0;
        results[i][2] = 0;
        results[i][3] = 0;
        
        size_t nRows = maxNRows, nCols = maxNCols;
        float fraction = ((float)(i+1))/((float)nSampleSize);
        nRows *= fraction;
        nCols *= fraction;
        
        cout << "Benchmark for size: " << nRows << " x " << nCols << " ... ";

        SubmatrixQueriesTest::multiBenchmarksPositionQueries(nRows,nCols, nSamplePerSize, results[i], results[i]+1, results[i]+2, results[i]+3);
        
        cout << " done\n";
    }
    return results;
}
