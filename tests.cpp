//
//  tests.cpp
//  SubmatrixQueries
//
//  Created by Raphael Bost on 08/02/13.
//  Copyright (c) 2013 Raphael Bost. All rights reserved.
//

#include "tests.h"

#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <climits>
#include <set>
#include <cfloat>
#include <time.h>
#include <fstream>

#include "debug_assert.h"

extern "C"
{
#include <pthread.h>
}

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}


#ifdef __MACH__
#define ZeroTime 0
#else
#define ZeroTime {0}
#endif

timespec diffTS(timespec start, timespec end)
{
	timespec temp;
	if ((end.tv_nsec-start.tv_nsec)<0) {
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	} else {
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp;
}

timespec addTS(timespec t1, timespec t2)
{
	timespec temp;
	temp.tv_sec = t1.tv_sec + t2.tv_sec;
	temp.tv_nsec = t1.tv_nsec + t2.tv_nsec;
	if(temp.tv_nsec >= 1000000000){
		temp.tv_sec++;
		temp.tv_nsec -= 1000000000;
	}
	return temp;
}

bench_time_t diff(bench_time_t t1, bench_time_t t2)
{
#ifdef __MACH__
    return t1 - t2;
#else
    return diffTS(t2,t1);
#endif
}

bench_time_t add(bench_time_t t1, bench_time_t t2)
{
#ifdef __MACH__
    return t1 + t2;
#else
    return addTS(t2,t1);
#endif
}

bench_time_t now()
{
#ifdef __MACH__
    return clock();
#else
    timespec t;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t);
    return t;
#endif
}

double timespecAsSeconds(timespec ts)
{
    return ts.tv_sec + (((double)ts.tv_nsec)/((double)1000000000));
}

double timespecAsMiliSeconds(timespec ts)
{
    return 1000*ts.tv_sec + (((double)ts.tv_nsec)/((double)1000000));
}

double benchTimeAsMiliSeconds(bench_time_t t)
{
#ifdef __MACH__
    return 1000*((double)t)/CLOCKS_PER_SEC;
#else
    return timespecAsMiliSeconds(t);
#endif
}

#define PRINT_TEST_MATRIX false
#define BENCHMARK false

#define MULTITHREAD_GENERATION true
#define GENERATION_THREAD_COUNT 15

SubmatrixQueriesTest::SubmatrixQueriesTest(Matrix<double> *m)
{
    DEBUG_ASSERT(m->isInverseMonge());
    // copy the matrix
    _testMatrix = new ComplexMatrix<double>(m);
    
#if BENCHMARK
    bench_time_t time = now();
#endif
    _queryDS = new SubmatrixQueriesDataStructure<double>(*_testMatrix);
#if BENCHMARK
    time = diff(now(),time);
    cout << "Building Data Structure: " << benchTimeAsMiliSeconds(time) << " ms" << endl;
#endif
}

SubmatrixQueriesTest::SubmatrixQueriesTest(size_t rows, size_t cols)
{
#if BENCHMARK
    bench_time_t time = now();
#endif

#if MULTITHREAD_GENERATION
    _testMatrix = generateInverseMongeMatrixSlopeMultithread(rows, cols,GENERATION_THREAD_COUNT);
#else
    _testMatrix = generateInverseMongeMatrixSlope(rows, cols);
#endif

#if BENCHMARK
    time = diff(now(),time);
    cout << "Building Matrix: " << benchTimeAsMiliSeconds(time) << " ms" << endl;
    time = now();
#endif
//    DEBUG_ASSERT(_testMatrix->isInverseMonge());
    
    _queryDS = new SubmatrixQueriesDataStructure<double>(*_testMatrix);
#if BENCHMARK
    time = diff(now(),time);
    cout << "Building Data Structure: " << benchTimeAsMiliSeconds(time) << " ms" << endl;
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

bool SubmatrixQueriesTest::testColumnQuery(Range rowRange, size_t col, bench_time_t *naiveTime, bench_time_t *queryTime)
{
    double naiveMax, queryMax;
    bench_time_t clock1, clock2, clock3;

    clock1 = now();
    queryMax = _queryDS->rowsTree()->maxForColumnInRange(col, rowRange.min, rowRange.max);
    clock2 = now();
    naiveMax = SubmatrixQueriesTest::naiveMaximumInColumn(_testMatrix, rowRange, col);
    
    clock3 = now();
    
    if (naiveTime) {
        *naiveTime = add(*naiveTime,diff(clock3,clock2));
    }
    if (queryTime) {
        *queryTime = add(*queryTime,diff(clock2,clock1));
    }
    
    return queryMax == naiveMax;
}

bool SubmatrixQueriesTest::testColumnQuery(bench_time_t *naiveTime, bench_time_t *queryTime)
{
    size_t col = rand() % (_testMatrix->cols());
    
    size_t r1, r2;
    r1 = rand() % (_testMatrix->rows());
    r2 = rand() % (_testMatrix->rows());
    
    Range r = Range(min(r1,r2),max(r1,r2));
    
    return testColumnQuery(r, col, naiveTime, queryTime);
}

bool SubmatrixQueriesTest::testCascadingColQuery(Range rowRange, size_t col, bench_time_t *queryTime, bench_time_t *cascadingTime, bench_time_t *simpleCascadingTime)
{
    double cascadingMax, simpleCascadingMax, queryMax;
    bench_time_t clock1, clock2, clock3, clock4;
    
    clock1 = now();
    cascadingMax = _queryDS->rowsTree()->cascadingMaxInRange(col, rowRange);
    
    clock2 = now();
    simpleCascadingMax = _queryDS->rowsTree()->simpleCascadingMaxInRange(col, rowRange);
    
    clock3 = now();
    queryMax = _queryDS->rowsTree()->maxForColumnInRange(col, rowRange.min, rowRange.max);
    
    clock4 = now();
    
    if (queryTime) {
        *queryTime = add(*queryTime,diff(clock4, clock3));
    }
    if (cascadingTime) {
        *cascadingTime = add(*cascadingTime,diff(clock2, clock1));
    }
    if (simpleCascadingTime) {
        *simpleCascadingTime = add(*simpleCascadingTime,diff(clock3, clock2));
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

bool SubmatrixQueriesTest::testCascadingColQuery(bench_time_t *queryTime, bench_time_t *cascadingTime, bench_time_t *simpleCascadingTime)
{
    size_t col = rand() % (_testMatrix->cols());
    
    size_t r1,r2;
    r1 = rand() % (_testMatrix->rows());
    r2 = rand() % (_testMatrix->rows());
    
    Range r = Range(min(r1,r2),max(r1,r2));
    
    return testCascadingColQuery(r, col,queryTime,cascadingTime,simpleCascadingTime);
}


bool SubmatrixQueriesTest::testRowQuery(Range colRange, size_t row, bench_time_t *naiveTime, bench_time_t *queryTime)
{
    double naiveMax, queryMax;
    bench_time_t clock1, clock2, clock3;
    
    clock1 = now();
    queryMax = _queryDS->columnTree()->maxForRowInRange(row, colRange.min, colRange.max);
    
    clock2 = now();
    naiveMax = SubmatrixQueriesTest::naiveMaximumInRow(_testMatrix, colRange, row);
    
    clock3 = now();

    if (naiveTime) {
        *naiveTime = add(*naiveTime,diff(clock3,clock2));
    }
    if (queryTime) {
        *queryTime = add(*queryTime,diff(clock2,clock1));
    }

    return queryMax == naiveMax;
}

bool SubmatrixQueriesTest::testRowQuery(bench_time_t *naiveTime, bench_time_t *queryTime)
{
    size_t row = rand() % (_testMatrix->rows());
    
    size_t c1,c2;
    c1 = rand() % (_testMatrix->cols());
    c2 = rand() % (_testMatrix->cols());
    
    Range c = Range(min(c1,c2),max(c1,c2));
    
    return testRowQuery(c, row,naiveTime,queryTime);
}

bool SubmatrixQueriesTest::testCascadingRowQuery(Range colRange, size_t row, bench_time_t *queryTime, bench_time_t *cascadingTime, bench_time_t *simpleCascadingTime)
{
    double cascadingMax, queryMax, simpleCascadingMax;
    bench_time_t clock1, clock2, clock3, clock4;
    
    clock1 = now();
    cascadingMax = _queryDS->columnTree()->cascadingMaxInRange(row, colRange);
    
    clock2 = now();
    simpleCascadingMax = _queryDS->columnTree()->simpleCascadingMaxInRange(row, colRange);

    clock3 = now();
    queryMax = _queryDS->columnTree()->maxForRowInRange(row, colRange.min, colRange.max);

    clock4 = now();

    if (queryTime) {
        *queryTime = add(*queryTime,diff(clock4, clock3));
    }
    if (cascadingTime) {
        *cascadingTime = add(*cascadingTime,diff(clock2, clock1));
    }
    if (simpleCascadingTime) {
        *simpleCascadingTime = add(*simpleCascadingTime,diff(clock3, clock2));
    }

    return queryMax == cascadingMax;
}

bool SubmatrixQueriesTest::testCascadingRowQuery(bench_time_t *queryTime, bench_time_t *cascadingTime, bench_time_t *simpleCascadingTime)
{
    size_t row = rand() % (_testMatrix->rows());
    
    size_t c1,c2;
    c1 = rand() % (_testMatrix->cols());
    c2 = rand() % (_testMatrix->cols());
    
    Range c = Range(min(c1,c2),max(c1,c2));
    
    return testCascadingRowQuery(c, row,queryTime,cascadingTime,simpleCascadingTime);
}

bool SubmatrixQueriesTest::testSubmatrixQuery(Range rowRange, Range colRange, bench_time_t *naiveTime, bench_time_t *queryTime)
{
    double naiveMax, queryMax;
    
    bench_time_t clock1, clock2, clock3;
    
    clock1 = now();
    queryMax = _queryDS->maxInRange(rowRange,colRange);
    clock2 = now();
    naiveMax = SubmatrixQueriesTest::naiveMaximumInSubmatrix(_testMatrix, rowRange, colRange);
    clock3 = now();
    
    if (naiveTime) {
        *naiveTime = add(*naiveTime,diff(clock3, clock2));
    }
    if (queryTime) {
        *queryTime = add(*queryTime,diff(clock2, clock1));
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

bool SubmatrixQueriesTest::testSubmatrixQuery(bench_time_t *naiveTime, bench_time_t *queryTime)
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

void SubmatrixQueriesTest::benchmarkAllRowQueries(Range colRange, size_t row, bench_time_t *naiveTime, bench_time_t *queryTime, bench_time_t *cascadingTime, bench_time_t *simpleCascadingTime)
{
    bench_time_t clock1, clock2, clock3, clock4, clock5;
    
    clock1 = now();
    _queryDS->columnTree()->cascadingMaxInRange(row, colRange);
    
    clock2 = now();
    _queryDS->columnTree()->simpleCascadingMaxInRange(row, colRange);

    clock3 = now();
     _queryDS->columnTree()->maxForRowInRange(row, colRange.min, colRange.max);
    
    clock4 = now();
    SubmatrixQueriesTest::naiveMaximumInRow(_testMatrix, colRange, row);

    clock5 = now();

    if (queryTime) {
        *queryTime = add(*queryTime,diff(clock4, clock3));
    }
    if (cascadingTime) {
        *cascadingTime = add(*cascadingTime,diff(clock2, clock1));
    }
    if (simpleCascadingTime) {
        *simpleCascadingTime = add(*simpleCascadingTime,diff(clock3, clock2));
    }
    if (naiveTime) {
        *naiveTime = add(*naiveTime,diff(clock5, clock4));
    }
}

void SubmatrixQueriesTest::benchmarkAllRowQueries(bench_time_t *naiveTime, bench_time_t *queryTime, bench_time_t *cascadingTime, bench_time_t *simpleCascadingTime)
{
    size_t row = rand() % (_testMatrix->rows());
    
    size_t c1,c2;
    c1 = rand() % (_testMatrix->cols());
    c2 = rand() % (_testMatrix->cols());
    
    Range c = Range(min(c1,c2),max(c1,c2));
    
    benchmarkAllRowQueries(c, row,naiveTime,queryTime,cascadingTime,simpleCascadingTime);
}

void SubmatrixQueriesTest::benchmarkAllColQueries(Range rowRange, size_t col, bench_time_t *naiveTime, bench_time_t *queryTime, bench_time_t *cascadingTime, bench_time_t *simpleCascadingTime)
{
    bench_time_t clock1, clock2, clock3, clock4, clock5;
    
    clock1 = now();
    _queryDS->rowsTree()->cascadingMaxInRange(col, rowRange);
    
    clock2 = now();
    _queryDS->rowsTree()->simpleCascadingMaxInRange(col, rowRange);
    
    clock3 = now();
    _queryDS->rowsTree()->maxForColumnInRange(col, rowRange.min, rowRange.max);
    
    clock4 = now();
    SubmatrixQueriesTest::naiveMaximumInColumn(_testMatrix, rowRange, col);
    
    clock5 = now();
    
    if (queryTime) {
        *queryTime = add(*queryTime,diff(clock4, clock3));
    }
    if (cascadingTime) {
        *cascadingTime = add(*cascadingTime,diff(clock2, clock1));
    }
    if (simpleCascadingTime) {
        *simpleCascadingTime = add(*simpleCascadingTime,diff(clock3, clock2));
    }
    if (naiveTime) {
        *naiveTime = add(*naiveTime,diff(clock5, clock4));
    }

}


void SubmatrixQueriesTest::benchmarkAllColQueries(bench_time_t *naiveTime, bench_time_t *queryTime, bench_time_t *cascadingTime, bench_time_t *simpleCascadingTime)
{
    size_t col = rand() % (_testMatrix->cols());
    
    size_t r1,r2;
    r1 = rand() % (_testMatrix->rows());
    r2 = rand() % (_testMatrix->rows());
    
    Range r = Range(min(r1,r2),max(r1,r2));
    
    benchmarkAllColQueries(r, col,naiveTime,queryTime,cascadingTime,simpleCascadingTime);
    
}

void SubmatrixQueriesTest::benchmarkAllSubmatrixQueries(Range rowRange, Range colRange, bench_time_t *naiveTime, bench_time_t *explicitNodesTime, bench_time_t *implicitNodesTime)
{
    bench_time_t clock1, clock2, clock3, clock4;
    
    
    clock1 = now();
    _queryDS->maxInRange(rowRange,colRange);
    
    clock2 = now();
    _queryDS->maxInSubmatrix(rowRange,colRange);
    
    clock3 = now();
    SubmatrixQueriesTest::naiveMaximumInSubmatrix(_testMatrix, rowRange, colRange);
    
    clock4 = now();
    
    if (naiveTime) {
        *naiveTime = add(*naiveTime,diff(clock4, clock3));
    }
    if (explicitNodesTime) {
        *explicitNodesTime = add(*explicitNodesTime,diff(clock2, clock1));
    }
    if (implicitNodesTime) {
        *implicitNodesTime = add(*implicitNodesTime,diff(clock3, clock2));
    }
    
}


void SubmatrixQueriesTest::benchmarkAllSubmatrixQueries(bench_time_t *naiveTime, bench_time_t *explicitNodesTime, bench_time_t *implicitNodesTime)
{
    size_t r1, r2;
    r1 = rand() % (_testMatrix->rows());
    r2 = rand() % (_testMatrix->rows());
    
    size_t c1,c2;
    c1 = rand() % (_testMatrix->cols());
    c2 = rand() % (_testMatrix->cols());
    
    
    Range r = Range(min(r1,r2),max(r1,r2));
    Range c = Range(min(c1,c2),max(c1,c2));
    
    benchmarkAllSubmatrixQueries(r, c,naiveTime,explicitNodesTime,implicitNodesTime);
    
}

bool SubmatrixQueriesTest::multipleColumnQueryTest(size_t n)
{
    bool result = true;
    bench_time_t naiveTime = ZeroTime, queryTime = ZeroTime;
    
    for (size_t i = 0; i < n && result; i++) {
        result = result && testColumnQuery(&naiveTime, &queryTime);
    }
#if BENCHMARK
    cout << "Benchmark for " << n << " column queries:" <<endl;
    cout << "Naive algorithm: " << benchTimeAsMiliSeconds(naiveTime) << " ms" << endl;
    cout << "Submatrix queries: " << benchTimeAsMiliSeconds(queryTime) << " ms" << endl;
#endif
    return result;
}

bool SubmatrixQueriesTest::multipleRowQueryTest(size_t n)
{
    bool result = true;
    bench_time_t naiveTime = ZeroTime, queryTime = ZeroTime;
    
    for (size_t i = 0; i < n && result; i++) {
        result = result && testRowQuery(&naiveTime, &queryTime);
    }

#if BENCHMARK
    cout << "Benchmark for " << n << " row queries:" <<endl;
    cout << "Naive algorithm: " << benchTimeAsMiliSeconds(naiveTime) << " ms" << endl;
    cout << "Submatrix queries: " << benchTimeAsMiliSeconds(queryTime) << " ms" << endl;
#endif
    return result;
}

bool SubmatrixQueriesTest::multipleRowQueryTestVsCascading(size_t n)
{
    bool result = true;
    bench_time_t queryTime = ZeroTime, cascadingTime = ZeroTime, simpleCascadingTime = ZeroTime;
    
    for (size_t i = 0; i < n && result; i++) {
        result = result && testCascadingRowQuery(&queryTime, &cascadingTime, &simpleCascadingTime);
    }
    
#if BENCHMARK
    cout << "Benchmark for " << n << " row queries:" <<endl;
    cout << "Submatrix queries: " << benchTimeAsMiliSeconds(queryTime) << " ms" << endl;
    cout << "Simple Cascading queries: " << benchTimeAsMiliSeconds(simpleCascadingTime) << " ms" << endl;
    cout << "Cascading queries: " << benchTimeAsMiliSeconds(cascadingTime) << " ms" << endl;
#endif
    return result;
}


bool SubmatrixQueriesTest::multipleColQueryTestVsCascading(size_t n)
{
    bool result = true;
    bench_time_t queryTime = ZeroTime, cascadingTime = ZeroTime, simpleCascadingTime = ZeroTime;
    
   
    for (size_t i = 0; i < n && result; i++) {
        result = result && testCascadingColQuery(&queryTime, &cascadingTime,  &simpleCascadingTime);
    }
    
#if BENCHMARK
    cout << "Cascading Benchmark for " << n << " col queries:" <<endl;
    cout << "Submatrix queries: " << benchTimeAsMiliSeconds(queryTime) << " ms" << endl;
    cout << "Simple Cascading queries: " << benchTimeAsMiliSeconds(simpleCascadingTime) << " ms" << endl;
    cout << "Cascading queries: " << benchTimeAsMiliSeconds(cascadingTime) << " ms" << endl;
#endif
    return result;
}

bool SubmatrixQueriesTest::multipleSubmatrixQueryTest(size_t n)
{
    bool result = true;
    bench_time_t naiveTime = ZeroTime, queryTime = ZeroTime;
    
    for (size_t i = 0; i < n && result; i++) {
        result = result && testSubmatrixQuery(&naiveTime, &queryTime);
    }
    
#if BENCHMARK
    cout << "Cascading Benchmark for " << n << " submatrix queries:" <<endl;
    cout << "Naive algorithm: " << benchTimeAsMiliSeconds(naiveTime) << " ms" << endl;
    cout << "Submatrix queries: " << benchTimeAsMiliSeconds(queryTime) << " ms" << endl;
#endif
    return result;
}

void SubmatrixQueriesTest::multipleBenchmarksRowQueries(size_t n)
{
    bench_time_t naiveTime = ZeroTime, queryTime = ZeroTime, cascadingTime = ZeroTime, simpleCascadingTime = ZeroTime;
    
    multipleBenchmarksRowQueries(n,&naiveTime, &queryTime, &cascadingTime,  &simpleCascadingTime);

    cout << "Benchmark for " << n << " row queries:" <<endl;
    cout << "Naive algorithm: " << benchTimeAsMiliSeconds(naiveTime) << " ms" << endl;
    cout << "Canonical nodes: " << benchTimeAsMiliSeconds(queryTime) << " ms" << endl;
    cout << "Cascading: " << benchTimeAsMiliSeconds(cascadingTime) << " ms" << endl;
    cout << "Simple Cascading queries: " << benchTimeAsMiliSeconds(simpleCascadingTime) << " ms" << endl;
}

void SubmatrixQueriesTest::multipleBenchmarksRowQueries(size_t n, bench_time_t *naiveTime, bench_time_t *queryTime, bench_time_t *cascadingTime, bench_time_t *simpleCascadingTime)
{
    for (size_t i = 0; i < n; i++) {
        benchmarkAllRowQueries(naiveTime, queryTime, cascadingTime,  simpleCascadingTime);
    }    
}

void SubmatrixQueriesTest::multipleBenchmarksColQueries(size_t n)
{
    bench_time_t naiveTime = ZeroTime, queryTime = ZeroTime, cascadingTime = ZeroTime, simpleCascadingTime = ZeroTime;
    
    multipleBenchmarksColQueries(n,&naiveTime, &queryTime, &cascadingTime,  &simpleCascadingTime);

    cout << "Benchmark for " << n << " row queries:" <<endl;
    cout << "Naive algorithm: " << benchTimeAsMiliSeconds(naiveTime) << " ms" << endl;
    cout << "Canonical nodes: " << benchTimeAsMiliSeconds(queryTime) << " ms" << endl;
    cout << "Cascading: " << benchTimeAsMiliSeconds(cascadingTime) << " ms" << endl;
    cout << "Simple Cascading queries: " << benchTimeAsMiliSeconds(simpleCascadingTime) << " ms" << endl;
    
}

void SubmatrixQueriesTest::multipleBenchmarksColQueries(size_t n, bench_time_t *naiveTime, bench_time_t *queryTime, bench_time_t *cascadingTime, bench_time_t *simpleCascadingTime)
{
    for (size_t i = 0; i < n; i++) {
        benchmarkAllColQueries(naiveTime, queryTime, cascadingTime,  simpleCascadingTime);
    }
}


void SubmatrixQueriesTest::multipleBenchmarksSubmatrixQueries(size_t n)
{
    bench_time_t naiveTime = ZeroTime, explicitNodesTime = ZeroTime, implicitNodesTime = ZeroTime;
    
    multipleBenchmarksSubmatrixQueries(n,&naiveTime, &explicitNodesTime, &implicitNodesTime);
    
    cout << "Benchmark for " << n << " submatrix queries:" <<endl;
    cout << "Naive algorithm: " << benchTimeAsMiliSeconds(naiveTime) << " ms" << endl;
    cout << "Explicit canonical nodes: " << benchTimeAsMiliSeconds(explicitNodesTime) << " ms" << endl;
    cout << "Implicit  canonical nodes: " << benchTimeAsMiliSeconds(implicitNodesTime) << " ms" << endl;
    
}

void SubmatrixQueriesTest::multipleBenchmarksSubmatrixQueries(size_t n,bench_time_t *naiveTime, bench_time_t *explicitNodesTime, bench_time_t *implicitNodesTime)
{
    for (size_t i = 0; i < n; i++) {
        benchmarkAllSubmatrixQueries(naiveTime, explicitNodesTime, implicitNodesTime);
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
    
    DEBUG_ASSERT(rowInterval > 0);
    DEBUG_ASSERT(colInterval > 0);
    
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
    
    DEBUG_ASSERT(rowInterval > 0);
    DEBUG_ASSERT(colInterval > 0);
    
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
    cout << "Fill the Inverse Monge Matrix ... ";
    bench_time_t t = now();
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

#if BENCHMARK
    t = diff(now(),t);
    cout << " done in " << benchTimeAsMiliSeconds(t) << " ms" << endl;
#endif
    
    return m;
}

struct _thread_arg_t
{
    Matrix<double> *matrix;
    vector<double> *slopes;
    size_t min;
    size_t max;
    
    _thread_arg_t(Matrix<double> *m, vector<double> *s, size_t min, size_t max) : matrix(m), slopes(s), min(min), max(max){}
};

void* _fillRows(void *args)
{
    _thread_arg_t input_args = *((_thread_arg_t *)args);
    
    size_t cols = input_args.matrix->cols();
    size_t rows = input_args.matrix->rows();

    for (size_t i = input_args.min; i < input_args.max; i++) {
        for (size_t j = 0; j < cols; j++) {
            double v;
            v = tan((*input_args.slopes)[i])*j + rows -1 -i;
            
            (*input_args.matrix)(i,j) = v;
        }
    }

    
    return NULL;
}

void* _fillCols(void *args)
{
    _thread_arg_t input_args = *((_thread_arg_t *)args);
    
    size_t cols = input_args.matrix->cols();
    size_t rows = input_args.matrix->rows();
    
    for (size_t i = input_args.min; i < input_args.max; i++) {
        for (size_t j = 0; j < rows; j++) {
            double v;
            v = tan((*input_args.slopes)[i])*j + rows -1 -i;
            
            (*input_args.matrix)(j,i) += v;
        }
    }

    return NULL;
}

Matrix<double>* SubmatrixQueriesTest::generateInverseMongeMatrixSlopeMultithread(size_t rows, size_t cols,size_t threadCount)
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
    cout << "Fill the Inverse Monge Matrix ... ";
    bench_time_t t = now();
#endif
    
    vector<double> slope(rows);
    
    srand ( time(NULL) );
    
    for (size_t i = 0; i < rows; i++) {
        slope[i] = fRand(-0.5*M_PI, 0.5*M_PI);
    }
    std::sort(slope.begin(), slope.end());
    
    pthread_t *threads = (pthread_t*)malloc(threadCount*sizeof(pthread_t));
    _thread_arg_t *args = (_thread_arg_t*)malloc(threadCount*sizeof(_thread_arg_t));
    size_t rowRange = rows/threadCount;
    int rc;
    
    for (size_t k = 0; k  < threadCount; k++) {
        args[k] = _thread_arg_t(m,&slope,k*rowRange,min((k+1)*rowRange,rows));
        rc = pthread_create(&threads[k], NULL, _fillRows, &args[k]);
        assert(0 == rc);
    }

    for (size_t k=0; k < threadCount; k++) {
        rc = pthread_join(threads[k], NULL);
        assert(0 == rc);
    }

    free(threads); free(args);

    slope = vector<double>(cols);
    for (size_t i = 0; i < cols; i++) {
        slope[i] = fRand(-0.5*M_PI, 0.5*M_PI);
    }
    std::sort(slope.begin(), slope.end());
    
    threads = (pthread_t*)malloc(threadCount*sizeof(pthread_t));
    args = (_thread_arg_t*)malloc(threadCount*sizeof(_thread_arg_t));
    size_t colRange = cols/threadCount;

    for (size_t k = 0; k  < threadCount; k++) {
        args[k] = _thread_arg_t(m,&slope,k*colRange,min((k+1)*colRange,cols));
        rc = pthread_create(&threads[k], NULL, _fillCols, &args[k]);
        assert(0 == rc);
    }
    
    for (size_t k=0; k < threadCount; k++) {
        rc = pthread_join(threads[k], NULL);
        assert(0 == rc);
    }

    free(threads); free(args);

#if BENCHMARK
    t = diff(now(),t);
    cout << " done in " << benchTimeAsMiliSeconds(t) << " ms" << endl;
#endif
    
    return m;
}

typedef bench_time_t* clock_ptr;

void SubmatrixQueriesTest::multiBenchmarksPositionQueries(size_t nRows, size_t nCols, size_t nSamples, bench_time_t *naiveTime, bench_time_t *queryTime, bench_time_t *cascadingTime, bench_time_t *simpleCascadingTime)
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

bench_time_t** SubmatrixQueriesTest::multiSizeBenchmarksPositionQueries(size_t maxNRows, size_t maxNCols, size_t nSampleSize, size_t nSamplePerSize)
{
    bench_time_t **results = new clock_ptr [nSampleSize];
    for (size_t i = 0; i < nSampleSize; i++) {
        results[i] = new bench_time_t [4];
        
#ifdef __MACH__
	results[i][0] = 0;
        results[i][1] = 0;
        results[i][2] = 0;
        results[i][3] = 0;
#else
        results[i][0].tv_sec = 0; results[i][0].tv_nsec = 0;
        results[i][1].tv_sec = 0; results[i][1].tv_nsec = 0;
        results[i][2].tv_sec = 0; results[i][2].tv_nsec = 0;
        results[i][3].tv_sec = 0; results[i][3].tv_nsec = 0;
#endif
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

void SubmatrixQueriesTest::multiBenchmarksSubmatrixQueries(size_t nRows, size_t nCols, size_t nSamples, bench_time_t *naiveTime, bench_time_t *explicitNodesTime, bench_time_t *implicitNodesTime)
{
    cout << "\n";
    for (size_t i = 0; i < nSamples; i++) {
        SubmatrixQueriesTest *t = new SubmatrixQueriesTest(nRows,nCols);
        cout<< "|";
        cout.flush();
        t->multipleBenchmarksSubmatrixQueries(100, naiveTime, explicitNodesTime, implicitNodesTime);
        t->multipleBenchmarksSubmatrixQueries(100, naiveTime, explicitNodesTime, implicitNodesTime);
        
        delete t;
    }
}

bench_time_t** SubmatrixQueriesTest::multiSizeBenchmarksSubmatrixQueries(size_t maxNRows, size_t maxNCols, size_t nSampleSize, size_t nSamplePerSize)
{
    bench_time_t **results = new clock_ptr [nSampleSize];
    for (size_t i = 0; i < nSampleSize; i++) {
        results[i] = new bench_time_t [3];
        
#ifdef __MACH__
        results[i][0] = 0;
        results[i][1] = 0;
        results[i][2] = 0;
#else
        results[i][0].tv_sec = 0; results[i][0].tv_nsec = 0;
        results[i][1].tv_sec = 0; results[i][1].tv_nsec = 0;
        results[i][2].tv_sec = 0; results[i][2].tv_nsec = 0;
#endif
        size_t nRows = maxNRows, nCols = maxNCols;
        float fraction = ((float)(i+1))/((float)nSampleSize);
        nRows *= fraction;
        nCols *= fraction;
        
        cout << "Benchmark for size: " << nRows << " x " << nCols << " ... ";
        
        SubmatrixQueriesTest::multiBenchmarksSubmatrixQueries(nRows,nCols, nSamplePerSize, results[i], results[i]+1, results[i]+2);
        
        cout << " done\n";
    }
    return results;
}

void SubmatrixQueriesTest::multiSizeBenchmarksSubmatrixQueries(size_t maxNRows, size_t maxNCols, size_t nSampleSize, size_t nSamplePerSize, ofstream &outputStream)
{
    for (size_t i = 0; i < nSampleSize; i++) {
        bench_time_t *benchmarks = new bench_time_t [3];
        
#ifdef __MACH__
        benchmarks[0] = 0;
        benchmarks[1] = 0;
        benchmarks[2] = 0;
#else
        benchmarks[0].tv_sec = 0; benchmarks[0].tv_nsec = 0;
        benchmarks[1].tv_sec = 0; benchmarks[1].tv_nsec = 0;
        benchmarks[2].tv_sec = 0; benchmarks[2].tv_nsec = 0;
#endif
        size_t nRows = maxNRows, nCols = maxNCols;
        float fraction = ((float)(i+1))/((float)nSampleSize);
        nRows *= fraction;
        nCols *= fraction;
        
        cout << "Benchmark for size: " << nRows << " x " << nCols << " ... ";
        
        SubmatrixQueriesTest::multiBenchmarksSubmatrixQueries(nRows,nCols, nSamplePerSize, benchmarks, benchmarks+1, benchmarks+2);
        
        cout << " done\n";
        
        outputStream << ((int)nRows) << " ; " << benchTimeAsMiliSeconds(benchmarks[0]) << " ; " << benchTimeAsMiliSeconds(benchmarks[1]) << " ; " << benchTimeAsMiliSeconds(benchmarks[2]) << "\n";
        outputStream.flush();
        
        delete [] benchmarks;
    }    
}
