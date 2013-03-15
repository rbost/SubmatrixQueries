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
#include <random>

#define PRINT_TEST_MATRIX false
#define BENCHMARK true

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
    cout << "\n";  
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

bool SubmatrixQueriesTest::testCascadingColQuery(Range rowRange, size_t col, clock_t *queryTime, clock_t *cascadingTime)
{
    double cascadingMax, queryMax;
    clock_t clock1, clock2, clock3;
    
    clock1 = clock();
    cascadingMax = _queryDS->rowsTree()->cascadingMaxInRange(col, rowRange);
    
    clock2 = clock();
    queryMax = _queryDS->rowsTree()->maxForColumnInRange(col, rowRange.min, rowRange.max);
    
    clock3 = clock();
    
    if (queryTime) {
        *queryTime += clock3 - clock2;
    }
    if (cascadingTime) {
        *cascadingTime += clock2 - clock1;
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

bool SubmatrixQueriesTest::testCascadingColQuery(clock_t *queryTime, clock_t *cascadingTime)
{
    size_t col = rand() % (_testMatrix->cols());
    
    size_t r1,r2;
    r1 = rand() % (_testMatrix->rows());
    r2 = rand() % (_testMatrix->rows());
    
    Range r = Range(min(r1,r2),max(r1,r2));
    
    return testCascadingColQuery(r, col,queryTime,cascadingTime);
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

bool SubmatrixQueriesTest::testCascadingRowQuery(Range colRange, size_t row, clock_t *queryTime, clock_t *cascadingTime)
{
    double cascadingMax, queryMax;
    clock_t clock1, clock2, clock3;
    
    clock1 = clock();
    cascadingMax = _queryDS->columnTree()->cascadingMaxInRange(row, colRange);
    
    clock2 = clock();
    queryMax = _queryDS->columnTree()->maxForRowInRange(row, colRange.min, colRange.max);

    clock3 = clock();

    if (queryTime) {
        *queryTime += clock3 - clock2;
    }
    if (cascadingTime) {
        *cascadingTime += clock2 - clock1;
    }

    return queryMax == cascadingMax;
}

bool SubmatrixQueriesTest::testCascadingRowQuery(clock_t *queryTime, clock_t *cascadingTime)
{
    size_t row = rand() % (_testMatrix->rows());
    
    size_t c1,c2;
    c1 = rand() % (_testMatrix->cols());
    c2 = rand() % (_testMatrix->cols());
    
    Range c = Range(min(c1,c2),max(c1,c2));
    
    return testCascadingRowQuery(c, row,queryTime,cascadingTime);    
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
    clock_t queryTime = 0, cascadingTime = 0;
    
    for (size_t i = 0; i < n && result; i++) {
        result = result && testCascadingRowQuery(&queryTime, &cascadingTime);
    }
    
#if BENCHMARK
    cout << "Benchmark for " << n << " row queries:" <<endl;
    cout << "Submatrix queries: " << 1000*((double)queryTime)/CLOCKS_PER_SEC << " ms" << endl;
    cout << "Cascading queries: " << 1000*((double)cascadingTime)/CLOCKS_PER_SEC << " ms" << endl;
#endif
    return result;
}


bool SubmatrixQueriesTest::multipleColQueryTestVsCascading(size_t n)
{
    bool result = true;
    clock_t queryTime = 0, cascadingTime = 0;
    
    for (size_t i = 0; i < n && result; i++) {
        result = result && testCascadingColQuery(&queryTime, &cascadingTime);
    }
    
#if BENCHMARK
    cout << "Cascading Benchmark for " << n << " row queries:" <<endl;
    cout << "Submatrix queries: " << 1000*((double)queryTime)/CLOCKS_PER_SEC << " ms" << endl;
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

#define LINE_DISTANCE2 0.25*max_abscissa*max_abscissa

Matrix<double>* SubmatrixQueriesTest::generateInverseMongeMatrixStripPerturbated(size_t rows, size_t cols, double pointsRange, double stripDistance)
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
    
    vector<double> rowsAbscissa, colsAbscissa;
    
    rowsAbscissa = vector<double>(rows);
    colsAbscissa = vector<double>(cols);
    std::pair<std::set<double>::iterator,bool> ret;
    
    if (stripDistance <= 0 ) {
        stripDistance = 2*pointsRange;
    }
    
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    normal_distribution<double> distribution(0, pointsRange);
    
    for (size_t i = 0; i < rows; i++) {
        rowsAbscissa[i] = distribution(generator);
    }

    for (size_t i = 0; i < rows; i++) {
        colsAbscissa[i] = distribution(generator);
    }
    
    std::set<double>::iterator rowIt, colIt;
    
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            double abscissaDiff = rowsAbscissa[i]-colsAbscissa[cols-1-j];
            double diffSquare = abscissaDiff*abscissaDiff;
            double distSquare = stripDistance*stripDistance + diffSquare;
            double v = sqrt(distSquare);

            (*m)(i,j) = v;
        }
    }
    return m;
}
#define LINE_DISTANCE3 300
Matrix<double>* SubmatrixQueriesTest::generateInverseMongeMatrixStripPerturbated2(size_t rows, size_t cols)
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
    
    double range = 50;
    
    vector<double> rowsAbscissa, colsAbscissa;
    
    rowsAbscissa = vector<double>(rows);
    colsAbscissa = vector<double>(cols);
    std::pair<std::set<double>::iterator,bool> ret;
    
     
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    normal_distribution<double> distribution(0, 3);
    
    for (size_t i = 0; i < rows; i++) {
        rowsAbscissa[i] = distribution(generator);
    }
    
    for (size_t i = 0; i < rows; i++) {
        colsAbscissa[i] = distribution(generator);
    }
    
    std::set<double>::iterator rowIt, colIt;
    
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            double rx,ry,cx,cy;
            rx = rowsAbscissa[i];
            cx = rowsAbscissa[cols-1-j] + LINE_DISTANCE3;
            ry = ((((double)i)/(double)rows)-0.5)*range;
            cy = (-(((double)j)/(double)cols)+0.5)*range;
            
            double xDiff = rx - cx;
            double yDiff = ry - cy;
            double v = sqrt(xDiff*xDiff + yDiff*yDiff);
            
            (*m)(i,j) = v;
        }
    }
    return m;
}

Matrix<double>* SubmatrixQueriesTest::generateInverseMongeMatrixSlope(size_t rows, size_t cols)
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
    
    vector<double> slope(rows);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    uniform_real_distribution<double> distribution(-0.5*M_PI, 0.5*M_PI);

    for (size_t i = 0; i < rows; i++) {
        slope[i] = distribution(generator);
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
        slope[i] = distribution(generator);
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