//
//  main.cpp
//  SubmatrixQueries
//
//  Created by Raphael Bost on 02/01/13.
//  Copyright (c) 2013 Raphael Bost. All rights reserved.
//

#include <iostream>
#include <cstdio>

#include "matrix.h"
#include "range.h"
#include "envelope.h"
#include "envelope_tree.h"
#include "tests.h"

using namespace std;
using namespace matrix;
using namespace envelope;

void printBreakpointList(const vector< Breakpoint > *bp)
{
    for (size_t i = 0; i < bp->size(); i++) {
        cout << '(' << (*bp)[i].row << ',' << (*bp)[i].col << ')';
    }
}

void testEnvelope()
{
    valarray<int> foo (35);
    foo[0] =  23;    foo[1] = 28;    foo[2] = 13;    foo[3] = 17;    foo[4] = 10;
    foo[5] =  23;    foo[6] = 29;    foo[7] = 16;    foo[8] = 22;    foo[9] = 17;
    foo[10] = 24;   foo[11] = 34;   foo[12] = 22;   foo[13] = 28;   foo[14] = 24;
    foo[15] = 7;    foo[16] = 17;   foo[17] = 6;    foo[18] = 13;   foo[19] = 11;
    foo[20] = 23;   foo[21] = 37;   foo[22] = 32;   foo[23] = 44;   foo[24] = 45;
    foo[25] = 6;    foo[26] = 21;   foo[27] = 19;   foo[28] = 33;   foo[29] = 36;
    foo[30] = 34;   foo[31] = 53;   foo[32] = 51;   foo[33] = 66;   foo[34] = 75;

    ComplexMatrix<int> m = ComplexMatrix<int>(7,5,foo);
    
    if (m.isMonge()) {
        cout << "Is Monge" << endl;
    }
    
    if (m.isInverseMonge()) {
        cout << "Is Inverse Monge" << endl;
    }
    
    RowEnvelope<int> *r3 = new RowEnvelope<int>(m,(size_t)3);
    RowEnvelope<int> *r5 = new RowEnvelope<int>(m,(size_t)5);
    RowEnvelope<int> * r35 = mergeRowEnvelopes(r3, r5);
    cout << "r35";
    printBreakpointList(r35->breakpoints());
    cout << endl;

    ColumnEnvelope<int> *e0 = new ColumnEnvelope<int>(m,(size_t)0);
    ColumnEnvelope<int> *e1 = new ColumnEnvelope<int>(m,(size_t)1);
    ColumnEnvelope<int> *e2 = new ColumnEnvelope<int>(m,(size_t)2);
    ColumnEnvelope<int> *e3 = new ColumnEnvelope<int>(m,(size_t)3);
    
    cout << "e0";
    printBreakpointList(e0->breakpoints());
    cout << endl;

    cout << "e1";
    printBreakpointList(e1->breakpoints());
    cout << endl;

    cout << "e2";
    printBreakpointList(e2->breakpoints());
    cout << endl;

    ColumnEnvelope<int> * e03 = mergeColumnEnvelopes(e0, e3);
    cout << "e03";
    printBreakpointList(e03->breakpoints());
    cout << endl;

    ColumnEnvelope<int> * e23 = mergeColumnEnvelopes(e2, e3);
    cout << "e23";
    printBreakpointList(e23->breakpoints());
    cout << endl;
}

void testRowTree()
{
    cout << "Test Row Tree" << endl;
    valarray<int> foo (35);
    foo[0] =  23;    foo[1] = 28;    foo[2] = 13;    foo[3] = 17;    foo[4] = 10;
    foo[5] =  23;    foo[6] = 29;    foo[7] = 16;    foo[8] = 22;    foo[9] = 17;
    foo[10] = 24;   foo[11] = 34;   foo[12] = 22;   foo[13] = 28;   foo[14] = 24;
    foo[15] = 7;    foo[16] = 17;   foo[17] = 6;    foo[18] = 13;   foo[19] = 11;
    foo[20] = 23;   foo[21] = 37;   foo[22] = 32;   foo[23] = 44;   foo[24] = 45;
    foo[25] = 6;    foo[26] = 21;   foo[27] = 19;   foo[28] = 33;   foo[29] = 36;
    foo[30] = 34;   foo[31] = 53;   foo[32] = 51;   foo[33] = 66;   foo[34] = 75;
    
    ComplexMatrix<int> m = ComplexMatrix<int>(7,5,foo);
    
    if (m.isMonge()) {
        cout << "Is Monge" << endl;
    }
    
    if (m.isInverseMonge()) {
        cout << "Is Inverse Monge" << endl;
    }
    
    ExtendedRowNode<int> root = ExtendedRowNode<int>(m);
    
    int max = root.maxForColumnInRange(2,2,2);
    
    cout << "max: " << max << endl;
    
}

void testColTree()
{
    cout << "Test Col Tree" << endl;
    valarray<int> foo (35);
    foo[0] =  23;    foo[1] = 28;    foo[2] = 13;    foo[3] = 17;    foo[4] = 10;
    foo[5] =  23;    foo[6] = 29;    foo[7] = 16;    foo[8] = 22;    foo[9] = 17;
    foo[10] = 24;   foo[11] = 34;   foo[12] = 22;   foo[13] = 28;   foo[14] = 24;
    foo[15] = 7;    foo[16] = 17;   foo[17] = 6;    foo[18] = 13;   foo[19] = 11;
    foo[20] = 23;   foo[21] = 37;   foo[22] = 32;   foo[23] = 44;   foo[24] = 45;
    foo[25] = 6;    foo[26] = 21;   foo[27] = 19;   foo[28] = 33;   foo[29] = 36;
    foo[30] = 34;   foo[31] = 53;   foo[32] = 51;   foo[33] = 66;   foo[34] = 75;
    
    ComplexMatrix<int> m = ComplexMatrix<int>(7,5,foo);
    
    if (m.isMonge()) {
        cout << "Is Monge" << endl;
    }
    
    if (m.isInverseMonge()) {
        cout << "Is Inverse Monge" << endl;
    }
    
    ColNode<int> root = ColNode<int>(m);
    
    int max = root.maxForRowInRange(1,0,4);
    
    cout << "max: " << max << endl;
    
}
void testMatrix()
{
    valarray<int> foo (12);
    for (int i=0; i<12; ++i) foo[i]=i;
    
    ComplexMatrix<int> m = ComplexMatrix<int>(3,4,foo);
    
    valarray<int> bar = m.row(1,0,2);
    
    cout << "row(1,0,2): ";
    for (size_t n=0; n<bar.size(); n++)
        cout << bar[n] << ' ';
    cout << endl;

    valarray<int> sub = m.submatrix(0, 2, 0, 3);
    
    cout << "submatrix(0, 2, 0, 3): ";
    for (size_t n=0; n<sub.size(); n++)
        cout << sub[n] << ' ';
    cout << endl;
}

void testMonge()
{
    valarray<int> foo (4);
    foo[0] = 11; foo[1] = 7;
    foo[2] = 17; foo[3] = 23;
    
    ComplexMatrix<int> m = ComplexMatrix<int>(2,2,foo);
    if (m.isMonge()) {
        cout << "Is Monge" << endl;
    }
    
    if (m.isInverseMonge()) {
        cout << "Is Inverse Monge" << endl;
    }
}

void  testSubmatrixQueries()
{
    cout << "Test Submatrix Queries" << endl;
    valarray<int> foo (25);

    foo[0] = 16054;     foo[1] = 11809; foo[2] = 7292;  foo[3] = 6225;  foo[4] = 1517;
    foo[5] = 14438;     foo[6] = 10193; foo[7] = 5676;  foo[8] = 4609;  foo[9] = 99;
    foo[10] = 11197;    foo[11] = 6952; foo[12] = 2435; foo[13] = 1368; foo[14] = 3340;
    foo[15] = 9786;     foo[16] = 5541; foo[17] = 1024; foo[18] = 43;   foo[19] = 4751;
    foo[20] = 9085; 	foo[21] = 4840; foo[22] = 323;  foo[23] = 744;  foo[24] = 5452;
    
    
    ComplexMatrix<int> m = ComplexMatrix<int>(5,5,foo);
    
    if (m.isMonge()) {
        cout << "Is Monge" << endl;
    }
    
    if (m.isInverseMonge()) {
        cout << "Is Inverse Monge" << endl;
    }
    
    SubmatrixQueriesDataStructure<int> structure = SubmatrixQueriesDataStructure<int>(m);
    
    int max = structure.maxInRange(0,4,1,3);
    
    cout << "max: " << max << endl;

}

void testTest( size_t nRows, size_t nCols)
{
    cout << "Test for " << nRows << " rows and " << nCols << " columns \n\n";
    
    SubmatrixQueriesTest test = SubmatrixQueriesTest(nRows, nCols);
    
    cout << "\n========================================";
    cout << "\nBeginning col queries tests ..." << endl;
    if (test.multipleColumnQueryTest(100)) {
        cout << "Col queries tests passed" << endl;
    }else{
        cout << "Tests failed" << endl;
    }
    cout << "\n";
    if (test.multipleColQueryTestVsCascading(1000)) {
        cout << "Cascading col queries tests passed" << endl;
    }else{
        cout << "Tests failed" << endl;
    }

    cout << "\n========================================";
    cout  << "\nBeginning row queries tests ..." << endl;
    if (test.multipleRowQueryTest(100)) {
        cout << "Row queries tests passed" << endl;
    }else{
        cout << "Tests failed" << endl;
    }
    cout << "\n";
    if (test.multipleRowQueryTestVsCascading(1000)) {
        cout << "Cascading row queries tests passed" << endl;
    }else{
        cout << "Tests failed" << endl;
    }

    cout << "\n========================================";
    cout << "\n\nBeginning submatrix queries tests ..." << endl;
    if (test.multipleSubmatrixQueryTest(50)) {
        cout << "Submatrix queries tests passed" << endl;
    }else{
        cout << "Tests failed" << endl;
    }
}

void benchmarks(size_t nRows, size_t nCols)
{
    cout << "Benchmarks for " << nRows << " rows and " << nCols << " columns \n\n";
    
    SubmatrixQueriesTest test = SubmatrixQueriesTest(nRows, nCols);
    
    cout << "\n========================================";
    cout << "\nBeginning column queries benchmarks ..." << endl;
    test.multipleBenchmarksColQueries(1000);
    
    cout << "\n========================================";
    cout << "\nBeginning row queries benchmarks ..." << endl;
    test.multipleBenchmarksRowQueries(1000);

}

void multiBenchmarks(size_t maxNRows, size_t maxNCols, size_t nSampleSize, size_t nSamplePerSize)
{
    cout << "Benchmarks on " << nSampleSize << " samples up to size " << maxNRows << " x " << maxNCols << ", "<< nSamplePerSize << " samples for eache size\n\n";
    
    clock_t **benchmarks;
    
    benchmarks = SubmatrixQueriesTest::multiSizeBenchmarksPositionQueries(maxNRows, maxNCols, nSampleSize, nSamplePerSize);
    
    for (size_t i = 0; i < nSampleSize; i++) {
        float fraction = ((float)(i+1))/((float)nSampleSize);
        cout << (size_t)(maxNRows*fraction) << " ; " << benchmarks[i][0] << " ; " << benchmarks[i][1] << " ; " << benchmarks[i][2] << " ; " << benchmarks[i][3] << "\n";
        delete [] benchmarks[i];
    }
    
    delete [] benchmarks;
}

int main(int argc, const char * argv[])
{
    size_t nRows = 10000; // default values for the number of columns and rows
    size_t nCols = 10000;
    
    if (argc >= 3) {
        sscanf(argv[1], "%lu", &nRows);
        sscanf(argv[2], "%lu", &nCols);
    }
    
//    testMatrix();
//    testMonge();
//    testEnvelope();
//    testRowTree();
//    testColTree();
//    testSubmatrixQueries();
    
//    testTest(nRows,nCols);
//    benchmarks(nRows, nCols);

    multiBenchmarks(nRows, nCols, ((float)nRows)/((float) 20),50);
    
    return 0;
}
