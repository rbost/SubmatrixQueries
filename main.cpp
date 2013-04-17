//
//  main.cpp
//  SubmatrixQueries
//
//  Created by Raphael Bost on 02/01/13.
//  Copyright (c) 2013 Raphael Bost. All rights reserved.
//

#include <iostream>
#include <cstdio>
#include <fstream>
#include <cstring>

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
    
    cout << "\n========================================";
    cout << "\nBeginning submatrix queries benchmarks ..." << endl;
    test.multipleBenchmarksSubmatrixQueries(1000);
    
}

void multiSizePositionQueriesBenchmarks(size_t maxNRows, size_t maxNCols, size_t nSampleSize, size_t nSamplePerSize)
{
    cout << "Benchmarks on " << nSampleSize << " samples up to size " << maxNRows << " x " << maxNCols << ", "<< nSamplePerSize << " samples for each size\n\n";
    
    bench_time_t **benchmarks;
    
    benchmarks = SubmatrixQueriesTest::multiSizeBenchmarksPositionQueries(maxNRows, maxNCols, nSampleSize, nSamplePerSize);
    
    for (size_t i = 0; i < nSampleSize; i++) {
        float fraction = ((float)(i+1))/((float)nSampleSize);
        cout << (size_t)(maxNRows*fraction) << " ; " << benchTimeAsMiliSeconds(benchmarks[i][0]) << " ; " << benchTimeAsMiliSeconds(benchmarks[i][1]) << " ; " << benchTimeAsMiliSeconds(benchmarks[i][2]) << " ; " << benchTimeAsMiliSeconds(benchmarks[i][3]) << "\n";
        delete [] benchmarks[i];
    }
    
    delete [] benchmarks;
}

void multiSizeSubmatrixQueriesBenchmarks(size_t maxNRows, size_t maxNCols, size_t minNRows, size_t minNCols,  size_t stepSize, size_t nSamplePerSize)
{
    if(minNRows == 0) minNRows += stepSize;
    if(minNCols == 0) minNCols += stepSize;
    
    cout << "Benchmarks on samples from size " << minNRows << " x " << minNCols << " up to size " << maxNRows << " x " << maxNCols << ", "<< nSamplePerSize << " samples for each size\n\n";
    bench_time_t **benchmarks;
    
    size_t sampleSize;
    benchmarks = SubmatrixQueriesTest::multiSizeBenchmarksSubmatrixQueries(maxNRows, maxNCols,minNRows, minNCols, stepSize, nSamplePerSize, &sampleSize);
    
    for (size_t i = 0; i < sampleSize; i++) {
        float fraction = ((float)(i+1))/((float)sampleSize);
        cout << (size_t)(maxNRows*fraction) << " ; " << benchTimeAsMiliSeconds(benchmarks[i][0]) << " ; " << benchTimeAsMiliSeconds(benchmarks[i][1]) << " ; " << benchTimeAsMiliSeconds(benchmarks[i][2]) << "\n";
        delete [] benchmarks[i];
    }
    
    delete [] benchmarks;

}

void multiSizeSubmatrixQueriesBenchmarks(size_t maxNRows, size_t maxNCols, size_t minNRows, size_t minNCols, size_t stepSize, size_t nSamplePerSize, const char* outputFilename)
{
    if(minNRows == 0) minNRows += stepSize;
    if(minNCols == 0) minNCols += stepSize;

    cout << "Benchmarks on samples from size " << minNRows << " x " << minNCols << " up to size " << maxNRows << " x " << maxNCols << ", "<< nSamplePerSize << " samples for each size\n\n";
    ofstream output (outputFilename,ios::out | ios::trunc);
    
    if (output.is_open()) {
        SubmatrixQueriesTest::multiSizeBenchmarksSubmatrixQueries(maxNRows, maxNCols, minNRows, minNCols, stepSize, nSamplePerSize,output);
    }else{
        cout << "Fail to open the benchmarking output file.\n";
    }
    output.close();
}

void testInitalization(size_t nRows, size_t nCols)
{
    cout << "Initialization test for " << nRows << " rows and " << nCols << " columns \n\n";
    
    SubmatrixQueriesTest test = SubmatrixQueriesTest(nRows, nCols);
}

void multiSizeFastestQueriesBenchmarks(size_t maxNRows, size_t maxNCols, size_t minNRows, size_t minNCols, size_t stepSize, size_t nSamplePerSize, size_t nPosQueries, size_t nSMQueries, const char* outputFilename)
{
    if(minNRows == 0) minNRows += stepSize;
    if(minNCols == 0) minNCols += stepSize;
    
    cout << "Benchmarks on samples from size " << minNRows << " x " << minNCols << " up to size " << maxNRows << " x " << maxNCols << ", "<< nSamplePerSize << " samples for each size\n";
    cout << nSMQueries << "submatrix queries, " << nPosQueries << "position queries for each sample\n\n" << flush;
    
    ofstream output (outputFilename,ios::out | ios::trunc);
    
    if (outputFilename == NULL) {
        SubmatrixQueriesTest::multiSizeBenchmarkBestPositionAndSubmatrixQueries(maxNRows, maxNCols, minNRows, minNCols, stepSize, nSamplePerSize, nPosQueries, nSMQueries,cout);
    }else if (output.is_open()) {
        SubmatrixQueriesTest::multiSizeBenchmarkBestPositionAndSubmatrixQueries(maxNRows, maxNCols, minNRows, minNCols, stepSize, nSamplePerSize, nPosQueries, nSMQueries,output);
        output.close();
    }else{
        cout << "Fail to open the benchmarking output file.\n";
    }
}

void envelopeSizeStatistics(size_t maxN, size_t minN, size_t stepSize, bool logSteps, size_t nSamplePerSize, const char* outputFilename)
{
    if(minN == 0){
        if(logSteps) minN = 2;
        else minN += stepSize;
    }
    
    ofstream output (outputFilename,ios::out | ios::trunc);
    
    if (outputFilename == NULL) {
        SubmatrixQueriesTest::envelopeSizesStats(maxN, minN, stepSize, logSteps, nSamplePerSize, cout);
    }else if (output.is_open()) {
        SubmatrixQueriesTest::envelopeSizesStats(maxN, minN, stepSize, logSteps, nSamplePerSize,output);
        output.close();
    }else{
        cout << "Fail to open the benchmarking output file.\n";
    }

}

int main(int argc, const char * argv[])
{
    size_t nRows = 10000; // default values for the number of columns and rows
    size_t nCols = 10000;
    size_t minRows = 0, minCols = 0;
    size_t stepSize = 50;
    size_t samplesPerSize = 50;
    
    int mode = -1; // for default behavior
    
    const char* filename = NULL;
    bool externalOutput = false;

    if (argc >= 3) {
        sscanf(argv[1], "%lu", &nRows);
        sscanf(argv[2], "%lu", &nCols);
    }

    int i = 3;

    while (i < argc) {
        if (strcmp("-o", argv[i]) == 0) {
            assert(i+1 < argc);
            filename = argv[i+1];
            externalOutput = true;
            i = i+2;
        }else if (strcmp("--minSize", argv[i]) == 0){
            assert(i+2 < argc);
            sscanf(argv[i+1], "%lu", &minRows);
            sscanf(argv[i+2], "%lu", &minCols);
            i = i+3;
        }else if (strcmp("--step", argv[i]) == 0){
            assert(i+1 < argc);
            sscanf(argv[i+1], "%lu", &stepSize);
            i = i+2;
        }else if (strcmp("--samples-per-size", argv[i]) == 0){
            assert(i+1 < argc);
            sscanf(argv[i+1], "%lu", &samplesPerSize);
            i = i+2;
        }else if (strcmp("--mode", argv[i]) == 0){
            assert(i+1 < argc);
            sscanf(argv[i+1], "%d", &mode);
            i = i+2;
        }
    }

    
    SubmatrixQueriesTest::benchmarkNaiveQueries = false;
    SubmatrixQueriesTest::verboseBenchmarks = false;
    SubmatrixQueriesTest::verboseMatrixGeneration = false;
    
    switch (mode) {
        case 0:
            SubmatrixQueriesTest::verboseMatrixGeneration = true;
            SubmatrixQueriesTest::verboseBenchmarks = true;
            testInitalization(nRows,nCols);
            break;
            
        case 1:
            SubmatrixQueriesTest::verboseBenchmarks = true;
            SubmatrixQueriesTest::verboseMatrixGeneration = true;
            testTest(nRows,nCols);
            break;

        case 2:
            multiSizePositionQueriesBenchmarks(nRows, nCols, ((float)nRows)/((float) stepSize),samplesPerSize);
            break;

        case 4:
            multiSizeFastestQueriesBenchmarks(nRows, nCols, minRows, minCols,stepSize,samplesPerSize, 1000,100, filename);
            break;
            
        case 5:
            envelopeSizeStatistics(max(nRows,nCols),min(minRows,minCols),stepSize,true,samplesPerSize,filename);
            break;
            
        case 3:
        default:
            if (externalOutput) {
                multiSizeSubmatrixQueriesBenchmarks(nRows, nCols, minRows, minCols,stepSize,samplesPerSize,filename);
            }else{
                multiSizeSubmatrixQueriesBenchmarks(nRows, nCols, minRows, minCols,stepSize,samplesPerSize);
            }
            
            break;
    }
    

    return 0;
}
