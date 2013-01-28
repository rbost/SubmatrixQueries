//
//  enveloppe_tree.h
//  KMNS
//
//  Created by Raphael Bost on 08/01/13.
//  Copyright (c) 2013 Raphael Bost. All rights reserved.
//

#ifndef __KMNS__envelope_tree__
#define __KMNS__envelope_tree__

#include <vector>
#include <cassert>
#include <iostream>

#include "envelope.h"
#include "matrix.h"
#include "max_value.h"

using namespace envelope;
using namespace matrix;

template <typename T>
class RowNode {
private:
    size_t _minRow, _maxRow;
    
    RowNode<T> *_lowIndicesNode, *_highIndicesNode;
    bool _isLeaf;

    RowEnvelope<T> *_envelope;
    
protected:
    void setEnvelope(RowEnvelope< T > *newEnvelope){
        _envelope = newEnvelope;
    }
    
    void setLowIndicesNode(RowNode<T> *newLowNode)
    {
        _lowIndicesNode = newLowNode;
    }
    
    void setHighIndicesNode(RowNode<T> *newHighNode)
    {
        _highIndicesNode = newHighNode;
    }
    
public:
    RowNode(size_t minRow, size_t maxRow,RowNode<T> *lowIndices, RowNode<T> *highIndices ,Matrix<T> const& matrix):
    _minRow(minRow),_maxRow(maxRow)
    {
        assert(minRow <= maxRow);
        if (minRow == maxRow) { // it is a leaf
            _isLeaf = true;
            _envelope = new RowEnvelope<T>(matrix, minRow);
        }else{
            _isLeaf = false;
            
            _lowIndicesNode = lowIndices;
            _highIndicesNode = highIndices;
        }
    }

    RowNode(size_t minRow, size_t maxRow, Matrix<T> const& matrix):
    _minRow(minRow),_maxRow(maxRow)
    {
        assert(minRow <= maxRow);
        if (minRow == maxRow) { // it is a leaf
            _isLeaf = true;
            _envelope = new RowEnvelope<T>(matrix, minRow);
        }else{
            _isLeaf = false;
            
            size_t midRow = minRow + ((maxRow - minRow)/2);
            
            _lowIndicesNode = new RowNode<T>(minRow,midRow,matrix);
            _highIndicesNode= new RowNode<T>(midRow+1,maxRow,matrix);
            
            this->setEnvelope(mergeRowEnvelopes(this->lowIndicesNode()->envelope(), this->highIndicesNode()->envelope()));
        }
    }
    

    RowNode(Matrix<T> const& matrix) : _minRow(0), _maxRow(matrix.rows()-1)
    {
        _isLeaf = false;
        
        size_t minRow = this-> minRow();
        size_t maxRow = this->maxRow();
        size_t midRow = minRow + ((maxRow - minRow)/2);
        
        _lowIndicesNode = new RowNode<T>(minRow,midRow,matrix);
        _highIndicesNode= new RowNode<T>(midRow+1,maxRow,matrix);
        
        this->setEnvelope(mergeRowEnvelopes(this->lowIndicesNode()->envelope(), this->highIndicesNode()->envelope()));
    }
    
    ~RowNode()
    {
        if(!_isLeaf )
        {
            delete _lowIndicesNode;
            delete _highIndicesNode;
        }
        delete _envelope;
    }
    
    RowEnvelope<T>* envelope() const
    {
        return _envelope;
    }
    
    bool isLeaf() const { return _isLeaf; }
    virtual RowNode<T>* lowIndicesNode() const { return _lowIndicesNode; }
    virtual RowNode<T>* highIndicesNode() const { return _highIndicesNode; }
    size_t minRow() const { return _minRow; }
    size_t maxRow() const { return _maxRow; }
    
    // Returns the canonical nodes (cf. the article) for the specified indices
    vector<const RowNode<T> *> canonicalNodes(size_t minRow, size_t maxRow) const
    {
        std::vector<const RowNode<T> *> buffer;
        this->getCanonicalNodes(buffer, minRow, maxRow);
        
        return buffer;
    }
    
    virtual void getCanonicalNodes(std::vector<const RowNode<T> *> & buffer, size_t minRow, size_t maxRow) const
    {
        assert(minRow <= maxRow);
        if (minRow > this->maxRow() || maxRow < this->minRow()) { // check if the interval intersects the node's rows
            return; // if not, exit
        }
        
        if (minRow <= this->minRow() && maxRow >= this->maxRow()) { // check if the interval includes the node's rows
            buffer.push_back(this); // in this the case, add the entire node to the buffer
            return;
        }
        
        if(!this->isLeaf()){
            this->lowIndicesNode()->getCanonicalNodes(buffer,minRow,maxRow);
            this->highIndicesNode()->getCanonicalNodes(buffer,minRow,maxRow);
        }
    }
    
    // Returns the maximum value of the matrix in the specified column and in the specified row range
    T maxForColumnInRange(size_t col, size_t minRow, size_t maxRow) const{
        std::vector<const RowNode<T> *> cNodes = this->canonicalNodes(minRow,maxRow);
        
        assert(cNodes.size() > 0);
        
        T max = cNodes[0]->envelope()->valueForColumn(col);
        
        for (size_t i = 1; i < cNodes.size(); i++) {
            T value = cNodes[i]->envelope()->valueForColumn(col);
            
            if (value > max) {
                max = value;
            }
        }
        
        return max;
    }
};
/*
 * class ColNode
 *
 * This class provides the same data structure implementation as for RowNode but for queries on rows and column ranges instead of queries on columns and row ranges
 */
template <typename T>
class ColNode {

    size_t _minCol, _maxCol;
    
    ColNode<T> *_lowIndicesNode, *_highIndicesNode;
    bool _isLeaf;
    
    ColumnEnvelope<T> *_envelope;
    
public:
    
    ColNode(size_t minCol, size_t maxCol, Matrix<T> const& matrix):
    _minCol(minCol),_maxCol(maxCol)
    {
        assert(minCol <= maxCol);
        if (minCol == maxCol) { // it is a leaf
            _isLeaf = true;
            _envelope = new ColumnEnvelope<T>(matrix, minCol);
        }else{
            _isLeaf = false;
            size_t midRow = minCol + ((maxCol - minCol)/2);
            
            _lowIndicesNode = new ColNode<T>(minCol,midRow,matrix);
            _highIndicesNode= new ColNode<T>(midRow+1,maxCol,matrix);
            
            _envelope = mergeColumnEnvelopes(this->lowIndicesNode()->envelope(), this->highIndicesNode()->envelope());
        }
    }
    
    
    ColNode(Matrix<T> const& matrix) : _minCol(0), _maxCol(matrix.cols()-1)
    {
        size_t midRow = _minCol + ((_maxCol - _minCol)/2);
        
        _lowIndicesNode = new ColNode<T>(_minCol,midRow,matrix);
        _highIndicesNode= new ColNode<T>(midRow+1,_maxCol,matrix);
        
        _envelope = mergeColumnEnvelopes(this->lowIndicesNode()->envelope(), this->highIndicesNode()->envelope());
    }
    
    ~ColNode()
    {
        if(!this->isLeaf() )
        {
            delete _lowIndicesNode;
            delete _highIndicesNode;
        }
        delete _envelope;
    }
    
    bool isLeaf() const { return _isLeaf; }
    size_t minCol() const { return _minCol; }
    size_t maxCol() const { return _maxCol; }
    
    ColumnEnvelope<T>* envelope() const
    {
        return _envelope;
    }
    ColNode<T>* lowIndicesNode() const { return _lowIndicesNode; }
    ColNode<T>* highIndicesNode() const { return _highIndicesNode; }

    // Returns the canonical nodes (cf. the article) for the specified indices
    std::vector<const ColNode<T> *> canonicalNodes(size_t minCol, size_t maxCol) const
    {
        std::vector<const ColNode<T> *> buffer;
        getCanonicalNodes(buffer, minCol, maxCol);
        
        return buffer;
    }
    
    void getCanonicalNodes(std::vector<const ColNode<T> *> & buffer, size_t minCol, size_t maxCol) const
    {
        assert(minCol <= maxCol);
        if (minCol > this->maxCol() || maxCol <  this->minCol()) { // check if the interval intersects the node's rows
            return;
        }
        
        if (minCol <=  this->minCol() && maxCol >= this->maxCol()) { // check if the interval include the node's rows
            buffer.push_back(this); // in this the case, add the entire node to the buffer
            return;
        }
        
        if(!this->isLeaf()){
            this->lowIndicesNode()->getCanonicalNodes(buffer,minCol,maxCol);
            this->highIndicesNode()->getCanonicalNodes(buffer,minCol,maxCol);
        }
    }
    
    T maxForRowInRange(size_t row, size_t minCol, size_t maxCol) const{
        std::vector<const ColNode<T> *> cNodes = canonicalNodes(minCol,maxCol);
        
        assert(cNodes.size() > 0);
        
        T max = cNodes[0]->envelope()->valueForRow(row);
        
        for (size_t i = 1; i < cNodes.size(); i++) {
            T value = cNodes[i]->envelope()->valueForRow(row);
            
            if (value > max) {
                max = value;
            }
        }
        
        return max;
    }
};

/*
 * class ExtendedRowNode
 *
 * This subclass of RowNode adds a _maxima field to store the breakpoints' intervals maximum values.
 *
 */
template <typename T>
class ExtendedRowNode : public RowNode<T> {
    vector< T > *_maxima; // the _maxima vector stores the maxima of breakpoints intervals
    
public:
    ExtendedRowNode(size_t minRow, size_t maxRow, Matrix<T> const& matrix): RowNode<T>(minRow,maxRow,NULL,NULL,matrix)
    {
        if (minRow < maxRow){
            size_t midRow = minRow + ((maxRow - minRow)/2);
            this->setLowIndicesNode(new ExtendedRowNode<T>(minRow,midRow,matrix));
            this->setHighIndicesNode(new ExtendedRowNode<T>(midRow+1,maxRow,matrix));
            this->setEnvelope(mergeRowEnvelopes(this->lowIndicesNode()->envelope(), this->highIndicesNode()->envelope()));
        }
    }
    
    ExtendedRowNode(Matrix<T> const& matrix): RowNode<T>(0,matrix.rows()-1,NULL,NULL,matrix)
    {
        size_t minRow = this-> minRow();
        size_t maxRow = this->maxRow();
        size_t midRow = minRow + ((maxRow - minRow)/2);
        
        this->setLowIndicesNode(new ExtendedRowNode<T>(minRow,midRow,matrix));
        this->setHighIndicesNode(new ExtendedRowNode<T>(midRow+1,maxRow,matrix));
        this->setEnvelope(mergeRowEnvelopes(this->lowIndicesNode()->envelope(), this->highIndicesNode()->envelope()));
    }
    
    ~ExtendedRowNode<T>()
    {
        delete _maxima;
    }
    
    virtual ExtendedRowNode<T>* lowIndicesNode() const { return (ExtendedRowNode<T>*) RowNode<T>::lowIndicesNode(); }
    virtual ExtendedRowNode<T>* highIndicesNode() const { return (ExtendedRowNode<T>*) RowNode<T>::highIndicesNode(); }

    const vector< T > *maxima() const
    {
        return _maxima;
    }
    
    // This method computes the _maxima vector
    void computeIntervalMaxima(const ColNode<T> *flippedTree)
    {
        const vector< Breakpoint > *breakpoints = this->envelope()->breakpoints();
        size_t n = breakpoints->size();
        
        _maxima = new vector< T >(n); // create a new vector of the same size than the breakpoints one
    
        size_t minCol, maxCol;
        size_t row;
        size_t i;
        
        for (i = 0; i < n-1; i++) { // for every interval ...
            
            minCol = (*breakpoints)[i].col; // ... get the interval first index ...
            maxCol = (*breakpoints)[i+1].col-1; // ... its last index ...
            row = (*breakpoints)[i].row; // ... and the corresponding row ...

            (*_maxima)[i] = flippedTree->maxForRowInRange(row,minCol,maxCol); // ... and finally compute the maximum for the column range
        }
        
        // Do not forget the last interval!
        minCol = (*breakpoints)[i].col;
        maxCol = this->envelope()->maxPosition()-1;
        row = (*breakpoints)[i].row;
        (*_maxima)[i] = flippedTree->maxForRowInRange(row,minCol,maxCol);

    }
    
    // Computes the interval maxima and tells the node's children to do the same
    void recursivelyComputeIntervalMaxima(const ColNode<T> *flippedTree)
    {
        computeIntervalMaxima(flippedTree);
        
        if (!this->isLeaf()) {
            this->lowIndicesNode()->computeIntervalMaxima(flippedTree);
            this->highIndicesNode()->computeIntervalMaxima(flippedTree);
        }
    }
    
    // We override this to avoid compile-time errors (thank you C++) 
    vector<const ExtendedRowNode<T> *> canonicalNodes(size_t minRow, size_t maxRow) const
    {
        std::vector<const ExtendedRowNode<T> *> buffer;
        getCanonicalNodes(buffer, minRow, maxRow);
        
        return buffer;
    }
    
    void getCanonicalNodes(std::vector<const ExtendedRowNode<T> *> & buffer, size_t minRow, size_t maxRow) const
    {
        assert(minRow <= maxRow);
        if (minRow > this->maxRow() || maxRow < this->minRow()) { // check if the interval intersects the node's rows
            return; // if not, exit
        }
        
        if (minRow <= this->minRow() && maxRow >= this->maxRow()) { // check if the interval includes the node's rows
            buffer.push_back(this); // in this the case, add the entire node to the buffer
            return;
        }
        
        if(!this->isLeaf()){
            this->lowIndicesNode()->getCanonicalNodes(buffer,minRow,maxRow);
            this->highIndicesNode()->getCanonicalNodes(buffer,minRow,maxRow);
        }
    }

};

/*
 * class SubmatrixQueriesDataStructure
 *
 * This is THE most important class as it combines the different classes and query method for the implementation of the fast data structure for maximum queries on submatrices.
 *
 */

template <typename T>
class SubmatrixQueriesDataStructure {
    ExtendedRowNode<T> *_rowsTree; // denoted T_h in the KMNS article
    ColNode<T> *_columnTree; // denoted \mathcal{B} in the KMNS article
    
public:
    SubmatrixQueriesDataStructure(Matrix<T> const& matrix)
    {
        _rowsTree = new ExtendedRowNode<T>(matrix);
        _columnTree = new ColNode<T>(matrix);
        
        _rowsTree->recursivelyComputeIntervalMaxima(_columnTree);
    }
    
    ~SubmatrixQueriesDataStructure()
    {
        delete _rowsTree;
        delete _columnTree;
    }
    
    T maxInRange(size_t minRow, size_t maxRow, size_t minCol, size_t maxCol) const
    {
        return maxInRange(Range(minRow,maxRow), Range(minCol,maxCol));
    }
    
    T maxInRange(Range rowRange, Range colRanges) const
    {
        vector<const ExtendedRowNode<T> *> rowNodes = _rowsTree->canonicalNodes(rowRange.min,rowRange.max);
        
        MaxValue<T> max;
        for (typename vector<const ExtendedRowNode<T> *>::iterator nodesIterator = rowNodes.begin(); nodesIterator != rowNodes.end(); ++nodesIterator) {
            

            RowEnvelope<T> *envelope = (*nodesIterator)->envelope();
            const vector< Breakpoint > *breakpoints = envelope->breakpoints();
            
            size_t i, numberOfBp = breakpoints->size();

            for (i = 0; i < numberOfBp; i++) {
                Breakpoint bp = (*breakpoints)[i];
                if (bp.col >= colRanges.min) {
                    break;
                }
            }
            
            
            // first degenerated case: the range is after the last breakpoint
            if (i == numberOfBp) {
                size_t row = (*breakpoints)[i-1].row;
                max.updateMax( _columnTree->maxForRowInRange(row,colRanges.min,colRanges.max));
                continue;
            }
            
            // at this point, we have the prefix: it is the interval [colRanges.min,(*breakpoints)[i].col-1]
            // we first have to check if it is empty or not ...
            if ((*breakpoints)[i].col > colRanges.min) {
                // it is not empty, go on ...
                size_t row = (*breakpoints)[i-1].row;
                max.updateMax(_columnTree->maxForRowInRange(row,colRanges.min, (*breakpoints)[i].col-1));
            }
            
            // now, we check for the fully contained intervals
            // in the case we still have some breakpoints to explore ...
            if(i < numberOfBp){
                
                for (i = i+1; i < numberOfBp; i++) {
                    // we are considering the interval between (*breakpoints)[i-1].col and (*breakpoints)[i].col-1
                    // we know that (*breakpoints)[i-1].col is in the column range so we only have to check for (*breakpoints)[i].col-1 
                    
                    Breakpoint bp = (*breakpoints)[i];
                    if (bp.col-1 <= colRanges.max) {
                        // alright, the interval is fully contained in the column range
                        const vector< T > *maxima = (*nodesIterator)->maxima();
                        // the maximum of interval [(*breakpoints)[i-1].col, (*breakpoints)[i].col-1] is in maxima[i-1]
                        max.updateMax((*maxima)[i-1]);
                    }else{
                        break;
                    }
                }
                
            }
                // check for the rest of the range
                // the remaining of the interval is between the last breakpoint column and the column range maximum :
                // the suffix is of range [(*breakpoints)[i-1].col,colRanges.max]
                // ( even if there are two ways to exit the loop, in the end we have to do the same processing for the remaining )
                
                size_t row = (*breakpoints)[i-1].row;
                max.updateMax(_columnTree->maxForRowInRange(row,(*breakpoints)[i-1].col,colRanges.max));
        }
        
        return max.value();
    }
};
#endif /* defined(__KMNS__enveloppe_tree__) */
