//
//  enveloppe_tree.h
//  SubmatrixQueries
//
//  Created by Raphael Bost on 08/01/13.
//  Copyright (c) 2013 Raphael Bost. All rights reserved.
//

#ifndef __SubmatrixQueries__envelope_tree__
#define __SubmatrixQueries__envelope_tree__

#include <vector>
#include <iostream>

#include "envelope.h"
#include "range.h"
#include "matrix.h"
#include "max_value.h"
#include "range_query.h"
#include "debug_assert.h"

using namespace envelope;
using namespace matrix;


template <typename T>
class EnvTreeNode {
private:
    Range _range;
    
    EnvTreeNode<T> *_lowIndicesNode, *_highIndicesNode;
 
    Envelope<T> *_envelope;
    size_t _crossingBpIndex; // the index of the breakpoint inserted when merging the children envelopes

protected:
    void setEnvelope(Envelope< T > *newEnvelope){
        _envelope = newEnvelope;
    }
    
    void setLowIndicesNode(EnvTreeNode<T> *newLowNode)
    {
        _lowIndicesNode = newLowNode;
    }
    
    void setHighIndicesNode(EnvTreeNode<T> *newHighNode)
    {
        _highIndicesNode = newHighNode;
    }
    
    // _crossingBpIndex is mutable only for the subclasses
    size_t& crossingBreakpointIndex() { return _crossingBpIndex; }

public:

    // This constructor creates a new RowNode with the specified children.
    // If they are not NULL, it will also compute the merged envelope.
    EnvTreeNode(Range r): EnvTreeNode(r.min,r.max)
    {
    }
    
    EnvTreeNode(size_t minRow, size_t maxRow): _range(minRow,maxRow)
    {
        assert(minRow <= maxRow);
    }

    virtual ~EnvTreeNode()
    {
        if(!(this->isLeaf()) )
        {
            delete _lowIndicesNode;
            delete _highIndicesNode;
        }
        delete _envelope;
    }
    
    inline Envelope<T>* envelope() const
    {
        return _envelope;
    }
    
    inline bool isLeaf() const {
        Range r = this->range();
        return (r.max - r.min) <= 0;
    }
    inline virtual EnvTreeNode<T>* lowIndicesNode() const { return _lowIndicesNode; }
    inline virtual EnvTreeNode<T>* highIndicesNode() const { return _highIndicesNode; }
    inline Range range() const { return _range; }
    inline size_t crossingBreakpointIndex() const { return _crossingBpIndex; }
    
    // Returns the canonical nodes (cf. the article) for the specified indices
    // COMPLEXITY : O(log(number of rows)
    vector<const EnvTreeNode<T> *> canonicalNodes(size_t minIndex, size_t maxIndex) const
    {
        std::vector<const EnvTreeNode<T> *> buffer;
        this->getCanonicalNodes(buffer, minIndex, maxIndex);
        
        return buffer;
    }
    
    // Auxiliary method for the previous one.
    // It adds itself to the buffer if the query range contains the row range the node represents.
    // Otherwise, it recursively calls its children.
    void getCanonicalNodes(std::vector<const EnvTreeNode<T> *> & buffer, size_t minIndex, size_t maxIndex) const
    {
        DEBUG_ASSERT(minIndex <= maxIndex);
        Range r = this->range();
        
        if (minIndex > r.max || maxIndex < r.min) { // check if the interval intersects the node's rows
            return; // if not, exit
        }
        
        if (minIndex <= r.min && maxIndex >= r.max) { // check if the interval includes the node's rows
            buffer.push_back(this); // in this the case, add the entire node to the buffer
            return;
        }
        
        if(!this->isLeaf()){
            this->lowIndicesNode()->getCanonicalNodes(buffer,minIndex,maxIndex);
            this->highIndicesNode()->getCanonicalNodes(buffer,minIndex,maxIndex);
        }
    }
    
    // Returns the maximum value of the matrix in the specified column and in the specified row range
    T maxInRange(size_t position, Range r) const{        
        // First of all, we get the canonical nodes
        std::vector<const EnvTreeNode<T> *> cNodes = this->canonicalNodes(r.min,r.max);
        
        // If the set of canonical nodes is empty, it means that the query range is empty
        DEBUG_ASSERT(cNodes.size() > 0);
        
        // Compute the maximum over the canonical nodes
        T max = cNodes[0]->envelope()->valueForPosition(position);
        
        for (size_t i = 1; i < cNodes.size(); i++) {
            T value = cNodes[i]->envelope()->valueForPosition(position);
            
            if (value > max) {
                max = value;
            }
        }
        
        return max;
    }
    
    size_t maxEnvelopeSize() const
    {
        if(this->isLeaf()){
            return this->envelope()->numberOfBreakpoints();
        }
        
        return max( max(this->lowIndicesNode()->maxEnvelopeSize(),this->highIndicesNode()->maxEnvelopeSize()),
                   this->envelope()->numberOfBreakpoints());
    }
    
    size_t minEnvelopeSize() const
    {
        if(this->isLeaf()){
            return this->envelope()->numberOfBreakpoints();
        }
        
        return min( min(this->lowIndicesNode()->minEnvelopeSize(),this->highIndicesNode()->minEnvelopeSize()),
                   this->envelope()->numberOfBreakpoints());
    }
    
    T cascadingMaxInRange(size_t position, Range r) const{
        MaxValue<T> max;
        
        this->updateRecursiveMaxInRange(position, r, &max, 0, this->envelope()->numberOfBreakpoints()-1);
        
        return max.value();
    }
    
    void updateRecursiveMaxInRange(size_t position, Range r, MaxValue<T>* max, size_t iMin, size_t iMax) const{
        if (!r.intersects(this->range())) {
            return;
        }
            
        size_t bpIndex;
        Breakpoint bp;
        
        if (iMin == iMax) {
            bpIndex = iMin;
            bp = (*this->envelope()->breakpoints())[bpIndex];
        }else{
            bp = this->envelope()->breakpointBeforePosition(position,iMin,iMax,&bpIndex);
        }
        
        if (r.isInRange(this->envelope()->mappedPositionForBreakpoint(bp))) {
            // If the breakpoint mapped position is in the range, we have the max for this node
            // Just update it an return
            T value = this->envelope()->valueForPositionAfterBreakpoint(position,bp);
            
            max->updateMax(value);
            return;
        }
        
        // Here, the returned breakpoint is out of range
        // We will have to make a recursive call on the node's children
        
        if (!this->isLeaf()) {
            
            size_t crossingIndex = this->crossingBreakpointIndex();

            // first of all, we treat some degenerate cases
            if (crossingIndex <= 0) {
                // the high indices node envelope is on the top
                this->lowIndicesNode()->updateRecursiveMaxInRange(position,r,max,0,this->lowIndicesNode()->envelope()->numberOfBreakpoints() -1);
                this->highIndicesNode()->updateRecursiveMaxInRange(position,r,max,bpIndex,bpIndex);
                return;
            }else if (crossingIndex >= this->envelope()->numberOfBreakpoints()){
                // the low indices node envelope is on the top
                this->lowIndicesNode()->updateRecursiveMaxInRange(position,r,max,bpIndex,bpIndex);
                this->highIndicesNode()->updateRecursiveMaxInRange(position,r,max,0,this->highIndicesNode()->envelope()->numberOfBreakpoints() -1);
                return;
            }

            if (bpIndex < crossingIndex) {
                // bp belongs to the part of the envelope that comes from the low indices node
                // That means, we do not have to search for that breakpoint in this child
                
                size_t reverseIndex = this->envelope()->numberOfBreakpoints() -1 - crossingIndex;
                size_t indexInHIN = this->highIndicesNode()->envelope()->numberOfBreakpoints() -1 - reverseIndex;

                this->lowIndicesNode()->updateRecursiveMaxInRange(position,r,max,bpIndex,bpIndex);
                this->highIndicesNode()->updateRecursiveMaxInRange(position,r,max,0,indexInHIN);
                return;
            }
            if (bpIndex > crossingIndex){
                size_t reverseIndex = this->envelope()->numberOfBreakpoints() -1 - bpIndex;
                size_t indexInHIN = this->highIndicesNode()->envelope()->numberOfBreakpoints() -1 - reverseIndex;

                this->lowIndicesNode()->updateRecursiveMaxInRange(position,r,max,crossingIndex-1,this->lowIndicesNode()->envelope()->numberOfBreakpoints() -1);
                this->highIndicesNode()->updateRecursiveMaxInRange(position,r,max,indexInHIN,indexInHIN);
                
                return;
            }
            if (bpIndex == crossingIndex) {
                size_t reverseIndex = this->envelope()->numberOfBreakpoints() -1 - crossingIndex;
                size_t indexInHIN = this->highIndicesNode()->envelope()->numberOfBreakpoints() -1 - reverseIndex;
                
                this->lowIndicesNode()->updateRecursiveMaxInRange(position,r,max,crossingIndex-1,this->lowIndicesNode()->envelope()->numberOfBreakpoints() -1);
                if (indexInHIN == 0) {
                    this->highIndicesNode()->updateRecursiveMaxInRange(position,r,max,0,indexInHIN);
                }else{
                    this->highIndicesNode()->updateRecursiveMaxInRange(position,r,max,indexInHIN-1,indexInHIN);
                }
                return;
            }

            this->lowIndicesNode()->updateRecursiveMaxInRange(position,r,max,0,this->lowIndicesNode()->envelope()->numberOfBreakpoints() -1);
            this->highIndicesNode()->updateRecursiveMaxInRange(position,r,max,0,this->highIndicesNode()->envelope()->numberOfBreakpoints() -1);
        }
    }
    
    T simpleCascadingMaxInRange(size_t position, Range r) const{
        MaxValue<T> max;
        
        this->updateRecursiveMaxInRange(position, r, &max);
        
        return max.value();
    }
    
    void updateRecursiveMaxInRange(size_t position, Range r, MaxValue<T>* max) const{
        if (!r.intersects(this->range())) {
            return;
        }
        
        Breakpoint bp = this->envelope()->breakpointBeforePosition(position,NULL);
        
        if (r.isInRange(this->envelope()->mappedPositionForBreakpoint(bp))) {
            // If the breakpoint mapped position is in the range, we have the max for this node
            // Just update it an return
            T value = this->envelope()->valueForPositionAfterBreakpoint(position,bp);
            
            max->updateMax(value);
            return;
        }
        
        // Here, the returned breakpoint is out of range
        // We will have to make a recursive call on the node's children
        
        if (!this->isLeaf()) {
            this->lowIndicesNode()->updateRecursiveMaxInRange(position,r,max);
            this->highIndicesNode()->updateRecursiveMaxInRange(position,r,max);
        }
    }
    
    inline T fastestMaxInRange(size_t position, Range r) const{
        return cascadingMaxInRange(position, r);
    }
};


/*
 * class RowNode
 *
 * This is the base class for the envelope trees on rows.
 * It provides in particular the efficient queries for maximum value of a column in a row range.
 *
 */
template <typename T>
class RowNode : public EnvTreeNode<T>{
private:

public:
    
    // This constructor creates a new RowNode with the specified children.
    // If they are not NULL, it will also compute the merged envelope.
    RowNode(Range r,RowNode<T> *lowIndices, RowNode<T> *highIndices ,Matrix<T> const* matrix): RowNode(r.min,r.max,lowIndices,highIndices,matrix)
    {
    }
    
    RowNode(size_t minRow, size_t maxRow,RowNode<T> *lowIndices, RowNode<T> *highIndices ,Matrix<T> const* matrix):EnvTreeNode<T>(minRow,maxRow)
    {
        if (minRow == maxRow) { // it is a leaf
            this->setEnvelope(new RowEnvelope<T>(matrix, minRow));
        }else{            
            this->setLowIndicesNode(lowIndices);
            this->setHighIndicesNode(highIndices);

            if (this->lowIndicesNode() && this->highIndicesNode()) { // if we can merge the envelopes of the children, do it immediately
                this->setEnvelope(mergeRowEnvelopes(this->lowIndicesNode()->envelope(), this->highIndicesNode()->envelope(), &(this->crossingBreakpointIndex()) ));
            }
        }
    }

    // Creates a row envelope binary tree node for the specified interval.
    // This also will creates its children and merge their envelopes.

    RowNode(size_t minRow, size_t maxRow, Matrix<T> const* matrix): EnvTreeNode<T>(minRow,maxRow)
    {
        if (minRow == maxRow) { // it is a leaf
            this->setEnvelope(new RowEnvelope<T>(matrix, minRow));
        }else{
            
            size_t midRow = minRow + ((maxRow - minRow)/2);
            
            RowNode<T> *lowIndices;
            RowNode<T> *highIndices;
            lowIndices = new RowNode<T>(minRow,midRow,matrix);
            highIndices= new RowNode<T>(midRow+1,maxRow,matrix);
            this->setLowIndicesNode(lowIndices);
            this->setHighIndicesNode(highIndices);

            this->setEnvelope(mergeRowEnvelopes(this->lowIndicesNode()->envelope(), this->highIndicesNode()->envelope(), &(this->crossingBreakpointIndex()) ));
        }
    }
    
    // Use this constructor to build the root of the row envelope binary tree for the specified matrix
    // COMPLEXITY: O(number_of_rows*( log(number_of_cols) + log(number_of_rows) ))
    // SIZE: O(number_of_rows* log(number_of_rows))
    RowNode(Matrix<T> const* matrix) : EnvTreeNode<T>(0,matrix->rows()-1)
    {        
        size_t minRow = this->minRow();
        size_t maxRow = this->maxRow();
        size_t midRow = minRow + ((maxRow - minRow)/2);
        
        RowNode<T> *lowIndices;
        RowNode<T> *highIndices;
        lowIndices = new RowNode<T>(minRow,midRow,matrix);
        highIndices= new RowNode<T>(midRow+1,maxRow,matrix);
        this->setLowIndicesNode(lowIndices);
        this->setHighIndicesNode(highIndices);

        this->setEnvelope(mergeRowEnvelopes(this->lowIndicesNode()->envelope(), this->highIndicesNode()->envelope()));
    }
    
    RowEnvelope<T>* envelope() const
    {
        return (RowEnvelope<T>*) EnvTreeNode<T>::envelope();
    }
    
    virtual RowNode<T>* lowIndicesNode() const { return (RowNode<T>*) EnvTreeNode<T>::lowIndicesNode(); }
    virtual RowNode<T>* highIndicesNode() const { return (RowNode<T>*) EnvTreeNode<T>::highIndicesNode(); }
    size_t minRow() const { return (this->range()).min; }
    size_t maxRow() const { return (this->range()).max; }
    
    // Returns the maximum value of the matrix in the specified column and in the specified row range
    T maxForColumnInRange(size_t col, size_t minRow, size_t maxRow) const{
        return this->maxInRange(col,Range(minRow,maxRow));
    }
};

/*
 * class ColNode
 *
 * This class provides the same data structure implementation as for RowNode
 * but for queries on rows and column ranges instead of queries on columns and row ranges.
 *
 *
 * For a better documentation of the methods, just refer to their equivalents in the RowNode class.
 */
template <typename T>
class ColNode : public EnvTreeNode<T>{

public:
    
    ColNode(size_t minCol, size_t maxCol, Matrix<T> const* matrix): EnvTreeNode<T>(minCol,maxCol)
    {
        if (minCol == maxCol) { // it is a leaf
            ColumnEnvelope<T>* envelope = new ColumnEnvelope<T>(matrix, minCol);
            this->setEnvelope(envelope);
        }else{
            size_t midCol = minCol + ((maxCol - minCol)/2);
            
            
            ColNode<T> *lowIndices;
            ColNode<T> *highIndices;
            lowIndices = new ColNode<T>(minCol,midCol,matrix);
            highIndices= new ColNode<T>(midCol+1,maxCol,matrix);
            this->setLowIndicesNode(lowIndices);
            this->setHighIndicesNode(highIndices);
            
            ColumnEnvelope<T>*  envelope = mergeColumnEnvelopes(this->lowIndicesNode()->envelope(), this->highIndicesNode()->envelope(), &(this->crossingBreakpointIndex()) );
            this->setEnvelope(envelope);
        }
    }
    
    
    ColNode(Matrix<T> const* matrix) : EnvTreeNode<T>(0,matrix->cols()-1)
    {
        size_t minCol = this->minCol();
        size_t maxCol = this->maxCol();
        size_t midCol = minCol + ((maxCol - minCol)/2);
        
        
        ColNode<T> *lowIndices;
        ColNode<T> *highIndices;
        lowIndices = new ColNode<T>(minCol,midCol,matrix);
        highIndices= new ColNode<T>(midCol+1,maxCol,matrix);
        this->setLowIndicesNode(lowIndices);
        this->setHighIndicesNode(highIndices);
        
        ColumnEnvelope<T>* envelope = mergeColumnEnvelopes(this->lowIndicesNode()->envelope(), this->highIndicesNode()->envelope(), &(this->crossingBreakpointIndex()) );
        
        this->setEnvelope(envelope);
    }
    
    size_t minCol() const  { return (this->range()).min; }
    size_t maxCol() const { return (this->range()).max; }

    ColumnEnvelope<T>* envelope() const
    {
        return (ColumnEnvelope<T>*)EnvTreeNode<T>::envelope();
    }
    ColNode<T>* lowIndicesNode() const { return (ColNode<T>* ) EnvTreeNode<T>::lowIndicesNode(); }
    ColNode<T>* highIndicesNode() const { return (ColNode<T>* ) EnvTreeNode<T>::highIndicesNode(); }

    
    T maxForRowInRange(size_t row, size_t minCol, size_t maxCol) const{
        return this->maxInRange(row,Range(minCol,maxCol));
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
    BasicRQNode< T > *_rangeMaxima;
    
public:
    ExtendedRowNode(size_t minRow, size_t maxRow, Matrix<T> const* matrix): RowNode<T>(minRow,maxRow,NULL,NULL,matrix)
    {
        if (minRow < maxRow){
            size_t midRow = minRow + ((maxRow - minRow)/2);
            this->setLowIndicesNode(new ExtendedRowNode<T>(minRow,midRow,matrix));
            this->setHighIndicesNode(new ExtendedRowNode<T>(midRow+1,maxRow,matrix));
            this->setEnvelope(mergeRowEnvelopes(this->lowIndicesNode()->envelope(), this->highIndicesNode()->envelope(),&(this->crossingBreakpointIndex()) ));
        }
    }
    
    ExtendedRowNode(Matrix<T> const* matrix): RowNode<T>(0,matrix->rows()-1,NULL,NULL,matrix)
    {
        size_t minRow = this-> minRow();
        size_t maxRow = this->maxRow();
        size_t midRow = minRow + ((maxRow - minRow)/2);
        
        this->setLowIndicesNode(new ExtendedRowNode<T>(minRow,midRow,matrix));
        this->setHighIndicesNode(new ExtendedRowNode<T>(midRow+1,maxRow,matrix));
        this->setEnvelope(mergeRowEnvelopes(this->lowIndicesNode()->envelope(), this->highIndicesNode()->envelope(),&(this->crossingBreakpointIndex()) ));
    }
    
    ~ExtendedRowNode<T>()
    {
        delete _rangeMaxima;
        delete _maxima;
    }
    
    virtual ExtendedRowNode<T>* lowIndicesNode() const { return (ExtendedRowNode<T>*) RowNode<T>::lowIndicesNode(); }
    virtual ExtendedRowNode<T>* highIndicesNode() const { return (ExtendedRowNode<T>*) RowNode<T>::highIndicesNode(); }
    
    inline const vector< T > *maxima() const
    {
        return _maxima;
    }
    
    inline const BasicRQNode<T> *rangeMaxima() const
    {
        return _rangeMaxima;
    }
    
    // This method computes the _maxima vector
    // COMPLEXITY: O(number_of_breakpoints * log(number_of_columns) ) = O(number_of_rows_in_the_envelope * log(number_of_columns) )
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
        maxCol = this->envelope()->maxPosition();
        row = (*breakpoints)[i].row;
        (*_maxima)[i] = flippedTree->maxForRowInRange(row,minCol,maxCol);
        
        _rangeMaxima = new BasicRQNode<T>(_maxima,0,_maxima->size()-1,&std::max<T>);

    }
    
    // Computes the interval maxima and tells the node's children to do the same
    //
    // COMPLEXITY: for the root node (i.e. computing the maximum in the entire tree), it should be
    // O(number_of_rows * log(number_of_rows) * log(number_of_columns))
    // At this point, I think it is rather O(number_of_rows^2 * log(number_of_columns)) (I know, this is bad!)
    void recursivelyComputeIntervalMaxima(const ColNode<T> *flippedTree)
    {
        computeIntervalMaxima(flippedTree);
        
        if (!this->isLeaf()) {
            this->lowIndicesNode()->recursivelyComputeIntervalMaxima(flippedTree);
            this->highIndicesNode()->recursivelyComputeIntervalMaxima(flippedTree);
        }
    }
    void recursivelyComputeIntervalMaxima_fast(const ColNode<T> *flippedTree)
    {
        if (this->isLeaf()) {
            // if we are at a leaf, we only have one breakpoint and it is easy to compute the maximum
            size_t minCol = 0, maxCol = this->envelope()->maxPosition();
            size_t row;
            
            row = (*this->envelope()->breakpoints())[0].row;
            T max = flippedTree->fastestMaxInRange(row,Range(minCol, maxCol));
            _rangeMaxima = new BasicRQNode<T>(max,&std::max<T>);
            _maxima = new vector<T>(1,max);
        }else{
            this->lowIndicesNode()->recursivelyComputeIntervalMaxima_fast(flippedTree);
            this->highIndicesNode()->recursivelyComputeIntervalMaxima_fast(flippedTree);
            
            if (this->crossingBreakpointIndex() <= 0) {
                
                // only the envelope of the highIndicesNode has been kept when merging
                // we duplicate the RMQ DS of this node
                _maxima = new vector<T>(*this->highIndicesNode()->maxima());
                _rangeMaxima = new BasicRQNode<T>(_maxima,0,_maxima->size()-1,&std::max<T>);

            }else if (this->crossingBreakpointIndex() == this->envelope()->numberOfBreakpoints()){
                // only the envelope of the lowIndicesNode has been kept when merging
                // we duplicate the RMQ DS of this node
                _maxima = new vector<T>(*this->lowIndicesNode()->maxima());
                _rangeMaxima = new BasicRQNode<T>(_maxima,0,_maxima->size()-1,&std::max<T>);

            }else{
                
                const vector< Breakpoint > *breakpoints = this->envelope()->breakpoints();
                size_t n = breakpoints->size();
                
                _maxima = new vector< T >(n); // create a new vector of the same size than the breakpoints one
                
                // for the first breakpoints, just copy the maxima table from the lowIndicesNode
                for (size_t i = 0; i < this->crossingBreakpointIndex() - 1; i++) {
                    (*_maxima)[i] = (*this->lowIndicesNode()->maxima())[i];
                }
                
                // for the intervals on both sides of the crossing breakpoint, we have to recompute the maximum using the flipped tree
                size_t minCol, maxCol;
                size_t row;

                // on the left side
                minCol = (*breakpoints)[this->crossingBreakpointIndex() - 1].col; // ... get the interval first index ...
                maxCol = (*breakpoints)[this->crossingBreakpointIndex()].col-1; // ... its last index ...
                row = (*breakpoints)[this->crossingBreakpointIndex() - 1].row; // ... and the corresponding row ...

                (*_maxima)[this->crossingBreakpointIndex() - 1] = flippedTree->fastestMaxInRange(row,Range(minCol,maxCol));

                // and on the right side 
                minCol = (*breakpoints)[this->crossingBreakpointIndex()].col; // ... get the interval first index ...
                // for the last index, be sure that we are not out of bounds
                if(this->crossingBreakpointIndex()+1 == breakpoints->size()){
                    maxCol = this->envelope()->maxPosition();
                }else{
                    maxCol = (*breakpoints)[this->crossingBreakpointIndex()+1].col-1;
                }
                row = (*breakpoints)[this->crossingBreakpointIndex()].row; // ... and the corresponding row ...
                
               (*_maxima)[this->crossingBreakpointIndex()] = flippedTree->fastestMaxInRange(row,Range(minCol,maxCol));
                
                // for the last part of the breakpoints, we again have to copy the maxima table for the highIndicesNode
                // to avoid computing the beginning index of the copy in the child max table, we do the copy backward
                size_t m = this->highIndicesNode()->envelope()->numberOfBreakpoints();
                
                for (size_t i = 1; n-i > this->crossingBreakpointIndex(); i++) {
                    (*_maxima)[n-i] = (*this->highIndicesNode()->maxima())[m-i];
                }

                // and we end by creating the RMQ data structure
                _rangeMaxima = new BasicRQNode<T>(_maxima,0,_maxima->size()-1,&std::max<T>);
            }
            this->lowIndicesNode()->deleteMaximaVector();
            this->highIndicesNode()->deleteMaximaVector();
        } //endif (this->isLeaf())
        
        Range range = this->range();
        if(range.min == 0 && range.max == this->envelope()->values()->rows()-1) // this is the root
        {
            this->deleteMaximaVector();
        }
    }
    
    // We override this to avoid compile-time errors (thank you C++) 
    vector<const ExtendedRowNode<T> *> canonicalNodes(size_t minRow, size_t maxRow) const
    {
        std::vector<const EnvTreeNode<T> *> buffer;
        this->getCanonicalNodes(buffer, minRow, maxRow);
        
        vector<const ExtendedRowNode<T> *> castedBuffer(buffer.size());
        
        for (size_t i = 0 ; i < buffer.size(); i++) {
            castedBuffer[i] = (const ExtendedRowNode<T> *)buffer[i];
        }
        return castedBuffer;
    }
    
protected:
    void deleteMaximaVector()
    {
        delete _maxima;
        _maxima = NULL;
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
    ExtendedRowNode<T> *_rowsTree; // denoted T_h in the SubmatrixQueries article
    ColNode<T> *_columnTree; // denoted \mathcal{B} in the SubmatrixQueries article
    
public:
    // Constructs a new submatrix query datastructure for the given inverse Monge matrix.
    // COMPLEXITY (expected): if m = number_of_rows and n = number_of_columns, O(m log(m) log(n) + n(log m + log n)) )
    // The current complexity if rather O(m^2 * log(n) + n(log m + log n))
    // SIZE: O(m log m)
    SubmatrixQueriesDataStructure(Matrix<T> const* matrix)
    {
        _rowsTree = new ExtendedRowNode<T>(matrix);
        _columnTree = new ColNode<T>(matrix);
        
        _rowsTree->recursivelyComputeIntervalMaxima_fast(_columnTree);
    }
    
    ~SubmatrixQueriesDataStructure()
    {
        delete _rowsTree;
        delete _columnTree;
    }
    
    inline const ExtendedRowNode<T>* rowsTree() const
    {
        return _rowsTree;
    }
    
    inline const ColNode<T>* columnTree() const
    {
        return _columnTree;
    }
    
    
    void updateMaxForRowNodeOverColumnRange(ExtendedRowNode<T> const *rowNode, Range colRange, MaxValue<T> *max) const
    {
        RowEnvelope<T> *envelope = rowNode->envelope();
        const vector< Breakpoint > *breakpoints = envelope->breakpoints();
        
        Breakpoint startBP, endBP;
        size_t startBPIndex, endBPIndex;
        
        // we get the last breakpoint before colRanges.min and its index ...
        startBP = envelope->breakpointBeforePosition(colRange.min, &startBPIndex);
        // ... and the last breakpoint before colRanges.max and its index
        endBP = envelope->breakpointBeforePosition(colRange.max, &endBPIndex);
        
        // if the range does not contain any interval, the returned breakpoints will be the sames
        // as a consequence, we have to check for that case to avoid undefined behavior
        if (endBPIndex == startBPIndex) {
            size_t row = startBP.row;
            max->updateMax( _columnTree->fastestMaxInRange(row,colRange));
            return;
        }
        
        // at this point, we have the prefix: it is the interval [colRanges.min,(*breakpoints)[startBPIndex+1].col-1]
        // we first have to check if it is empty or not ...
        if ((*breakpoints)[startBPIndex].col < colRange.min) {
            // it is not empty, go on ...
            size_t row = (*breakpoints)[startBPIndex].row;
            max->updateMax(_columnTree->fastestMaxInRange(row,Range(colRange.min, (*breakpoints)[startBPIndex+1].col-1)));
        }else{
            // it is empty, we just move  the start index so it is consistent with the call on the RMQ data structure
            startBPIndex--;
        }
        
        // now, we check for the fully contained intervals
        // in the case we still have some breakpoints to explore ...
        
        if (endBPIndex - startBPIndex > 1) { // to have at least one interval, we need at least two breakpoints ...
            // the range of the set of intervals is then [(*breakpoints)[startBPIndex+1].col,(*breakpoints)[endBPIndex].col-1]
            const BasicRQNode<T> *rangeMaxima = rowNode->rangeMaxima();
            max->updateMax(rangeMaxima->query(startBPIndex+1,endBPIndex-1));
        }
        
        // check for the rest of the range
        // the remaining of the interval is between the last breakpoint column and the column range maximum :
        // the suffix is of range [(*breakpoints)[i-1].col,colRanges.max]
        // ( even if there are two ways to exit the loop, in the end we have to do the same processing for the remaining )
        
        size_t row = endBP.row;
        max->updateMax(_columnTree->fastestMaxInRange(row,Range(endBP.col,colRange.max)));
    }
    
    T maxInRange(size_t minRow, size_t maxRow, size_t minCol, size_t maxCol) const
    {
        return maxInRange(Range(minRow,maxRow), Range(minCol,maxCol));
    }
    
    // This the query method: it returns the maximum of the submatrix of the Monge inverse matrix within the specified row and column ranges.
    // COMPLEXITY (expected): O(log(number_of_rows) * (log(number_of_rows) + log(number_of_cols)) )
    // In this slow implementation, we have instead O(log(number_of_rows) * (number_of_rows + log(number_of_cols)) ) (see the WARNING)
    T maxInRangeSlow(Range rowRange, Range colRanges) const
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
                
                // !!WARNING!!
                // HERE, INSTEAD OF NAIVELY COMPUTING THE MAXIMUM BY SUCCESSIVELY QUERYING THE INTERVALS,
                // WE SHOULD USE A RMQ DATA STRUCTURE TO REDUCE THE RUNNING TIME TO A LOG FACTOR (WE CAN EVEN HAVE O(1) TIME QUERY !)
                //
                // THIS IS DONE IN THE "QUICK" VERSION OF THE QUERY (SEE THE METHOD maxInRange )
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
    
    // This the query method: it returns the maximum of the submatrix of the Monge inverse matrix within the specified row and column ranges.
    // COMPLEXITY: O(log(number_of_rows) * (log(number_of_rows) + log(number_of_cols)) )
    T maxInRange(Range rowRange, Range colRanges) const
    {
        vector<const ExtendedRowNode<T> *> rowNodes = _rowsTree->canonicalNodes(rowRange.min,rowRange.max);
        
        MaxValue<T> max;
        for (typename vector<const ExtendedRowNode<T> *>::iterator nodesIterator = rowNodes.begin(); nodesIterator != rowNodes.end(); ++nodesIterator) {
            updateMaxForRowNodeOverColumnRange((*nodesIterator), colRanges, &max);
        }
        
        return max.value();
    }
    
    T maxInSubmatrix(Range rowRange, Range colRange) const
    {
        MaxValue<T> max;
        
        updateRecursivelyMaxInRange(rowRange,colRange,rowsTree(),&max);
        
        return max.value();
    }
    
    void updateRecursivelyMaxInRange(Range rowRange, Range colRange, ExtendedRowNode<T> const *rowNode, MaxValue<T> *max) const
    {
        if (!rowRange.intersects(rowNode->range())) {
            return;
        }
        
        if (!rowRange.contains(rowNode->range())) {
            // make a recursive call
            updateRecursivelyMaxInRange(rowRange, colRange, rowNode->lowIndicesNode(), max);
            updateRecursivelyMaxInRange(rowRange, colRange, rowNode->highIndicesNode(), max);
            return;
        }
        
        updateMaxForRowNodeOverColumnRange(rowNode, colRange, max);
    }
    
    inline T fastestMaxInSubmatrix(Range rowRange, Range colRange) const
    {
        return maxInSubmatrix(rowRange, colRange);
    }
};
#endif /* defined(__SubmatrixQueries__enveloppe_tree__) */
