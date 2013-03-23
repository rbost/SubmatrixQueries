//
//  envelope.h
//  SubmatrixQueries
//
//  Created by Raphael Bost on 03/01/13.
//  Copyright (c) 2013 Raphael Bost. All rights reserved.
//

#ifndef __SubmatrixQueries__envelope__
#define __SubmatrixQueries__envelope__

#include <vector>
#include "matrix.h"
#include "debug_assert.h"

using namespace std;
using namespace matrix;

namespace envelope {
    
    /*
     *  class Breakpoint
     *  This class represents the Breakpoint data structure used to implicitly represent envelopes.
     *  It just stores the position of a breakpoint in the matrix
     *
     */
    struct Breakpoint{
        size_t row;
        size_t col;

        inline Breakpoint() : row(0), col(0) {}
        inline Breakpoint(size_t r, size_t c) : row(r), col(c) {}
        inline bool operator==(Breakpoint bp) const{ return (row == bp.row) && (col == bp.col);}
        inline bool operator!=(Breakpoint bp) const{ return (row != bp.row) || (col != bp.col);}
    };
    
    /*
     *  class Envelope
     *  This represents the envelope data structure of the SubmatrixQueries's article.
     *  We implicitly represent it by using a breakpoint vector.
     *
     *  The matrix used to build the envelope is referenced by the _values reference.
     */
    
    template <typename T>
    class Envelope {
    protected:
        Matrix<T> const& _values;
        const vector< Breakpoint > *_breakpoints;
        
    private:
        Envelope(Matrix<T> const& values) : _values(values)
        {
            
        }
    public:
        // Initializes a new envelope with the specified breakpoints.
        Envelope(Matrix<T> const& values, const vector< Breakpoint > *breakPoints): _values(values), _breakpoints(breakPoints)
        {
        }
        
        virtual ~Envelope()
        {
            delete _breakpoints;
        }
        
        inline const vector< Breakpoint >* breakpoints() const
        {
            return this->_breakpoints;
        }
        
        inline size_t numberOfBreakpoints() const {
            return this->breakpoints()->size();
        }
        
        inline Matrix<T> const& values() const
        {
            return this->_values;
        }

        
        // Finds the last breakpoint before col in the breakpoint list using binary search
        // In the first method, we specify the index in between we have to search
        // COMPLEXITY: O( log(number_of_breakpoints) )
        
        Breakpoint breakpointBeforePosition(size_t pos, size_t iMin, size_t iMax, size_t *foundPosition) const
        {
            DEBUG_ASSERT(iMin >= 0);
            DEBUG_ASSERT(iMax < this->numberOfBreakpoints());
            DEBUG_ASSERT(iMin <= iMax);
            
            size_t iMid;
            size_t i;
            i = iMin;
            
            while (iMax > iMin) {
                iMid = iMin + ((iMax - iMin)/2); // to avoid overflows
                
                size_t selectedPosition = this->positionForBreakpointAtIndex(iMid);// (*this->breakpoints())[iMid].col;
                
                size_t nextPosition = this->positionForBreakpointAtIndex(iMid+1);
                
                if (selectedPosition <= pos && pos < nextPosition) {
                    i = iMid;
                    break; // we are in the interval between two breakpoints!
                }else if (selectedPosition > pos){
                    iMax = iMid - 1;
                }else{
                    iMin = iMid + 1;
                }
                i = iMin;
            }

            if (foundPosition) {
                *foundPosition = i;
            }
            return (*this->breakpoints())[i];
        }
        
        inline Breakpoint breakpointBeforePosition(size_t pos, size_t *foundPosition) const
        {
            return breakpointBeforePosition(pos,0,this->numberOfBreakpoints() - 1,foundPosition);
        }
        
        inline Breakpoint breakpointBeforePosition(size_t pos) const
        {
            return breakpointBeforePosition(pos, NULL);
        }
        
        virtual size_t mappedPosition(size_t index) const = 0;
        
        // ABSTRACT METHOD
        // Implementation must return a breakpoint that will be on the envelope at the given position
        virtual Breakpoint newBreakPointAtPosition(size_t index) const = 0;
        
        // Returns the value of the envelope for the given position
        virtual T valueForPosition(size_t position) const
        {
            Breakpoint bp = this->breakpointBeforePosition(position);
            return this->valueForPositionAfterBreakpoint(position,bp);
        }
        
        // ABSTRACT METHOD
        // Implementations must return the maximum position for a breakpoint.
        virtual inline size_t maxPosition() const = 0;
        
        // ABSTRACT METHOD
        // Implementations must return the position of the breakpoint given as an argument.
        virtual inline size_t positionForBreakpoint(Breakpoint bp) const = 0;
        virtual inline size_t mappedPositionForBreakpoint(Breakpoint bp) const = 0;

        // ABSTRACT METHOD
        // Implementations must return the value at position given that bp is the last breakpoint before position 
        virtual T valueForPositionAfterBreakpoint(size_t position, Breakpoint bp) const = 0;
        
        virtual inline size_t positionForBreakpointAtIndex(size_t i) const
        {
            Breakpoint bp = (*this->breakpoints())[i];
            return this->positionForBreakpoint(bp);
        }

        virtual T firstValue() const = 0;
        virtual T lastValue() const = 0;
        

        Breakpoint firstBreakpoint() const
        {
            return this->breakpoints()->front();
        }
        
        Breakpoint lastBreakpoint() const
        {
            return this->breakpoints()->back();
        }
        

        inline T valueForBreakpoint(Breakpoint bp) const
        {
            return this->_values(bp.row,bp.col);
        }
        
    };
    
    template <typename T>
    class RowEnvelope : public Envelope<T>
    {
    public:
        // Initializes a new envelope corresponding to an envelope of one row (considered as a pseudo-line).
        RowEnvelope(Matrix<T> const& values, size_t row) : Envelope<T>(values,new vector< Breakpoint >(1,Breakpoint(row,0)))
        {
        }
        
        // Initializes a new envelope withe the specified breakpoints.
        RowEnvelope(Matrix<T> const& values, const vector< Breakpoint > *breakPoints): Envelope<T>(values,breakPoints)
        {
        }

        // Finds the row corresponding to the upper-line of the envelope for the column col using binary search
        // COMPLEXITY: O( log(number_of_breakpoints) )
        size_t rowForColumn(size_t col) const
        {
            
            return this->breakpointBeforePosition(col).row;
        }
        
        // COMPLEXITY: O( log(number_of_breakpoints) )
        T valueForColumn(size_t col) const
        {
            size_t row = this->rowForColumn(col);
            
            return  this->_values(row,col);
        }
        
        size_t numberOfColumns() const
        {
            return this->_values.cols();
        }
        
        size_t mappedPosition(size_t index) const
        {
            return rowForColumn(index);
        }
        
        size_t maxPosition() const
        {
            return numberOfColumns()-1;
        }
        
        Breakpoint newBreakPointAtPosition(size_t index) const
        {
            return Breakpoint(mappedPosition(index),index);
        }
        
        inline size_t positionForBreakpoint(Breakpoint bp) const
        {
            return bp.col;
        }
        
        inline size_t mappedPositionForBreakpoint(Breakpoint bp) const
        {
            return bp.row;
        }
        
        T valueForPositionAfterBreakpoint(size_t position, Breakpoint bp) const
        {
            return (this->values())(bp.row,position);
        }
        
        T firstValue() const
        {
            return (this->values())((this->breakpoints())->front().row,0);
        }
        T lastValue() const
        {
            size_t lastCol = this->values().cols()-1;
            return (this->values())((this->breakpoints())->back().row,lastCol);
        }
    };

    template <typename T>
    class ColumnEnvelope : public Envelope<T>
    {
    public:
        ColumnEnvelope(Matrix<T> const& values, size_t col) : Envelope<T>(values,new vector< Breakpoint >(1,Breakpoint(0,col)))
        {
        }
        
        // Initializes a new envelope withe the specified breakpoints.
        ColumnEnvelope(Matrix<T> const& values, const vector< Breakpoint > *breakPoints): Envelope<T>(values,breakPoints)
        {
        }
        
        // Finds the row corresponding to the upper-line of the envelope for the column col using binary search
        // COMPLEXITY: O( log(number_of_breakpoints) )
        size_t columnForRow(size_t row) const
        {
            return this->breakpointBeforePosition(row).col;
        }
        
        // COMPLEXITY: O( log(number_of_breakpoints) )
        T valueForRow(size_t row) const
        {
            size_t col = this->columnForRow(row);
            
            return  this->_values(row,col);
        }
        
        size_t numberOfRows() const
        {
            return this->_values.rows();
        }
        
        size_t mappedPosition(size_t index) const
        {
            return columnForRow(index);
        }
        
        size_t maxPosition() const
        {
            return numberOfRows()-1;
        }
        
        Breakpoint newBreakPointAtPosition(size_t index) const
        {
            return Breakpoint(index,mappedPosition(index));
        }

        inline size_t positionForBreakpoint(Breakpoint bp) const
        {
            return bp.row;
        }
        
        inline size_t mappedPositionForBreakpoint(Breakpoint bp) const
        {
            return bp.col;
        }

        T valueForPositionAfterBreakpoint(size_t position, Breakpoint bp) const
        {
            return (this->values())(position,bp.col);
        }

        T firstValue() const
        {
            return (this->values())(0,(this->breakpoints())->front().col);
        }
        T lastValue() const
        {
            size_t lastRow = this->values().rows()-1;
            return (this->values())(lastRow,(this->breakpoints())->back().col);
        }

    };

#pragma mark Envelope Merging
    
    // we suppose that the envelope e1 is the envelope for lower indices than e2
    // this function assumes that there is a breakpoint.
    // COMPLEXITY: O( log(number_of_breakpoints) * log(number_of_columns) )
    template <typename T> Breakpoint mergeBreakPoint(Envelope<T> * e1, Envelope<T> * e2, size_t *indexInE1, size_t *indexInE2) {
        size_t minPos, maxPos, pos;
        size_t minIndex1, maxIndex1, index1;
        size_t minIndex2, maxIndex2, index2;
        Breakpoint bp1, bp2;
        
        minPos = 0; maxPos = e1->maxPosition();
        minIndex1 = 0; maxIndex1 = e1->numberOfBreakpoints()-1;
        minIndex2 = 0; maxIndex2 = e2->numberOfBreakpoints()-1;
        
        while (maxPos > minPos) {
            pos = minPos + ((maxPos - minPos)/2); // to avoid overflows

            bp1 = e1->breakpointBeforePosition(pos,minIndex1,maxIndex1,&index1);
            bp2 = e2->breakpointBeforePosition(pos,minIndex2,maxIndex2,&index2);

            T val1 = e1->valueForPositionAfterBreakpoint(pos,bp1);
            T val2 = e2->valueForPositionAfterBreakpoint(pos,bp2);
            
            if (val1 > val2) {
                // the envelope e1 is this over e2 for this column, so search on the right part of the the envelopes
                minPos = pos + 1;
                
                minIndex1 = index1;
                minIndex2 = index2;
            }else{
                maxPos = pos;
                maxIndex1 = index1;
                maxIndex2 = index2;
            }
        }
        if (indexInE1) {
            *indexInE1 = index1+1;
        }
        if (indexInE2) {
            *indexInE2 = index2+1;
        }
        return e2->newBreakPointAtPosition(maxPos);
    }

    // we suppose that the envelope e1 is the envelope for lower indices than e2
    // if a valid pointer crossingBpIndex is provided. In the new enveloppe, every breakpoint whose index is strictly less than crossingBpIndex
    // comes from e1, and if striclty more than crossingBpIndex comes from e2. If the returned envelope is not a "merge" of e1 and e2,
    // crossingBpIndex is set to 0 if a copy of e2 is returned, and is set to number_of_breakpoints of e1 if a copy of e1 is returned. 
    
    // COMPLEXITY: O( log(number_of_breakpoints) * log(number_of_columns) + number_of_breakpoints)
    
    // WARNING: THIS RETURNS A NEWLY ALLOCATED ENVELOPE. IT MEANS YOU ARE RESPONSIBLE FOR DESTROYING IT WHEN YOU ARE DONE WITH IT!
    template <typename T> vector<Breakpoint>* mergeEnvelopes(Envelope<T> *e1, Envelope<T> *e2, size_t *crossingBpIndex)
    {
        T firstValue1 = e1->firstValue();
        T firstValue2 = e2->firstValue();
        
        T lastValue1 = e1->lastValue(); 
        T lastValue2 = e2->lastValue();
        
        if ((firstValue1 - firstValue2)*(lastValue1 - lastValue2) >= 0) {
            // The two expressions have the same sign, we do not have to compute a breakpoint.
            // The pseudo lines corresponding to the envelopes are not crossing, so we have to return the "upper" one.
            if (firstValue1 > firstValue2 || lastValue1 >= lastValue2) {
                // e1 is "over" e2
                if (crossingBpIndex) {
                    *crossingBpIndex = e1->breakpoints()->size();
                }
                return new vector<Breakpoint>(*e1->breakpoints());
            }else{
                // e2 is "over" e1
                if (crossingBpIndex) {
                    *crossingBpIndex = 0;
                }
                return new vector<Breakpoint>(*e2->breakpoints());
            }
        }

        // first, we find the new breakpoint ...
        size_t indexInE1, indexInE2;
        Breakpoint newBreakpoint = mergeBreakPoint(e1, e2, &indexInE1, &indexInE2);
        size_t newBpPosition = e1->positionForBreakpoint(newBreakpoint);

        // ... and then we concatenate the two envelopes and the new breakpoint:
        
        // first of all, compute the size of the envelope and create the breakpoint vector
        size_t eLength;
        eLength = indexInE1 + (e2->numberOfBreakpoints() - indexInE2 + 1);
        size_t e1MaxIndex = indexInE1;

        if (newBpPosition == e2->positionForBreakpointAtIndex(indexInE2)) {
            eLength--;
        }
        if (newBpPosition == e1->positionForBreakpointAtIndex(indexInE1-1)) {
            eLength--;
            e1MaxIndex--;
        }
        
        vector<Breakpoint> *newBpList = new vector<Breakpoint>(eLength);
        
        // then, insert the breakpoints of the first envelope
        size_t i = 0;
        
        while (i < e1->breakpoints()->size() && i < e1MaxIndex) {
            (*newBpList)[i] = (*(e1->breakpoints()))[i];
            i++;
        }
        
        // add the new breakpoint
        (*newBpList)[i] = newBreakpoint;

        // save the position of the crossing breakpoint in the new enveloppe
        if (crossingBpIndex) {
            *crossingBpIndex = i;
        }
        i++;
        
        // and copy the breakpoint of the second envelope
        // To improve performances, here, we just just jump to the first breakpoint to insert using a
        // logarithmic complexity method instead of just going through all the useless breakpoints.
        size_t beginningIndex;
        e2->breakpointBeforePosition(newBpPosition,&beginningIndex);
        
        for(size_t j = beginningIndex+1; j < e2->breakpoints()->size(); j++, i++){
            (*newBpList)[i] = (*(e2->breakpoints()))[j];
        }
        
        return newBpList;
    }
    
    template <typename T> RowEnvelope<T> * mergeRowEnvelopes(RowEnvelope<T> * e1, RowEnvelope<T> * e2, size_t *crossingBpIndex) {
        vector<Breakpoint>* newBreakpoints = mergeEnvelopes(e1, e2, crossingBpIndex);
        
        return new RowEnvelope<T>(e1->values(),newBreakpoints);
        
    }
    
    template <typename T> RowEnvelope<T> * mergeRowEnvelopes(RowEnvelope<T> * e1, RowEnvelope<T> * e2)
    {
        return mergeRowEnvelopes(e1,e2,NULL);
    }

    template <typename T> ColumnEnvelope<T> * mergeColumnEnvelopes(ColumnEnvelope<T> * e1, ColumnEnvelope<T> * e2, size_t *crossingBpIndex) {
        vector<Breakpoint>* newBreakpoints = mergeEnvelopes(e1, e2, crossingBpIndex);
        
        return new ColumnEnvelope<T>(e1->values(),newBreakpoints);    
    }

    template <typename T> ColumnEnvelope<T> * mergeColumnEnvelopes(ColumnEnvelope<T> * e1, ColumnEnvelope<T> * e2)
    {
        return mergeColumnEnvelopes(e1,e2,NULL);
    }

}

#endif /* defined(__SubmatrixQueries__envelope__) */
