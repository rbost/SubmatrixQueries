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
        
        const vector< Breakpoint >* breakpoints() const
        {
            return this->_breakpoints;
        }
        
        size_t numberOfBreakpoints() const {
            return this->breakpoints()->size();
        }
        
        Matrix<T> const& values() const
        {
            return this->_values;
        }

        
        // Finds the last breakpoint before col in the breakpoint list using binary search
        // COMPLEXITY: O( log(number_of_breakpoints) )

        Breakpoint breakpointBeforePosition(size_t pos, size_t *foundPosition) const
        {
            size_t iMin, iMax, iMid;
            size_t i;
            iMin = 0; iMax = this->numberOfBreakpoints() - 1;
            i = 0;
            
            while (iMax > iMin) {
                iMid = iMin + ((iMax - iMin)/2); // to avoid overflows
                
                size_t selectedPosition = this->positionForBreakpointAtIndex(iMid);// (*this->breakpoints())[iMid].col;
                
                if (iMid == numberOfBreakpoints()-1) {
                    // in case we selected the last breakpoint
                    if (selectedPosition > pos){ // we have to select a lower breakpoint
                        iMax = iMid - 1;
                    }else{
                        i = iMid;
                        break;
                    }
                }else{
                    size_t nextPosition = this->positionForBreakpointAtIndex(iMid+1);
                    
                    if (selectedPosition <= pos && pos < nextPosition) {
                        i = iMid;
                        break; // we are in the interval between two breakpoints!
                    }else if (selectedPosition > pos){
                        iMax = iMid - 1;
                    }else{
                        iMin = iMid + 1;
                    }
                }
                i = iMin;
            }
            
            if (foundPosition) {
                *foundPosition = i;
            }
            return (*this->breakpoints())[i];
        }
        
        Breakpoint breakpointBeforePosition(size_t pos) const
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
        virtual size_t maxPosition() const = 0;
        
        // ABSTRACT METHOD
        // Implementations must return the position of the breakpoint given as an argument.
        virtual inline size_t positionForBreakpoint(Breakpoint bp) const = 0;
        virtual inline size_t mappedPositionForBreakpoint(Breakpoint bp) const = 0;

        // ABSTRACT METHOD
        // Implementations must return the value at position given that bp is the last breakpoint before position 
        virtual T valueForPositionAfterBreakpoint(size_t position, Breakpoint bp) const = 0;
        
        virtual size_t positionForBreakpointAtIndex(size_t i) const
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
        

        T valueForBreakpoint(Breakpoint bp) const
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
    template <typename T> Breakpoint mergeBreakPoint(Envelope<T> * e1, Envelope<T> * e2) {
        size_t iMin, iMax, index;
        iMin = 0; iMax = e1->maxPosition();
        
        while (iMax > iMin) {
            index = iMin + ((iMax - iMin)/2); // to avoid overflows

            T val1 = e1->valueForPosition(index);
            T val2 = e2->valueForPosition(index);
            

            if (val1 > val2) {
                // the envelope e1 is this over e2 for this column, so search on the right part of the the envelopes
                iMin = index + 1;
            }else{
                iMax = index;
            }
        }
        
        return e2->newBreakPointAtPosition(iMax);
    }

    // we suppose that the envelope e1 is the envelope for lower indices than e2
    // if a valid pointer crossingBpIndex is provided. In the new enveloppe, every breakpoint whose index is strictly less than crossingBpIndex
    // comes from e1, and if striclty more than crossingBpIndex comes from e2. If the returned envelope is not a "merge" of e1 and e2,
    // crossingBpIndex is set to -1 if a copy of e2 is returned, and is set to number_of_breakpoints of e1 if a copy of e1 is returned. 
    
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
                    *crossingBpIndex =-1;
                }
                return new vector<Breakpoint>(*e2->breakpoints());
            }
        }

        // first, we find the new breakpoint ...
        Breakpoint newBreakpoint = mergeBreakPoint(e1, e2);
        
        // ... and then we concatenate the two envelopes and the new breakpoint:
        
        // first, insert the breakpoints of the first envelope
        vector<Breakpoint> *newBpList = new vector<Breakpoint>();
        size_t i = 0;
        size_t newBpPosition = e1->positionForBreakpoint(newBreakpoint);
        
        while (i < e1->breakpoints()->size() && e1->positionForBreakpointAtIndex(i) < newBpPosition) {
            newBpList->push_back((*(e1->breakpoints()))[i]);
            i++;
        }

        
        // add the new breakpoint
        newBpList->push_back(newBreakpoint);
        
        // save the position of the crossing breakpoint in the new enveloppe
        if (crossingBpIndex) {
            *crossingBpIndex = i;
        }
        
        // and copy the breakpoint of the second envelope
        // To improve performances, here, we just just jump to the first breakpoint to insert using a
        // logarithmic complexity method instead of just going through all the useless breakpoints.
        size_t beginningIndex;
        e2->breakpointBeforePosition(newBpPosition,&beginningIndex);
        i = beginningIndex+1;
        
        while (i < e2->breakpoints()->size()){
            newBpList->push_back((*(e2->breakpoints()))[i]);
            i++;
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
