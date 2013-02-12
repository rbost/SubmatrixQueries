//
//  rmq.h
//  SubmatrixQueries
//
//  Created by Raphael Bost on 30/01/13.
//  Copyright (c) 2013 Raphael Bost. All rights reserved.
//

#ifndef __SubmatrixQueries__rmq__
#define __SubmatrixQueries__rmq__

#include <iostream>
#include <vector>
#include <valarray>
#include <cassert>
#include <algorithm>


template <typename T>
class BasicRQNode {
    T _value;
    
    size_t _minIndex, _maxIndex;
    bool _isLeaf;
    const BasicRQNode<T> *_highIndicesNode;
    const BasicRQNode<T> *_lowIndicesNode;
    
    typedef const T& (*compareFunctionPtr)(T const&, T const&);
    compareFunctionPtr _compare;
    
public:
    
    T value() const { return _value; }
    size_t minIndex() const { return _minIndex; }
    size_t maxIndex() const { return _maxIndex; }
    bool isLeaf() const { return _isLeaf; }
    const BasicRQNode<T> *highIndicesNode() const { return _highIndicesNode; }
    const BasicRQNode<T> *lowIndicesNode() const { return _lowIndicesNode; }
    compareFunctionPtr compareFunction() const { return _compare; }
    
    BasicRQNode(const std::vector<T> *values, size_t minIndex, size_t maxIndex, const T& (*compareFunc)(T const&, T const&) ) :
        _minIndex(minIndex),_maxIndex(maxIndex),_compare(compareFunc)
    {
        assert(minIndex <= maxIndex);
        
        if (minIndex == maxIndex) {
            _isLeaf = true;
            _value = (*values)[minIndex];
        }else{
            _isLeaf = false;
            
            size_t midIndex = minIndex + (maxIndex - minIndex)/2;
            _lowIndicesNode  = new BasicRQNode<T>(values, minIndex, midIndex, compareFunc);
            _highIndicesNode = new BasicRQNode<T>(values, midIndex+1, maxIndex, compareFunc);
            
            _value = (*_compare)(_lowIndicesNode->value(),_highIndicesNode->value());
        }
    }
    
    BasicRQNode(const std::valarray<T> *values, size_t minIndex, size_t maxIndex, const T& (*compareFunc)(T const&, T const&) ) :
    _minIndex(minIndex),_maxIndex(maxIndex),_compare(compareFunc)
    {
        assert(minIndex <= maxIndex);
        
        if (minIndex == maxIndex) {
            _isLeaf = true;
            _value = (*values)[minIndex];
        }else{
            _isLeaf = false;
            
            size_t midIndex = minIndex + (maxIndex - minIndex)/2;
            _lowIndicesNode  = new BasicRQNode<T>(values, minIndex, midIndex, compareFunc);
            _highIndicesNode = new BasicRQNode<T>(values, midIndex+1, maxIndex, compareFunc);
            
            _value = (*_compare)(_lowIndicesNode->value(),_highIndicesNode->value());
        }
    }
    
    BasicRQNode(T val,const T& (*compareFunc)(T const&, T const&)) : _minIndex(0), _maxIndex(0), _isLeaf(true), _value(val),_compare(compareFunc)
    {}
    
    BasicRQNode(BasicRQNode<T> *nodePtr) : _minIndex(nodePtr->minIndex()), _maxIndex(nodePtr->maxIndex()), _isLeaf(nodePtr->isLeaf()), _value(nodePtr->value()), _compare(nodePtr->compareFunction())
    {
        if (!this->isLeaf()) {
            _lowIndicesNode = new BasicRQNode<T>(nodePtr->lowIndicesNode());
            _highIndicesNode = new BasicRQNode<T>(nodePtr->highIndicesNode());
        }
    }
    ~BasicRQNode()
    {
        if (_isLeaf) {
            delete _highIndicesNode;
            delete _lowIndicesNode;
        }
    }
    
    T query(size_t startIndex, size_t endIndex) const
    {
        assert(startIndex <= endIndex);

        size_t maxIndex = this->maxIndex();
        size_t minIndex = this->minIndex();
        
        assert(startIndex >= minIndex);
        assert(endIndex <= maxIndex);

        if (minIndex >= startIndex && maxIndex <= endIndex) {
            return this->value();
        }

        size_t midIndex = minIndex + (maxIndex - minIndex)/2;
        
        return (*_compare)(this->lowIndicesNode()->query(startIndex,midIndex),this->lowIndicesNode()->query(midIndex+1,endIndex));
    }
};

#endif /* defined(__SubmatrixQueries__rmq__) */
