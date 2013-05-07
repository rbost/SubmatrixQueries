//
//  max_value.h
//  SubmatrixQueries
//
//  Created by Raphael Bost on 16/01/13.
//  Copyright (c) 2013 Raphael Bost. All rights reserved.
//

#ifndef __SubmatrixQueries__max_value__
#define __SubmatrixQueries__max_value__

#include <algorithm>
#include <cassert>

/* class MaxValue
 *
 * This class is a wrapper class for maximum values.
 * The problem is that for a template implementation, you just cannot predefine the minimum value for a specific type:
 * for example FLT_MAX is different of DBL_MAX or INT_MAX ...
 * So, to easily implement maxima, we just create this wrapper that encaspulates a value and a flag that is set since a value has been put in the wrapper.
 */
 
template <typename T>
class MaxValue {
    T _value;
    bool _hasBeenSet;
    
public:
    MaxValue() : _hasBeenSet(false) {};
    MaxValue(T v) : _value(v), _hasBeenSet(true) {};
    
    // Replaces the wrapper's value by the input value if it is bigger 
    T updateMax(T v)
    {
        if (_hasBeenSet) {
            _value = std::max(_value,v);
        }else{
            _value = v;
            _hasBeenSet = true;
        }
        
        return _value;
    }
    
    T updateMax(MaxValue<T> m) {
      updateMax(m.value());
    }
    
    inline bool hasBeenSet() const { return _hasBeenSet;}
    inline T  value() const { assert(_hasBeenSet); return _value; }
    
    bool operator==(const MaxValue<T> &other) const {
      return value() == other.value();
    }
    bool operator!=(const MaxValue<T> &other) const {
      return !(*this == other);
    }
    bool operator<(const MaxValue<T> &other) const {
      return value() < other.value();
    }
    bool operator<=(const MaxValue<T> &other) const {
      return value() <= other.value();
    }
    bool operator>(const MaxValue<T> &other) const {
      return value() > other.value();
    }
    bool operator>=(const MaxValue<T> &other) const {
      return value() >= other.value();
    }
};

template <typename T>
class MaxInMatrix : public MaxValue<T> {
    size_t _row;
    size_t _col;
    
public:
    MaxInMatrix() : MaxValue<T>() {};
    MaxInMatrix(T v, size_t r, size_t c) : MaxValue<T>(v), _row(r), _col(c){};
    
    T updateMax(T v, size_t r, size_t c)
    {
        if (!this->hasBeenSet()) {
            _row = r;
            _col = c;
        }else if(v > MaxValue<T>::value()){
            _row = r;
            _col = c;
        }

        return  MaxValue<T>::updateMax(v);
    }
    
    T updateMax(MaxInMatrix<T> m) {
      size_t r, c;
      T v = m.value(&r, &c);
      
      if (!this->hasBeenSet()) {
          _row = r;
          _col = c;
      } else if(v > MaxValue<T>::value()){
          _row = r;
          _col = c;
      }
      
      return MaxValue<T>::updateMax(v);
    }
    
    inline void getMaxPosition(size_t *r, size_t *c)
    {
        assert(this->hasBeenSet());
        if (r != NULL) {
            *r = _row;
        }
        if (c != NULL) {
            *c = _col;
        }
    }
    
    inline T value(size_t *r, size_t *c)
    {
        getMaxPosition(r,c);
        return MaxValue<T>::value();
    }
};
#endif /* defined(__SubmatrixQueries__max_value__) */
