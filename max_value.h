//
//  max_value.h
//  KMNS
//
//  Created by Raphael Bost on 16/01/13.
//  Copyright (c) 2013 Raphael Bost. All rights reserved.
//

#ifndef __KMNS__max_value__
#define __KMNS__max_value__

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
    
    bool hasBeenSet() const { return _hasBeenSet;}
    T  value() const { assert(_hasBeenSet); return _value; }
};
#endif /* defined(__KMNS__max_value__) */
