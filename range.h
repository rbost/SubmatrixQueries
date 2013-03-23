//
//  range.h
//  SubmatrixQueries
//
//  Created by Raphael Bost on 21/02/13.
//  Copyright (c) 2013 Raphael Bost. All rights reserved.
//

#ifndef __SubmatrixQueries__range__
#define __SubmatrixQueries__range__

#include <algorithm>
#include <cassert>

/*
 *  class Breakpoint
 *  This structure represents a range of indices
 *
 */

struct Range {
    const size_t min, max;
    
    inline Range(size_t minimun, size_t maximum) : min(minimun), max(maximum){
        assert(minimun <= maximum);
    }
    
    inline Range(Range const &r) : min(r.min), max(r.max)
    {
    }
    
    inline bool isInRange(size_t i)
    {
        return (i >= min) && (i <= max);
    }
    
    inline bool intersects(Range r)
    {
        if (r.max > max) {
            return r.min <= max;
        }else if (r.max >= min){
            return true;
        }
        return false;
    }
    
    inline Range intersection(Range r)
    {
        assert(intersects(r));
        
        return Range(std::min<size_t>(min,r.min),std::max<size_t>(max,r.max));
    }
    
    inline bool contains(Range r)
    {
        return min <= r.min && max >= r.max;
    }
};



#endif /* defined(__SubmatrixQueries__range__) */
