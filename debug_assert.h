//
//  debug_assert.h
//  SubmatrixQueries
//
//  Created by Raphael Bost on 22/03/13.
//  Copyright (c) 2013 Raphael Bost. All rights reserved.
//

#ifndef SubmatrixQueries_debug_assert_h
#define SubmatrixQueries_debug_assert_h

#include <cassert>

#ifdef DEBUG
#define DEBUG_ASSERT(e) assert(e)
#else
#define DEBUG_ASSERT(e) 
#endif

#endif
