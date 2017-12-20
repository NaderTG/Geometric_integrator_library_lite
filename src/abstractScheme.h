//
//  abstractScheme.h
//  vorticity_lite
//
//  Created by Nader on 20/12/2017.
//  Copyright © 2017 Nader Ganaba. All rights reserved.
//

#ifndef abstractScheme_h
#define abstractScheme_h

#include "field_traits.h"


template<
class FieldClass,
class Scheme,
class Time = double ,
class Traits = field_traits< FieldClass >
>
class AbstractSolver{
public:
    typedef Time TimeType;
    typedef Traits TraitsType;
    
    typedef typename TraitsType::StateType    StateType;
    typedef typename TraitsType::Iterator      Iterator;
    //Add section type
private:
    //        double tol;
    //        int max_iter;
    
public:
    
    virtual void solve( FieldClass& _state, Scheme& _scheme) = 0;
    //void setTolerance(double _tol){tol = _tol;}
    
};

#endif /* abstractScheme_h */
