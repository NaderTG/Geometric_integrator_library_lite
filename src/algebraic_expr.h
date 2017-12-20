//
//  algebraic_expr.h
//  vorticity_lite
//
//  Created by Nader on 20/12/2017.
//  Copyright © 2017 Nader Ganaba. All rights reserved.
//

#ifndef algebraic_expr_h
#define algebraic_expr_h

#include "field_traits.h"



template< class FieldClass,
class Traits = field_traits< FieldClass >
>
//Abstract class to help defining algebraic expressions
class Stencil{
public:
    typedef Traits TraitsType;
    typedef typename TraitsType::StateType StateType;
    typedef typename TraitsType::Iterator Iterator;
    
    virtual double iteration( Iterator, double, double, int) = 0;
    //virtual void iteration(Iterator ,  Iterator, Iterator, double, double, int) = 0;
    virtual double timeIteration( Iterator, Iterator, double, double, double, int) = 0;
    virtual double residual(Iterator , Iterator, double, double, int)= 0;
    virtual void  boundaryConditions(Iterator, Iterator, int, int, double, double) = 0;
    
    
}; //Scheme abstract class



#endif /* algebraic_expr_h */
