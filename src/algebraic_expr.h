/*
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Copyright (C) 2017 Nader Ganaba
 * Created by Nader on 19/12/2017.
 * algebraic_expr.h
 */

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
