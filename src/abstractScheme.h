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
 * abstractScheme.h
 */

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
