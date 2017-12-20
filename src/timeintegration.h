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
 * timeintegration.h
 */

#ifndef timeIntegration_h
#define timeIntegration_h

#include "field_traits.h"
template<
class FieldClass,
class Scheme,
class Time = double,
class Traits = field_traits< FieldClass >
>
class timeAdvance{
    
public:
    typedef Time TimeType;
    typedef Traits TraitsType;
    typedef typename TraitsType::StateType    StateType;
    typedef typename TraitsType::Iterator      Iterator;
    
    timeAdvance(Time dt){dt_ = dt;}
    
    void nextTimeStep(FieldClass& _field, Scheme& _scheme){
        
        
        
        //Defining the needed iterators
        Iterator it_temp, it_w, it_v, it_temp1, it_w1;
        
        
        
        //Needed integers and coefficients
        int x_size = _field.getNx();
        int y_size = _field.getNy();
        double dx = _field.getDx();
        double dy = _field.getDy();
        double visc = _field.getViscosity();
        double temp = 0.0;
        
        //Setting the boundary conditions
        //_field.setBoundaryCondW();
        
        it_w = _field.getW(0,0);
        it_v = _field.getV(0,0);
        _scheme.boundaryConditions(it_w,  it_v,  x_size,  y_size ,  dx,  dy);
        
        int idx;
        for(int i = 1; i < x_size+1; i++){
            for(int j = 1; j < y_size+1; j++){
                it_temp = _field.getWTemp(i,j);
                it_w = _field.getW(i,j);
                it_v = _field.getV(i,j);
                idx = i*(y_size + 2) + j;
                temp = _scheme.timeIteration(it_w,it_v, dx, dy,  visc,   y_size +2  );
                
                (*it_temp) =  dt_*temp;
            }
        }
        
        //Time stepping
        _field.fieldSAXPY();
        
    }
    
private:
    Time dt_;
    
};
#endif /* timeIntegration_h */
