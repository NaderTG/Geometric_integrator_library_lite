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
 * pde_def.h
 */

#ifndef pde_def_h
#define pde_def_h
#include <string>
#include "algebraic_expr.h"


template< class FieldClass,
class Traits = field_traits< FieldClass >
>

//Definition of the stencil used for vorticity formulation of the Navier-Stokes equation
class VNS_EQ : public Stencil<FieldClass> {
private:
    std::string pde_name_;
public:
    
    typedef Traits TraitsType;
    typedef typename TraitsType::StateType StateType;
    typedef typename TraitsType::Iterator Iterator;
    
    explicit
    
    VNS_EQ(std::string pde_name){
        pde_name_ = pde_name;
    }
    void printName(){
        std::cout << pde_name_ << std::endl;
    }
    //Laplacian part of the vorticity-stream Laplace equation
    //Ap = \frac{p_{i,j+1} -2p_{i,j} + p_{i,j-1}}{dy^2} + \frac{p_{i+1,j} -2p_{i,j} + p_{i-1,j}}{dx^2}
    double iteration(Iterator iter1,  double alpha, double beta, int x_size){
        double result = alpha*(-2*(*iter1) + (*(iter1 + 1)) + (*(iter1 - 1))) + beta*(-2*(*iter1) + (*(iter1 + x_size)) + (*(iter1 - x_size)));
        return result;
    }
    
    using Stencil<FieldClass>::iteration;
    //\frac{v_{i,j+1} -2v_{i,j} + v_{i,j-1}}{dy^2} + \frac{v_{i+1,j} -2v_{i,j} + v_{i-1,j}}{dx^2} = -w_{ij}
    double iteration(Iterator iter1, Iterator iter2,  double alpha, double beta, int x_size){
        double result = (-(*iter2) - alpha*( (*(iter1 + 1)) + (*(iter1 - 1))) - beta*((*(iter1 + x_size)) + (*(iter1 - x_size))) ) / (-2*(alpha + beta));
        return result;
        
    }
    
    //Residual vector
    //r_{i,j} = -\frac{v_{i,j+1} -2v_{i,j} + v_{i,j-1}}{dy^2} + \frac{v_{i+1,j} -2v_{i,j} + v_{i-1,j}}{dx^2} = -w_{ij}
    double residual(Iterator iter1, Iterator iter2, double alpha, double beta, int x_size){
        double result;
        result = -alpha*(-2*(*iter1) + (*(iter1 + 1)) + (*(iter1 - 1))) - beta*(-2*(*iter1) + (*(iter1 + x_size)) + (*(iter1 - x_size))) - (*iter2);
        return result;
    }
    
    //Time step
    //w_{i,j}^{n+1} = w_{i,j}^{n} + dt ( \nu \frac{w_{i,j+1} -2w_{i,j} + w_{i,j-1}}{dy^2} + \frac{w_{i+1,j} -2w_{i,j} + w_{i-1,j}}{dx^2}) + dt(\frac{(w_{i,j+1}-w_{i,j-1})(v_{i+1,j} + v_{i-1,j})}{4 dy dx}) - dt(\frac{(v_{i,j+1}-v_{i,j-1})(w_{i+1,j} + w_{i-1,j})}{4 dy dx})
    double timeIteration(Iterator iter1, Iterator iter2, double alpha, double beta, double gamma, int x_size){
        double  visc =gamma*((1.0/(alpha*alpha))*(-2*(*iter1) + (*(iter1 + 1)) + (*(iter1 - 1))) + (1.0/(beta*beta))*(-2*(*iter1) + (*(iter1 + x_size)) + (*(iter1 - x_size))));
        
        double adv1 = (0.25/(beta*alpha))*(  (*(iter2 + 1)) - (*(iter2 - 1)))*( (*(iter1 + x_size)) - (*(iter1 - x_size)));
        double adv2 =-(0.25/(beta*alpha))*((*(iter2 + x_size)) - (*(iter2 - x_size)))*(  (*(iter1 + 1)) - (*(iter1 - 1)));
        
        return   (adv1 + adv2 +visc );
        
    }
    
    //Boundary condition
    void boundaryConditions(Iterator iter1,  Iterator iter2, int nx, int ny, double dx, double dy){
        
        int idx,idx_1;
        
        for(int i = 1; i < nx+1;i++){
            //w_{i,0} = -2\frac{v_{i, 1}}{dy^2}
            idx = i*(ny+2) ;
            idx_1 = i*(ny+2) + 1;
            (*(iter1 + idx)) = -2.0*((*(iter2 + idx_1))  )/(dy*dy);
            
            //Moving wall
            //w_{i,ny+1} = -2\frac{v_{i, ny}}{dy^2} - 2 \frac{U_{wall}}{dy}
            idx = i*(ny+2) + ny +1;
            idx_1 = i*(ny +2) + ny ;
            (*(iter1 + idx))  = (-2.0*((*(iter2 + idx_1)))/(dy*dy) )- (2.0/dy);
            
        }
        
        
        for(int j = 1; j < ny+1;j++){
            //w_{0,j} = -2\frac{v_{1, j}}{dy^2}
            idx =  j;
            idx_1 = 1*(ny+2) + j;
            (*(iter1 + idx))  = -2.0*((*(iter2 + idx_1)) )/(dx*dx);
            
            //w_{nx+1,j} = -2\frac{v_{nx, j}}{dy^2}
            idx = (nx+1)*(ny+2) + j;
            idx_1 = (nx)*(nx+2) + j;
            (*(iter1 + idx))  = -2.0*((*(iter2 + idx_1)) )/(dx*dx);
        }
    }
    
}; //VNS_EQ class

#endif /* pde_def_h */
