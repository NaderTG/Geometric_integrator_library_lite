//
//  scheme.h
//  Vorticity_simple
//
//  Created by Nader on 19/12/2017.
//  Copyright © 2017 Nader Ganaba. All rights reserved.
//

#ifndef scheme_h
#define scheme_h

#include <iostream>
#include <cmath>
#include "field_traits.h"
#include "abstractScheme.h"
//
template<
class FieldClass,
class Scheme,
class Time = double,
class Traits = field_traits< FieldClass >
>
class SORMethod: public AbstractSolver<FieldClass, Scheme>{
    
public:
    typedef Time TimeType;
    typedef Traits TraitsType;
    typedef typename TraitsType::StateType    StateType;
    typedef typename TraitsType::Iterator      Iterator;
    
    //template<class Scheme>
    SORMethod(double tol, int max_iter, double omg){max_iter_ = max_iter; tol_ = tol; omega_ = omg;}
    
    void solve( FieldClass& _field, Scheme& _scheme){
        //The way it should work is that it should take the section as the state and
        double res = 0.0; //change it one
        
        
        int num_iters = 0;
        
        int x_size = _field.getNx();
        int y_size = _field.getNy();
        double dx = _field.getDx();
        double dy = _field.getDy();
        double dx2 = 1.0/(dx*dx);
        double dy2 = 1.0/(dy*dy);
        double temp;
        Iterator  it_w, it_v;
        
        while(num_iters < max_iter_ && res < tol_){
            res = 0.0;
            
            for(int i = 1; i < x_size+1; i++){
                for(int j = 1; j < y_size+1; j++){
                    
                    it_w = _field.getW(i, j);
                    it_v = _field.getV(i, j);
                    
                    temp = _scheme.iteration(it_v, it_w,dx2, dy2, x_size +2);
                    (*it_v) = omega_*temp + (1.0 - omega_)*(*it_v);
                    
                }
            }
            
            
            
            for(int i = 1; i < x_size+1; i++){
                for(int j = 1; j < y_size+1; j++){
                    it_w = _field.getW(i, j);
                    it_v = _field.getV(i, j);
                    //res +=  residual(it_v,it_w, dx2, dy2, x_size +2);
                    temp = _scheme.residual(it_v,it_w,dx2, dy2, x_size +2);
                    res += temp*temp;
                }
            } //Residual
            // res = sqrt(res);
            
            num_iters++;
        } //While main loop
        
        
        
    }
private:
    //StateType *_x0; //initial condition
    int max_iter_;
    double tol_;
    double omega_;
};


template<
class FieldClass,
class Scheme,
class Time = double,
class Traits = field_traits< FieldClass >
>
class ConjugateGrad : public AbstractSolver<FieldClass, Scheme> {
    
public:
    typedef Time TimeType;
    typedef Traits TraitsType;
    typedef typename TraitsType::StateType    StateType;
    typedef typename TraitsType::Iterator      Iterator;
    
    //template<class Scheme>
    ConjugateGrad(double tol, int max_iter){max_iter_ = max_iter; tol_ = tol;}
    
    void solve( FieldClass& _field, Scheme& _scheme){
        //The way it should work is that it should take the section as the state and
        
        int total_size = _field.getTotalSize();
        StateType residual_field(total_size);
        StateType A_times_p(total_size);
        StateType p_dir(total_size);
        
        Iterator it_temp, it_w, it_v, it_p1,it_r,it_v1;
        
        int num_iters = 0;
        double delta0, delta1, alpha, beta, temp;
        int x_size = _field.getNx();
        int y_size = _field.getNy();
        double dx = _field.getDx();
        double dy = _field.getDy();
        
        double dx2 = 1.0/(dx*dx);
        double dy2 = 1.0/(dy*dy);
        
        
        delta0 = 0.0;
        
        _field.initV();
        for(int i = 1; i < x_size+1; i++){
            for(int j = 1; j < y_size+1; j++){
                it_temp = residual_field.begin() + (i*(y_size+2) +j);
                it_w = _field.getW(i, j);
                it_v = _field.getV(i, j);
                (*it_temp)= _scheme.residual(it_v, it_w,dx2, dy2,  y_size +2 );
                p_dir.at(i*(y_size+2) +j) = (*it_temp);
                delta0 += (*it_temp)*(*it_temp);
            }
        }
       
        if(delta0 < tol_){return;}
        
        do{
            alpha = 0.0;
            
            for(int i = 1; i < x_size+1; i++){
                for(int j = 1; j < y_size+1; j++){
                    it_v1 = A_times_p.begin() + (i*(y_size+2) +j);
                    it_p1 = p_dir.begin() + (i*(y_size+2) +j);
                    
                    (*it_v1) = _scheme.iteration(it_p1, dx2, dy2, y_size +2 );
                    
                    alpha += (*it_p1)*(*it_v1);
                }
            }
            
            alpha = delta0 / alpha;
            
            delta1 = 0.0;
            
            for(int i = 1; i < x_size+1; i++){
                for(int j = 1; j < y_size+1; j++){
                    it_v = _field.getV(i, j);
                    temp = residual_field[i*(y_size+2) +j] - alpha*A_times_p[i*(y_size+2) +j];
                    
                    (*it_v) += alpha*(p_dir[i*(y_size+2) +j]) ;
                    residual_field[i*(y_size+2) +j] = temp;
                    delta1 += temp*temp;
                }
            }
            if(sqrt(delta1) < tol_){break; }
            
            beta = delta1 / delta0;
            for(int i = 0; i < x_size+1; i++){
                for(int j = 0; j < y_size+1; j++){
                    
                    p_dir[i*(y_size+2) +j] = residual_field[i*(y_size+2) +j] + beta*p_dir[i*(y_size+2) +j];
                }
            }
            delta0 = delta1;
            
            
            num_iters++;
            
        }while(num_iters < max_iter_ && delta0 > tol_); //While main loop
        
         std::cout << "Delta = " << delta0  << "\t Iters = " << num_iters << std::endl;
    }
private:
    
    int max_iter_;
    double tol_;
};

#endif /* scheme_h */
