//
//  field.h
//  Vorticity_simple
//
//  Created by Nader on 18/12/2017.
//  Copyright © 2017 Nader Ganaba. All rights reserved.
//

#ifndef field_h
#define field_h

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>

//Macro for Lexicographic order of R^2
#define idx_func(a,b,c,d,i) ((a*(nx_+1) + b*i)*(ny_+2) + c*(ny_+1) + d*i )
#define min(a, b)  (((a) < (b)) ? (a) : (b))
class Field{
public:
    typedef std::vector<double> StateType;
    typedef std::vector<double>::iterator Iterator;
private:
    std::vector<double> w_field_; //Scalar field to store values of vorticity function
    std::vector<double> v_field_; //Scalar field to store values of stream function
    std::vector<double> w_field_temp_; //Temporary scalar field to store values of stream function
    int nx_, ny_, total_size_;
    double visc_, dx_, dy_;
    
public:
    //Constructor for Field class
    Field(int nx, int ny):nx_(nx), ny_(ny), w_field_((nx + 2)*(ny + 2), 0.0f), v_field_((nx + 2)*(ny + 2), 0.0f), w_field_temp_((nx + 2)*(ny + 2), 0.0f){
        total_size_ = (nx_ + 2)*(ny_ + 2);
        visc_ = 0.1;
        dx_ = 1.0/(nx_+1);
        dy_ = 1.0/(ny_+1);
        
        std::cout << "Field is [" <<nx_ << ", " << ny_ << "] \t spacing [" << dx_ << ", " << dy_ << "]\t w_wall = " << (2.0/dy_) << std::endl;
    }
    
    //Destructor
    //~Field();
    
    //Getters and setters
    
    void setViscosity(double _val){ visc_ = _val; }
    double getViscosity(){ return visc_;}
    
    Iterator getW(int i, int j) {return w_field_.begin() + (i*(ny_+2) + j);}
    Iterator getWTemp(int i, int j) {return w_field_temp_.begin() + (i*(ny_+2) + j);}
    Iterator getV(int i, int j) {return v_field_.begin() +(i*(ny_+2) + j);}
    
    int getNx(){return nx_;}
    int getNy(){return ny_;}
    int getTotalSize(){return total_size_;}
    
    double getDx(){return dx_;}
    double getDy(){return dy_;}
    
    //The time step should satisfy the following condition \frac{\nu dt}{h^2} < \frac{1}{4}
    double estimateDt(){
        double result;
        double h_min = min(dx_, dy_);
        result = 0.25*(h_min*h_min)/visc_;
        result *= 0.75;
        return result;
    }
    
    void fieldSAXPY(){ std::transform(w_field_.begin( ), w_field_.end( ),  w_field_temp_.begin( ), w_field_.begin( ),std::plus<double>( ));}
    
    //Functions for initialising the scalar fields
    void initField(double val){
        for(int i = 0; i < total_size_; i++){
            w_field_.at(i) = val;
            v_field_.at(i) = val;
            w_field_temp_.at(i) = val;
        }
    }
    void initField(){
        initField(0.0);
    }
    
    void initV(){
        for(int i = 0; i < total_size_; i++){
            v_field_.at(i) = 0.0;
        }
    }
    
//
//    //Used for setting the boundary condition
//    void setBoundaryCondW(){
//        //vt(2:nx_-1,1)=-2.0*sf(2:nx_-1,2)/(h*h);
//        //vt(2:nx_-1,ny_)=-2.0*sf(2:nx_-1,ny_-1)/(h*h)-2.0/h;
//
//        int idx,idx_1;
//
//        for(int i = 1; i < nx_+1;i++){
//            idx = i*(ny_+2) ; //i*(ny_+2) + j
//            idx_1 = i*(ny_+2) + 1; //i*(ny_+2) + 1
//            w_field_.at(idx) = -2.0*(v_field_.at(idx_1)  )/(dy_*dy_);
//
//            idx = i*(ny_+2) + ny_+1; //i*(ny_+2) + ny_+1
//            idx_1 = i*(ny_+2) + ny_; //-2.0*sf(2:nx_-1,ny_-1)/(h*h)-2.0/h;
//            w_field_.at(idx) = (-2.0*(v_field_.at(idx_1)  )/(dy_*dy_) )- (2.0/dy_);
//
//        }
//
//        //        vt(1,2:ny_-1)=-2.0*sf(2,2:ny_-1)/(h*h); % vorticity on right wall
//        //        vt(nx_,2:ny_-1)=-2.0*sf(nx_-1,2:ny_-1)/(h*h); % vorticity on left wall
//        for(int j = 1; j < ny_+1;j++){
//            idx =  j; //i*(ny_+2) + j
//            idx_1 = 1*(ny_+2) + j;
//            w_field_.at(idx) = -2.0*(v_field_.at(idx_1) )/(dx_*dx_);
//
//            idx = (nx_+1)*(ny_+2) + j;
//            idx_1 = (nx_)*(nx_+2) + j;
//            w_field_.at(idx) = -2.0*(v_field_.at(idx_1)  )/(dx_*dx_);
//        }
//
//    }
//
    void printW(){
        Iterator w_it;
        for(int i = 0; i < nx_ + 2; i++){
            for(int j = 0; j < ny_ + 2; j++){
                w_it = w_field_.begin() + i*(ny_+2) + j;
                std::cout << (*w_it) << " ";
            }
            std::cout << std::endl;
        }
    }
    void printWT(){
        Iterator wt_it;
        for(int i = 0; i < nx_ + 2; i++){
            for(int j = 0; j < ny_ + 2; j++){
                wt_it = w_field_temp_.begin() + i*(ny_+2) + j;
                std::cout << (*wt_it) << " ";
            }
            std::cout << std::endl;
        }
    }
    
    void printV(){
        Iterator v_it;
        
        for(int i = 0; i < nx_ + 2; i++){
            for(int j = 0; j < ny_ + 2; j++){
                v_it = v_field_.begin() + i*(ny_+2) + j;
                std::cout << (*v_it) << " ";
            }
            std::cout << std::endl;
        }
    }
    
    void update(double dt){
        Iterator w_it;
        Iterator wt_it;
        int idx;
        for(int i = 1; i < nx_ + 1; i++){
            for(int j = 1; j < ny_ + 1; j++){
                idx = (i)*(ny_+2) + j;
                w_field_.at(idx) += dt*w_field_temp_.at(idx);
            }
            
        }
    }
    
    
    //Function used for Parallel implementation of the algorithm. This function takes the ghost fields (where the vorticity is stored in [0:N+1] elements and the stream field is stored in [N+2:2N+1] elements
    //corner_flag indicates which boarder is exchanged North = 0, South = 1, East = 2 and West = 3
    
    void setGhost(double slab[], int corner_flag){ //Slab contains the data exchanged between processes
        
        int limit  = 0;
        int params[4];
        switch (corner_flag) {
            case 0:             //North
                params[0] = 1; params[1] = 0; params[2] = 0; params[3] = 1;
                limit = ny_ + 2;
                break;
            case 1:             //South
                params[0] = 0; params[1] = 0; params[2] = 0; params[3] = 1;
                limit = ny_ + 2;
                break;
            case 2:             //East
                params[0] = 0; params[1] = 1; params[2] = 0; params[3] = 0;
                limit = nx_ + 2;
                break;
            case 3:             //West
                params[0] = 0; params[1] = 1; params[2] = 1; params[3] = 0;
                limit = nx_ + 2;
                break;
            default:
                break;
        }
        int idx_temp;
        for(int i = 0; i < limit; i++){
            idx_temp = idx_func(params[0],params[1],params[2],params[3],i);
            w_field_.at(idx_temp)  = slab[2*i];
            v_field_.at(idx_temp)  = slab[2*i+1];
        }
        
    }
    
    //Function used for Parallel implementation of the algorithm. This function returns the ghost fields (where the vorticity is stored in [0:N+1] elements and the stream field is stored in [N+2:2N+1] elements
    //corner_flag indicates which boarder is exchanged North = 0, South = 1, East = 2 and West = 3
    double* saveSlab( int corner_flag){
        
        int limit  = 0;
        int params[4];
        switch (corner_flag) {
            case 0:             //North
                params[0] = 1; params[1] = 0; params[2] = 0; params[3] = 1;
                limit = ny_ + 2;
                break;
            case 1:             //South
                params[0] = 0; params[1] = 0; params[2] = 0; params[3] = 1;
                limit = ny_ + 2;
                break;
            case 2:             //East
                params[0] = 0; params[1] = 1; params[2] = 0; params[3] = 0;
                limit = nx_ + 2;
                break;
            case 3:             //West
                params[0] = 0; params[1] = 1; params[2] = 1; params[3] = 0;
                limit = nx_ + 2;
                break;
            default:
                break;
        }
        
        double *slab = (double*) calloc (2*(limit),sizeof(double));
        if(slab == NULL){
            std::cerr << "Allocation failure in saveSlab() method";
        }
        
        int idx_temp;
        for(int i = 0; i < limit; i++){
            idx_temp = idx_func(params[0],params[1],params[2],params[3],i);
            slab[2*i] = w_field_.at(idx_temp);
            slab[2*i+1] = v_field_.at(idx_temp) ;
        }
        
        return slab;
    }
    
    
}; //End Field class
#endif /* field_h */
