//
//  main.cpp
//  vorticity_lite
//
//  Copyright © 2017 Nader Ganaba. All rights reserved.
//  License:       Apache License 2.0
//  Created by Nader on 19/12/2017.
//

 

#include <iostream>
#include "field.h"
#include "scheme.h"
#include "timeintegration.h"
#include "pde_def.h"
#include <vector>
#include "visualise.h"

int main(int argc, const char * argv[]) {
    
    int Nx = 32;
    int Ny = 32;
    double time = 0.0;
    double dt ;
    double T_end = 1.2;
    int N_iter;
    
    
    Field _scalarField(Nx,Ny);
    _scalarField.initField();
    dt = _scalarField.estimateDt();
    N_iter = (int)(T_end/dt) + 1;
    
    std::cout << "T_end = " << T_end << "\t dt = " << dt << "\t N_iter = " << N_iter << std::endl;
    std::cout << "T_end est = " << (dt*N_iter) << std::endl;
    std::cout << std::endl;
    typedef VNS_EQ<Field> SchemeType;
    VNS_EQ<Field> scheme("vorticity");
    
    //SORMethod<Field> _ellipSolver(0.001, 100, 1.5);
    ConjugateGrad<Field, SchemeType> _ellipSolver(0.001, 1000);
    timeAdvance<Field, SchemeType> _timeStepper(dt);
    
    
    
    //  int N_iter = 60;
    for(int i = 0; i <   N_iter ; i++){
        time += dt;
        std::cout << "time = " << time << std::endl;
        
        _ellipSolver.solve(_scalarField, scheme);
        _timeStepper.nextTimeStep(_scalarField, scheme);
        
    }
    
    std::cout << std::endl;
    _scalarField.printW();
    
    std::cout << std::endl;
    _scalarField.printV();
    

    pgmFile<Field> _test("VNS_test");
    _test.saveData(_scalarField);
    //    //_scalarField.saveData("tst2.pgm");
    
    
    return 0;
}

