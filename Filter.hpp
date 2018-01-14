//
//  Filter.hpp
//  FSI_2
//
//  Created by Aidan Hamilton on 1/13/18.
//  Copyright © 2018 Aidan Hamilton. All rights reserved.
//

#ifndef Filter_hpp
#define Filter_hpp

#include <stdio.h>
#include <Eigen/Core>
#include <string>
#include <iostream>
#include <complex>
#include <cmath>
#include <tuple>
#include <Eigen/Eigenvalues>
#include <vector>
#include "Typedefs_statics.hpp"



//implements the filters for filtered subspace iteration in its methods. This is seperated from the class Spectralproj for readability
// This class will also implement all of the 'discretization' of the eigenvalues search set. Ie. it will decide how many seperate fsi algorithms to run


//all methods must write to the internal complex weights, complex translations and number of points in the filter (n) variables.

class Filter {
    //add more to here as you add new filter methods. 
    enum filter_type {
        circle_trapezoid_snug
    };

    
private:
    
    //FIELDS
    unsigned int n = 0;
    //NOTE: see if you can figure out a way to change this to a static Matrix. Would be preferable for efficiency
    VectorXcd weights,translation;
    // internal variable for storing the filter type of the current Filter object. Used to determine which version of the inside method to call.
    filter_type currentfilter;
    
    //METHODS
    
    // the inside methods, to determine if some given points are inside of the filter. Must take as arguments VectorXcd points and double fudge and return a pointer to an unsigned array containing the indices of VectorXcd whose values are not inside the filter. 
    unsigned int* circ_trapez_snug_inside(const VectorXcd &points, const double &fudge);
    
    
public:
    //FIELDS
    bool verbose = 0; double radius; std::complex<double> center;
    
    //METHODS
    
    // filter methods
    void circ_trapez_snug();
    
    
    //general methods
    void setn(unsigned int n);
    
    unsigned int* inside(const VectorXcd &points, const double &fudge);
    
    // quality of life methods
    void prnt(std::string str);
    void string_call(std::string filter);
    
    //return private fields methods
    VectorXcd getweights(){return weights;}
    VectorXcd gettranslation(){return translation;}
    
    //Constructors
    Filter(unsigned int n);
    
    Filter(unsigned int n, std::string filter);
        
};



#endif /* Filter_hpp */
