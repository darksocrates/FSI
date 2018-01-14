//
//  Filter.cpp
//  FSI_2
//
//  Created by Aidan Hamilton on 1/13/18.
//  Copyright © 2018 Aidan Hamilton. All rights reserved.
//

#include "Filter.hpp"



// Filters

void Filter::circ_trapez_snug(){
    
    prnt("Setting snug trapezoidal filter on circular contour");
    double eta = radius/(pow(2.0,(1.0/double(n))));
    
    if(weights.rows()!=n){
        weights.resize(n,1);
    }
    if(translation.rows() !=  n){
        translation.resize(n,1);
    }
    
    for(int k =0; k<n;++k) {
        weights(k,0) = eta*exp(onei*2.0*pi/double(n)*double(k));
        translation(k,0) = weights(k,0)+center;
    }
    
    currentfilter = circle_trapezoid_snug;

}


// inside methods, methods to determine if the given values inside of an complex Eigen vector are inside or outside the filter. Return an array containing the indices of the values in the given complex Eigen Vector that are outside fo the filter. These methods are all private.
unsigned int* Filter::circ_trapez_snug_inside(const VectorXcd &points, const double &fudge){
    std::vector<unsigned int> indexholder;
    unsigned int sizehold = 0;
    //iterate through the inputted vector of points to check.
    for (int i = 0; i <points.size(); i++){
        //checks if the value specified at points[i] is outside of the circle of radius fudge+this->radius centered at this->center.
        if ((points[i].real()-center.real())*(points[i].real()-center.real())+(points[i].imag()-center.imag())*(points[i].imag()-center.imag()) <= ((fudge+radius)*(fudge+radius))){
            indexholder.push_back(i); sizehold++;
        }
    }
    //create an array that contains the indices in the vector indexholder.
    if(sizehold != 0){
    unsigned int indices[sizehold];
        //fill array
        for (int  i =0; i<sizehold; i++){
            indices[i] = indexholder[i];
        }
        return indices;
        
    } else {
        // if there were no points outside of the filter, return a null pointer.
        return 0;
    }
}



//general methods


unsigned int* Filter::inside(const VectorXcd &points, const double &fudge){
    switch (currentfilter){
        case circle_trapezoid_snug:
            auto output =  this->circ_trapez_snug_inside(points, fudge);
            return output;
            break;
    }
}

void Filter::setn(unsigned int n){
    //return next highest even number from initialized n.
    this->n = n + n%2;
}

//quality of life methods

void Filter::string_call(std::string filter){
    if (filter == "circ_trapez_snug"){
        Filter::circ_trapez_snug();
    } else{
        throw("The requested filter is not implemented");
    }
};


void Filter::prnt(std::string str) {
    if (this->verbose == true){
        std::cout << str <<"\n";
    }
}

//Filter Constructors

Filter::Filter(unsigned int n){
    //return next highest even number from initialized n.
    this->n = n + n%2;
}

Filter::Filter(unsigned int n, std::string filter){
    this->n = n+ n%2;
    Filter::string_call(filter);
}