//
//  SpectralProj.cpp
//  FSI
//
//  Created by Aidan Hamilton on 1/1/18.
//
//

#include "SpectralProj.hpp"


// Span methods

void Span::prnt(std::string str) {
    if (this->verbose ==true){
        std::cout << str <<"\n";
    }
    }

// while this function cannot be made purely virtual as it is templated, so I can't have it abstract, without overriding it it simply throws an error. 
template<typename Derived>
std::complex<double> Span::innerproduct(const Eigen::MatrixBase<Derived>& u, const Eigen::MatrixBase<Derived>& v){
    
    throw err;
}



// Filter methods

void Filter::prnt(std::string str) {
    if (this->verbose ==true){
        std::cout << str <<"\n";
    }
}

void Filter::circ_trapez_snug(){
    
    prnt("Setting snug trapezoidal filter on circular contour");
    double eta = radius/(pow(2.0,(1.0/double(n))));
    
    if(weights.rows()!=n){
        weights.resize(n,1);
    }
    if(translation.rows() !=  n){
        weights.resize(n,1);
    }
    
    for(int k =0; k<n;++k) {
        weights(k,1) = eta*exp(onei*2.0*pi/double(n)*double(k));
        translation(k,1) = weights(k,1)+center;
    }
}

void Filter::string_call(std::string filter){
    if (filter == "circ_trapez_snug"){
        Filter::circ_trapez_snug();
    } else{
        throw("The requested filter is not implemented");
    }
};

//Constructors

Filter::Filter(unsigned int n){
    //return next highest even number from initialized n.
    this->n = n + n%2;
}

Filter::Filter(unsigned int n, std::string filter){
    this->n = n+ n%2;
    Filter::string_call(filter);
}
