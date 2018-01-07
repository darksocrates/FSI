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
    if (this->verbose == true){
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
        translation.resize(n,1);
    }
    
    for(int k =0; k<n;++k) {
        weights(k,0) = eta*exp(onei*2.0*pi/double(n)*double(k));
        translation(k,0) = weights(k,0)+center;
    }
}

void Filter::string_call(std::string filter){
    if (filter == "circ_trapez_snug"){
        Filter::circ_trapez_snug();
    } else{
        throw("The requested filter is not implemented");
    }
};

//Filter Constructors

Filter::Filter(unsigned int n){
    //return next highest even number from initialized n.
    this->n = n + n%2;
}

Filter::Filter(unsigned int n, std::string filter){
    this->n = n+ n%2;
    Filter::string_call(filter);
}



//Spectralproj methods

void Spectralproj::prnt(std::string str) {
    if (this->verbose == true){
        std::cout << str <<"\n";
    }
}

VectorXcd Spectralproj::feast_step_hermitian(Span* q){
    
    //apply spectral projector to input span
    Spectralproj::mult(q);
    
    //call the rayliegh ritz procedure implemented in derived class
    MatrixXcd qAq;
    MatrixXcd qq;
    std::tie(qAq,qq) = rayleigh(q);
    
   //NOTE: Replace with try catch statement once the nature of Eigen's error system for failed convergence is known. For now the algorithm assumes that this works, which is risky.
    Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXcd> es(qAq,qq);
    
    //recover the span
    q->operator*(es.eigenvectors());
    
    //output the eigenvalues
    return es.eigenvalues();
}



//Spectralproj constructors
