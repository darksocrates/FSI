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


void Filter::setn(unsigned int n){
    //return next highest even number from initialized n.
    this->n = n + n%2;
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



//Spectralproj methods

void Spectralproj::prnt(std::string str) {
    if (this->verbose == true){
        std::cout << str <<"\n";
    }
}

VectorXcd Spectralproj::fsi_step_hermitian(Span* q){
    
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

FSI_data Spectralproj::fsi_hermitian(Span *q, double stop_tol,double eigen_fudge , unsigned int maxiterations, bool report_residuals){

    //NOTE: Shouldn't I just hold my eigenvalues as a vector<vector<std::complex>>, It fits quite well, and the point of having the data be a specific object was so how I store the data is up to my discretion, any operations that come up that require a Eigen vector of eigenvalues could have the data change happen behind the scenes. 
    MatrixXcd eigenvalues;
    
    
    Spectralproj::prnt("\n ============Staring FSI iterations==========" );
    Spectralproj::prnt("Trying with " + std::to_string(q->numcolumns())+ " columns in the span \n with filter with center " +std::to_string(filter->center.real())+ "+ " + std::to_string(filter->center.imag()) + "i and radius " + std::to_string(filter->radius) +"\n");
    
    for(int iter = 0; iter<maxiterations; iter++) {
        Spectralproj::prnt("Iteration " + std::to_string(iter)+"\n");
        VectorXcd eigenhold = Spectralproj::fsi_step_hermitian(q);
        
        
        
    
    }

}

//Spectralproj constructors



//FSI_data Methods

//FSI_data Constructors

FSI_data::FSI_data(MatrixXcd eigenvalues, Span *span, std::string spectralproj_type){
    
    this->span = span;
    this->eigenvalues = eigenvalues;
    this->spectralproj_type = spectralproj_type;
    
}

FSI_data::FSI_data(MatrixXcd eigenvalues, Span *span, std::string spectralproj_type, Eigen::MatrixXd residual){
    this->span = span;
    this->eigenvalues = eigenvalues;
    this->spectralproj_type = spectralproj_type;
    this->residual_calculated = true;
    this->residual = residual; 
}

