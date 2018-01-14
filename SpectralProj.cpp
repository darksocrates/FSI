//
//  SpectralProj.cpp
//  FSI
//
//  Created by Aidan Hamilton on 1/1/18.
//
//

#include "SpectralProj.hpp"


//function I needed for now
template<class T>
Eigen::Matrix<T, Eigen::Dynamic,1>  splice(Eigen::Matrix<T, Eigen::Dynamic,1> A, unsigned int *cut) {
    
    Eigen::Matrix<T, Eigen::Dynamic,1>  B(A.size());
    int index = 0; int l = sizeof(*cut)/sizeof(unsigned int);
    for (int i = 0; i<A.size(); i++){
        if (*(cut+index) == i){
            if (index < l){
            index++;
            }
        } else {
            B(i - index) = A(i);
        }
    }
    return B;
}





// Span methods

void Span::prnt(std::string str) {
    if (this->verbose ==true){
        std::cout << str <<"\n";
    }
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

    //create vector of complex Eigen vectors to store eigenvalues. reserve the maximum amount of memeory that could possibly be required.
    std::vector<VectorXcd> eigenvalues;
    eigenvalues.reserve(maxiterations);
    
    // initialize necessary variable, this is initialized to be larger than stop_tol to allow me to place an if then statement in a certain spot to make the rest of the code cleaner.
    double maxerr = 2*stop_tol;
    
    //required to initialize this variable even if  I don't use it. However I don't reserve any memory unless the residuals are actually going to be calculated. 
    std::vector<Eigen::VectorXd> resi;
    if (report_residuals == true) {
    resi.reserve(maxiterations);
    }
    
    
    Spectralproj::prnt("\n ============Staring FSI iterations==========" );
    Spectralproj::prnt("Trying with " + std::to_string(q->numcolumns())+ " columns in the span \n with filter with center " +std::to_string(filter->center.real())+ "+ " + std::to_string(filter->center.imag()) + "i and radius " + std::to_string(filter->radius) +"\n");
    
    for(int iter = 0; iter<maxiterations; iter++) {
        Spectralproj::prnt("Iteration " + std::to_string(iter)+"\n");
        auto eigenhold = Spectralproj::fsi_step_hermitian(q);
        
        //Checks if the calculated eigenvalues are outside of the given filter. If so removes the corresponding elements from the span. Then calculates the error between succesive iterations for the remaining eigenvalues. If the max of these errors is less than stop_tol, the FSI iterations stop.
        if (iter>0){
            Eigen::VectorXd errV = (eigenhold - eigenvalues[iter-1]).cwiseAbs();
            auto indexs = filter->inside(eigenvalues[iter], eigen_fudge);
            if (indexs != nullptr) {
                q->remove(indexs);
                eigenvalues.push_back(splice(eigenhold,indexs));
                maxerr = splice(errV,indexs).maxCoeff();
            }else{
                eigenvalues.push_back(eigenhold);
                maxerr = errV.maxCoeff();
            }
            
        } else {
        //stuff that is done during the first iteration only, when iteration errors can't be calculated or the eigenvalues, at least not really.
            eigenvalues.push_back(eigenhold);
        }
        
        if (report_residuals == true){
            resi.push_back(this->residual(eigenvalues[iter], q));
        }
        
        
        if (stop_tol>=maxerr){
            prnt("Stop tolerance reached. Endin FSI");
            break;
        }
        
    }
    //remove excess memory from eigenvalues, useful when instance stop_tol was reached before maxiterations.
    eigenvalues.shrink_to_fit();
    
    
    //create FSI_data object to return, a different constructor is called depending on whether residuals were asked or or not.
    if (report_residuals == false){
        FSI_data data(eigenvalues, q,this->whatareyou());
        return data;
    } else {
        resi.shrink_to_fit();
        FSI_data data(eigenvalues, q,this-> whatareyou(), resi);
        return data;
    }
}

//Spectralproj constructors



//FSI_data Methods

//FSI_data Constructors

FSI_data::FSI_data( std::vector<VectorXcd> eigenvalues, Span *span, std::string spectralproj_type){
    
    this->span = span;
    this->eigenvalues = eigenvalues;
    this->spectralproj_type = spectralproj_type;
    
}

FSI_data::FSI_data( std::vector<VectorXcd>  eigenvalues, Span *span, std::string spectralproj_type, std::vector<Eigen::VectorXd> residual){
    this->span = span;
    this->eigenvalues = eigenvalues;
    this->spectralproj_type = spectralproj_type;
    this->residual_calculated = true;
    this->residual = residual; 
}

