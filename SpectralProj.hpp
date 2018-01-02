//
//  SpectralProj.hpp
//  FSI
//
//  Created by Aidan Hamilton on 1/1/18.
//
//


/* Creates an abstract class for defining the spectral projector used in filtered subspace iterations (FSI). Defines the methods to be implemented by a derived class and explictly declares how FSI (or FEAST) will operate.
 
 Two Classes are defined
 
 Span - a span of abstract vectors, contains implementation of inner product, orthagonaliation, removing elements.
 
 Spectralproj - Implementation of FSI (FEAST), depends on Span. */

#ifndef SpectralProj_hpp
#define SpectralProj_hpp

#include <stdio.h>
#include <Eigen/Core>
#include <string>
#include <iostream>
#include <complex>
#include <cmath>
#include <tuple>


const static double pi =3.14159265358979323846;
const static std::complex<double> onei(0,1);

class Span {
    
private:
//internal values used in the methods. m corresponds to the number of rows in the Span, n is the number of columns. Verbose is designater used to determine when to print various things when methods are called. err is a string that holds the standard error message to display for functions that are optional to override in a derived class (they are not pure virtual functions), but do not have any meaning here.
    
    unsigned int m,n;
    std::string err = "Function not implemented in derived class, base Span function called";

    
    
public:
    
    bool verbose =0;
//The following functions have a sensible definition so do not require overriding in a derived class.
    
     void prnt(std::string str);
    
//The following functions must be implemented by derived classes
    
    virtual void orthagonalize() = 0;
    /* replace the current span by an orthagonal span in a suitable H space */
    
    virtual void remove(unsigned int *cut[]) =0;
    /* remove vectors with indices given by the array cut*/
    
    
    template<typename Derived>
    std::complex<double> innerproduct(const Eigen::MatrixBase<Derived>& u, const Eigen::MatrixBase<Derived>& v);
    /* required implementation of inner product for two memberes in the underlying space of the Span, output is assumed to complex. This function assumes that they underlying space of the span is represented by an Matrix (Vector) from the Eigen library.*/
    
};



//implements the filters for filtered subspace iteration in its methods. This is seperated from the class Spectralproj for readability

//all methods must write to the internal complex weights, complex translations and number of points in the filter (n) variables.

class Filter {
private:
    unsigned int n = 0;
    //NOTE: see if you can figure out a way to change this to a static Matrix. Would be preferable for efficiency
    Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> weights,translation;
public:
    //fields
    bool verbose = 0; double radius; std::complex<double> center;
    
    //methods
    
    // filter methods
    void circ_trapez_snug();
    
    // quality of life methods
    void prnt(std::string str);
    void string_call(std::string filter);
    
    //return private fields methods
    Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> getweights(){return weights;}
    Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> gettranslation(){return translation;}
    
    //Constructors
    Filter(unsigned int n);
    
    Filter(unsigned int n, std::string filter);
};


class Spectralproj {
private:
public:
    //the following methods must be implemented by derived classes
    
    //overload * operator to corrrespond with the action of the particular operator in a derived class on a Span.
    
    virtual Span operator *(Span y) = 0;
    
    //compute (Aq,q) and (q,q) in the H-inner product for each vector in Span. NOTE: as i assume span implements inner product and the obove * operator is overloaded, I think I can just write this one down. Output is assumed to be two m by m Eigen matrices where m is the number of columns in span.

    virtual std::tuple<Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic>, Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic>> rayleigh (Span y) = 0;
    
    
    //return suitable norms of the residuals A y(:,i) - eigenvalue(i) * y(:,i).
    
    virtual double residual(Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> eigs, Span y) = 0;
    
    // Enrich the given span 'y' by more vectors or somehow with more eigenspace content. Called when FSI restarts after a convergence failure. Must work by altering the original span passed in.
    
    virtual void augment(Span& y) = 0;
    
/* ASSUMING THAT ABOVE METHODS ARE AVAILABLE AND SPAN IS PROPERLY BUILT THIS IS THE FSI (or FEAST) ALGORITHM */
    
    
    
};

#endif /* SpectralProj_hpp */
