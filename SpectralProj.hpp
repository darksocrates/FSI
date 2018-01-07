//
//  SpectralProj.hpp
//  FSI
//
//  Created by Aidan Hamilton on 1/1/18.
//
//


/* Creates an abstract class for defining the spectral projector used in filtered subspace iterations (FSI). Defines the methods to be implemented by a derived class and explictly declares how FSI (or FEAST) will operate.
 
 Three Classes are defined
 
 Span - a span of abstract vectors, contains implementation of inner product, orthagonaliation, removing elements. An abstract class.
 
 Filter - a class that contains filters for performing FSI within its methods. Seperated from Spectralproj to allow for more reusability and ease of 'growth' and readability. It is nice to have all the filter information contained within one place. This class does not require creating a derived class that overwrites abstract functions contained within.
 
 Spectralproj - Implementation of FSI (FEAST), depends on Span and Filter. An abstract class requiring a derived class to operate */



/* Dependencies
 Requires at least C++11 compilation, because of tuple usage.
 
 These non-standard libraries are required for the implementation to work.
 
 - Eigen (An opensource matrix library with optimized dense and sparse matrix operations/solvers)
 
 
 */

#ifndef SpectralProj_hpp
#define SpectralProj_hpp

#include <stdio.h>
#include <Eigen/Core>
#include <string>
#include <iostream>
#include <complex>
#include <cmath>
#include <tuple>
#include <Eigen/Eigenvalues>


const static double pi =3.14159265358979323846;
const static std::complex<double> onei(0,1);

typedef Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> MatrixXcd;
typedef Eigen::Matrix<std::complex<double>,Eigen::Dynamic, 1> VectorXcd;

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
    
    virtual void operator * (MatrixXcd) = 0;
    /* implement the actions of the span object times an Eigen matrix of complex doubles */

    //NOTE couldn't I just have this function have a return type of a Span object, therefore no assumptions are made on how the underlying span is represented? YEAH DUMBASS YOU CAN!!!! - LATE NIGHT AIDAN

    template<typename Derived>
    std::complex<double> innerproduct(const Eigen::MatrixBase<Derived>& u, const Eigen::MatrixBase<Derived>& v);
    /* required implementation of inner product for two members in the underlying space of the Span, output is assumed to be complex. This function assumes that they underlying space of the span is represented by an Matrix (Vector) from the Eigen library.*/
    
};



//implements the filters for filtered subspace iteration in its methods. This is seperated from the class Spectralproj for readability

//all methods must write to the internal complex weights, complex translations and number of points in the filter (n) variables.

class Filter {
private:
    unsigned int n = 0;
    //NOTE: see if you can figure out a way to change this to a static Matrix. Would be preferable for efficiency
    VectorXcd weights,translation;
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
    VectorXcd getweights(){return weights;}
    VectorXcd gettranslation(){return translation;}
    
    //Constructors
    Filter(unsigned int n);
    
    Filter(unsigned int n, std::string filter);
};


class Spectralproj {
private:
    
public:
    //fields
    bool verbose = 0; Filter filter;
    
    //quality of life methods
    void prnt(std::string str); 
    
    
    //the following methods must be implemented by derived classes
    
    //mult corrresponds with the action of the approximate spectral projection operator on a pointer to a Span. Alters the given span object.
    
    virtual void mult (Span *y) = 0;
    
    //compute (Aq,q) and (q,q) in the H-inner product for each vector in Span. NOTE: as i assume span implements inner product and the above * operator is overloaded, I think I can just write this one down. Output is assumed to be two m by m Eigen matrices where m is the number of columns in span.

    virtual std::tuple<MatrixXcd, MatrixXcd> rayleigh (Span *y) = 0;
    
    
    //return suitable norms of the residuals A y(:,i) - eigenvalue(i) * y(:,i).
    
    virtual double residual(MatrixXcd eigs, Span *y) = 0;
    
    // Enrich the given span 'y' by more vectors or somehow with more eigenspace content. Called when FSI restarts after a convergence failure. Must work by altering the original span passed in.
    
    virtual void augment(Span *y) = 0;
    
/* ASSUMING THAT ABOVE METHODS ARE AVAILABLE AND SPAN IS PROPERLY BUILT THIS IS THE FSI (or FEAST) ALGORITHM */
    
   VectorXcd feast_step_hermitian(Span *q);
    /* Perform one step of the FSI algorithm  assuming that the operator is self-adjoint. Takes a pointer to a span as an input and returns a  Eigen vector than contains the ritzvalues and edits the given Span object to be the span after the spectral projector has been applied */
    
    
    
    
};

#endif /* SpectralProj_hpp */
