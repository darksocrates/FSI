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
protected:
    
    //n is the number of 'columns' in the Span
 unsigned int n;
    
private:
//err is a string that holds the standard error message to display for functions that are optional to override in a derived class (they are not pure virtual functions), but do not have any meaning here.
    
    std::string err = "Function not implemented in derived class, base Span function called";

    
    
public:
    
    bool verbose =0;
//The following functions have a sensible definition so do not require overriding in a derived class.
    
     void prnt(std::string str);
    
     unsigned int numcolumns(){return n;}
    
//The following functions must be implemented by derived classes
    
    virtual void orthagonalize() = 0;
    /* replace the current span by an orthagonal span in a suitable H space */
    
    virtual void remove(unsigned int *cut[]) =0;
    /* remove vectors with indices given by the array cut*/
    
    virtual void operator * (MatrixXcd) = 0;
    /* implement the actions of the span object times an Eigen matrix of complex doubles */

    virtual std::complex<double> operator *(Span *u) =0;
    /* required implementation of inner product for two members in the underlying space of the Span, output is assumed to be complex.*/
    
    virtual std::string whatareyou() = 0;
    /*Derived classes must implement this functions, which expects to return a string saying what the derived class is. This is to allow some way for a user to from a Span pointer, which can be pointing to some unknown derived class of Span, query what derived class of Span they are working with (ie an L2 span, H1 span C1 span etc etc) */
    
};



//implements the filters for filtered subspace iteration in its methods. This is seperated from the class Spectralproj for readability
// This class will also implement all of the 'discretization' of the eigenvalues search set. Ie. it will decide how many seperate fsi algorithms to run


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
    void setn(unsigned int n); 
    
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


/* Class that holds the data returned by FSI, implemented as a class and not a structure to allow for the addition of methods for plotting and etc, to ease routine tasks that would be done with the data. Also lets me screw with how the data is actually stored internally without effecting any current code using this. Fields can only be set via a constructor or via the read method.*/

class FSI_data {
private:
    /* fields, all are private so that the user can't accidently screw with them. */
    
    // pointer to a span object, will be the same pointer that is inputted into the fsi algorithm
    Span *span;
    //Eigen Matrix containing the ritzvalues for each iteration indexed from top to bottom.
    MatrixXcd eigenvalues;
    //Eigen Matrix containing the residuals (if calculated), along with boolean to track if the residual is calculated, defaults to not.
    Eigen::MatrixXd residual; bool residual_calculated = false;
    
    //String that contains the result of the whatareyou() function the of particular Spectral projector that called this particular instance of FSI data in its fsi function.
    std::string spectralproj_type;
    
public:
    /*access methods to the private data*/
    //NOTE: add user protections to check if the various fields have actually been declared or not.
    
    MatrixXcd geteigenvalues() {return eigenvalues;}
    
    Eigen::MatrixXd getresidual(){if (residual_calculated == false){
        throw "residual was not calculated, there is no residual data here";
    } else {
        return residual;
    } }
    
    //labeled explictly as Span pointer to make clear that you are getting a pointer to a Span object
    Span* getSpan_pointer(){return span;}
    
    std::string getSpectralproj_name(){return spectralproj_type;}
    
    /* Constructors*/
    
    FSI_data(MatrixXcd eigenvalues, Span *span, std::string spectralproj_type);
    
    FSI_data(MatrixXcd eigenvalues, Span *span, std::string spectralproj_type, Eigen::MatrixXd residual);
    
    
    
    //Quality of life methods to implement!!!!!
    // -  plot function(s), generalized one and short hands for the common requests.
    //though with the below implemented you could easily just import the saved data into matlab and have a set script created for plotting, as a tempory measure
    
    //NOTE: to implement, save function, writes the data to a text file in a sensible labeled format.
    //NOTE: to implement, read function, reads the data in the above sensible format to and instance of this class.
    
    
};


class Spectralproj {
protected:
    
public:
    //fields
    bool verbose = 0; Filter *filter;
    
    //quality of life methods
    void prnt(std::string str);
    
    
    //the following methods must be implemented by derived classes
    
    
    virtual std::string whatareyou() = 0;
    /*Method to on query return string saying what type of spectral projector you are, used to allow more polymorphism behavior in applications*/
    
    //mult corrresponds with the action of the approximate spectral projection operator on a pointer to a Span. Alters the given span object.
    
    virtual void mult (Span *y) = 0;
    
    //compute (Aq,q) and (q,q) in the H-inner product for each vector in Span. NOTE: as i assume span implements inner product and the above * operator is overloaded, I think I can just write this one down. Output is assumed to be two m by m Eigen matrices where m is the number of columns in span.

    virtual std::tuple<MatrixXcd, MatrixXcd> rayleigh (Span *y) = 0;
    
    
    //return suitable norms of the residuals A y(:,i) - eigenvalue(i) * y(:,i).
    
    virtual double residual(MatrixXcd eigs, Span *y) = 0;
    
    // Enrich the given span 'y' by more vectors or somehow with more eigenspace content. Called when FSI restarts after a convergence failure. Must work by altering the original span passed in.
    
    virtual void augment(Span *y) = 0;
    
/* ASSUMING THAT ABOVE METHODS ARE AVAILABLE AND SPAN IS PROPERLY BUILT THIS IS THE FSI (or FEAST) ALGORITHM */
    
VectorXcd fsi_step_hermitian(Span *q);
    /* Perform one step of the FSI algorithm  assuming that the operator is self-adjoint. Takes a pointer to a span as an input and returns a  Eigen vector than contains the ritzvalues and edits the given Span object to be the span after the spectral projector has been applied */
    
    //NOTE: currently implemented with restart
    FSI_data fsi_hermitian(Span *q, double stop_tol = 1e-9,double eigen_fudge = 1e-9, unsigned int maxiterations = 10, bool report_residuals = false);
    
};


#endif /* SpectralProj_hpp */
