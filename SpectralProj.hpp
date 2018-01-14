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
#include <vector>
#include "Filter.hpp"
#include"Typedefs_statics.hpp"



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
    
    //the following two methods are used heavily in the FSI_data class
        virtual void store(std::string path) = 0;
        /*method to store the contents of the span in a user specified format at the location specified by path */
    
        virtual void load(std::string path) = 0;
        /*method to load the given span object with the contents dumped to the file given by path. Assumed that the file given by path stores the contents of a span object in the same format as teh store method would*/
    
    virtual void orthagonalize() = 0;
    /* replace the current span by an orthagonal span in a suitable H space */
    
    virtual void remove(unsigned int *cut) =0;
    /* remove vectors with indices given by the array cut */
    
    virtual void operator * (MatrixXcd) = 0;
    /* implement the actions of the span object times an Eigen matrix of complex doubles */

    virtual std::complex<double> operator *(Span *u) =0;
    /* required implementation of inner product for two members in the underlying space of the Span, output is assumed to be complex.*/
    
    virtual std::string whatareyou() = 0;
    /*Derived classes must implement this functions, which expects to return a string saying what the derived class is. This is to allow some way for a user to from a Span pointer, which can be pointing to some unknown derived class of Span, query what derived class of Span they are working with (ie an L2 span, H1 span C1 span etc etc) */
    
};




/* Class that holds the data returned by FSI, implemented as a class and not a structure to allow for the addition of methods for plotting and etc, to ease routine tasks that would be done with the data. Also lets me screw with how the data is actually stored internally without effecting any current code using this. Fields can only be set via a constructor or via the read method.*/

class FSI_data {
private:
    /* fields, all are private so that the user can't accidently screw with them. */
    
    // pointer to a span object, will be the same pointer that is inputted into the fsi algorithm
    Span *span;
    //Eigen Matrix containing the ritzvalues for each iteration indexed from top to bottom.
    std::vector<VectorXcd> eigenvalues;
    //std::vector<Eigen::VectorXd> containing the residuals (if calculated), along with boolean to track if the residual is calculated, defaults to not.
    std::vector<Eigen::VectorXd> residual; bool residual_calculated = false;
    
    //String that contains the result of the whatareyou() function the of particular Spectral projector that called this particular instance of FSI data in its fsi function.
    std::string spectralproj_type;
    
public:
    /*access methods to the private data*/
    //NOTE: add user protections to check if the various fields have actually been declared or not.
    
    //returns the eigenvalues calculated in the last iteration of the FSI as a complex Eigen vector
    VectorXcd geteigenvalues() {return eigenvalues[eigenvalues.size()-1];}
    
    std::vector<VectorXcd>  getalleigenvalues() {return eigenvalues;}
    
    //get all eigenvalues in the form of an Eigen Matrix with NAN in the spaces where nothing was calculated
    MatrixXcd getalleigenvaluesM();
    
    std::vector<Eigen::VectorXd> getresiduals(){if (residual_calculated == false){
        throw "residual was not calculated, there is no residual data here";
    } else {
        return residual;
    } }
    
    //labeled explictly as Span pointer to make clear that you are getting a pointer to a Span object
    Span* getSpan_pointer(){return span;}
    
    std::string getSpectralproj_name(){return spectralproj_type;}
    
    
    /* various quality of life methods in the usage of the class */
    
    void store(std::string path, std::string name);
    /* stores the contents of the current FSI_data object in folder specified by the path, path as a collection of two files in a folder in the path folder with name, name. The two files are named FSI_data.txt and another file with format given by the particular derived class of Span used*/
    
    void load(std::string path, std::string name);
    /* loads the contents given by the path, (path\name) to the current FSI_data object. Assumes the file structure within the given directory specified is of the format that would be created by using the store method*/
    
    
    /* Constructors*/
    
    FSI_data(std::vector<VectorXcd>  eigenvalues, Span *span, std::string spectralproj_type);
    
    FSI_data(std::vector<VectorXcd> eigenvalues, Span *span, std::string spectralproj_type, std::vector<Eigen::VectorXd> residual);
    
    
    
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
    
    virtual Eigen::VectorXd residual(VectorXcd eigs, const Span *y) = 0;
    
    // Enrich the given span 'y' by more vectors or somehow with more eigenspace content. Called when FSI restarts after a convergence failure. Must work by altering the original span passed in.
    
    virtual void augment(Span *y) = 0;
    
/* ASSUMING THAT ABOVE METHODS ARE AVAILABLE AND SPAN IS PROPERLY BUILT THIS IS THE FSI (or FEAST) ALGORITHM */
    
VectorXcd fsi_step_hermitian(Span *q);
    /* Perform one step of the FSI algorithm  assuming that the operator is self-adjoint. Takes a pointer to a span as an input and returns a  Eigen vector than contains the ritzvalues and edits the given Span object to be the span after the spectral projector has been applied */
    
    //NOTE: currently implemented without restart capability
    FSI_data fsi_hermitian(Span *q, double stop_tol = 1e-9,double eigen_fudge = 1e-9, unsigned int maxiterations = 10, bool report_residuals = false);
    
};


#endif /* SpectralProj_hpp */
