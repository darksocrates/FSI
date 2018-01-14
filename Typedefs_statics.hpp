//
//  Typedefs_statics.hpp
//  FSI_2
//
//  Created by Aidan Hamilton on 1/13/18.
//  Copyright © 2018 Aidan Hamilton. All rights reserved.
//

#ifndef Typedefs_statics_hpp
#define Typedefs_statics_hpp

#include <stdio.h>
#include <complex>
#include<Eigen/Core>

const static double pi =3.14159265358979323846;
const static std::complex<double> onei(0,1);

typedef Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic> MatrixXcd;
typedef Eigen::Matrix<std::complex<double>,Eigen::Dynamic, 1> VectorXcd;


#endif /* Typedefs_statics_hpp */
