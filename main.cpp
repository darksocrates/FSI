//
//  main.cpp
//  FSI_2
//
//  Created by Aidan Hamilton on 1/1/18.
//  Copyright © 2018 Aidan Hamilton. All rights reserved.
//

#include <iostream>
#include "SpectralProj.hpp"

int main(int argc, const char * argv[]) {

    Filter test(12);
    test.center = std::complex<double>(1,1); test.radius = 2.0;
    
    test.circ_trapez_snug();
    std::cout<<test.gettranslation()<<"\n" <<test.getweights()<<"\n";
    
    
    return 0;
}
