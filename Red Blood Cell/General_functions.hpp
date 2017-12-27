//
//  General_functions.hpp
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 27/08/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//

#ifndef General_functions_hpp
#define General_functions_hpp

void crossvector( double c[3],double d[3],double b[3] ); // cross porduct
double innerproduct(double n1[3],double n2[3]);
double vectorlength(double v[3]); // calculate length of vector
double sign_function(double x);


#include <stdio.h>

#endif /* General_functions_hpp */
