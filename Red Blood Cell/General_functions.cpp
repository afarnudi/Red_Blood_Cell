//
//  General_functions.cpp
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 27/08/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//

#include "General_functions.hpp"
#include <math.h>
void crossvector( double c[3],double d[3],double b[3] ) // cross porduct
{
    c[0]=d[1]*b[2]-d[2]*b[1];    // normal vector to plane (not unitary length)
    c[1]=d[2]*b[0]-d[0]*b[2];
    c[2]=d[0]*b[1]-d[1]*b[0];
}

double innerproduct(double n1[3],double n2[3])
{
    return n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];
}

double vectorlength(double v[3]) // calculate length of vector
{
    return  sqrt( v[0]*v[0]+v[1]*v[1]+v[2]*v[2] );
}
double sign_function(double x)
{
    if(x<0){
        return -1.0;
    }
    
    return 1.0;
}
