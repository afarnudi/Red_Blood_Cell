//
//  ECM_functions.hpp
//  Cell_Project
//
//  Created by Ali Farnudi on 22/10/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//

#ifndef ECM_functions_hpp
#define ECM_functions_hpp

#include <stdio.h>
#include <vector>
#include <fstream>
#include <iostream>
#include "General_ECM.h"
#include "General_functions.hpp"

//  Breif discription:
//      -This function reads the Gmesh generated file 'ECM' that contains the node position of the actin network (triangular pyramids).
//      -It should be also noted that we manually delete some of the lines that start with '$' in the actual Gmesh file. So the file is modified before handing it down to the programme.

//  In this function:
//      1-We read through the file, ignorong the node positions and triangle (since we have will read them farther into the programme) until we reach the pyramid object list.
//      2-check if we have already read them from another pyramid or not
//      3-build the list.
//      4-return the size of the list.

//  Suggested improvments:
//        1-Nope
int  ECM_Node_Pair_Identifier( );

void ECM_constructor(double  ECM_Node_Position [][3],double  ECM_Node_Velocity [][3],double  ECM_Node_Force [][3],int ECM_surface_triangle_list[ECM_Surface_num_of_Triangles][3],double ECM_Node_Pair_List[][3],double ECM_upper_surface_Node_Pairs[], int ECM_num_of_Bonds);

void ECM_triangle_normal_direction_justifier(double  ECM_Node_Position [][3],int ECM_surface_triangle_list[ECM_Surface_num_of_Triangles][3]);// n=AB x AC should point out side (A,B,C=ECM_surface_triangle_list[][0,1,2])  for every triangle!



#endif /* ECM_functions_hpp */
