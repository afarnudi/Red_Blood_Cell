//
//  Actin_functions.hpp
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 27/09/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//

#ifndef Actin_functions_hpp
#define Actin_functions_hpp

#include "General_Actin.h"
#include "General_Membrane.h"
#include "General_functions.hpp"
#include <stdio.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <math.h>

//  Breif discription:
//      -This function reads the Gmesh generated file 'actin' that contains the node position of the actin network (triangular pyramids).
//      -It should be also noted that we manually delete some of the lines that start with '$' in the actual Gmesh file. So the file is modified before handing it down to the programme.

//  In this function:
//      1-We read 4 nodes (trianguar pyramid)
//      2-check if we have already read them from another pyramid or not
//      3-build the list.
//      4-return the size of the list.

//  Suggested improvments:
//        1-Nope
int  Actin_Node_Pair_Identifier( );  // gets input file and the output is network


//  Breif discription:
//      -This function reads the Gmesh generated file 'actin' that contains the node position of the actin network (triangular pyramids).
//      -It should be also noted that we manually delete some of the lines that start with '$' in the actual Gmesh file. So the file is modified before handing it down to the programme.

//  In this function:
//      1-The 'Actin_Node_Velocity', 'Actin_Node_Force',  and the 'Actin_Node_Position' lists are initialised.
//      2-4 nodes (trianguar pyramid) are read from the Gmesh file and the 'Actin_Node_Pair_List' is constructed
//      3-At the very end we set theActin_Node_Pair_List[][2] element for all Actin node pairs. This value is the distance between the actin_node pairs. We use this element to store a 'one step' history for the actin pair distances where we use to calculate the viscus force betweeen the actin node pairs in the dash-pot model.

//  Suggested improvments:
//        1-Nope
void Actin_constructor(double Actin_Node_Position[][3], double Actin_Node_Velocity[][3], double Actin_Node_Force [][3], double Actin_Node_Pair_List[][3], int Actin_num_of_Bonds, double Actin_Node_Pair_List_2[][5]);

#endif /* Actin_functions_hpp */
