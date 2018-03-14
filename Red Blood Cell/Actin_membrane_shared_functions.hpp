//
//  Actin_membrane_shared_functions.hpp
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 11/10/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//

#ifndef Actin_membrane_shared_functions_hpp
#define Actin_membrane_shared_functions_hpp

#include <stdio.h>
#include <iostream>
#include "General_Membrane.h"
#include "General_Actin_Membrane_shared.h"
#include "General_Actin.h"

//  Breif discription:
//      -A very simple programme that identifies the Actin nodes on the Membrane

//  In this function:
//      1-Since they are built from the same Gmesh file, The Actin Nodes on the Membrane are located at the exact coordinates as the Membrane nodes they are supposed to be connected with (The connection is established via a spring force in another function)

//  Suggested improvments:
//        1-Nope
void Membrane_Actin_shared_Node_Identifier(int Membrane_Actin_shared_Node_list[Actin_Membrane_shared_num_of_Nodes][2],double Membrane_Node_Position[Membrane_num_of_Nodes][3],double  Actin_Node_Position [][3]);


#endif /* Actin_membrane_shared_functions_hpp */

