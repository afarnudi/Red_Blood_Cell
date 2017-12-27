//
//  Chromatin_functions.hpp
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 11/10/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//

#ifndef Chromatin_functions_hpp
#define Chromatin_functions_hpp

#include <stdio.h>
#include "General_Chromatin.h"
#include <fstream>

//  Breif discription:
//      -This function reads the 'xvchromatin.txt' file that should be generated prior to the programme initiation. The file contains the position of the chromatin beads. This file also contains the list of chromatins that make up chains, which we also need for our calculations.

//  In this function:
//      1-We read the chromatin bead initial position and velocity from a file.
//  Suggested improvments:
//        1-Nope.
void Chromatin_constructor(double Chromatin_Bead_Position[Chromatin_num_of_Beads][3], double Chromatin_Bead_Velocity[Chromatin_num_of_Beads][3], double Chromatin_Bead_Force[Chromatin_num_of_Beads][3],  bool chromatin_contant_matrix_calculation_flag, float CmLaststeps[Chromatin_num_of_Beads][Chromatin_num_of_Beads]);



#endif /* Chromatin_functions_hpp */

