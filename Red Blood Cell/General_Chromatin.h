//
//  General_Chromatin.h
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 11/10/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//

#ifndef General_Chromatin_h
#define General_Chromatin_h

#define Chromatin_num_of_Beads    500//200  // number of beads  per chromatin
#define Chromatin_num_of_chains 10   ///  Chromatin_num_of_Beads should be devidable to Chromatin_num_of_chains!
#define Chromatin_Bead_Mass    1.0 // mass of each bead
#define Chromatin_Scaling_Factor  0.7 //This is the scaling length. It (usually) shrinks the chromatines that are the result of the chromatine confinement code, 'xvchromatin.txt'.
#define Chromatin_spring_coefficient     400.0    // streching constant of chromatin
#define Chromatin_bending_coefficient 0.0
#define sigmachromatin  1.0 // L-j parameter for FENE
#define Rmaxchromatin  1.5  // FENE parameter
#define Chromatin_Membrane_Radius_of_Hard_Sphere_Interaction   0.5
#define ThermostatOnChromatin 100
#define chromatin_force_cut_off 3000







#endif /* General_Chromatin_h */
