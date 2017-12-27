//
//  General_Actin.h
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 27/09/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//

#ifndef General_Actin_h
#define General_Actin_h

#define Actin_num_of_Nodes 587//566// //This is the number of nodes in the actin volume. It is the first number that appears in the 'actin' file.
#define Actin_spring_coefficient      5.0//10.0//5.0//0.3 // force constant of springs
#define Actin_kelvin_damping_coefficient        1000.0 // viscouse constant for avtin Maxwell fluids: maxwell model >> IT RESISTS AGAINST CHANGING BOND ENGTHS
#define CytoskeletonNetworkType  0 // 0=normal hookian network    1=passive cable network  2=active cable network
//incase:  CytoskeletonNetworkType=1  - - - - - - -
#define Actin_Passive_Cable_Network_Coefficient  2.0   //more detail in the paper: Contractile network models for adherent cells

//incase:  CytoskeletonNetworkType=2  - - - - - - -
#define KActinACN_EA  2.0
#define ACN_TL0       0.50
#define ACN_LC        0.1



#define  Actin_damping_Coefficient    0.1   // viscouse constant for actin Maxwell fluids: voijdt model >> IT RESISTS AGAINST MOVING  BEADS
#define  Actin_Node_Mass 1.0//0.30
#define minimumlengthActin 1.5
#define maximumlengthActin 2.5


#endif /* General_Actin_h */
