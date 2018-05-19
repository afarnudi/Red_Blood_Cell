//
//  General_Membrane.h
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 05/09/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//

#ifndef General_Membrane_h
#define General_Membrane_h

//#define Membrane_num_of_Nodes 1678//1648 // //This is the number of nodes on the membrane (Both the outer membrane and the Nucleus). This is the first number that appears in the 'membrane' file (once opend with a text editor)
//#define Membrane_num_of_Triangles_temp 3488//3288//  //This is the number of triangles on the membrane (Both the outer membrane and the Nucleus). This is the number that appears in the 'membrane' file after the node position list is finished and before Gmesh lists the nodes that make a triangle.
#define Membrane_Radius 10.0 //Is the outer membrane radius. The radius of the membrane should be the same as the one used to generate the Gmesh membrane file otherwise the programme crashes.
//#define Nucleus_Membrane_radius   4.0  //Is the nucleus membrane radius. The radius of the nucleus should be the same as the one used to generate the Gmesh membrane file otherwise the programme crashes.
#define kvolume    0.0

//#define Membrane_spring_coefficient     2.0  //1.0   // streching constant
//#define Membrane_bending_coefficient   0.0//30.0//15.0     // bending constant
#define membrane_damping_coefficient 0.0   // Viscosity of the Mmmbrane. It is applied in Force calculation for the Membrane Node pairs. I have commented out these parts in the 'Membrane_Force_Calculator' because I think the current code does not need it (some energy consuming array calculations were invloved).
#define Membrane_Node_Mass 1.0
#define K_surfaceConstant_local 100
#define Membrane_spring_coefficient     2.0  //1.0   // streching constant
#define Membrane_bending_coefficient   2.0//30.0//15.0     // bending constant
#define Membrane_spring_force_cutt_off  10000.0
#define membrane_damping_coefficient  0.0   // Viscosity of the Mmmbrane. It is applied in Force calculation for the Membrane Node pairs. I have commented out these parts in the 'Membrane_Force_Calculator' because I think the current code does not need it (some energy consuming array calculations were invloved).
#define Membrane_barrier_calculation_rate     5
#define Node_radius 1.0 // width of interacion range of nodes: use to be 1.0 !

#endif /* General_Membrane_h */

