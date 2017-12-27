//
//  General_ECM.h
//  Cell_Project
//
//  Created by Ali Farnudi on 22/10/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//

#ifndef General_ECM_h
#define General_ECM_h

#define ECM_num_of_Nodes 4079//20458//1595//5743//4079//1859//1595 //Number of Nodes in the ECM Volume. It is the first number in the 'ECM' file (manually).
#define ECM_Surface_num_of_Triangles 1052//4252//520//6548//1052//584//520 //Number of triangles on the ECM. It is read from the 'ECM' file (manually), and it is the index of the last 7 member number series.
#define ECM_Thickness  10.0//15.0//20.0 //The thickness of the ECM (throughout this code the y axis will difine the Vertical distances)
#define Membrane_Centre_distance_from_ECM 11.30  // distance between extra cellular matrix and center of the cell
#define ECM_x_scaling_factor   1.7
#define ECM_y_scaling_factor   1.7
#define ECM_z_scaling_factor   1.7

#endif /* General_ECM_h */
