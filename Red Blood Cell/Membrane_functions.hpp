//
//  Membrane_functions.hpp
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 05/09/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//

#ifndef Membrane_functions_hpp
#define Membrane_functions_hpp

#include <stdio.h>
#include "General_Membrane.h"


//  Breif discription:
//      -This function reads the Gmesh generated file 'membrane' that contains the node position of both the outer membrane and nucleus            nodes. This file also contains the list of nodes that make up tringles, which we also need for our calculations.
//      -It should be also noted that we manually delete some of the lines that start with '$' in the actual Gmesh file. So the file is modified before handing it down to the programme.

//  In this function:
//      1-We set the initial Node velocities and forces to zero.
//      2-Read membrane node positions from the Gmesh file.
//      3-Read the list of nodes that make triangles.

//  Suggested improvments:
//        1-We should probably modify this function a bit so that we would not need to manually modify the Gmesh file before handing it down to the programme.
void Membrane_constructor(double Membrane_Node_Position[Membrane_num_of_Nodes][3], double Membrane_Node_Velocity [Membrane_num_of_Nodes][3], double Membrane_Node_Force[Membrane_num_of_Nodes][3], int Membrane_triangle_list [Membrane_num_of_Triangles][3]);


//  Breif discription:
//      -The triangles on the membrane and nucleus are identified.
//      -The Membrane_Normal_direction is built so that all ABxAC product of triangles on these membranes will point out of the cell.

//  In this function:
//      -It reads the position of all membrane nodes and by comparing them to the user input radii it identifies the nodes on the membrane (nucleus), and saves them in the 'Membrane_Normal_direction[i][0]'= +1 (-1) list for each triangle index 'i';
//      -Throughout the code the ABC vertexes of the membrane triangles are taken as A=Membrane_triangle_list[][0], B=Membrane_triangle_list[][1], C=Membrane_triangle_list[][2]. Also we often use the ABxAC cross product and we want the triangles on the membrane to point out of the cell. At the beginning of the cell construction, the centre of the cell is the same as the origin. So for the ABxAC product to point outwards, the inner product of the position of the triangle and the ABxAC should be positive. We put +/- 1 in the 'Membrane_Normal_direction[][1]' list for each triangle and define the normal direction of each triangle as Membrane_Normal_direction[i][1]*ABxAC that will always be positive, hence pointing out of the cell.

//  Suggested improvments:
//        1-I do not think that we need to set a list of +/- to use as a multiplier for ABxAC throughout the code for us to have ABxAC pointing out of the cell. Similar to the method I used for the ECM surface triangles, we simply need to swap B <-> C whenever the the inner product of the position of the triangle and the ABxAC vector is negative. This will also releave us of many calculations, confusion, and potential bugs in the code.

void Membrane_Normal_direction_Identifier( double  Membrane_Node_Position [Membrane_num_of_Nodes][3], int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2], int  &Outer_Membrane_num_of_triangles, int &Nucleus_Membrane_num_of_triangles, bool cell_has_nucleus);


//  Breif discription:
//      -Sorts the 'Membrane_triangle_list' by putting the triangles belonging to the outer membrane at the beginning of the list followed by the nucleus triangles.
//      -Then it goes through the 'Membrane_triangle_list' and counts the number of nodes on it.

//  In this function:
//      -First off it creats a temporary copy of the 'Membrane_triangle_list'.
//      -Then by using the Membrane_Normal_direction[i][0] identifier, it locates the triangles on the outer membrane.
//      -After which by combing through the triangles on the outer membrane, we count the nodes.
//  Suggested improvments:
//        1-Nothing particular.

void Outer_Membrane_Identifier(int Membrane_Normal_direction[Membrane_num_of_Triangles][2], int Membrane_triangle_list[Membrane_num_of_Triangles][3], int  Outer_Membrane_num_of_triangles, int  &Outer_Membrane_num_of_Nodes);


//  Breif discription:
//      -The Nodes on the outer Membrane and Nucleus are sorted in the 'Outer_Membrane_list_of_Nodes' and 'Nucleus_Membrane_list_of_Nodes' lists respectfuly. Remember that the triangle list of the membrane was sorted to contain the triangles of the membrane at the beginning og the array in the 'Outer_Membrane_Identifier' function.
//  Suggested improvments:
//        1-Nothing particular.

void Membrane_and_Nucleus_Node_list_builder(double Membrane_Node_Position [Membrane_num_of_Nodes][3],int Nucleus_Membrane_list_of_Nodes[],int Outer_Membrane_list_of_Nodes[], int Membrane_triangle_list[Membrane_num_of_Triangles][3], int  Outer_Membrane_num_of_triangles);


//  Breif discription:
//      -Counts the number of neighbouring triangle pairs on the membrane.
//  Suggested improvments:
//        1-If we actually don't need the sorting procedure in the next function ('Membrane_Triangle_Pair_Identifier') we can altogether remove this function and integrate it into the next function.

int Membrane_triangle_pair_counter( int Membrane_triangle_list[Membrane_num_of_Triangles][3]);

//  Breif discription:
//      -Neighbouring triangle pairs on the membrane are identified and written to the 'Membrane_Triangle_Pair_Nodes' list.
//      -The 'Membrane_Triangle_Pair_Nodes' is sorted (ascending)[don't know the purpose of this sorting]
//  Suggested improvments:
//        1-Nothing particular.

void Membrane_Triangle_Pair_Identifier(int Membrane_reiangle_list[Membrane_num_of_Triangles][3], int Membrane_Triangle_Pair_Nodes[][4], int Membrane_num_of_Triangle_Pairs);

//  Breif discription:
//        1-Counts the 'Outer_Membrane_num_of_Node_Pairs'.
//        2-Counts the total number of node pairs on the membrane.

//  In this function:
//      //      -creats a large temporary list of node pairs, 'bondslist'
//      -Then it goes through the 'Membrane_triangle_list' and counts the number of nodes on it.
//      -If the node pair was identified on triangles belonging to the outer membrane, we add a +1 to the  'Outer_Membrane_num_of_Node_Pairs'.
//  Suggested improvments:
//        1-Well we can modify the 'Membrane_Node_Pair_list' to a vector so that we can merge this function with 'Membrane_num_of_Node_Pair_Counter_2'. But it will take too much of my time at the moment.
int Membrane_num_of_Node_Pair_Counter(int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Outer_Membrane_num_of_triangles, int &Outer_Membrane_num_of_Node_Pairs);

//  Breif discription:
//        1-Builds the 'Membrane_Node_Pair_list'.

//  Suggested improvments:
//        1-Well we can modify the 'Membrane_Node_Pair_list' to a vector so that we can merge this function with 'Membrane_num_of_Node_Pair_Counter'. But it will take too much of my time at the moment.

void Membrane_num_of_Node_Pair_Counter_2(int Membrane_Node_Pair_list[][2], int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Outer_Membrane_num_of_triangles, int Membrane_num_of_Node_Pairs);





#endif /* Membrane_functions_hpp */
