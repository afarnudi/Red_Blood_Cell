//
//  Monte_Carlo_Bond_Flip.hpp
//  Test_2
//
//  Created by Ali Farnudi on 12/03/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#ifndef Monte_Carlo_Bond_Flip_hpp
#define Monte_Carlo_Bond_Flip_hpp

#include <stdio.h>
#include <vector>
#include <math.h>
#include "General_Membrane.h"
#include "General_functions.hpp"

using namespace std;

void Monte_carlo_bond_flip (double Membrane_Node_Position[][3], vector<vector<int> > &Membrane_new_triangle_list, vector<vector<int> > &membrane_triangle_pair_list, int Membrane_Node_Pair_list[][2], int Membrane_num_of_Node_Pairs, int Membrane_num_of_Nodes);

void trinalge_pair_node_sorter(int triangle_1, int triangle_2, vector<vector<int> > Membrane_new_triangle_list, int sorted_nodes[4]);

void new_trinalge_pair_node_sorter(int triangle_A_node_A, int triangle_A_node_B, int triangle_A_node_C, int triangle_B_node_D, int triangle_A_neighbours[2], int triangle_B_neighbours[2], vector<vector<int> > Membrane_new_triangle_list, int new_sorted_nodes[4][4], int triangle_A_new_neighbours[2], int triangle_B_new_neighbours[2]);



#endif /* Monte_Carlo_Bond_Flip_hpp */

