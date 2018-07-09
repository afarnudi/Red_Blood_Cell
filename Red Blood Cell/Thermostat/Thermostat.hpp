//
//  Thermostat.hpp
//  Red Blood Cell
//
//  Created by Ali Farnudi on 08/07/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#ifndef Thermostat_hpp
#define Thermostat_hpp

#include <stdio.h>
#include "General_Membrane.h"
#include "General_functions.hpp"
#include "General_Actin.h"
#include <math.h>

void Thermostat_2(double Membrane_Node_Velocity[][3], double Actin_Node_Velocity[][3], int Membrane_num_of_Nodes);
void Thermostat_N6(double Membrane_Node_Position[][3], double Actin_Node_Position[Actin_num_of_Nodes][3], double Membrane_Node_Velocity[][3], double Actin_Node_Velocity[][3], int Membrane_num_of_Nodes);

void two_vector_COM_calculator(double vec_COM[3], double vec_1[][3], double vec_2[][3], double vec_1_mass, double vec_2_mass, int vec_1_size, int vec_2_size);

void vector_COM_calculator(double vec_COM[3], double vec[][3], int vec_size);

void COM_omega_calculator(double Omega_com[3], double Membrane_Node_Velocity[][3], double Actin_Node_Velocity[][3], int Membrane_num_of_Nodes,  double Membrane_Node_Position[][3], double Actin_Node_Position[Actin_num_of_Nodes][3], double Position_COM[3], double V_com[3]);

void alpha_regulator(double Omega_com[3], double Membrane_Node_Velocity[][3], double Actin_Node_Velocity[][3], int Membrane_num_of_Nodes,  double Membrane_Node_Position[][3], double Actin_Node_Position[Actin_num_of_Nodes][3], double Position_COM[3], double V_com[3], double alpha);

#endif /* Thermostat_hpp */
