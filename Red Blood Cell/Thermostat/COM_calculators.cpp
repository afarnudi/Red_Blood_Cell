//
//  COM_calculators.cpp
//  Red Blood Cell
//
//  Created by Ali Farnudi on 09/07/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include <stdio.h>
#include "Thermostat.hpp"

void two_vector_COM_calculator(double vec_COM[3], double vec_1[][3], double vec_2[][3], double vec_1_mass, double vec_2_mass, int vec_1_size, int vec_2_size){
    double vec_2_COM[3]={0,0,0};
    
    for (int i=0; i<3; i++) {
        vec_COM[i]=0;
    }
    
    for (int i=0; i<vec_1_size; i++) {
        vec_COM[0]+=vec_1[i][0];
        vec_COM[1]+=vec_1[i][1];
        vec_COM[2]+=vec_1[i][2];
    }
    
    vec_COM[0]*=vec_1_mass;
    vec_COM[1]*=vec_1_mass;
    vec_COM[2]*=vec_1_mass;
    
    for (int i=0; i<vec_2_size; i++) {
        vec_2_COM[0]+=vec_2[i][0];
        vec_2_COM[1]+=vec_2[i][1];
        vec_2_COM[2]+=vec_2[i][2];
    }
    
    vec_2_COM[0]*=vec_2_mass;
    vec_2_COM[1]*=vec_2_mass;
    vec_2_COM[2]*=vec_2_mass;
    
    vec_COM[0]=(vec_COM[0]+vec_2_COM[0])/(vec_1_size*vec_1_mass+vec_2_size*vec_2_mass);
    vec_COM[1]=(vec_COM[1]+vec_2_COM[1])/(vec_1_size*vec_1_mass+vec_2_size*vec_2_mass);
    vec_COM[2]=(vec_COM[2]+vec_2_COM[2])/(vec_1_size*vec_1_mass+vec_2_size*vec_2_mass);
}

void vector_COM_calculator(double vec_COM[3], double vec[][3], int vec_size){
    for (int i=0; i<3; i++) {
        vec_COM[i]=0;
    }
    
    for (int i=0; i<vec_size; i++) {
        vec_COM[0]+=vec[i][0];
        vec_COM[1]+=vec[i][1];
        vec_COM[2]+=vec[i][2];
    }
    
    vec_COM[0]/=vec_size;
    vec_COM[1]/=vec_size;
    vec_COM[2]/=vec_size;
}

void COM_omega_calculator(double Omega_com[3], double Membrane_Node_Velocity[][3], double Actin_Node_Velocity[][3], int Membrane_num_of_Nodes, double Membrane_Node_Position[][3], double Actin_Node_Position[Actin_num_of_Nodes][3], double Position_COM[3], double V_com[3]){
    double angular_momentum[3]={0,0,0}, actin_angular_momentum[3]={0,0,0}, moment_of_inertia=0, actin_moment_of_inertia=0, temp_cross_product[3], temp_position[3], temp_velocity[3];
    
    for (int i=0; i<Membrane_num_of_Nodes; i++) {
        for (int j=0; j<3; j++) {
            temp_position[j]=Membrane_Node_Position[i][j]-Position_COM[j];
            temp_velocity[j]=Membrane_Node_Velocity[i][j]-V_com[j];
        }
        
        crossvector(temp_cross_product, temp_position, temp_velocity);
        moment_of_inertia+=innerproduct(temp_position, temp_position);
        
        for (int j=0; j<3; j++) {
            angular_momentum[j]+=temp_cross_product[j];
        }
        
    }
    moment_of_inertia*=Membrane_Node_Mass;
    for (int j=0; j<3; j++) {
        angular_momentum[j]*=Membrane_Node_Mass;
    }
    
    //========================== Actin ==========================
    for (int i=0; i<Actin_num_of_Nodes; i++) {
        for (int j=0; j<3; j++) {
            temp_position[j]=Actin_Node_Position[i][j]-Position_COM[j];
            temp_velocity[j]=Actin_Node_Velocity[i][j]-V_com[j];
        }
        crossvector(temp_cross_product, temp_position, temp_velocity);
        actin_moment_of_inertia+=innerproduct(temp_position, temp_position);
        for (int j=0; j<3; j++) {
            actin_angular_momentum[j]+=temp_cross_product[j];
        }
        
    }
    moment_of_inertia+=actin_moment_of_inertia*Actin_Node_Mass;
    for (int j=0; j<3; j++) {
        angular_momentum[j]+=actin_angular_momentum[j]*Actin_Node_Mass;
        Omega_com[j]=angular_momentum[j]/moment_of_inertia;
    }
    
    
}
