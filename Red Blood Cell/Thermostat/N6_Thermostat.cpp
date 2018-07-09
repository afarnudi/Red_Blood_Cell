//
//  N6_Thermostat.cpp
//  Red Blood Cell
//
//  Created by Ali Farnudi on 09/07/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include <stdio.h>
#include "Thermostat.hpp"

void Thermostat_N6(double Membrane_Node_Position[][3], double Actin_Node_Position[Actin_num_of_Nodes][3], double Membrane_Node_Velocity[][3], double Actin_Node_Velocity[][3], int Membrane_num_of_Nodes)
{
    double Position_COM[3]={0,0,0}, Membrane_Position_COM[3]={0,0,0}, Actin_Position_COM[3]={0,0,0}, Velocity_COM[3]={0,0,0}, Membrane_Velocity_COM[3]={0,0,0}, Actin_Velocity_COM[3]={0,0,0}, Omega_com[3]={0,0,0};
    
    vector_COM_calculator(Membrane_Velocity_COM, Membrane_Node_Velocity, Membrane_num_of_Nodes);
    vector_COM_calculator(Actin_Velocity_COM, Actin_Node_Velocity, Actin_num_of_Nodes);
    
    for (int i=0; i<3; i++) {
        Velocity_COM[i]=(Membrane_num_of_Nodes*Membrane_Node_Mass*Membrane_Velocity_COM[i]+Actin_num_of_Nodes*Actin_Node_Mass*Actin_Velocity_COM[i])/(Membrane_num_of_Nodes*Membrane_Node_Mass+Actin_num_of_Nodes*Actin_Node_Mass);
    }
    
    vector_COM_calculator(Membrane_Position_COM, Membrane_Node_Position, Membrane_num_of_Nodes);
    vector_COM_calculator(Actin_Position_COM, Actin_Node_Position, Actin_num_of_Nodes);
    
    for (int i=0; i<3; i++) {
        Position_COM[i]=(Membrane_num_of_Nodes*Membrane_Node_Mass*Membrane_Position_COM[i]+Actin_num_of_Nodes*Actin_Node_Mass*Actin_Position_COM[i])/(Membrane_num_of_Nodes*Membrane_Node_Mass+Actin_num_of_Nodes*Actin_Node_Mass);
    }
    
    //    two_vector_COM_calculator(Velocity_COM, Membrane_Node_Velocity, Actin_Node_Velocity, Membrane_Node_Mass, Actin_Node_Mass, Membrane_num_of_Nodes, Actin_num_of_Nodes);
    
    //    two_vector_COM_calculator(Position_COM, Membrane_Node_Position, Actin_Node_Position, Membrane_Node_Mass, Actin_Node_Mass, Membrane_num_of_Nodes, Actin_num_of_Nodes);
    
    COM_omega_calculator(Omega_com, Membrane_Node_Velocity, Actin_Node_Velocity, Membrane_num_of_Nodes, Membrane_Node_Position, Actin_Node_Position, Position_COM, Velocity_COM);
    //    COM_velocity_calculator(V_com, Membrane_Node_Velocity, Actin_Node_Velocity, Membrane_num_of_Nodes);
    
    
    double alpha;
    //----------------------membrane---------------------
    double Kinetic_energy=0, temp_Kinetic_energy=0, temp_cross_product[3]={0,0,0}, temp_vector[3]={0,0,0};
    for (int i=0; i<Membrane_num_of_Nodes; i++) {
        
        for (int j=0; j<3; j++) {
            temp_vector[j]=Membrane_Node_Position[i][j]-Position_COM[j];
        }
        crossvector(temp_cross_product, Omega_com, temp_vector);
        for (int j=0; j<3; j++) {
            temp_vector[j]=Membrane_Node_Velocity[i][j]-Velocity_COM[j]-temp_cross_product[j];
        }
        Kinetic_energy+=vectorlength_squared(temp_vector);
    }
    Kinetic_energy*=Membrane_Node_Mass;
    
    for (int i=0; i<Actin_num_of_Nodes; i++) {
        for (int j=0; j<3; j++) {
            temp_vector[j]=Actin_Node_Position[i][j]-Position_COM[j];
        }
        crossvector(temp_cross_product, Omega_com, temp_vector);
        for (int j=0; j<3; j++) {
            temp_vector[j]=Actin_Node_Velocity[i][j]-Velocity_COM[j]-temp_cross_product[j];
        }
        temp_Kinetic_energy+=vectorlength_squared(temp_vector);
    }
    temp_Kinetic_energy*=Actin_Node_Mass;
    Kinetic_energy+=temp_Kinetic_energy;
    
    alpha=sqrt(6*(Membrane_num_of_Nodes+Actin_num_of_Nodes)*KT)/Kinetic_energy;
    
    alpha_regulator(Omega_com, Membrane_Node_Velocity, Actin_Node_Velocity, Membrane_num_of_Nodes, Membrane_Node_Position, Actin_Node_Position, Position_COM, Velocity_COM, alpha);
}

void alpha_regulator(double Omega_com[3], double Membrane_Node_Velocity[][3], double Actin_Node_Velocity[][3], int Membrane_num_of_Nodes, double Membrane_Node_Position[][3], double Actin_Node_Position[Actin_num_of_Nodes][3], double Position_COM[3], double V_com[3], double alpha){
    double temp_cross_product[3], temp_position[3], temp_velocity[3], temp_constant;
    
    for (int i=0; i<Membrane_num_of_Nodes; i++) {
        temp_position[0]=Membrane_Node_Position[i][0]-Position_COM[0];
        temp_position[1]=Membrane_Node_Position[i][1]-Position_COM[1];
        temp_position[2]=Membrane_Node_Position[i][2]-Position_COM[2];
        crossvector(temp_cross_product, Omega_com, temp_position);
        for (int j=0; j<3; j++) {
            Membrane_Node_Velocity[i][j]+=-V_com[j]-temp_cross_product[j];
            Membrane_Node_Velocity[i][j]*=alpha;
            Membrane_Node_Velocity[i][j]+=+V_com[j]+temp_cross_product[j];
        }
        
    }
    
    //========================== Actin ==========================
    for (int i=0; i<Actin_num_of_Nodes; i++) {
        temp_position[0]=Actin_Node_Position[i][0]-Position_COM[0];
        temp_position[1]=Actin_Node_Position[i][1]-Position_COM[1];
        temp_position[2]=Actin_Node_Position[i][2]-Position_COM[2];
        crossvector(temp_cross_product, Omega_com, temp_position);
        for (int j=0; j<3; j++) {
            Actin_Node_Velocity[i][j]+=-V_com[j]-temp_cross_product[j];
            Actin_Node_Velocity[i][j]*=alpha;
            Actin_Node_Velocity[i][j]+=+V_com[j]+temp_cross_product[j];
        }
    }
}
