//
//  Thermostat.cpp
//  Red Blood Cell
//
//  Created by Ali Farnudi on 08/07/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Thermostat.hpp"

void Thermostat_2(double Membrane_Node_Velocity[][3], double Actin_Node_Velocity[][3], int Membrane_num_of_Nodes)
{
    double V_com[3], temp_V[3];
    
    temp_V[0]=0;
    temp_V[1]=0;
    temp_V[2]=0;
    
    V_com[0]=0;
    V_com[1]=0;
    V_com[2]=0;
    for (int i=0; i<Membrane_num_of_Nodes; i++) {
        V_com[0]+=Membrane_Node_Velocity[i][0];
        V_com[1]+=Membrane_Node_Velocity[i][1];
        V_com[2]+=Membrane_Node_Velocity[i][2];
    }
    
    V_com[0]*=Membrane_Node_Mass;
    V_com[1]*=Membrane_Node_Mass;
    V_com[2]*=Membrane_Node_Mass;
    
    for (int i=0; i<Actin_num_of_Nodes; i++) {
        temp_V[0]+=Actin_Node_Velocity[i][0];
        temp_V[1]+=Actin_Node_Velocity[i][1];
        temp_V[2]+=Actin_Node_Velocity[i][2];
    }
    
    temp_V[0]*=Actin_Node_Mass;
    temp_V[1]*=Actin_Node_Mass;
    temp_V[2]*=Actin_Node_Mass;
    
    V_com[0]=(V_com[0]+temp_V[0])/(Membrane_num_of_Nodes*Membrane_Node_Mass+Actin_num_of_Nodes*Actin_Node_Mass);
    V_com[1]=(V_com[1]+temp_V[1])/(Membrane_num_of_Nodes*Membrane_Node_Mass+Actin_num_of_Nodes*Actin_Node_Mass);
    V_com[2]=(V_com[2]+temp_V[2])/(Membrane_num_of_Nodes*Membrane_Node_Mass+Actin_num_of_Nodes*Actin_Node_Mass);
    
    //    cout<<"COM before= "<<sqrt(V_com[0]*V_com[0]+V_com[1]*V_com[1]+V_com[2]*V_com[2])<<"\nV_X="<<V_com[0]<<"\tV_Y="<<V_com[1]<<"\tV_Z="<<V_com[2]<<endl;
    
    
    double alpha;
    //----------------------membrane---------------------
    
    //    alpha=sqrt(  (3*Membrane_num_of_Nodes*KT) / kineticenergymembrane( Membrane_Node_Velocity )         );/// NOTE THAT THERMOSTATE IS FOR MEMBRANE YET. IN ORDER TO
    /// UPDATE IT FOR BOTH THE MEMBRANE AND NUCLEI WE HAVE TO
    /// WRITE  alpha=sqrt(      (2*3*Membrane_num_of_Nodes*KT) / kineticenergy( Membrane_Node_Velocity,vnuclei
    double Kinetic_energy=0, temp_Kinetic_energy=0;
    for (int i=0; i<Membrane_num_of_Nodes; i++) {
        Membrane_Node_Velocity[i][0]-=V_com[0];
        Membrane_Node_Velocity[i][1]-=V_com[1];
        Membrane_Node_Velocity[i][2]-=V_com[2];
    }
    
    for (int i=0; i<Actin_num_of_Nodes; i++) {
        Actin_Node_Velocity[i][0]-=V_com[0];
        Actin_Node_Velocity[i][1]-=V_com[1];
        Actin_Node_Velocity[i][2]-=V_com[2];
    }
    for (int i=0; i<Membrane_num_of_Nodes; i++) {
        Kinetic_energy+=Membrane_Node_Velocity[i][0]*Membrane_Node_Velocity[i][0]+Membrane_Node_Velocity[i][1]*Membrane_Node_Velocity[i][1]+Membrane_Node_Velocity[i][2]*Membrane_Node_Velocity[i][2];
    }
    Kinetic_energy*=Membrane_Node_Mass;
    
    for (int i=0; i<Actin_num_of_Nodes; i++) {
        temp_Kinetic_energy+=Actin_Node_Velocity[i][0]*Actin_Node_Velocity[i][0]+Actin_Node_Velocity[i][1]*Actin_Node_Velocity[i][1]+Actin_Node_Velocity[i][2]*Actin_Node_Velocity[i][2];
    }
    temp_Kinetic_energy*=Actin_Node_Mass;
    Kinetic_energy+=temp_Kinetic_energy;
    
    alpha=sqrt(3*(Membrane_num_of_Nodes+Actin_num_of_Nodes)*KT)/Kinetic_energy;
    
    for(int i=0;i<Membrane_num_of_Nodes;i++)
    {
        Membrane_Node_Velocity [i][0]*=alpha;
        Membrane_Node_Velocity [i][1]*=alpha;
        Membrane_Node_Velocity [i][2]*=alpha;
        Membrane_Node_Velocity[i][0]+=V_com[0];
        Membrane_Node_Velocity[i][1]+=V_com[1];
        Membrane_Node_Velocity[i][2]+=V_com[2];
    }
    
    for(int i=0;i<Actin_num_of_Nodes;i++)
    {
        Actin_Node_Velocity [i] [0]*=alpha;
        Actin_Node_Velocity [i] [1]*=alpha;
        Actin_Node_Velocity [i] [2]*=alpha;
        Actin_Node_Velocity[i][0]+=V_com[0];
        Actin_Node_Velocity[i][1]+=V_com[1];
        Actin_Node_Velocity[i][2]+=V_com[2];
    }
    
}






