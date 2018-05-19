//
//  Actin_membrane_shared_functions.cpp
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 11/10/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//


#include "Actin_membrane_shared_functions.hpp"

using namespace std;


void Membrane_Actin_shared_Node_Identifier(vector<vector<int> > &Membrane_Actin_shared_Node_list, double Membrane_Node_Position[][3], double Actin_Node_Position[][3], int Membrane_num_of_Nodes){
    
    //    int actin_membrane_num_of_shared_nodes=0;
    vector<int> push;
    push.resize(2);
    for (int mem_node_counter=0; mem_node_counter<Membrane_num_of_Nodes; mem_node_counter++) {
        for (int act_node_counter=0; act_node_counter<Actin_num_of_Nodes; act_node_counter++) {
            
            if( (Membrane_Node_Position[mem_node_counter][0]==Actin_Node_Position[act_node_counter][0]) &&  (Membrane_Node_Position[mem_node_counter][1]==Actin_Node_Position[act_node_counter][1]) &&  (Membrane_Node_Position[mem_node_counter][2]==Actin_Node_Position[act_node_counter][2])   )
            {
                push[0]=mem_node_counter;
                push[1]=act_node_counter;
                Membrane_Actin_shared_Node_list.push_back(push);
                //                actin_membrane_num_of_shared_nodes++;
                break;
            }
        }
    }
    //    if(actin_membrane_num_of_shared_nodes!=Actin_Membrane_shared_num_of_Nodes)
    //    {
    //        cout<<"Problem in the 'Membrane_Actin_shared_Node_Identifier' function. The number of 'Actin_Membrane_shared_num_of_Nodes' is different from the number of Nodes identified on the Membrane"<<endl;
    //        exit (EXIT_FAILURE);
    //    }
    
}

