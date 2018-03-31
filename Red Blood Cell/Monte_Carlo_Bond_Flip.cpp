//
//  Monte_Carlo_Bond_Flip.cpp
//  Test_2
//
//  Created by Ali Farnudi on 12/03/2018.
//  Copyright © 2018 Ali Farnudi. All rights reserved.
//

#include "Monte_Carlo_Bond_Flip.hpp"

void Monte_carlo_bond_flip(double Membrane_Node_Position[][3], vector<vector<int> > &Membrane_new_triangle_list, vector<vector<int> > &membrane_triangle_pair_list, int Membrane_Node_Pair_list[][2], int Membrane_num_of_Node_Pairs, int Membrane_num_of_Nodes){
    //This function is designed to localy calculate the difference in the energy status of the membrane triangle bonds and by using a metropolise weight flip bonds between nodes but conserve the total number of triangles.
    
    int triangle_A, triangle_B;
    int triangle_A_node_A, triangle_A_node_B, triangle_A_node_C, triangle_B_node_D;
    int sorted_nodes[4];
    double energy_state_initial=0, energy_state_final=0;
    
    //    We will first pick out a random triangle from the membrane and identify it as triangle_A.
    triangle_A=rand()%(Membrane_new_triangle_list.size()-1);
    //    Next we will choose one of its three neighbours as the candidate for the bond flipping process, triangle_B.
    triangle_B=membrane_triangle_pair_list[triangle_A][rand()%3];
    
    //    Here we will sort the nodes. In this function the first node is the uncommon node of triangle A, next we have the two common nodes, and finally triangle B's uncommon node.
    trinalge_pair_node_sorter(triangle_A, triangle_B, Membrane_new_triangle_list, sorted_nodes);
    triangle_A_node_A=sorted_nodes[0];
    triangle_A_node_B=sorted_nodes[1];
    triangle_A_node_C=sorted_nodes[2];
    triangle_B_node_D=sorted_nodes[3];
    
    
    //Calculating the bond energy
    
    double le0,le1,lmax,lmin;
    double deltax_i, deltay_i, deltaz_i, temp_Node_distance_i;
    double deltax_f, deltay_f, deltaz_f, temp_Node_distance_f;
    int pos1, pos2, pos3, pos4;  // to making calculation of surface force easier
    double temp_p1[3], temp_p2[3], temp_p3[3], temp_p4[3];
    double  N1[3], N2[3], N3[3], p3p1[3], p3p2[3], p4p2[3], p4p1[3], sinus;// for exmple p3p1 is p3-p1 and ....
    
    double temp_AB[3], temp_AC[3], temp_ABxAC[3];
    
    le0=1.15000*Node_radius;
    lmax=1.33000*Node_radius;
    le1=0.85000*Node_radius;
    lmin=0.67000*Node_radius;
    
    deltax_i=Membrane_Node_Position[triangle_A_node_C][0]-Membrane_Node_Position[triangle_A_node_B][0];
    deltay_i=Membrane_Node_Position[triangle_A_node_C][1]-Membrane_Node_Position[triangle_A_node_B][1];
    deltaz_i=Membrane_Node_Position[triangle_A_node_C][2]-Membrane_Node_Position[triangle_A_node_B][2];
    
    deltax_f=Membrane_Node_Position[triangle_B_node_D][0]-Membrane_Node_Position[triangle_A_node_A][0];
    deltay_f=Membrane_Node_Position[triangle_B_node_D][1]-Membrane_Node_Position[triangle_A_node_A][1];
    deltaz_f=Membrane_Node_Position[triangle_B_node_D][2]-Membrane_Node_Position[triangle_A_node_A][2];
    
    temp_Node_distance_i=sqrt(deltax_i*deltax_i+deltay_i*deltay_i+deltaz_i*deltaz_i);
    temp_Node_distance_f=sqrt(deltax_f*deltax_f+deltay_f*deltay_f+deltaz_f*deltaz_f);
    
    double temp_exp_le0_i=exp(1.0/(le0-temp_Node_distance_i));
    double temp_exp_le1_i=exp(1.0/(temp_Node_distance_i-le1));
    double temp_exp_le0_f=exp(1.0/(le0-temp_Node_distance_f));
    double temp_exp_le1_f=exp(1.0/(temp_Node_distance_f-le1));
    
    
    // Initial state bond energy
    if(temp_Node_distance_i >le1  & temp_Node_distance_i < le0 )  //free zone
    {
        energy_state_initial=0 ; // free zone
    } else if(temp_Node_distance_i > le0  & temp_Node_distance_i <lmax )  //bondforce
    {
        energy_state_initial= Membrane_spring_coefficient*temp_exp_le0_i/(lmax-temp_Node_distance_i);
        
    } else if(temp_Node_distance_i < le1   &  temp_Node_distance_i > lmin  )  // repulsive force
    {
        energy_state_initial= Membrane_spring_coefficient*temp_exp_le1_i/(temp_Node_distance_i-lmin);
    } else if(temp_Node_distance_i > lmax )
    {
        energy_state_initial=   1.81599  + 965.31 * ( temp_Node_distance_i - 1.3280*Node_radius )+0.5*Membrane_spring_force_cutt_off * ( temp_Node_distance_i - 1.3280*Node_radius ) * ( temp_Node_distance_i - 1.3280*Node_radius );
    } else if(temp_Node_distance_i < lmin )
    {
        energy_state_initial = 1.85038 + 1005.05 * ( 0.671965*Node_radius - temp_Node_distance_i )+0.5*Membrane_spring_force_cutt_off*( 0.671965*Node_radius - temp_Node_distance_i )*( 0.671965*Node_radius - temp_Node_distance_i );
    }
    //    cout<<"energy_state_initial:\tbond="<<energy_state_initial<<endl;
    
    // Final state bond energy
    if(temp_Node_distance_f >le1  & temp_Node_distance_f < le0 )  //free zone
    {
        energy_state_final=0 ; // free zone
    } else if(temp_Node_distance_f > le0  & temp_Node_distance_f <lmax )  //bondforce
    {
        energy_state_final= Membrane_spring_coefficient*temp_exp_le0_f/(lmax-temp_Node_distance_f);
        
    } else if(temp_Node_distance_f < le1   &  temp_Node_distance_f > lmin  )  // repulsive force
    {
        energy_state_final= Membrane_spring_coefficient*temp_exp_le1_f/(temp_Node_distance_f-lmin);
    } else if(temp_Node_distance_f > lmax )
    {
        energy_state_final=   1.81599  + 965.31 * ( temp_Node_distance_f - 1.3280*Node_radius )+0.5*Membrane_spring_force_cutt_off * ( temp_Node_distance_f - 1.3280*Node_radius ) * ( temp_Node_distance_f - 1.3280*Node_radius );
    } else if(temp_Node_distance_f < lmin )
    {
        energy_state_final = 1.85038 + 1005.05 * ( 0.671965*Node_radius - temp_Node_distance_f )+0.5*Membrane_spring_force_cutt_off*( 0.671965*Node_radius - temp_Node_distance_f )*( 0.671965*Node_radius - temp_Node_distance_f );
    }
    
    //    cout<<"energy_state_final:\tbond="<<energy_state_final<<endl<<endl;
    
    
    
    //    cout<<"Triangle A two node neighbours before flip: ";
    // Determining triangle A and B's other two neighbours in the initial state.
    int triangle_A_neighbours[2], triangle_B_neighbours[2], neighbour_number=0;
    for (int i=0; i<3; i++) {
        if (triangle_B != membrane_triangle_pair_list[triangle_A][i]) {
            triangle_A_neighbours[neighbour_number]=membrane_triangle_pair_list[triangle_A][i];
            //            cout<<"\t A"<<neighbour_number<<"="<<membrane_triangle_pair_list[triangle_A][i];
            neighbour_number++;
        }
    }
    //    cout<<"\nTriangle B two node neighbours before flip: ";
    neighbour_number=0;
    for (int i=0; i<3; i++) {
        if (triangle_A != membrane_triangle_pair_list[triangle_B][i]) {
            triangle_B_neighbours[neighbour_number]=membrane_triangle_pair_list[triangle_B][i];
            //            cout<<"\t B"<<neighbour_number<<"="<<membrane_triangle_pair_list[triangle_B][i];
            neighbour_number++;
        }
    }
    
    //    cout<<"Triangle A="<<triangle_A<<"\tTriangle A0="<<triangle_A_neighbours[0]<<"\tTriangle A1="<<triangle_A_neighbours[1]<<"\nTriangle B="<<triangle_B<<"\tTriangle B0="<<triangle_B_neighbours[0]<<"\tTriangle B1="<<triangle_B_neighbours[1]<<endl;
    
    //    I need the com of the membrane so I can determin if the cross product of the triangle edges are pointing in or out of the membrane.
    double membrane_com[3];
    membrane_com[0]=0;
    membrane_com[1]=0;
    membrane_com[2]=0;
    //----------------------membrane---------------------
    for(int i=0;i<Membrane_num_of_Nodes;i++)
    {
        membrane_com [0] += Membrane_Node_Position [i][0];
        membrane_com [1] += Membrane_Node_Position [i][1];
        membrane_com [2] += Membrane_Node_Position [i][2];
        
    }
    membrane_com[0]/=Membrane_num_of_Nodes;
    membrane_com[1]/=Membrane_num_of_Nodes;
    membrane_com[2]/=Membrane_num_of_Nodes;
    
    double outward_vector[3];
    
    //calculating triangle A and B's bending energy in the initial state
    
    //Because I have used Tiam's old code for the energy calculations I wil re sort the nodes so that the common nodes are the first two nodes (pos1 and pos2).
    pos1=triangle_A_node_B;
    pos2=triangle_A_node_C;
    pos3=triangle_A_node_A;
    pos4=triangle_B_node_D;
    
    for (int index=0; index<3; index++) {
        temp_p1[index]=Membrane_Node_Position[pos1][index];
        temp_p2[index]=Membrane_Node_Position[pos2][index];
        temp_p3[index]=Membrane_Node_Position[pos3][index];
        temp_p4[index]=Membrane_Node_Position[pos4][index];
        
        p3p1[index]=temp_p3[index]-temp_p1[index];
        p3p2[index]=temp_p3[index]-temp_p2[index];
        p4p2[index]=temp_p4[index]-temp_p2[index];
        p4p1[index]=temp_p4[index]-temp_p1[index];
        outward_vector[index]=temp_p1[index]-membrane_com[index];
    }
    
    crossvector(N1, p3p2, p3p1);
    crossvector(N2, p4p1, p4p2);
    
    //Making sure N1 and N2 are pointing out of the membrane
    if (innerproduct(N1, outward_vector)<0) {
        N1[0]=-N1[0];
        N1[1]=-N1[1];
        N1[2]=-N1[2];
    }
    if (innerproduct(N2, outward_vector)<0) {
        N2[0]=-N2[0];
        N2[1]=-N2[1];
        N2[2]=-N2[2];
    }
    
    crossvector(N3, N2, N1);
    sinus=vectorlength(N3)/(vectorlength(N2)*vectorlength(N1));
    
    energy_state_initial += Membrane_bending_coefficient*(1.0 -  ( innerproduct(N1,N2)/(vectorlength(N1)*vectorlength(N2)) ));
    //    cout<<"A and B bending energy (initial state)="<<Membrane_bending_coefficient*(1.0 -  ( innerproduct(N1,N2)/(vectorlength(N1)*vectorlength(N2)) ))<<endl;
    
    //bending energy of triangle_A and its neighbours
    //    double temp_energy=0;
    for(int j=0; j<2 ;j++)
    {
        trinalge_pair_node_sorter(triangle_A, triangle_A_neighbours[j], Membrane_new_triangle_list, sorted_nodes);
        
        pos1=sorted_nodes[1];
        pos2=sorted_nodes[2];
        pos3=sorted_nodes[0];
        pos4=sorted_nodes[3];
        
        
        for (int index=0; index<3; index++) {
            temp_p1[index]=Membrane_Node_Position[pos1][index];
            temp_p2[index]=Membrane_Node_Position[pos2][index];
            temp_p3[index]=Membrane_Node_Position[pos3][index];
            temp_p4[index]=Membrane_Node_Position[pos4][index];
            
            p3p1[index]=temp_p3[index]-temp_p1[index];
            p3p2[index]=temp_p3[index]-temp_p2[index];
            p4p2[index]=temp_p4[index]-temp_p2[index];
            p4p1[index]=temp_p4[index]-temp_p1[index];
            outward_vector[index]=temp_p1[index]-membrane_com[index];
        }
        
        
        crossvector(N1, p3p2, p3p1);
        crossvector(N2, p4p1, p4p2);
        
        if (innerproduct(N1, outward_vector)<0) {
            N1[0]=-N1[0];
            N1[1]=-N1[1];
            N1[2]=-N1[2];
        }
        if (innerproduct(N2, outward_vector)<0) {
            N2[0]=-N2[0];
            N2[1]=-N2[1];
            N2[2]=-N2[2];
        }
        
        crossvector(N3, N2, N1);
        sinus=vectorlength(N3)/(vectorlength(N2)*vectorlength(N1));
        
        
        energy_state_initial += Membrane_bending_coefficient*(1.0 -  ( innerproduct(N1,N2)/(vectorlength(N1)*vectorlength(N2)) ));
        //        temp_energy+=Membrane_bending_coefficient*(1.0 -  ( innerproduct(N1,N2)/(vectorlength(N1)*vectorlength(N2)) ));
    }
    
    //bending energy of triangle_B and its neighbours
    for(int j=0; j<2 ;j++)
    {
        trinalge_pair_node_sorter(triangle_B, triangle_B_neighbours[j], Membrane_new_triangle_list, sorted_nodes);
        
        pos1=sorted_nodes[1];
        pos2=sorted_nodes[2];
        pos3=sorted_nodes[0];
        pos4=sorted_nodes[3];
        
        
        for (int index=0; index<3; index++) {
            temp_p1[index]=Membrane_Node_Position[pos1][index];
            temp_p2[index]=Membrane_Node_Position[pos2][index];
            temp_p3[index]=Membrane_Node_Position[pos3][index];
            temp_p4[index]=Membrane_Node_Position[pos4][index];
            
            p3p1[index]=temp_p3[index]-temp_p1[index];
            p3p2[index]=temp_p3[index]-temp_p2[index];
            p4p2[index]=temp_p4[index]-temp_p2[index];
            p4p1[index]=temp_p4[index]-temp_p1[index];
            outward_vector[index]=temp_p1[index]-membrane_com[index];
        }
        
        
        crossvector(N1, p3p2, p3p1);
        crossvector(N2, p4p1, p4p2);
        
        if (innerproduct(N1, outward_vector)<0) {
            N1[0]=-N1[0];
            N1[1]=-N1[1];
            N1[2]=-N1[2];
        }
        if (innerproduct(N2, outward_vector)<0) {
            N2[0]=-N2[0];
            N2[1]=-N2[1];
            N2[2]=-N2[2];
        }
        
        crossvector(N3, N2, N1);
        sinus=vectorlength(N3)/(vectorlength(N2)*vectorlength(N1));
        
        energy_state_initial += Membrane_bending_coefficient*(1.0 -  ( innerproduct(N1,N2)/(vectorlength(N1)*vectorlength(N2)) ));
        //        temp_energy+=Membrane_bending_coefficient*(1.0 -  ( innerproduct(N1,N2)/(vectorlength(N1)*vectorlength(N2)) ));
    }
    //    cout<<"A & B neihbours bending energy (initial state)="<<temp_energy<<endl;
    //    temp_energy=0;
    
    
    //Final_state:
    
    //Triangle_A and Triangle_B bending energy in the final state
    pos1=triangle_A_node_A;
    pos2=triangle_B_node_D;
    pos3=triangle_A_node_B;
    pos4=triangle_A_node_C;
    
    for (int index=0; index<3; index++) {
        temp_p1[index]=Membrane_Node_Position[pos1][index];
        temp_p2[index]=Membrane_Node_Position[pos2][index];
        temp_p3[index]=Membrane_Node_Position[pos3][index];
        temp_p4[index]=Membrane_Node_Position[pos4][index];
        
        p3p1[index]=temp_p3[index]-temp_p1[index];
        p3p2[index]=temp_p3[index]-temp_p2[index];
        p4p2[index]=temp_p4[index]-temp_p2[index];
        p4p1[index]=temp_p4[index]-temp_p1[index];
        outward_vector[index]=temp_p1[index]-membrane_com[index];
    }
    
    crossvector(N1, p3p2, p3p1);
    crossvector(N2, p4p1, p4p2);
    
    if (innerproduct(N1, outward_vector)<0) {
        N1[0]=-N1[0];
        N1[1]=-N1[1];
        N1[2]=-N1[2];
    }
    if (innerproduct(N2, outward_vector)<0) {
        N2[0]=-N2[0];
        N2[1]=-N2[1];
        N2[2]=-N2[2];
    }
    
    crossvector(N3, N2, N1);
    sinus=vectorlength(N3)/(vectorlength(N2)*vectorlength(N1));
    
    energy_state_final += Membrane_bending_coefficient*(1.0 -  ( innerproduct(N1,N2)/(vectorlength(N1)*vectorlength(N2)) ));
    //    cout<<"\nA and B bending energy (final state)="<<Membrane_bending_coefficient*(1.0 -  ( innerproduct(N1,N2)/(vectorlength(N1)*vectorlength(N2)) ))<<endl;
    
    int new_sorted_nodes[4][4];
    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
            new_sorted_nodes[i][j]=-1;
        }
    }
    //Calculating the bending energy of triangle A and B and their neighbours in the final state
    int triangle_A_new_neighbours[2], triangle_B_new_neighbours[2];
    
    new_trinalge_pair_node_sorter(triangle_A_node_A, triangle_A_node_B, triangle_A_node_C, triangle_B_node_D, triangle_A_neighbours, triangle_B_neighbours, Membrane_new_triangle_list, new_sorted_nodes, triangle_A_new_neighbours, triangle_B_new_neighbours);
    //    cout<<"\nTriangle A two node neighbours after flip: ";
    //    cout<<"\t A0="<<triangle_A_new_neighbours[0];
    //    cout<<"\t A1="<<triangle_A_new_neighbours[1];
    //    cout<<"\nTriangle B two node neighbours after flip: ";
    //    cout<<"\t B0="<<triangle_B_new_neighbours[0];
    //    cout<<"\t B1="<<triangle_B_new_neighbours[1]<<endl;
    //    cout<<"============\n";
    for(int j=0; j<4 ;j++)
    {
        pos1=new_sorted_nodes[j][1];
        pos2=new_sorted_nodes[j][2];
        pos3=new_sorted_nodes[j][0];
        pos4=new_sorted_nodes[j][3];
        
        //        cout<<"\nnew_sorted_nodes["<<j<<"][1]="<<pos1<<"\nnew_sorted_nodes["<<j<<"][2]="<<pos2<<"\nnew_sorted_nodes["<<j<<"][0]="<<pos3<<"\nnew_sorted_nodes["<<j<<"][3]="<<pos4<<endl;
        
        for (int index=0; index<3; index++) {
            temp_p1[index]=Membrane_Node_Position[pos1][index];
            temp_p2[index]=Membrane_Node_Position[pos2][index];
            temp_p3[index]=Membrane_Node_Position[pos3][index];
            temp_p4[index]=Membrane_Node_Position[pos4][index];
            
            p3p1[index]=temp_p3[index]-temp_p1[index];
            p3p2[index]=temp_p3[index]-temp_p2[index];
            p4p2[index]=temp_p4[index]-temp_p2[index];
            p4p1[index]=temp_p4[index]-temp_p1[index];
            outward_vector[index]=temp_p1[index]-membrane_com[index];
        }
        
        
        crossvector(N1, p3p2, p3p1);
        crossvector(N2, p4p1, p4p2);
        
        if (innerproduct(N1, outward_vector)<0) {
            N1[0]=-N1[0];
            N1[1]=-N1[1];
            N1[2]=-N1[2];
        }
        if (innerproduct(N2, outward_vector)<0) {
            N2[0]=-N2[0];
            N2[1]=-N2[1];
            N2[2]=-N2[2];
        }
        
        crossvector(N3, N2, N1);
        sinus=vectorlength(N3)/(vectorlength(N2)*vectorlength(N1));
        
        energy_state_final += Membrane_bending_coefficient*(1.0 -  ( innerproduct(N1,N2)/(vectorlength(N1)*vectorlength(N2)) ));
        //        temp_energy+=Membrane_bending_coefficient*(1.0 -  ( innerproduct(N1,N2)/(vectorlength(N1)*vectorlength(N2)) ));
    }
    //    cout<<"A & B neihbours bending energy (final state)="<<temp_energy<<endl;
    //    temp_energy=0;
    //    cout<<"\nenergy_state_initial="<<energy_state_initial<<"\n"<<"energy_state_final="<<energy_state_final<<endl<<endl;
    //    exit(EXIT_FAILURE);
    
    if ( ((double) rand() / (RAND_MAX))<exp((energy_state_initial-energy_state_final)/KT) ){
        //        node pair list update
        //        cout<<"Monte Carlo\n";
        for (int i=0; i<Membrane_num_of_Node_Pairs; i++) {
            if ((Membrane_Node_Pair_list[i][0]==triangle_A_node_A && Membrane_Node_Pair_list[i][1]==triangle_B_node_D) || (Membrane_Node_Pair_list[i][0]==triangle_B_node_D && Membrane_Node_Pair_list[i][1]==triangle_A_node_A) ) {
                
                Membrane_Node_Pair_list[i][0]=triangle_A_node_B;
                Membrane_Node_Pair_list[i][1]=triangle_A_node_C;
            }
        }
        
        //        triangle node list update
        
        Membrane_new_triangle_list[triangle_A][0]=triangle_A_node_A;
        Membrane_new_triangle_list[triangle_A][1]=triangle_A_node_B;
        Membrane_new_triangle_list[triangle_A][2]=triangle_B_node_D;
        
        for (int i=0; i<3; i++) {
            temp_AB[i]=Membrane_Node_Position[triangle_A_node_B][i]-Membrane_Node_Position[triangle_A_node_A][i];
            temp_AC[i]=Membrane_Node_Position[triangle_B_node_D][i]-Membrane_Node_Position[triangle_A_node_A][i];
            outward_vector[i]=Membrane_Node_Position[triangle_A_node_A][i]-membrane_com[i];
        }
        crossvector(temp_ABxAC, temp_AB, temp_AC);
        if (innerproduct(temp_ABxAC, outward_vector)<0) {
            Membrane_new_triangle_list[triangle_A][1]=triangle_B_node_D;
            Membrane_new_triangle_list[triangle_A][2]=triangle_A_node_B;
        }
        
        
        Membrane_new_triangle_list[triangle_B][0]=triangle_A_node_A;
        Membrane_new_triangle_list[triangle_B][1]=triangle_A_node_C;
        Membrane_new_triangle_list[triangle_B][2]=triangle_B_node_D;
        
        for (int i=0; i<3; i++) {
            temp_AB[i]=Membrane_Node_Position[triangle_A_node_C][i]-Membrane_Node_Position[triangle_A_node_A][i];
            temp_AC[i]=Membrane_Node_Position[triangle_B_node_D][i]-Membrane_Node_Position[triangle_A_node_A][i];
        }
        crossvector(temp_ABxAC, temp_AB, temp_AC);
        if (innerproduct(temp_ABxAC, outward_vector)<0) {
            Membrane_new_triangle_list[triangle_A][1]=triangle_B_node_D;
            Membrane_new_triangle_list[triangle_A][2]=triangle_A_node_C;
        }
        
        //        triangle triangle neighbour list update
        membrane_triangle_pair_list[triangle_A][0]=triangle_A_neighbours[0];
        membrane_triangle_pair_list[triangle_A][1]=triangle_A_neighbours[1];
        membrane_triangle_pair_list[triangle_A][2]=triangle_B;
        
        membrane_triangle_pair_list[triangle_B][0]=triangle_B_neighbours[0];
        membrane_triangle_pair_list[triangle_B][1]=triangle_B_neighbours[1];
        membrane_triangle_pair_list[triangle_B][2]=triangle_A;
        
        for (int i=0; i<3; i++) {
            if (membrane_triangle_pair_list[triangle_A_neighbours[0]][i]==triangle_B) {
                membrane_triangle_pair_list[triangle_A_neighbours[0]][i]=triangle_A;
            }
            if (membrane_triangle_pair_list[triangle_A_neighbours[1]][i]==triangle_B) {
                membrane_triangle_pair_list[triangle_A_neighbours[1]][i]=triangle_A;
            }
            if (membrane_triangle_pair_list[triangle_B_neighbours[0]][i]==triangle_A) {
                membrane_triangle_pair_list[triangle_B_neighbours[0]][i]=triangle_B;
            }
            if (membrane_triangle_pair_list[triangle_B_neighbours[1]][i]==triangle_A) {
                membrane_triangle_pair_list[triangle_B_neighbours[1]][i]=triangle_B;
            }
        }
        
        //need to remember to resort the triangle nodes for the normal direction stuff!!
    }
}

void trinalge_pair_node_sorter(int triangle_1, int triangle_2, vector<vector<int> > Membrane_new_triangle_list, int sorted_nodes[4]){
    
    if (Membrane_new_triangle_list[triangle_1][0]!=Membrane_new_triangle_list[triangle_2][0] && Membrane_new_triangle_list[triangle_1][0]!=Membrane_new_triangle_list[triangle_2][1] &&
        Membrane_new_triangle_list[triangle_1][0]!=Membrane_new_triangle_list[triangle_2][2])  {
        sorted_nodes[0]=Membrane_new_triangle_list[triangle_1][0];
        sorted_nodes[1]=Membrane_new_triangle_list[triangle_1][1];
        sorted_nodes[2]=Membrane_new_triangle_list[triangle_1][2];
        if (Membrane_new_triangle_list[triangle_2][0]!=Membrane_new_triangle_list[triangle_1][1] && Membrane_new_triangle_list[triangle_2][0]!=Membrane_new_triangle_list[triangle_1][2]) {
            sorted_nodes[3]=Membrane_new_triangle_list[triangle_2][0];
        } else if (Membrane_new_triangle_list[triangle_2][1]!=Membrane_new_triangle_list[triangle_1][1] && Membrane_new_triangle_list[triangle_2][1]!=Membrane_new_triangle_list[triangle_1][2]){
            sorted_nodes[3]=Membrane_new_triangle_list[triangle_2][1];
        } else {
            sorted_nodes[3]=Membrane_new_triangle_list[triangle_2][2];
        }
    } else if (Membrane_new_triangle_list[triangle_1][1]!=Membrane_new_triangle_list[triangle_2][0] && Membrane_new_triangle_list[triangle_1][1]!=Membrane_new_triangle_list[triangle_2][1] &&
               Membrane_new_triangle_list[triangle_1][1]!=Membrane_new_triangle_list[triangle_2][2])  {
        sorted_nodes[0]=Membrane_new_triangle_list[triangle_1][1];
        sorted_nodes[1]=Membrane_new_triangle_list[triangle_1][0];
        sorted_nodes[2]=Membrane_new_triangle_list[triangle_1][2];
        if (Membrane_new_triangle_list[triangle_2][0]!=Membrane_new_triangle_list[triangle_1][0] && Membrane_new_triangle_list[triangle_2][0]!=Membrane_new_triangle_list[triangle_1][2]) {
            sorted_nodes[3]=Membrane_new_triangle_list[triangle_2][0];
        } else if (Membrane_new_triangle_list[triangle_2][1]!=Membrane_new_triangle_list[triangle_1][0] && Membrane_new_triangle_list[triangle_2][1]!=Membrane_new_triangle_list[triangle_1][2]){
            sorted_nodes[3]=Membrane_new_triangle_list[triangle_2][1];
        } else {
            sorted_nodes[3]=Membrane_new_triangle_list[triangle_2][2];
        }
    } else {
        sorted_nodes[0]=Membrane_new_triangle_list[triangle_1][2];
        sorted_nodes[1]=Membrane_new_triangle_list[triangle_1][0];
        sorted_nodes[2]=Membrane_new_triangle_list[triangle_1][1];
        if (Membrane_new_triangle_list[triangle_2][0]!=Membrane_new_triangle_list[triangle_1][0] && Membrane_new_triangle_list[triangle_2][0]!=Membrane_new_triangle_list[triangle_1][1]) {
            sorted_nodes[3]=Membrane_new_triangle_list[triangle_2][0];
        } else if (Membrane_new_triangle_list[triangle_2][1]!=Membrane_new_triangle_list[triangle_1][0] && Membrane_new_triangle_list[triangle_2][1]!=Membrane_new_triangle_list[triangle_1][1]){
            sorted_nodes[3]=Membrane_new_triangle_list[triangle_2][1];
        } else {
            sorted_nodes[3]=Membrane_new_triangle_list[triangle_2][2];
        }
    }
}


void new_trinalge_pair_node_sorter(int triangle_A_node_A, int triangle_A_node_B, int triangle_A_node_C, int triangle_B_node_D, int triangle_A_neighbours[2], int triangle_B_neighbours[2], vector<vector<int> > Membrane_new_triangle_list, int new_sorted_nodes[4][4], int triangle_A_new_neighbours[2], int triangle_B_new_neighbours[2]){
    
    //    cout<<"\n\nTriangle A neighbours: "<<"\t A0="<<triangle_A_neighbours[0]<<"\t A1="<<triangle_A_neighbours[1];
    //    cout<<"\nTriangle B neighbours: "<<"\t B0="<<triangle_B_neighbours[0]<<"\t B1="<<triangle_B_neighbours[1];
    
    for (int i=0; i<2; i++) {
        for (int neighbours_browser=0; neighbours_browser<3; neighbours_browser++) {
            
            if ((Membrane_new_triangle_list[triangle_A_neighbours[i]][neighbours_browser]!=triangle_A_node_A && Membrane_new_triangle_list[triangle_A_neighbours[i]][neighbours_browser]!=triangle_A_node_B) &&
                ((Membrane_new_triangle_list[triangle_A_neighbours[i]][(neighbours_browser+1)%3]==triangle_A_node_B && Membrane_new_triangle_list[triangle_A_neighbours[i]][(neighbours_browser+2)%3]==triangle_A_node_A) || (Membrane_new_triangle_list[triangle_A_neighbours[i]][(neighbours_browser+2)%3]==triangle_A_node_B && Membrane_new_triangle_list[triangle_A_neighbours[i]][(neighbours_browser+1)%3]==triangle_A_node_A))
                )  {
                new_sorted_nodes[0][0]=Membrane_new_triangle_list[triangle_A_neighbours[i]][neighbours_browser];
                new_sorted_nodes[0][1]=triangle_A_node_A;
                new_sorted_nodes[0][2]=triangle_A_node_B;
                new_sorted_nodes[0][3]=triangle_B_node_D;
                triangle_A_new_neighbours[0]=triangle_A_neighbours[i];
                
                for (int j=0; j<3; j++) {
                    if ((Membrane_new_triangle_list[triangle_A_neighbours[(i+1)%2]][j]!=triangle_A_node_A && Membrane_new_triangle_list[triangle_A_neighbours[(i+1)%2]][j]!=triangle_A_node_C) &&
                        ((Membrane_new_triangle_list[triangle_A_neighbours[(i+1)%2]][(j+1)%3]==triangle_A_node_A && Membrane_new_triangle_list[triangle_A_neighbours[(i+1)%2]][(j+2)%3]==triangle_A_node_C) || (Membrane_new_triangle_list[triangle_A_neighbours[(i+1)%2]][(j+2)%3]==triangle_A_node_A && Membrane_new_triangle_list[triangle_A_neighbours[(i+1)%2]][(j+1)%3]==triangle_A_node_C))
                        ){
                        new_sorted_nodes[2][0]=Membrane_new_triangle_list[triangle_A_neighbours[(i+1)%2]][j];
                        new_sorted_nodes[2][1]=triangle_A_node_A;
                        new_sorted_nodes[2][2]=triangle_A_node_C;
                        new_sorted_nodes[2][3]=triangle_B_node_D;
                        triangle_B_new_neighbours[0]=triangle_A_neighbours[(i+1)%2];
                        break;
                    }
                }
                break;
            } // end of if (Membrane_triangle_list[triangle_A_neighbours[i]][neighbours_browser]!=triangle_A_node_A && Membrane_triangle_list[triangle_A_neighbours[i]][neighbours_browser]!=triangle_A_node_B && (Membrane_triangle_list[triangle_A_neighbours[i]][(neighbours_browser+1)%3]==triangle_A_node_B && Membrane_triangle_list[triangle_A_neighbours[i]][(neighbours_browser+2)%3]==triangle_A_node_A) || (Membrane_triangle_list[triangle_A_neighbours[i]][(neighbours_browser+2)%3]==triangle_A_node_B && Membrane_triangle_list[triangle_A_neighbours[i]][(neighbours_browser+1)%3]==triangle_A_node_A) )
            
        }// end of for (int neighbours_browser=0; neighbours_browser<3; neighbours_browser++) {
        
    }// end of for (int i=0; i<2; i++) {
    
    for (int i=0; i<2; i++) {
        for (int neighbours_browser=0; neighbours_browser<3; neighbours_browser++) {
            
            if ((Membrane_new_triangle_list[triangle_B_neighbours[i]][neighbours_browser]!=triangle_A_node_B && Membrane_new_triangle_list[triangle_B_neighbours[i]][neighbours_browser]!=triangle_B_node_D) &&
                ((Membrane_new_triangle_list[triangle_B_neighbours[i]][(neighbours_browser+1)%3]==triangle_A_node_B && Membrane_new_triangle_list[triangle_B_neighbours[i]][(neighbours_browser+2)%3]==triangle_B_node_D) || (Membrane_new_triangle_list[triangle_B_neighbours[i]][(neighbours_browser+2)%3]==triangle_A_node_B && Membrane_new_triangle_list[triangle_B_neighbours[i]][(neighbours_browser+1)%3]==triangle_B_node_D))
                )  {
                new_sorted_nodes[1][0]=Membrane_new_triangle_list[triangle_B_neighbours[i]][neighbours_browser];
                new_sorted_nodes[1][1]=triangle_A_node_B;
                new_sorted_nodes[1][2]=triangle_B_node_D;
                new_sorted_nodes[1][3]=triangle_A_node_A;
                triangle_A_new_neighbours[1]=triangle_B_neighbours[i];
                for (int j=0; j<3; j++) {
                    if ((Membrane_new_triangle_list[triangle_B_neighbours[(i+1)%2]][j]!=triangle_B_node_D && Membrane_new_triangle_list[triangle_B_neighbours[(i+1)%2]][j]!=triangle_A_node_C) &&
                        ((Membrane_new_triangle_list[triangle_B_neighbours[(i+1)%2]][(j+1)%3]==triangle_B_node_D && Membrane_new_triangle_list[triangle_B_neighbours[(i+1)%2]][(j+2)%3]==triangle_A_node_C) || (Membrane_new_triangle_list[triangle_B_neighbours[(i+1)%2]][(j+2)%3]==triangle_B_node_D && Membrane_new_triangle_list[triangle_B_neighbours[(i+1)%2]][(j+1)%3]==triangle_A_node_C))
                        ){
                        new_sorted_nodes[3][0]=Membrane_new_triangle_list[triangle_B_neighbours[(i+1)%2]][j];
                        new_sorted_nodes[3][1]=triangle_B_node_D;
                        new_sorted_nodes[3][2]=triangle_A_node_C;
                        new_sorted_nodes[3][3]=triangle_A_node_A;
                        triangle_B_new_neighbours[1]=triangle_B_neighbours[(i+1)%2];
                        break;
                    }
                }
                break;
            } // end of if (Membrane_triangle_list[triangle_A_neighbours[i]][neighbours_browser]!=triangle_A_node_A && Membrane_triangle_list[triangle_A_neighbours[i]][neighbours_browser]!=triangle_A_node_B && (Membrane_triangle_list[triangle_A_neighbours[i]][(neighbours_browser+1)%3]==triangle_A_node_B && Membrane_triangle_list[triangle_A_neighbours[i]][(neighbours_browser+2)%3]==triangle_A_node_A) || (Membrane_triangle_list[triangle_A_neighbours[i]][(neighbours_browser+2)%3]==triangle_A_node_B && Membrane_triangle_list[triangle_A_neighbours[i]][(neighbours_browser+1)%3]==triangle_A_node_A) )
            
        }// end of for (int neighbours_browser=0; neighbours_browser<3; neighbours_browser++) {
        
    }// end of for (int i=0; i<2; i++) {
    
}









