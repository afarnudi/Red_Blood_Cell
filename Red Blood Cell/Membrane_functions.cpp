//
//  Membrane_functions.cpp
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 05/09/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//

#include "Membrane_functions.hpp"
#include <fstream>
#include <iostream>
#include <math.h>
#include "General_functions.hpp"

using namespace std;

void Membrane_constructor(double Membrane_Node_Position[Membrane_num_of_Nodes][3], double Membrane_Node_Velocity [Membrane_num_of_Nodes][3], double Membrane_Node_Force[Membrane_num_of_Nodes][3], int  Membrane_triangle_list [Membrane_num_of_Triangles][3])
{
    //    In this section we set all the Node forces and velocities to zero.
    for(int i=0;i<Membrane_num_of_Nodes;i++)
    {
        Membrane_Node_Velocity[i][0]= 0.0;
        Membrane_Node_Velocity[i][1]=0.0;
        Membrane_Node_Velocity[i][2]= 0.0;
        
        Membrane_Node_Force[i][0]= 0.0;
        Membrane_Node_Force[i][1]= 0.0;
        Membrane_Node_Force[i][2]= 0.0;
    }
    ifstream read; //This is the main ifstream that will read the Gmesh-Membrane generated file
    read.open("membrane"); //It should be noted that the name of the file should not contain '-'. I don't know why but the memory managnet of the arrays (at the very least) in the programme will collapse when we use '-' in the file name.
    int temp_int; // This is just a temp intiger charachter that we use to read unnecessary Gmesh generated intigers. We never use these intigers in the actual programme.
    
    read>> temp_int;
    // In this section the Node coordinates are read from the Gmesh membrane generated file. These include both the Nodes on the Membrane and on the nucleus membrane.
    for(int i=0;i<Membrane_num_of_Nodes;i++)
    {
        read>> temp_int;
        read>> Membrane_Node_Position [i][0];
        read>> Membrane_Node_Position [i][1];
        read>> Membrane_Node_Position [i][2];
    }
    
    // In this section the Node list that make up triangles on the outer membrane and nucleus are read from the Gmesh generated file.
    read>> temp_int;
    for(int i=0;i<Membrane_num_of_Triangles;i++)
    {
        read>>temp_int;
        read>>temp_int;
        read>>temp_int;
        read>>temp_int;
        read>>temp_int;
        
        read>>Membrane_triangle_list[i][0];
        read>>Membrane_triangle_list[i][1];
        read>>Membrane_triangle_list[i][2];
        //We have written the programme so that the Node indecies start from '0'. The node indecies in the Gmesh generated file start from '1', so we correct the 'Membrane_triangle_list'.
        Membrane_triangle_list[i][0]--;
        Membrane_triangle_list[i][1]--;
        Membrane_triangle_list[i][2]--;
    }
} //END OF: Membrane_constructor function



void Membrane_Normal_direction_Identifier( double  Membrane_Node_Position [Membrane_num_of_Nodes][3], int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2], int  &Outer_Membrane_num_of_triangles, int &Nucleus_Membrane_num_of_triangles, bool cell_has_nucleus)
{
    //    The 'Outer_Membrane_num_of_triangles' and 'Nucleus_Membrane_num_of_triangles' are passed by reference to this function and the values are set here. These values will never change throughout the code.
    Outer_Membrane_num_of_triangles=0;
    Nucleus_Membrane_num_of_triangles=0;
    
    //     Each triangle has three nodes A, B, and C. Throughout the programme, A=Membrane_triangle_list[][0], B=Membrane_triangle_list[][1], C=Membrane_triangle_list[][2].
    double AC[3], AB[3], ABxAC[3], xyz[3], temp_radius;
    
    for(int i=0; i<Membrane_num_of_Triangles; i++) // membrane or nucleus
    {
        //        temp_radius is used to store the distance of the triangle from the origin
        temp_radius=sqrt( Membrane_Node_Position[Membrane_triangle_list[i][0]][0]*Membrane_Node_Position[Membrane_triangle_list[i][0]][0] + Membrane_Node_Position[Membrane_triangle_list[i][0]][1]*Membrane_Node_Position[Membrane_triangle_list[i][0]][1] + Membrane_Node_Position[Membrane_triangle_list[i][0]][2]*Membrane_Node_Position[Membrane_triangle_list[i][0]][2]  );//cout<< r <<endl;
        
        //        The Membrane_Normal_direction[i][0]= +/- 1; is used to identify triangles on the membrane (+1) and nucleus (-1)
        if(temp_radius < Membrane_Radius+0.5 & temp_radius> Membrane_Radius - 0.5)
        {
            Membrane_Normal_direction[i][0]= +1;
            Outer_Membrane_num_of_triangles++;
        }
        
        else if((temp_radius < Nucleus_Membrane_radius+0.5 & temp_radius> Nucleus_Membrane_radius - 0.5) && cell_has_nucleus==true)
        {
            Membrane_Normal_direction[i][0]= -1;
            Nucleus_Membrane_num_of_triangles++;
        }
        else
        {
            cout<< " Error in Normal vector initialisation:  Double check the inputs" <<endl;
        }
    }// END OF: for(int i=0; i<Membrane_num_of_Triangles; i++)
    
    
    for(  int i=0;i<Membrane_num_of_Triangles;i++  )
    {
        
        AB[0]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][0]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][0];
        AB[1]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][1]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][1];
        AB[2]=Membrane_Node_Position[ Membrane_triangle_list[i][1]][2]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][2];
        
        AC[0]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][0]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][0];
        AC[1]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][1]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][1];
        AC[2]=Membrane_Node_Position[ Membrane_triangle_list[i][2]][2]-Membrane_Node_Position[ Membrane_triangle_list[i][0]][2];
        
        xyz[0]=Membrane_Node_Position[ Membrane_triangle_list[i][0]][0];
        xyz[1]=Membrane_Node_Position[ Membrane_triangle_list[i][0]][1];
        xyz[2]=Membrane_Node_Position[ Membrane_triangle_list[i][0]][2];
        
        crossvector(ABxAC,AB,AC);
        //        Throughout the code the ABC vertexes of the membrane triangles are taken as A=Membrane_triangle_list[][0], B=Membrane_triangle_list[][1], C=Membrane_triangle_list[][2]. Also we often use the ABxAC cross product and we want the triangles on the membrane to point out of the cell. At the beginning of the cell construction, the centre of the cell is the same as the origin. So for the ABxAC product to point outwards, the inner product of the position of the triangle and the ABxAC should be positive. We put +/- 1 in the 'Membrane_Normal_direction[][1]' list for each triangle and define the normal direction of each triangle as Membrane_Normal_direction[i][1]*ABxAC that will always be positive, hence pointing out of the cell.
        if(innerproduct(ABxAC,xyz)>0 )
        {
            Membrane_Normal_direction[i][1]=+1;
            //cout<<"pos"<<endl;
        }
        if(innerproduct(ABxAC,xyz)<0 )
        {
            Membrane_Normal_direction[i][1]=-1;
            //cout<<"min"<<endl;
        }
        if(innerproduct(ABxAC,xyz )*Membrane_Normal_direction[i][1]<0)
        {
            cout <<"error 01 \nThis error was generated in the 'Membrane_Normal_direction_Identifier' function."<<innerproduct(ABxAC,xyz )*Membrane_Normal_direction[i][1]<<endl;
        }
        
    } // END OF: for(  int i=0;i<Membrane_num_of_Triangles;i++  )
}// END OF: Membrane_Normal_direction_Identifier function




void Outer_Membrane_Identifier(int Membrane_Normal_direction[Membrane_num_of_Triangles][2] , int Membrane_triangle_list[Membrane_num_of_Triangles][3], int  Outer_Membrane_num_of_triangles, int &Outer_Membrane_num_of_Nodes)
{
    
    
    int temp_Membrane_normal_direction[Membrane_num_of_Triangles][2];
    int temp_Membrane_triangle_list[Membrane_num_of_Triangles][3];
    //    We first make a temporary copy of the 'Membrane_Normal_direction' and the 'Membrane_triangle_list'. we will sort them in the next loop
    for(  int i=0;i<Membrane_num_of_Triangles;i++  ) // membrane or nucleus
    {
        temp_Membrane_normal_direction[i][0]=Membrane_Normal_direction[i][0];
        temp_Membrane_normal_direction[i][1]=Membrane_Normal_direction[i][1];
        
        temp_Membrane_triangle_list[i][0]=Membrane_triangle_list[i][0];
        temp_Membrane_triangle_list[i][1]=Membrane_triangle_list[i][1];
        temp_Membrane_triangle_list[i][2]=Membrane_triangle_list[i][2];
        
    }
    
    int temp_index=0;
    //    we use this index to build the new sorted 'Membrane_triangle_list' and 'Membrane_Normal_direction' lists.
    for(  int i=0;i<Membrane_num_of_Triangles;i++  )
    {
        if( temp_Membrane_normal_direction[i][0]==+1 )
        {
            Membrane_Normal_direction[temp_index][0] = temp_Membrane_normal_direction[i][0];
            Membrane_Normal_direction[temp_index][1] = temp_Membrane_normal_direction[i][1];
            
            Membrane_triangle_list[temp_index][0]=temp_Membrane_triangle_list[i][0];
            Membrane_triangle_list[temp_index][1]=temp_Membrane_triangle_list[i][1];
            Membrane_triangle_list[temp_index][2]=temp_Membrane_triangle_list[i][2];
            
            temp_index++;
        }
    }
    //    After we are done with the triangles on the outer membrane, we will continue to the nucleus triangles.
    for(  int i=0;i<Membrane_num_of_Triangles;i++  )
    {
        if( temp_Membrane_normal_direction[i][0]==-1)
        {
            Membrane_Normal_direction[temp_index][0] = temp_Membrane_normal_direction[i][0];
            Membrane_Normal_direction[temp_index][1] = temp_Membrane_normal_direction[i][1];
            
            Membrane_triangle_list[temp_index][0]=temp_Membrane_triangle_list[i][0];
            Membrane_triangle_list[temp_index][1]=temp_Membrane_triangle_list[i][1];
            Membrane_triangle_list[temp_index][2]=temp_Membrane_triangle_list[i][2];
            
            temp_index++;
        }
    }
    
    //    Now we proceed to counting the number of nodes on the outer membrane by combing through the triangles on the outer membrane.
    int counter=0;
    int temp_node;
    bool found_new_node=true;
    int elemets[Outer_Membrane_num_of_triangles]; // maximum posssible size is   Outer_Membrane_num_of_triangles! butt its much less!
    
    elemets[0]=Membrane_triangle_list[0][0];
    counter++;
    
    for(  int i=0;i<Outer_Membrane_num_of_triangles;i++  ) // membrane or nucleus
    {
        for(int j=0;j<3;j++)
        {
            found_new_node=true;
            temp_node=Membrane_triangle_list[i][j];
            for(int k=0;k<counter;k++)
            {
                if( temp_node == elemets[k]  )
                {
                    found_new_node=false;
                }
            }
            if(found_new_node==true)
            {
                elemets[counter]=temp_node;
                counter++;
            }
        }
    }
    Outer_Membrane_num_of_Nodes=counter;
    cout << "Outer Membrane num of Nodes= \t"<<counter<<"\nMembrane + Nucleus # of triangles: "<< Outer_Membrane_num_of_triangles <<endl;
}



void Membrane_and_Nucleus_Node_list_builder(double Membrane_Node_Position [Membrane_num_of_Nodes][3],int Nucleus_Membrane_list_of_Nodes[],int Outer_Membrane_list_of_Nodes[], int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Outer_Membrane_num_of_triangles)
{
    int counter=0;
    int temp_node;
    bool found_new_node=true;
    
    //The first node in the first triangle in the 'Membrane_triangle_list' is obviously the first member of the 'Outer_Membrane_list_of_Nodes'.
    Outer_Membrane_list_of_Nodes[0]=Membrane_triangle_list[0][0];
    counter++;
    for(  int i=0;i<Outer_Membrane_num_of_triangles;i++  ) // membrane or nucleus
    {
        for(int j=0;j<3;j++)
        {
            found_new_node=true;
            temp_node=Membrane_triangle_list[i][j];
            for(int k=0;k<counter;k++)
            {
                if( temp_node == Outer_Membrane_list_of_Nodes[k]  )
                {
                    found_new_node=false;
                    break;
                }
            }
            if(found_new_node==true)
            {
                Outer_Membrane_list_of_Nodes[counter]=temp_node;
                counter++;
            }
        }
    }
    //Now that we have identified all of the outer membrane nodes, we will continue to the neucleus.
    counter=0;
    found_new_node=true;
    
    Nucleus_Membrane_list_of_Nodes[0]=Membrane_triangle_list[Outer_Membrane_num_of_triangles][0];
    counter++;
    
    for(  int i=Outer_Membrane_num_of_triangles;i<Membrane_num_of_Triangles;i++  ) // membrane or nucleus
    {
        for(int j=0;j<3;j++)
        {
            found_new_node=true;
            temp_node=Membrane_triangle_list[i][j];
            for(int k=0;k<counter;k++)
            {
                if( temp_node == Nucleus_Membrane_list_of_Nodes[k]  )
                {
                    found_new_node=false;
                    break;
                }
            }
            if(found_new_node==true)
            {
                Nucleus_Membrane_list_of_Nodes[counter]=temp_node;
                counter++;
            }
        }
    }
    
}

int Membrane_triangle_pair_counter( int Membrane_triangle_list[Membrane_num_of_Triangles][3])
{
    //In this function we count the total number of triangles that have a common edge (we count them twice, hence report half the number at the end).
    int temp_triangle_node_A, temp_triangle_node_B, temp_triangle_node_C;
    int triangle_pairs=0;  // This counts the number of triangle pairs that have an edge in common.
    for(int i=0 ;i<Membrane_num_of_Triangles;i++)  // who are neighbors??
    {
        temp_triangle_node_A=Membrane_triangle_list[i][0];  // read the tree lable number of nodes  of every triangle
        temp_triangle_node_B=Membrane_triangle_list[i][1];
        temp_triangle_node_C=Membrane_triangle_list[i][2];
        
        for(int j=0;j<Membrane_num_of_Triangles;j++)
        {
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_A  &  Membrane_triangle_list[j][1]==temp_triangle_node_B  & Membrane_triangle_list[j][2]!=temp_triangle_node_C ){
                triangle_pairs++;
            }
            if     ( Membrane_triangle_list[j][0]==temp_triangle_node_B  &  Membrane_triangle_list[j][1]==temp_triangle_node_A  & Membrane_triangle_list[j][2]!=temp_triangle_node_C ){
                triangle_pairs++;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_A  &  Membrane_triangle_list[j][1]!=temp_triangle_node_C  & Membrane_triangle_list[j][2]==temp_triangle_node_B ){
                triangle_pairs++;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_B  &  Membrane_triangle_list[j][1]!=temp_triangle_node_C  & Membrane_triangle_list[j][2]==temp_triangle_node_A ){
                triangle_pairs++;
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_C  &  Membrane_triangle_list[j][1]==temp_triangle_node_A  & Membrane_triangle_list[j][2]==temp_triangle_node_B ){
                triangle_pairs++;
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_C  &  Membrane_triangle_list[j][1]==temp_triangle_node_B  & Membrane_triangle_list[j][2]==temp_triangle_node_A ){
                triangle_pairs++;
            }
            // neibors of temp_triangle_node_B-temp_triangle_node_C :
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_B  &  Membrane_triangle_list[j][1]==temp_triangle_node_C  & Membrane_triangle_list[j][2]!=temp_triangle_node_A ){
                triangle_pairs++;
            }
            if     ( Membrane_triangle_list[j][0]==temp_triangle_node_C  &  Membrane_triangle_list[j][1]==temp_triangle_node_B  & Membrane_triangle_list[j][2]!=temp_triangle_node_A ){
                triangle_pairs++;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_B  &  Membrane_triangle_list[j][1]!=temp_triangle_node_A  & Membrane_triangle_list[j][2]==temp_triangle_node_C ){
                triangle_pairs++;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_C  &  Membrane_triangle_list[j][1]!=temp_triangle_node_A  & Membrane_triangle_list[j][2]==temp_triangle_node_B ){
                triangle_pairs++;
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_A  &  Membrane_triangle_list[j][1]==temp_triangle_node_B  & Membrane_triangle_list[j][2]==temp_triangle_node_C ){
                triangle_pairs++;
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_A  &  Membrane_triangle_list[j][1]==temp_triangle_node_C  & Membrane_triangle_list[j][2]==temp_triangle_node_B ){
                triangle_pairs++;
            }
            // neibors of temp_triangle_node_C-temp_triangle_node_A :
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_C  &  Membrane_triangle_list[j][1]==temp_triangle_node_A  & Membrane_triangle_list[j][2]!=temp_triangle_node_B ){
                triangle_pairs++;
            }
            if     ( Membrane_triangle_list[j][0]==temp_triangle_node_A  &  Membrane_triangle_list[j][1]==temp_triangle_node_C  & Membrane_triangle_list[j][2]!=temp_triangle_node_B ){
                triangle_pairs++;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_C  &  Membrane_triangle_list[j][1]!=temp_triangle_node_B  & Membrane_triangle_list[j][2]==temp_triangle_node_A ){
                triangle_pairs++;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_A  &  Membrane_triangle_list[j][1]!=temp_triangle_node_B  & Membrane_triangle_list[j][2]==temp_triangle_node_C ){
                triangle_pairs++;
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_B  &  Membrane_triangle_list[j][1]==temp_triangle_node_C  & Membrane_triangle_list[j][2]==temp_triangle_node_A ){
                triangle_pairs++;
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_B  &  Membrane_triangle_list[j][1]==temp_triangle_node_A  & Membrane_triangle_list[j][2]==temp_triangle_node_C ){
                triangle_pairs++;
            }
        }
    }
    return triangle_pairs/2;
}


void Membrane_Triangle_Pair_Identifier(int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Triangle_Pair_Nodes[][4], int Membrane_num_of_Triangle_Pairs){
    
    int temp_triangle_node_A, temp_triangle_node_B, temp_triangle_node_C, temp_triangle_node_D, neighbour=0, neighbour_indicator;
    int triangle_pairs=0;
    int temp[4][2*Membrane_num_of_Triangle_Pairs];
    int temp2[4][2*Membrane_num_of_Triangle_Pairs];
    
    
    for(int i=0 ;i<Membrane_num_of_Triangles;i++)
    {
        temp_triangle_node_A=Membrane_triangle_list[i][0];  // read the tree lable number of nodes  of every triangle
        temp_triangle_node_B=Membrane_triangle_list[i][1];
        temp_triangle_node_C=Membrane_triangle_list[i][2];
        neighbour_indicator=0;
        //        Indicates the existence of Node neighbour for a node pair (other than the membrane of the triangle):
        //        neighbour_indicator=0 No Node Pairs; neighbour_indicator=1, for temp_triangle_node_A-temp_triangle_node_B; neighbour_indicator=2, for temp_triangle_node_B-temp_triangle_node_C; And neighbour_indicator=3 for temp_triangle_node_C-temp_triangle_node_A.
        for(int j=0;j<Membrane_num_of_Triangles;j++)
        {
            //************************** finding neighbours **************************
            // neibours of temp_triangle_node_A-temp_triangle_node_B:
            if ( Membrane_triangle_list[j][0]==temp_triangle_node_A  &  Membrane_triangle_list[j][1]==temp_triangle_node_B  & Membrane_triangle_list[j][2]!=temp_triangle_node_C ){
                neighbour=Membrane_triangle_list[j][2];
                neighbour_indicator=1;
            }
            if     ( Membrane_triangle_list[j][0]==temp_triangle_node_B  &  Membrane_triangle_list[j][1]==temp_triangle_node_A  & Membrane_triangle_list[j][2]!=temp_triangle_node_C )
            {
                neighbour=Membrane_triangle_list[j][2];
                neighbour_indicator=1;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_A  &  Membrane_triangle_list[j][1]!=temp_triangle_node_C  & Membrane_triangle_list[j][2]==temp_triangle_node_B ){
                neighbour=Membrane_triangle_list[j][1];
                neighbour_indicator=1;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_B  &  Membrane_triangle_list[j][1]!=temp_triangle_node_C  & Membrane_triangle_list[j][2]==temp_triangle_node_A ){
                neighbour=Membrane_triangle_list[j][1];
                neighbour_indicator=1;
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_C  &  Membrane_triangle_list[j][1]==temp_triangle_node_A  & Membrane_triangle_list[j][2]==temp_triangle_node_B ){
                neighbour=Membrane_triangle_list[j][0];
                neighbour_indicator=1;
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_C  &  Membrane_triangle_list[j][1]==temp_triangle_node_B  & Membrane_triangle_list[j][2]==temp_triangle_node_A ){
                neighbour=Membrane_triangle_list[j][0];
                neighbour_indicator=1;
            }
            // neibors of temp_triangle_node_B-temp_triangle_node_C :
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_B  &  Membrane_triangle_list[j][1]==temp_triangle_node_C  & Membrane_triangle_list[j][2]!=temp_triangle_node_A ){
                neighbour=Membrane_triangle_list[j][2];
                neighbour_indicator=2;
            }
            if     ( Membrane_triangle_list[j][0]==temp_triangle_node_C  &  Membrane_triangle_list[j][1]==temp_triangle_node_B  & Membrane_triangle_list[j][2]!=temp_triangle_node_A ){
                neighbour=Membrane_triangle_list[j][2];
                neighbour_indicator=2;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_B  &  Membrane_triangle_list[j][1]!=temp_triangle_node_A  & Membrane_triangle_list[j][2]==temp_triangle_node_C ){
                neighbour=Membrane_triangle_list[j][1];
                neighbour_indicator=2;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_C  &  Membrane_triangle_list[j][1]!=temp_triangle_node_A  & Membrane_triangle_list[j][2]==temp_triangle_node_B )
            {
                neighbour=Membrane_triangle_list[j][1];
                neighbour_indicator=2;
                
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_A  &  Membrane_triangle_list[j][1]==temp_triangle_node_B  & Membrane_triangle_list[j][2]==temp_triangle_node_C ){
                neighbour=Membrane_triangle_list[j][0];
                neighbour_indicator=2;
                
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_A  &  Membrane_triangle_list[j][1]==temp_triangle_node_C  & Membrane_triangle_list[j][2]==temp_triangle_node_B ){
                neighbour=Membrane_triangle_list[j][0];
                neighbour_indicator=2;
            }
            // neibors of temp_triangle_node_C-temp_triangle_node_A :
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_C  &  Membrane_triangle_list[j][1]==temp_triangle_node_A  & Membrane_triangle_list[j][2]!=temp_triangle_node_B ){
                neighbour=Membrane_triangle_list[j][2];
                neighbour_indicator=3;
                
            }
            if     ( Membrane_triangle_list[j][0]==temp_triangle_node_A  &  Membrane_triangle_list[j][1]==temp_triangle_node_C  & Membrane_triangle_list[j][2]!=temp_triangle_node_B ){
                neighbour=Membrane_triangle_list[j][2];
                neighbour_indicator=3;
                
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_C  &  Membrane_triangle_list[j][1]!=temp_triangle_node_B  & Membrane_triangle_list[j][2]==temp_triangle_node_A ){
                neighbour=Membrane_triangle_list[j][1];
                neighbour_indicator=3;
            }
            if      ( Membrane_triangle_list[j][0]==temp_triangle_node_A  &  Membrane_triangle_list[j][1]!=temp_triangle_node_B  & Membrane_triangle_list[j][2]==temp_triangle_node_C ){
                neighbour=Membrane_triangle_list[j][1];
                neighbour_indicator=3;
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_B  &  Membrane_triangle_list[j][1]==temp_triangle_node_C  & Membrane_triangle_list[j][2]==temp_triangle_node_A ){
                neighbour=Membrane_triangle_list[j][0];
                neighbour_indicator=3;
            }
            if      ( Membrane_triangle_list[j][0]!=temp_triangle_node_B  &  Membrane_triangle_list[j][1]==temp_triangle_node_A  & Membrane_triangle_list[j][2]==temp_triangle_node_C ){
                neighbour=Membrane_triangle_list[j][0];
                neighbour_indicator=3;
            }
            
            if(neighbour_indicator!=0)  //  to speed up  the programme we first check if we have found a neighbour or not
            {
                // note that temp_triangle_node_A-temp_triangle_node_B-temp_triangle_node_C-neighbour  are 4 point of two triangle wich will interact
                if(neighbour_indicator==1)
                {
                    temp[0] [triangle_pairs]=temp_triangle_node_A;
                    temp[1] [triangle_pairs]=temp_triangle_node_B;
                    temp[2] [triangle_pairs]=temp_triangle_node_C;
                    temp[3] [triangle_pairs]=neighbour;
                    
                }
                if(neighbour_indicator==2)
                {
                    
                    temp[0] [triangle_pairs]=temp_triangle_node_B;
                    temp[1] [triangle_pairs]=temp_triangle_node_C;
                    temp[2] [triangle_pairs]=temp_triangle_node_A;
                    temp[3] [triangle_pairs]=neighbour;
                }
                if(neighbour_indicator==3)
                {
                    
                    temp[0] [triangle_pairs]=temp_triangle_node_C;
                    temp[1] [triangle_pairs]=temp_triangle_node_A;
                    temp[2] [triangle_pairs]=temp_triangle_node_B;
                    temp[3] [triangle_pairs]=neighbour;
                }
                triangle_pairs++;
            }
            neighbour_indicator=0;
        }
    }
    
    for(int i=0;i<2*Membrane_num_of_Triangle_Pairs;i++)//saving temp in temp2
    {
        for(int j=0;j<4;j++)
        {
            temp2[j][i]=temp[j][i];
        }
    }
    
    for(int abc=0;abc<2*Membrane_num_of_Triangle_Pairs;abc++)// sorting temp
    {
        if( temp[0] [abc] > temp[1] [abc]   )
        {
            temp_triangle_node_A=temp[0] [abc];
            temp_triangle_node_B=temp[1] [abc];
            temp[0] [abc]=temp_triangle_node_B;
            temp[1] [abc]=temp_triangle_node_A;
        }
        if( temp[0] [abc] > temp[2] [abc]   )
        {
            temp_triangle_node_A=temp[0] [abc];
            temp_triangle_node_B=temp[2] [abc];
            temp[0] [abc]=temp_triangle_node_B;
            temp[2] [abc]=temp_triangle_node_A;
        }
        if( temp[0] [abc] > temp[3] [abc]   )
        {
            temp_triangle_node_A=temp[0] [abc];
            temp_triangle_node_B=temp[3] [abc];
            temp[0] [abc]=temp_triangle_node_B;
            temp[3] [abc]=temp_triangle_node_A;
        }
        
        if( temp[1] [abc] > temp[2] [abc]   )
        {
            temp_triangle_node_A=temp[1] [abc];
            temp_triangle_node_B=temp[2] [abc];
            temp[1] [abc]=temp_triangle_node_B;
            temp[2] [abc]=temp_triangle_node_A;
        }
        
        if( temp[1] [abc] > temp[3] [abc]   )
        {
            temp_triangle_node_A=temp[1] [abc];
            temp_triangle_node_B=temp[3] [abc];
            temp[1] [abc]=temp_triangle_node_B;
            temp[3] [abc]=temp_triangle_node_A;
        }
        
        if( temp[2] [abc] > temp[3] [abc]   )
        {
            temp_triangle_node_A=temp[2] [abc];
            temp_triangle_node_B=temp[3] [abc];
            temp[2] [abc]=temp_triangle_node_B;
            temp[3] [abc]=temp_triangle_node_A;
        }
    }
    
    for(int abc=0;abc<2*Membrane_num_of_Triangle_Pairs;abc++)
    {
        if(temp[0] [abc]!=-1)
        {
            temp_triangle_node_A=temp[0] [abc];
            temp_triangle_node_B=temp[1] [abc];
            temp_triangle_node_C=temp[2] [abc];
            temp_triangle_node_D=temp[3] [abc];
            for(int cab=0;cab<2*Membrane_num_of_Triangle_Pairs;cab++)
            {
                if( temp_triangle_node_A==temp[0] [cab] &   temp_triangle_node_B==temp[1] [cab]   &   temp_triangle_node_C==temp[2] [cab]   &   temp_triangle_node_D==temp[3] [cab] & abc!=cab  )
                {
                    temp[0] [cab] =-1;
                }
                
            }
            
        }
        
    }
    
    
    int temp_int=0;
    for(int abc=0;abc<2*Membrane_num_of_Triangle_Pairs;abc++)
    {
        if( temp[0] [abc]!=-1)
        {
            Membrane_Triangle_Pair_Nodes[temp_int][0]=temp2[0] [abc];
            Membrane_Triangle_Pair_Nodes[temp_int][1]=temp2[1] [abc];
            Membrane_Triangle_Pair_Nodes[temp_int][2]=temp2[2] [abc];
            Membrane_Triangle_Pair_Nodes[temp_int][3]=temp2[3] [abc];
            
            temp_int++;
        }
    }
}

int Membrane_num_of_Node_Pair_Counter(int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Outer_Membrane_num_of_triangles, int &Outer_Membrane_num_of_Node_Pairs)
{
    Outer_Membrane_num_of_Node_Pairs=0;
    int bondslist[10*Membrane_num_of_Nodes][2];
    for (int j=0 ; j< Membrane_num_of_Nodes*10 ; j++)
    {
        bondslist[j][0]=-1;
        bondslist[j][1]=-1;
    }
    
    int temp_Membrane_num_of_Node_Pairs=0;
    int temp_Membrane_triangle_Node_A, temp_Membrane_triangle_Node_B, temp_Membrane_triangle_Node_C;
    
    int repeatednumber1=0;
    int repeatednumber2=0;
    int repeatednumber3=0;
    
    for(int i=0;i<Membrane_num_of_Triangles;i++)
    {
        temp_Membrane_triangle_Node_A= Membrane_triangle_list[i][0];
        temp_Membrane_triangle_Node_B= Membrane_triangle_list[i][1];
        temp_Membrane_triangle_Node_C= Membrane_triangle_list[i][2];
        
        for(int j=0;j<10*Membrane_num_of_Nodes;j++)
        {
            if(  ( bondslist[j][0]==temp_Membrane_triangle_Node_A &  bondslist[j][1]==temp_Membrane_triangle_Node_B )  || ( bondslist[j][0]==temp_Membrane_triangle_Node_B &  bondslist[j][1]==temp_Membrane_triangle_Node_A )    )
            {
                repeatednumber1=1;
            }
            
            if(  ( bondslist[j][0]==temp_Membrane_triangle_Node_B &  bondslist[j][1]==temp_Membrane_triangle_Node_C )  || ( bondslist[j][0]==temp_Membrane_triangle_Node_C &  bondslist[j][1]==temp_Membrane_triangle_Node_B )    )
            {
                repeatednumber2=1;
            }
            
            if(  ( bondslist[j][0]==temp_Membrane_triangle_Node_A &  bondslist[j][1]==temp_Membrane_triangle_Node_C )  || ( bondslist[j][0]==temp_Membrane_triangle_Node_C &  bondslist[j][1]==temp_Membrane_triangle_Node_A )    )
            {
                repeatednumber3=1;
            }
        }
        
        if(repeatednumber1==0)
        {
            bondslist[temp_Membrane_num_of_Node_Pairs][0]=temp_Membrane_triangle_Node_A;       //note that first node store in i and second in i+numofbonds  ---Membrane_Node_Pair_list[2*numofbonds]
            bondslist[temp_Membrane_num_of_Node_Pairs][1]=temp_Membrane_triangle_Node_B;
            temp_Membrane_num_of_Node_Pairs++;
            if(i<Outer_Membrane_num_of_triangles)
            {  Outer_Membrane_num_of_Node_Pairs++; }
        }
        
        if(repeatednumber2==0)
        {
            bondslist[temp_Membrane_num_of_Node_Pairs][0]=temp_Membrane_triangle_Node_B;       //note that first node store in i and second in i+numofbonds  ---Membrane_Node_Pair_list[2*numofbonds]
            bondslist[temp_Membrane_num_of_Node_Pairs][1]=temp_Membrane_triangle_Node_C;
            temp_Membrane_num_of_Node_Pairs++;
            if(i<Outer_Membrane_num_of_triangles)
            {  Outer_Membrane_num_of_Node_Pairs++; }
        }
        
        if(repeatednumber3==0)
        {
            bondslist[temp_Membrane_num_of_Node_Pairs][0]=temp_Membrane_triangle_Node_A;       //note that first node store in i and second in i+numofbonds  ---Membrane_Node_Pair_list[2*numofbonds]
            bondslist[temp_Membrane_num_of_Node_Pairs][1]=temp_Membrane_triangle_Node_C;
            temp_Membrane_num_of_Node_Pairs++;
            if(i<Outer_Membrane_num_of_triangles)
            {  Outer_Membrane_num_of_Node_Pairs++; }
        }
        
        repeatednumber1=0;
        repeatednumber2=0;
        repeatednumber3=0;
    }
    
    //    cout<<"Outer_Membrane_num_of_Node_Pairs_test=\t"<<Outer_Membrane_num_of_Node_Pairs<<endl;
    //    exit (EXIT_FAILURE);
    //*******************************************************************************************************
    /*BUG
     |---\   |    |  /---\
     |    |  |    |  |
     |---<   |    |  |  -\
     |    |  |    |  |   |
     |---/   \----/  \---/
     */
    //*******************************************************************************************************
    //***************** Potential BUG: This counter gives exactley the same number as the *******************
    //***************** 'Membrane_triangle_pair_counter'. I think we can use that number, *******************
    //***************** but we have to check it first *******************************************************
    //*******************************************************************************************************
    return temp_Membrane_num_of_Node_Pairs;
    
}
void Membrane_num_of_Node_Pair_Counter_2(int Membrane_Node_Pair_list[][2],int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Outer_Membrane_num_of_triangles, int Membrane_num_of_Node_Pairs)
{
    for (int j=0 ; j< Membrane_num_of_Node_Pairs ; j++)
    {
        Membrane_Node_Pair_list[j][0]=-1;
        Membrane_Node_Pair_list[j][1]=-1;
    }
    
    int temp_Membrane_num_of_Node_Pairs=0;
    int temp_Membrane_triangle_Node_A, temp_Membrane_triangle_Node_B, temp_Membrane_triangle_Node_C;
    
    int repeatednumber1=0;
    int repeatednumber2=0;
    int repeatednumber3=0;
    for(int i=0;i<Membrane_num_of_Triangles;i++)
    {
        temp_Membrane_triangle_Node_A= Membrane_triangle_list[i][0];
        temp_Membrane_triangle_Node_B= Membrane_triangle_list[i][1];
        temp_Membrane_triangle_Node_C= Membrane_triangle_list[i][2];
        
        for(int j=0;j<Membrane_num_of_Node_Pairs;j++)
        {
            if(  ( Membrane_Node_Pair_list[j][0]==temp_Membrane_triangle_Node_A &  Membrane_Node_Pair_list[j][1]==temp_Membrane_triangle_Node_B )  || ( Membrane_Node_Pair_list[j][0]==temp_Membrane_triangle_Node_B &  Membrane_Node_Pair_list[j][1]==temp_Membrane_triangle_Node_A )    )
            {
                repeatednumber1=1;
            }
            
            if(  ( Membrane_Node_Pair_list[j][0]==temp_Membrane_triangle_Node_B &  Membrane_Node_Pair_list[j][1]==temp_Membrane_triangle_Node_C )  || ( Membrane_Node_Pair_list[j][0]==temp_Membrane_triangle_Node_C &  Membrane_Node_Pair_list[j][1]==temp_Membrane_triangle_Node_B )    )
            {
                repeatednumber2=1;
            }
            
            if(  ( Membrane_Node_Pair_list[j][0]==temp_Membrane_triangle_Node_A &  Membrane_Node_Pair_list[j][1]==temp_Membrane_triangle_Node_C )  || ( Membrane_Node_Pair_list[j][0]==temp_Membrane_triangle_Node_C &  Membrane_Node_Pair_list[j][1]==temp_Membrane_triangle_Node_A )    )
            {
                repeatednumber3=1;
            }
        }
        
        if(repeatednumber1==0)
        {
            Membrane_Node_Pair_list[temp_Membrane_num_of_Node_Pairs][0]=temp_Membrane_triangle_Node_A;       //note that first node store in i and second in i+numofbonds  ---Membrane_Node_Pair_list[2*numofbonds]
            Membrane_Node_Pair_list[temp_Membrane_num_of_Node_Pairs][1]=temp_Membrane_triangle_Node_B;
            
            temp_Membrane_num_of_Node_Pairs++;
            
        }
        
        if(repeatednumber2==0)
        {
            Membrane_Node_Pair_list[temp_Membrane_num_of_Node_Pairs][0]=temp_Membrane_triangle_Node_B;       //note that first node store in i and second in i+numofbonds  ---Membrane_Node_Pair_list[2*numofbonds]
            Membrane_Node_Pair_list[temp_Membrane_num_of_Node_Pairs][1]=temp_Membrane_triangle_Node_C;
            
            temp_Membrane_num_of_Node_Pairs++;
            
        }
        
        if(repeatednumber3==0)
        {
            Membrane_Node_Pair_list[temp_Membrane_num_of_Node_Pairs][0]=temp_Membrane_triangle_Node_A;       //note that first node store in i and second in i+numofbonds  ---Membrane_Node_Pair_list[2*numofbonds]
            Membrane_Node_Pair_list[temp_Membrane_num_of_Node_Pairs][1]=temp_Membrane_triangle_Node_C;
            
            temp_Membrane_num_of_Node_Pairs++;
            
        }
        
        repeatednumber1=0;
        repeatednumber2=0;
        repeatednumber3=0;
    }
    
}


