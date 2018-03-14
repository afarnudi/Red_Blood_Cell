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

void Membrane_constructor(double Membrane_Node_Position[Membrane_num_of_Nodes][3], double Membrane_Node_Velocity [Membrane_num_of_Nodes][3], double Membrane_Node_Force[Membrane_num_of_Nodes][3], int Membrane_triangle_list[Membrane_num_of_Triangles][3])
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



void Membrane_Normal_direction_Identifier( double  Membrane_Node_Position [Membrane_num_of_Nodes][3], int Membrane_triangle_list[Membrane_num_of_Triangles][3])
{
    double AC[3], AB[3], ABxAC[3], xyz[3];
    int Point_A, Point_B, Point_C;
    
    for(  int i=0;i<Membrane_num_of_Triangles;i++  )
    {
        Point_A=Membrane_triangle_list[i][0];
        Point_B=Membrane_triangle_list[i][1];
        Point_C=Membrane_triangle_list[i][2];
        
        AB[0]=Membrane_Node_Position[Point_B][0]-Membrane_Node_Position[Point_A][0];
        AB[1]=Membrane_Node_Position[Point_B][1]-Membrane_Node_Position[Point_A][1];
        AB[2]=Membrane_Node_Position[Point_B][2]-Membrane_Node_Position[Point_A][2];
        
        AC[0]=Membrane_Node_Position[Point_C][0]-Membrane_Node_Position[Point_A][0];
        AC[1]=Membrane_Node_Position[Point_C][1]-Membrane_Node_Position[Point_A][1];
        AC[2]=Membrane_Node_Position[Point_C][2]-Membrane_Node_Position[Point_A][2];
        
        xyz[0]=Membrane_Node_Position[ Membrane_triangle_list[i][0]][0];
        xyz[1]=Membrane_Node_Position[ Membrane_triangle_list[i][0]][1];
        xyz[2]=Membrane_Node_Position[ Membrane_triangle_list[i][0]][2];
        
        crossvector(ABxAC, AB, AC);
        //        Throughout the code the ABC vertexes of the membrane triangles are taken as A=Membrane_triangle_list[][0], B=Membrane_triangle_list[][1], C=Membrane_triangle_list[][2]. Also we often use the ABxAC cross product and we want the triangles on the membrane to point out of the cell. At the beginning of the cell construction, the centre of the cell is the same as the origin. So for the ABxAC product to point outwards, the inner product of the position of the triangle and the ABxAC should be positive. We put +/- 1 in the 'Membrane_Normal_direction[][1]' list for each triangle and define the normal direction of each triangle as Membrane_Normal_direction[i][1]*ABxAC that will always be positive, hence pointing out of the cell.
        
        if(innerproduct(ABxAC, xyz)<0 )
        {
            Membrane_triangle_list[i][1]=Point_C;
            Membrane_triangle_list[i][2]=Point_B;
            //cout<<"min"<<endl;
        }
    } // END OF: for(  int i=0;i<Membrane_num_of_Triangles;i++  )
}// END OF: Membrane_Normal_direction_Identifier function


int Membrane_triangle_pair_counter(int Membrane_triangle_list[Membrane_num_of_Triangles][3])
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


void Membrane_Triangle_Pair_Identifier(int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Triangle_Pair_Nodes[][4], int Membrane_num_of_Triangle_Pairs, vector<vector<int> > &membrane_triangle_pair_list){
    
    int temp_triangle_node_A, temp_triangle_node_B, temp_triangle_node_C, temp_triangle_node_D, neighbour=0, neighbour_indicator;
    int triangle_pairs=0;
    int temp[2*Membrane_num_of_Triangle_Pairs][4];
    int temp2[2*Membrane_num_of_Triangle_Pairs][4];
    
    
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
                membrane_triangle_pair_list[i].push_back(j);
                
                if(neighbour_indicator==1)
                {
                    temp[triangle_pairs][0] =temp_triangle_node_A;
                    temp[triangle_pairs][1] =temp_triangle_node_B;
                    temp[triangle_pairs][2] =temp_triangle_node_C;
                    temp[triangle_pairs][3] =neighbour;
                    
                }
                if(neighbour_indicator==2)
                {
                    
                    temp[triangle_pairs][0] =temp_triangle_node_B;
                    temp[triangle_pairs][1] =temp_triangle_node_C;
                    temp[triangle_pairs][2] =temp_triangle_node_A;
                    temp[triangle_pairs][3] =neighbour;
                }
                if(neighbour_indicator==3)
                {
                    
                    temp[triangle_pairs][0] =temp_triangle_node_C;
                    temp[triangle_pairs][1] =temp_triangle_node_A;
                    temp[triangle_pairs][2] =temp_triangle_node_B;
                    temp[triangle_pairs][3] =neighbour;
                }
                triangle_pairs++;
            }
            neighbour_indicator=0;
        }
    }
    
    for (int i=0; i<Membrane_num_of_Triangles; i++) {
        if (membrane_triangle_pair_list[i].size()!=3) {
            cout<<"There is an error in the 'Membrane_Triangle_Pair_Identifier' function. This error indicates that there is a membrane triangle that has more/less than 3 triangle neighbours."<<endl;
        }
    }
    
    for(int i=0;i<2*Membrane_num_of_Triangle_Pairs;i++)//saving temp in temp2
    {
        for(int j=0;j<4;j++)
        {
            temp2[i][j]=temp[i][j];
        }
    }
    
    for(int abc=0;abc<2*Membrane_num_of_Triangle_Pairs;abc++)// sorting temp
    {
        if( temp[abc][0]  > temp[abc][1]    )
        {
            temp_triangle_node_A=temp[abc][0] ;
            temp_triangle_node_B=temp[abc][1] ;
            temp[abc][0] =temp_triangle_node_B;
            temp[abc][1] =temp_triangle_node_A;
        }
        if( temp[abc][0]  > temp[abc][2]    )
        {
            temp_triangle_node_A=temp[abc][0] ;
            temp_triangle_node_B=temp[abc][2] ;
            temp[abc][0] =temp_triangle_node_B;
            temp[abc][2] =temp_triangle_node_A;
        }
        if( temp[abc][0]  > temp[abc][3]    )
        {
            temp_triangle_node_A=temp[abc][0] ;
            temp_triangle_node_B=temp[abc][3] ;
            temp[abc][0] =temp_triangle_node_B;
            temp[abc][3] =temp_triangle_node_A;
        }
        
        if( temp[abc][1]  > temp[abc][2]    )
        {
            temp_triangle_node_A=temp[abc][1] ;
            temp_triangle_node_B=temp[abc][2] ;
            temp[abc][1] =temp_triangle_node_B;
            temp[abc][2] =temp_triangle_node_A;
        }
        
        if( temp[abc][1]  > temp[abc][3]    )
        {
            temp_triangle_node_A=temp[abc][1] ;
            temp_triangle_node_B=temp[abc][3] ;
            temp[abc][1] =temp_triangle_node_B;
            temp[abc][3] =temp_triangle_node_A;
        }
        
        if( temp[abc][2]  > temp[abc][3]    )
        {
            temp_triangle_node_A=temp[abc][2] ;
            temp_triangle_node_B=temp[abc][3] ;
            temp[abc][2] =temp_triangle_node_B;
            temp[abc][3] =temp_triangle_node_A;
        }
    }
    
    for(int abc=0;abc<2*Membrane_num_of_Triangle_Pairs;abc++)
    {
        if(temp[abc][0] !=-1)
        {
            temp_triangle_node_A=temp[abc][0] ;
            temp_triangle_node_B=temp[abc][1] ;
            temp_triangle_node_C=temp[abc][2] ;
            temp_triangle_node_D=temp[abc][3] ;
            for(int cab=0;cab<2*Membrane_num_of_Triangle_Pairs;cab++)
            {
                if( temp_triangle_node_A==temp[cab][0]  &   temp_triangle_node_B==temp[cab][1]    &   temp_triangle_node_C==temp[cab][2]    &   temp_triangle_node_D==temp[cab][3]  & abc!=cab  )
                {
                    temp[cab][0]  =-1;
                }
                
            }
            
        }
        
    }
    
    
    int temp_int=0;
    for(int abc=0;abc<2*Membrane_num_of_Triangle_Pairs;abc++)
    {
        if( temp[abc][0] !=-1)
        {
            Membrane_Triangle_Pair_Nodes[temp_int][0]=temp2[abc][0] ;
            Membrane_Triangle_Pair_Nodes[temp_int][1]=temp2[abc][1] ;
            Membrane_Triangle_Pair_Nodes[temp_int][2]=temp2[abc][2] ;
            Membrane_Triangle_Pair_Nodes[temp_int][3]=temp2[abc][3] ;
            
            temp_int++;
        }
    }
}

int Membrane_num_of_Node_Pair_Counter(int Membrane_triangle_list[Membrane_num_of_Triangles][3])
{
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
            
        }
        
        if(repeatednumber2==0)
        {
            bondslist[temp_Membrane_num_of_Node_Pairs][0]=temp_Membrane_triangle_Node_B;       //note that first node store in i and second in i+numofbonds  ---Membrane_Node_Pair_list[2*numofbonds]
            bondslist[temp_Membrane_num_of_Node_Pairs][1]=temp_Membrane_triangle_Node_C;
            temp_Membrane_num_of_Node_Pairs++;
            
        }
        
        if(repeatednumber3==0)
        {
            bondslist[temp_Membrane_num_of_Node_Pairs][0]=temp_Membrane_triangle_Node_A;       //note that first node store in i and second in i+numofbonds  ---Membrane_Node_Pair_list[2*numofbonds]
            bondslist[temp_Membrane_num_of_Node_Pairs][1]=temp_Membrane_triangle_Node_C;
            temp_Membrane_num_of_Node_Pairs++;
            
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

void Membrane_num_of_Node_Pair_Counter_2(int Membrane_Node_Pair_list[][2], int Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_num_of_Node_Pairs)
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


void Membrane_Force_Calculator (double Membrane_Node_Position[Membrane_num_of_Nodes][3],double Membrane_Node_Velocity[Membrane_num_of_Nodes][3],double Membrane_Node_Force [Membrane_num_of_Nodes][3],int Membrane_Node_Pair_list[][2],int Membrane_Triangle_Pair_Nodes[][4],double &Total_Potential_Energy, int Membrane_num_of_Triangle_Pairs, int Membrane_num_of_Node_Pairs)
{
    double le0,le1,lmax,lmin;
    double deltax,deltay,deltaz,temp_Node_distance,temp_force;
    int pos1,pos2,pos3,pos4;  // to making calculation of surface force easier
    double temp_potential_energy = 0.0;
    double temp_p1[3], temp_p2[3], temp_p3[3], temp_p4[3];
    double  N1[3], N2[3], N3[3], p3p1[3], p3p2[3], p4p2[3], p4p1[3], Ep2p1[3], sinus, F0, F1[3], F2[3], F3[3], F4[3];// for exmple p3p1 is p3-p1 and ....
    
    /// calculate network force:
    int temp_Node_A, temp_Node_B;
    le0=1.15000*Node_radius;
    lmax=1.33000*Node_radius;
    le1=0.85000*Node_radius;
    lmin=0.67000*Node_radius;
    Total_Potential_Energy=0.0;
    
    for (int k=0 ; k< Membrane_num_of_Node_Pairs ; k++)
    {
        temp_Node_B=Membrane_Node_Pair_list[k][0];
        temp_Node_A=Membrane_Node_Pair_list[k][1];
        
        deltax=Membrane_Node_Position[temp_Node_A][0]-Membrane_Node_Position[temp_Node_B][0];
        deltay=Membrane_Node_Position[temp_Node_A][1]-Membrane_Node_Position[temp_Node_B][1];
        deltaz=Membrane_Node_Position[temp_Node_A][2]-Membrane_Node_Position[temp_Node_B][2];
        
        
        temp_Node_distance=sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz);
        temp_force=0.0;
        double temp_exp_le0=exp(1.0/(le0-temp_Node_distance));
        double temp_exp_le1=exp(1.0/(temp_Node_distance-le1));
        //*******************************************************************************************************
        /*BUG
         |---\   |    |  /---\
         |    |  |    |  |
         |---<   |    |  |  -\
         |    |  |    |  |   |
         |---/   \----/  \---/
         */
        //*******************************************************************************************************
        //***************** Potential BUG: F=-dU/dr but in many cases I cannot determin wheather ****************
        //***************** the '-' has been implemented or not. Since the potential energy is   ****************
        //***************** never used in the code it does not a threat. ****************************************
        //*******************************************************************************************************
        
        if(temp_Node_distance >le1  & temp_Node_distance < le0 )  //free zone
        {
            temp_potential_energy=0 ; // free zone
        }
        
        if(temp_Node_distance > le0  & temp_Node_distance <lmax )  //bondforce
        {
            temp_force = (Membrane_spring_coefficient*temp_exp_le0/(lmax-temp_Node_distance))*( 1.0/(lmax-temp_Node_distance) +  1.0/((le0-temp_Node_distance)*(le0-temp_Node_distance)));
            temp_potential_energy= Membrane_spring_coefficient*temp_exp_le0/(lmax-temp_Node_distance);
            
        }
        
        if(temp_Node_distance < le1   &  temp_Node_distance > lmin  )  // repulsive force
        {
            temp_force= -(Membrane_spring_coefficient*temp_exp_le1/(temp_Node_distance-lmin))*( 1.0/(temp_Node_distance-lmin) + 1.0/((temp_Node_distance-le1)*(temp_Node_distance-le1)));                 // force on i th from j
            temp_potential_energy= Membrane_spring_coefficient*temp_exp_le1/(temp_Node_distance-lmin);
        }
        /// my cutoff for force amplitute and for avoiding leting particle scape from force trap
        if(temp_force>965.31  || temp_Node_distance>lmax )
        {
            temp_force = 965.31+Membrane_spring_force_cutt_off* ( temp_Node_distance - 1.3280*Node_radius );
            temp_potential_energy=   1.81599  + 965.31 * ( temp_Node_distance - 1.3280*Node_radius )+0.5*Membrane_spring_force_cutt_off * ( temp_Node_distance - 1.3280*Node_radius ) * ( temp_Node_distance - 1.3280*Node_radius );
        }
        
        
        if(temp_force<-1000.05   ||  temp_Node_distance<lmin )
        {
            temp_force =-1000.05-Membrane_spring_force_cutt_off* ( 0.671965*Node_radius - temp_Node_distance );
            temp_potential_energy = 1.85038 + 1005.05 * ( 0.671965*Node_radius - temp_Node_distance )+0.5*Membrane_spring_force_cutt_off*( 0.671965*Node_radius - temp_Node_distance )*( 0.671965*Node_radius - temp_Node_distance );
        }
        
        Total_Potential_Energy += temp_potential_energy;
        
        // implimentation of forces:
        Membrane_Node_Force[temp_Node_A][0] += temp_force*deltax/temp_Node_distance+membrane_damping_coefficient*(Membrane_Node_Velocity[temp_Node_A][0]-Membrane_Node_Velocity[temp_Node_B][0]);
        Membrane_Node_Force[temp_Node_A][1] += temp_force*deltay/temp_Node_distance+membrane_damping_coefficient*(Membrane_Node_Velocity[temp_Node_A][1]-Membrane_Node_Velocity[temp_Node_B][1]);
        Membrane_Node_Force[temp_Node_A][2] += temp_force*deltaz/temp_Node_distance+membrane_damping_coefficient*(Membrane_Node_Velocity[temp_Node_A][2]-Membrane_Node_Velocity[temp_Node_B][2]);
        
        Membrane_Node_Force[temp_Node_B][0] += -temp_force*deltax/temp_Node_distance-membrane_damping_coefficient*(Membrane_Node_Velocity[temp_Node_A][0]-Membrane_Node_Velocity[temp_Node_B][0]); //from j  to i
        Membrane_Node_Force[temp_Node_B][1] += -temp_force*deltay/temp_Node_distance-membrane_damping_coefficient*(Membrane_Node_Velocity[temp_Node_A][1]-Membrane_Node_Velocity[temp_Node_B][1]);
        Membrane_Node_Force[temp_Node_B][2] += -temp_force*deltaz/temp_Node_distance-membrane_damping_coefficient*(Membrane_Node_Velocity[temp_Node_A][2]-Membrane_Node_Velocity[temp_Node_B][2]);
    }
    // End of Membrane Node Pair forces
    
    // Beginning of the  triangle-triangle (bending) force calculations
    for(int i=0 ;i<Membrane_num_of_Triangle_Pairs;i++)  // who are neighbors?
    {
        
        pos1=Membrane_Triangle_Pair_Nodes[i][0];
        pos2=Membrane_Triangle_Pair_Nodes[i][1];
        pos3=Membrane_Triangle_Pair_Nodes[i][2];
        pos4=Membrane_Triangle_Pair_Nodes[i][3];
        
        for (int index=0; index<3; index++) {
            temp_p1[index]=Membrane_Node_Position[pos1][index];
            temp_p2[index]=Membrane_Node_Position[pos2][index];
            temp_p3[index]=Membrane_Node_Position[pos3][index];
            temp_p4[index]=Membrane_Node_Position[pos4][index];
            
            p3p1[index]=temp_p3[index]-temp_p1[index];
            p3p2[index]=temp_p3[index]-temp_p2[index];
            p4p2[index]=temp_p4[index]-temp_p2[index];
            p4p1[index]=temp_p4[index]-temp_p1[index];
            Ep2p1[index]=temp_p2[index]-temp_p1[index];
        }
        
        
        crossvector(N1, p3p1, p3p2);
        crossvector(N2, p4p2, p4p1);
        crossvector(N3, N2, N1);
        sinus=vectorlength(N3)/(vectorlength(N2)*vectorlength(N1));
        
        crossvector(N3, N1, N2);
        double temp_Ep2p1_length=vectorlength(Ep2p1);
        //        if ((sinus- temp_sinus)>0.0001) {
        //            cout<< "oops!"<<endl;
        //            cout<<"temp_sinus="<<temp_sinus<<"\nsinus="<<sinus<<endl;
        //            exit (EXIT_FAILURE);
        //        }
        
        F0 = -(1-2*signbit(innerproduct(N3, Ep2p1)))*Membrane_bending_coefficient*sinus;
        //        cout<<"\nF0="<<F0<<endl;
        //        if( parallelORantiparallel(xpos ) == +1 )
        //        {
        //            F0=-F0;
        //        }
        double temp_N1_length_squared=vectorlength(N1)*vectorlength(N1);
        double temp_N2_length_squared=vectorlength(N2)*vectorlength(N2);
        
        //**************************************************** force calculation
        for (int l=0; l<3; l++) {
            F3[l]= F0 * temp_Ep2p1_length* N1[l]/ temp_N1_length_squared;
            
            F4[l]= F0 * temp_Ep2p1_length* N2[l]/ temp_N2_length_squared;
            
            F1[l]= (F0/temp_Ep2p1_length)*( innerproduct(p3p2,Ep2p1)*N1[l]/temp_N1_length_squared + innerproduct(p4p2,Ep2p1)*N2[l]/temp_N2_length_squared );
            
            F2[l]= (-F0/temp_Ep2p1_length)*( innerproduct(p3p1,Ep2p1)*N1[l]/temp_N1_length_squared + innerproduct(p4p1,Ep2p1)*N2[l]/temp_N2_length_squared );
            
            
            Membrane_Node_Force[pos1][l] += F1[l];
            Membrane_Node_Force[pos2][l] += F2[l];
            Membrane_Node_Force[pos3][l] += F3[l];
            Membrane_Node_Force[pos4][l] += F4[l];
            
        }
        
        //*******************************************************************************************************
        /*BUG
         |---\   |    |  /---\
         |    |  |    |  |
         |---<   |    |  |  -\
         |    |  |    |  |   |
         |---/   \----/  \---/
         */
        //*******************************************************************************************************
        //***************** Potential BUG: Not very sure about the calculations *****************************
        //***************** for 'Total_Potential_Energy. Need to double check *********************************
        //*******************************************************************************************************
        
        Total_Potential_Energy += Membrane_bending_coefficient*(1.0 -  ( innerproduct(N1,N2)/(vectorlength(N1)*vectorlength(N2)  )   ) );
    }  // end of 'for-i'
    //    exit (EXIT_FAILURE);
}///end of function


void ConstantSurfaceForceLocalTriangles(double Membrane_Node_Position[Membrane_num_of_Nodes][3],double Membrane_Node_Force[Membrane_num_of_Nodes][3], int Membrane_triangle_list[Membrane_num_of_Triangles][3])
{
    
    int temp_node_A, temp_node_B, temp_node_C;
    double AB[3], AC[3], temp_AB_length_squared, temp_AC_length_squared, f0, temp_ABxAC[3];
    
    // cout <<surfacearea(x,tri )<<"   s0="<< s0<<endl;
    double s0_i=0.41*Node_radius*Node_radius; //1.732=3^0.5
    double s_i=0;
    
    
    for(  int i=0;i<Membrane_num_of_Triangles;i++)
    {
        temp_node_A=Membrane_triangle_list[i][0];
        temp_node_B=Membrane_triangle_list[i][1];
        temp_node_C=Membrane_triangle_list[i][2];
        
        AB[0]=Membrane_Node_Position[temp_node_B][0]-Membrane_Node_Position[temp_node_A][0];
        AB[1]=Membrane_Node_Position[temp_node_B][1]-Membrane_Node_Position[temp_node_A][1];
        AB[2]=Membrane_Node_Position[temp_node_B][2]-Membrane_Node_Position[temp_node_A][2];
        
        AC[0]=Membrane_Node_Position[temp_node_C][0]-Membrane_Node_Position[temp_node_A][0];
        AC[1]=Membrane_Node_Position[temp_node_C][1]-Membrane_Node_Position[temp_node_A][1];
        AC[2]=Membrane_Node_Position[temp_node_C][2]-Membrane_Node_Position[temp_node_A][2];
        
        crossvector(temp_ABxAC, AB, AC);
        s_i = vectorlength(temp_ABxAC)/2.0;
        f0= K_surfaceConstant_local*(s_i -  s0_i )/2.0*vectorlength(temp_ABxAC);
        
        for (int j=0; j<3; j++) {
            temp_node_A=Membrane_triangle_list[i][j%3];
            temp_node_B=Membrane_triangle_list[i][(j+1)%3];
            temp_node_C=Membrane_triangle_list[i][(j+2)%3];
            
            AB[0]=Membrane_Node_Position[temp_node_B][0]-Membrane_Node_Position[temp_node_A][0];
            AB[1]=Membrane_Node_Position[temp_node_B][1]-Membrane_Node_Position[temp_node_A][1];
            AB[2]=Membrane_Node_Position[temp_node_B][2]-Membrane_Node_Position[temp_node_A][2];
            
            AC[0]=Membrane_Node_Position[temp_node_C][0]-Membrane_Node_Position[temp_node_A][0];
            AC[1]=Membrane_Node_Position[temp_node_C][1]-Membrane_Node_Position[temp_node_A][1];
            AC[2]=Membrane_Node_Position[temp_node_C][2]-Membrane_Node_Position[temp_node_A][2];
            
            temp_AB_length_squared=vectorlength(AB)*vectorlength(AB);
            temp_AC_length_squared=vectorlength(AC)*vectorlength(AC);
            
            Membrane_Node_Force[temp_node_A][0] +=  f0 *2.0* ( -temp_AC_length_squared *  AB[0]  -temp_AB_length_squared *  AC[0] + innerproduct(AB,AC)*  ( AB[0]+AC[0] )   ) ;
            Membrane_Node_Force[temp_node_A][1] +=  f0 *2.0* ( -temp_AC_length_squared *  AB[1]  -temp_AB_length_squared *  AC[1] + innerproduct(AB,AC)*  ( AB[1]+AC[1] )   ) ;
            Membrane_Node_Force[temp_node_A][2] +=  f0 *2.0* ( -temp_AC_length_squared *  AB[2]  -temp_AB_length_squared *  AC[2] + innerproduct(AB,AC)*  ( AB[2]+AC[2] )   ) ;
        }
    }
    
}





