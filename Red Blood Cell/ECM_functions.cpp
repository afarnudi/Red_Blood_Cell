//
//  ECM_functions.cpp
//  Cell_Project
//
//  Created by Ali Farnudi on 22/10/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//

#include "ECM_functions.hpp"
using namespace std;

int  ECM_Node_Pair_Identifier ( )  // gets input file and the output is network
{
    vector<vector<int> > ECM_Node_Pair_list_temp;
    vector<int> Node_Pairs;
    
    int ECM_pyramid_Node_1, ECM_pyramid_Node_2, ECM_pyramid_Node_3, ECM_pyramid_Node_4;
    double temp_double;
    int Num_of_objects;
    
    Node_Pairs.resize(2,-1);
    ECM_Node_Pair_list_temp.push_back(Node_Pairs);
    
    int temp_int_repeated_pair_counter_1=0;
    int temp_int_repeated_pair_counter_2=0;
    int temp_int_repeated_pair_counter_3=0;
    int temp_int_repeated_pair_counter_4=0;
    int temp_int_repeated_pair_counter_5=0;
    int temp_int_repeated_pair_counter_6=0;
    
    ifstream read;
    read.open("ECM");
    
    read>> temp_double;//only sweep not needed elements
    for(int i=0;i<ECM_num_of_Nodes;i++)   // The first 'ECM_num_of_Nodes' lines are the Node coordinates
    {
        read>> temp_double;//only sweep not needed elements
        read>> temp_double;
        read>> temp_double;
        read>> temp_double;
    }
    
    read>> Num_of_objects;
    for(int i=0;i<ECM_Surface_num_of_Triangles;i++)   // After the Node coordinates, we find the triangle coordinate that we do not need at the moment.
    {
        read>> temp_double; //Ignoring these lines
        read>> temp_double;
        read>> temp_double;
        read>> temp_double;
        read>> temp_double;
        read>> temp_double;
        read>> temp_double;
        read>> temp_double;
    }
    
    for(int i=0;i<Num_of_objects-ECM_Surface_num_of_Triangles;i++)    // Now we finally arive at where the triangular pyramid Nodes are written.
    {
        read>>temp_double;//Ignoring first five numbers
        read>>temp_double;
        read>>temp_double;
        read>>temp_double;
        read>>temp_double;
        
        // 4 pyramid nodes
        read>>ECM_pyramid_Node_1;
        read>>ECM_pyramid_Node_2;
        read>>ECM_pyramid_Node_3;
        read>>ECM_pyramid_Node_4;
        
        ECM_pyramid_Node_1=ECM_pyramid_Node_1-1;   //   The lables begin from 0 .
        ECM_pyramid_Node_2=ECM_pyramid_Node_2-1;
        ECM_pyramid_Node_3=ECM_pyramid_Node_3-1;
        ECM_pyramid_Node_4=ECM_pyramid_Node_4-1;
        
        for(int j=0;j<ECM_Node_Pair_list_temp.size()-1;j++)
        {
            if(  ( ECM_Node_Pair_list_temp[j][0]==ECM_pyramid_Node_1 &  ECM_Node_Pair_list_temp[j][1]==ECM_pyramid_Node_2 )  || ( ECM_Node_Pair_list_temp[j][0]==ECM_pyramid_Node_2 &  ECM_Node_Pair_list_temp[j][1]==ECM_pyramid_Node_1 )    )
            {
                temp_int_repeated_pair_counter_1=1;
            }
            
            if(  ( ECM_Node_Pair_list_temp[j][0]==ECM_pyramid_Node_2 &  ECM_Node_Pair_list_temp[j][1]==ECM_pyramid_Node_3 )  || ( ECM_Node_Pair_list_temp[j][0]==ECM_pyramid_Node_3 &  ECM_Node_Pair_list_temp[j][1]==ECM_pyramid_Node_2 )    )
            {
                temp_int_repeated_pair_counter_2=1;
            }
            
            if(  ( ECM_Node_Pair_list_temp[j][0]==ECM_pyramid_Node_1 &  ECM_Node_Pair_list_temp[j][1]==ECM_pyramid_Node_3 )  || ( ECM_Node_Pair_list_temp[j][0]==ECM_pyramid_Node_3 &  ECM_Node_Pair_list_temp[j][1]==ECM_pyramid_Node_1 )    )
            {
                temp_int_repeated_pair_counter_3=1;
            }
            
            if(  ( ECM_Node_Pair_list_temp[j][0]==ECM_pyramid_Node_1 &  ECM_Node_Pair_list_temp[j][1]==ECM_pyramid_Node_4 )  || ( ECM_Node_Pair_list_temp[j][0]==ECM_pyramid_Node_4 &  ECM_Node_Pair_list_temp[j][1]==ECM_pyramid_Node_1 )    )
            {
                temp_int_repeated_pair_counter_4=1;
            }
            
            if(  ( ECM_Node_Pair_list_temp[j][0]==ECM_pyramid_Node_2 &  ECM_Node_Pair_list_temp[j][1]==ECM_pyramid_Node_4 )  || ( ECM_Node_Pair_list_temp[j][0]==ECM_pyramid_Node_4 &  ECM_Node_Pair_list_temp[j][1]==ECM_pyramid_Node_2 )    )
            {
                temp_int_repeated_pair_counter_5=1;
            }
            
            if(  ( ECM_Node_Pair_list_temp[j][0]==ECM_pyramid_Node_3 &  ECM_Node_Pair_list_temp[j][1]==ECM_pyramid_Node_4 )  || ( ECM_Node_Pair_list_temp[j][0]==ECM_pyramid_Node_4 &  ECM_Node_Pair_list_temp[j][1]==ECM_pyramid_Node_3 )    )
            {
                temp_int_repeated_pair_counter_6=1;
            }
        }
        
        
        if(temp_int_repeated_pair_counter_1==0)
        {
            ECM_Node_Pair_list_temp[ECM_Node_Pair_list_temp.size()-1][0]=ECM_pyramid_Node_1;       //note that first node store in i and second in i+numofbonds  ---Membrane_Node_Pair_list[2*numofbonds]
            ECM_Node_Pair_list_temp[ECM_Node_Pair_list_temp.size()-1][1]=ECM_pyramid_Node_2;
            ECM_Node_Pair_list_temp.push_back(Node_Pairs);
            
        }
        if(temp_int_repeated_pair_counter_2==0)
        {
            ECM_Node_Pair_list_temp[ECM_Node_Pair_list_temp.size()-1][0]=ECM_pyramid_Node_2;       //note that first node store in i and second in i+numofbonds  ---Membrane_Node_Pair_list[2*numofbonds]
            ECM_Node_Pair_list_temp[ECM_Node_Pair_list_temp.size()-1][1]=ECM_pyramid_Node_3;
            ECM_Node_Pair_list_temp.push_back(Node_Pairs);
            
        }
        if(temp_int_repeated_pair_counter_3==0)
        {
            ECM_Node_Pair_list_temp[ECM_Node_Pair_list_temp.size()-1][0]=ECM_pyramid_Node_1;       //note that first node store in i and second in i+numofbonds  ---Membrane_Node_Pair_list[2*numofbonds]
            ECM_Node_Pair_list_temp[ECM_Node_Pair_list_temp.size()-1][1]=ECM_pyramid_Node_3;
            ECM_Node_Pair_list_temp.push_back(Node_Pairs);
            
        }
        
        
        if(temp_int_repeated_pair_counter_4==0)
        {
            ECM_Node_Pair_list_temp[ECM_Node_Pair_list_temp.size()-1][0]=ECM_pyramid_Node_1;       //note that first node store in i and second in i+numofbonds  ---Membrane_Node_Pair_list[2*numofbonds]
            ECM_Node_Pair_list_temp[ECM_Node_Pair_list_temp.size()-1][1]=ECM_pyramid_Node_4;
            ECM_Node_Pair_list_temp.push_back(Node_Pairs);
            
        }
        
        if(temp_int_repeated_pair_counter_5==0)
        {
            ECM_Node_Pair_list_temp[ECM_Node_Pair_list_temp.size()-1][0]=ECM_pyramid_Node_2;       //note that first node store in i and second in i+numofbonds  ---Membrane_Node_Pair_list[2*numofbonds]
            ECM_Node_Pair_list_temp[ECM_Node_Pair_list_temp.size()-1][1]=ECM_pyramid_Node_4;
            ECM_Node_Pair_list_temp.push_back(Node_Pairs);
            
        }
        
        if(temp_int_repeated_pair_counter_6==0)
        {
            ECM_Node_Pair_list_temp[ECM_Node_Pair_list_temp.size()-1][0]=ECM_pyramid_Node_3;       //note that first node store in i and second in i+numofbonds  ---Membrane_Node_Pair_list[2*numofbonds]
            ECM_Node_Pair_list_temp[ECM_Node_Pair_list_temp.size()-1][1]=ECM_pyramid_Node_4;
            ECM_Node_Pair_list_temp.push_back(Node_Pairs);
        }
        temp_int_repeated_pair_counter_1=0;
        temp_int_repeated_pair_counter_2=0;
        temp_int_repeated_pair_counter_3=0;
        temp_int_repeated_pair_counter_4=0;
        temp_int_repeated_pair_counter_5=0;
        temp_int_repeated_pair_counter_6=0;
    }
    
    cout <<"ECM_num_of_Bonds="<<ECM_Node_Pair_list_temp.size()-1<<endl;
    return int(ECM_Node_Pair_list_temp.size()-1);
}

//  Breif discription:
//      -This function reads the Gmesh generated file 'ECM' that contains the node position of the actin network (triangular pyramids).
//      -It should be also noted that we manually delete some of the lines that start with '$' in the actual Gmesh file. So the file is modified before handing it down to the programme.

//  In this function:
//      1-The 'ECM_Node_Velocity', 'ECM_Node_Force',  and the 'ECM_Node_Position' lists are initialised.
//      2-4 nodes (trianguar pyramid) are read from the Gmesh file and the 'Actin_Node_Pair_List' is constructed
//      3-At the very end we set theActin_Node_Pair_List[][2] element for all Actin node pairs. This value is the distance between the actin_node pairs. We use this element to store a 'one step' history for the actin pair distances where we use to calculate the viscus force betweeen the actin node pairs in the dash-pot model.

//  Suggested improvments:
//        1-Nope
void ECM_constructor (double ECM_Node_Position [][3], double ECM_Node_Velocity[][3], double ECM_Node_Force[][3], int ECM_surface_triangle_list[ECM_Surface_num_of_Triangles][3], double ECM_Node_Pair_List[][3], double ECM_upper_surface_Node_Pairs[], int ECM_num_of_Bonds)
{

    int ECM_pyramid_Node_1, ECM_pyramid_Node_2, ECM_pyramid_Node_3, ECM_pyramid_Node_4;
    double temp_double;
    int Num_of_objects;
    
    for(int i=0;i<ECM_num_of_Nodes;i++)
    {
        ECM_Node_Velocity[i][0]= 0.0;
        ECM_Node_Velocity[i][1]= 0.0;
        ECM_Node_Velocity[i][2]= 0.0;
        ECM_Node_Force[i][0]= 0.0;
        ECM_Node_Force[i][1]= 0.0;
        ECM_Node_Force[i][2]= 0.0;
    }
    ifstream read;
    read.open("ECM");
    read>> temp_double;
    
    // When we read the Node positions from the 'ECM' file, we use a 'ECM xyz scaling factor' to rescale the ECM.
    // We also move the ECM, 'Membrane_Centre_distance_from_ECM' below the origin. So we will always put the cell in the origin and move the ECM down.
    for(int i=0;i<ECM_num_of_Nodes;i++)   // This loop reads the initial coordinates of the ECM Nodes
    {
        read>> temp_double;//only sweep not needed elements
        read>> ECM_Node_Position[i][0];
        read>> ECM_Node_Position[i][1];
        read>> ECM_Node_Position[i][2];
        
        ECM_Node_Position[i][0] = ECM_x_scaling_factor * ECM_Node_Position[i][0];
        ECM_Node_Position[i][1] = ECM_y_scaling_factor * ECM_Node_Position[i][1]-Membrane_Centre_distance_from_ECM;
        ECM_Node_Position[i][2] = ECM_z_scaling_factor * ECM_Node_Position[i][2];
    }
    
    read>> Num_of_objects;
    //    Num_of_objects=Num_of_objects-ECM_Surface_num_of_Triangles; // I put this in the loop where it was needed.
    
    for(int i=0;i<ECM_Surface_num_of_Triangles;i++)   //
    {
        read>> temp_double; //Reading not needed elements
        read>> temp_double;
        read>> temp_double;
        read>> temp_double;
        read>> temp_double;
        // Note that gmesh only considers triangles on a 2d surface as such, hence this list does not contain the coordinates of triangles inside the volume (which works for us).
        read>> ECM_surface_triangle_list[i][0];
        read>> ECM_surface_triangle_list[i][1];
        read>> ECM_surface_triangle_list[i][2];
        // Node labels start from 0.
        ECM_surface_triangle_list[i][0]--;
        ECM_surface_triangle_list[i][1]--;
        ECM_surface_triangle_list[i][2]--;
    }
    //We always use the triangles on the ECM surface by taking A=ECM_Node_Position[ ECM_surface_triangle_list[i][0]], B=ECM_Node_Position[ ECM_surface_triangle_list[i][1]], C=ECM_Node_Position[ ECM_surface_triangle_list[i][2]]. We will use this loop to rearrange A, B, and C in order for the ABC triangle normal direction, ABxAC, to point upwards (positive). This is especially usefull in the Membrane_ECM_interaction_2 function.
    
    double AC[3], AB[3], ABxAC[3], xyz[3];
    for(int i=0;i<ECM_Surface_num_of_Triangles;i++)   //
    {
        AB[0]=ECM_Node_Position[ECM_surface_triangle_list[i][1]][0]-ECM_Node_Position[ECM_surface_triangle_list[i][0]][0];
        AB[1]=ECM_Node_Position[ECM_surface_triangle_list[i][1]][1]-ECM_Node_Position[ECM_surface_triangle_list[i][0]][1];
        AB[2]=ECM_Node_Position[ECM_surface_triangle_list[i][1]][2]-ECM_Node_Position[ECM_surface_triangle_list[i][0]][2];
        
        AC[0]=ECM_Node_Position[ECM_surface_triangle_list[i][2]][0]-ECM_Node_Position[ECM_surface_triangle_list[i][0]][0];
        AC[1]=ECM_Node_Position[ECM_surface_triangle_list[i][2]][1]-ECM_Node_Position[ ECM_surface_triangle_list[i][0]][1];
        AC[2]=ECM_Node_Position[ECM_surface_triangle_list[i][2]][2]-ECM_Node_Position[ ECM_surface_triangle_list[i][0]][2];
        
        xyz[0]=ECM_Node_Position[ECM_surface_triangle_list[i][0]][0];
        xyz[1]=ECM_Node_Position[ECM_surface_triangle_list[i][0]][1] + ECM_Thickness/2.0;//I will add a slight shift in the y directio so the position of the ECM is always pointing outwards in respect to the origin. Note that the ECM starts from the y=0 zero plane and ends at y=-ECM_Thickness. We need the outer inner product of the ABxAC and the xyz to be positive to show that the ECM triangle is pointing out of the ECM for simple structures.
        
        xyz[2]=ECM_Node_Position[ECM_surface_triangle_list[i][0]][2];
        //*******************************************************************************************************
        /*BUG
         |---\   |    |  /---\
         |    |  |    |  |
         |---<   |    |  |  -\
         |    |  |    |  |   |
         |---/   \----/  \---/
         */
        //*******************************************************************************************************
        //***************** This method (rearranging the ABC triangles on the ECM to point outwards, will not wotk if the ECM has complicated structures.)
        //*******************************************************************************************************
        crossvector(ABxAC, AB, AC);
        if(innerproduct(ABxAC, xyz)<0 )
        {
            int temp_ECM_Node_B, temp_ECM_Node_C;
            temp_ECM_Node_B=ECM_surface_triangle_list[i][1];
            temp_ECM_Node_C=ECM_surface_triangle_list[i][2];
            ECM_surface_triangle_list[i][1]=temp_ECM_Node_C;
            ECM_surface_triangle_list[i][2]=temp_ECM_Node_B;
        }
    }
    
    int counter=0;
    int temp_int_repeated_pair_counter_1=0;
    int temp_int_repeated_pair_counter_2=0;
    int temp_int_repeated_pair_counter_3=0;
    int temp_int_repeated_pair_counter_4=0;
    int temp_int_repeated_pair_counter_5=0;
    int temp_int_repeated_pair_counter_6=0;
    
    for(int i=0;i<Num_of_objects-ECM_Surface_num_of_Triangles;i++)    // this loop read network between nodes
    {
        
        read>>temp_double; //Reading not needed elements
        read>>temp_double;
        read>>temp_double;
        read>>temp_double;
        read>>temp_double;
        
        // 4 Nodes for the triangular pyramids
        read>>ECM_pyramid_Node_1;
        read>>ECM_pyramid_Node_2;
        read>>ECM_pyramid_Node_3;
        read>>ECM_pyramid_Node_4;
        
        
        
        ECM_pyramid_Node_1--;   //   to lables begin from 0 .
        ECM_pyramid_Node_2--;
        ECM_pyramid_Node_3--;
        ECM_pyramid_Node_4--;
        
        
        
        for(int j=0;j<counter;j++)
        {
            if(  (  ECM_Node_Pair_List[j][0]==ECM_pyramid_Node_1 &   ECM_Node_Pair_List[j][1]==ECM_pyramid_Node_2 )  || (  ECM_Node_Pair_List[j][0]==ECM_pyramid_Node_2 &   ECM_Node_Pair_List[j][1]==ECM_pyramid_Node_1 )    )
            {
                temp_int_repeated_pair_counter_1=1;
            }
            
            if(  (  ECM_Node_Pair_List[j][0]==ECM_pyramid_Node_2 &   ECM_Node_Pair_List[j][1]==ECM_pyramid_Node_3 )  || (  ECM_Node_Pair_List[j][0]==ECM_pyramid_Node_3 &   ECM_Node_Pair_List[j][1]==ECM_pyramid_Node_2 )    )
            {
                temp_int_repeated_pair_counter_2=1;
            }
            
            if(  (  ECM_Node_Pair_List[j][0]==ECM_pyramid_Node_1 &   ECM_Node_Pair_List[j][1]==ECM_pyramid_Node_3 )  || (  ECM_Node_Pair_List[j][0]==ECM_pyramid_Node_3 &   ECM_Node_Pair_List[j][1]==ECM_pyramid_Node_1 )    )
            {
                temp_int_repeated_pair_counter_3=1;
            }
            
            if(  (  ECM_Node_Pair_List[j][0]==ECM_pyramid_Node_1 &   ECM_Node_Pair_List[j][1]==ECM_pyramid_Node_4 )  || (  ECM_Node_Pair_List[j][0]==ECM_pyramid_Node_4 &   ECM_Node_Pair_List[j][1]==ECM_pyramid_Node_1 )    )
            {
                temp_int_repeated_pair_counter_4=1;
            }
            
            if(  (  ECM_Node_Pair_List[j][0]==ECM_pyramid_Node_2 &   ECM_Node_Pair_List[j][1]==ECM_pyramid_Node_4 )  || (  ECM_Node_Pair_List[j][0]==ECM_pyramid_Node_4 &   ECM_Node_Pair_List[j][1]==ECM_pyramid_Node_2 )    )
            {
                temp_int_repeated_pair_counter_5=1;
            }
            
            if(  (  ECM_Node_Pair_List[j][0]==ECM_pyramid_Node_3 &   ECM_Node_Pair_List[j][1]==ECM_pyramid_Node_4 )  || (  ECM_Node_Pair_List[j][0]==ECM_pyramid_Node_4 &   ECM_Node_Pair_List[j][1]==ECM_pyramid_Node_3 )    )
            {
                temp_int_repeated_pair_counter_6=1;
            }
        }
        
        if(temp_int_repeated_pair_counter_1==0)
        {
            if( ECM_Node_Position[ECM_pyramid_Node_1][1]==-Membrane_Centre_distance_from_ECM & ECM_Node_Position[ECM_pyramid_Node_2][1]==-Membrane_Centre_distance_from_ECM  )
            {
                ECM_upper_surface_Node_Pairs[counter]=0.0; // 0.0 means the pair are located on the upper surface of the ECM
            }
            else
            {
                ECM_upper_surface_Node_Pairs[counter]=1.0; // 1.0 means the pair are not on the upper surface of the ECM, hence in the volume or edges and buttom surface
            }
            ECM_Node_Pair_List[counter][0]=ECM_pyramid_Node_1;
            ECM_Node_Pair_List[counter][1]=ECM_pyramid_Node_2;
            
            counter++;
            
        }
        
        if(temp_int_repeated_pair_counter_2==0)
        {
            if( ECM_Node_Position[ECM_pyramid_Node_3][1]==-Membrane_Centre_distance_from_ECM & ECM_Node_Position[ECM_pyramid_Node_2][1]==-Membrane_Centre_distance_from_ECM  )
            {
                ECM_upper_surface_Node_Pairs[counter]=0.0;
            }
            else
            {
                ECM_upper_surface_Node_Pairs[counter]=1.0;
            }
            ECM_Node_Pair_List[counter][0]=ECM_pyramid_Node_2;
            ECM_Node_Pair_List[counter][1]=ECM_pyramid_Node_3;
            counter++;
            
        }
        if(temp_int_repeated_pair_counter_3==0)
        {
            if( ECM_Node_Position[ECM_pyramid_Node_1][1]==-Membrane_Centre_distance_from_ECM & ECM_Node_Position[ECM_pyramid_Node_3][1]==-Membrane_Centre_distance_from_ECM  )
            {
                ECM_upper_surface_Node_Pairs[counter]=0.0;
            }
            else
            {
                ECM_upper_surface_Node_Pairs[counter]=1.0;
            }
            
            ECM_Node_Pair_List[counter][0]=ECM_pyramid_Node_1;
            ECM_Node_Pair_List[counter][1]=ECM_pyramid_Node_3;
            counter++;
            
        }
        
        
        if(temp_int_repeated_pair_counter_4==0)
        {
            if( ECM_Node_Position[ECM_pyramid_Node_1][1]==-Membrane_Centre_distance_from_ECM & ECM_Node_Position[ECM_pyramid_Node_4][1]==-Membrane_Centre_distance_from_ECM  )
            {
                ECM_upper_surface_Node_Pairs[counter]=0.0;
            }
            else
            {
                ECM_upper_surface_Node_Pairs[counter]=1.0;
            }
            ECM_Node_Pair_List[counter][0]=ECM_pyramid_Node_1;
            ECM_Node_Pair_List[counter][1]=ECM_pyramid_Node_4;
            counter++;
            
        }
        
        if(temp_int_repeated_pair_counter_5==0)
        {
            if( ECM_Node_Position[ECM_pyramid_Node_4][1]==-Membrane_Centre_distance_from_ECM & ECM_Node_Position[ECM_pyramid_Node_2][1]==-Membrane_Centre_distance_from_ECM  )
            {
                ECM_upper_surface_Node_Pairs[counter]=0.0;
            }
            else
            {
                ECM_upper_surface_Node_Pairs[counter]=1.0;
            }
            ECM_Node_Pair_List[counter][0]=ECM_pyramid_Node_2;
            ECM_Node_Pair_List[counter][1]=ECM_pyramid_Node_4;
            counter++;
            
        }
        
        if(temp_int_repeated_pair_counter_6==0)
        {
            if( ECM_Node_Position[ECM_pyramid_Node_4][1]==-Membrane_Centre_distance_from_ECM & ECM_Node_Position[ECM_pyramid_Node_3][1]==-Membrane_Centre_distance_from_ECM  )
            {
                ECM_upper_surface_Node_Pairs[counter]=0.0;
            }
            else
            {
                ECM_upper_surface_Node_Pairs[counter]=1.0;
            }
            
            ECM_Node_Pair_List[counter][0]=ECM_pyramid_Node_3;
            ECM_Node_Pair_List[counter][1]=ECM_pyramid_Node_4;
            counter++;
            
        }
        temp_int_repeated_pair_counter_1=0;
        temp_int_repeated_pair_counter_2=0;
        temp_int_repeated_pair_counter_3=0;
        temp_int_repeated_pair_counter_4=0;
        temp_int_repeated_pair_counter_5=0;
        temp_int_repeated_pair_counter_6=0;
    }
    
    double temp_Actin_node_1_position[3],temp_Actin_node_2_position[3],temp_node_pair_distance[3];
    int node_1,node_2;
    
    
    for(int i=0;i<ECM_num_of_Bonds;i++)
    {
        node_1=ECM_Node_Pair_List[i][0];
        node_2=ECM_Node_Pair_List[i][1];
        
        temp_Actin_node_1_position[0]=ECM_Node_Position[node_1][0];
        temp_Actin_node_1_position[1]=ECM_Node_Position[node_1][1];
        temp_Actin_node_1_position[2]=ECM_Node_Position[node_1][2];
        
        temp_Actin_node_2_position[0]=ECM_Node_Position[node_2][0];
        temp_Actin_node_2_position[1]=ECM_Node_Position[node_2][1];
        temp_Actin_node_2_position[2]=ECM_Node_Position[node_2][2];
        
        temp_node_pair_distance[0]=temp_Actin_node_2_position[0]-temp_Actin_node_1_position[0];
        temp_node_pair_distance[1]=temp_Actin_node_2_position[1]-temp_Actin_node_1_position[1];
        temp_node_pair_distance[2]=temp_Actin_node_2_position[2]-temp_Actin_node_1_position[2];
        
        ECM_Node_Pair_List[i][2]=vectorlength(temp_node_pair_distance);
    }
    
    ECM_triangle_normal_direction_justifier(ECM_Node_Position, ECM_surface_triangle_list);
}


void  ECM_triangle_normal_direction_justifier(double  ECM_Node_Position [][3],int ECM_surface_triangle_list[ECM_Surface_num_of_Triangles][3])// n=AB x AC should point out side (A,B,C=ECM_surface_triangle_list[][0,1,2])  for every triangle!
{
    
    int Point_A,Point_B,Point_C; // Think of a triangle having 3 points, A, B, and C.
    double AB[3],AC[3],ABxAC[3];
    //*******************************************************************************************************
    /*BUG
     |---\   |    |  /---\
     |    |  |    |  |
     |---<   |    |  |  -\
     |    |  |    |  |   |
     |---/   \----/  \---/
     */
    //*******************************************************************************************************
    //***************** This function goes through all the triangles on all the surfaces, we usually need only the
    // triangles on the upper surface. Also if we want to calculate the vector product for all the triangles
    // on the surface using this method, 1- we will incounter probloms for the triangles lying on the sides of
    // the ECM (assuming we usually deal with cubes), and 2- the bottum surface needs to point in the -y direction
    // in order to point out of the ECM.
    //*******************************************************************************************************
    
    for(int i=0;i<ECM_Surface_num_of_Triangles;i++) // This number counts all surface triangles not just the ones on the upper surface.
    {
        Point_A=ECM_surface_triangle_list[i][0];
        Point_B=ECM_surface_triangle_list[i][1];
        Point_C=ECM_surface_triangle_list[i][2];
        
        AB[0]=ECM_Node_Position[Point_B][0]-ECM_Node_Position[Point_A][0];
        AB[1]=ECM_Node_Position[Point_B][1]-ECM_Node_Position[Point_A][1];
        AB[2]=ECM_Node_Position[Point_B][2]-ECM_Node_Position[Point_A][2];
        
        AC[0]=ECM_Node_Position[Point_C][0]-ECM_Node_Position[Point_A][0];
        AC[1]=ECM_Node_Position[Point_C][1]-ECM_Node_Position[Point_A][1];
        AC[2]=ECM_Node_Position[Point_C][2]-ECM_Node_Position[Point_A][2];
        crossvector(ABxAC, AB, AC);
        
        if(ABxAC[1] <0) // makes  ABxAC  points aoutward of surcface ONLY WORK FOR N || Y-direction
        {
            // Switching 'Point_B' and 'Point_c' to make the product positive
            ECM_surface_triangle_list[i][1]=Point_C;
            ECM_surface_triangle_list[i][2]=Point_B;
        }
    }
    
}
