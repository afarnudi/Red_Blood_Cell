//
//  Actin_functions.cpp
//  Cell-Durotaxis
//
//  Created by Ali Farnudi on 27/09/2017.
//  Copyright © 2017 Ali Farnudi. All rights reserved.
//

#include "Actin_functions.hpp"

using namespace std;

int  Actin_Node_Pair_Identifier( )  // Reads the 'Actin' files and Identifies the Nodes that have a bond between them, hence are a pair.
{
    vector<vector<int> > Actin_Node_Pair_list_temp;
    vector<int> Node_Pairs;
    
    int Actin_pyramid_Node_1, Actin_pyramid_Node_2, Actin_pyramid_Node_3, Actin_pyramid_Node_4;
    double temp_double;
    int Num_of_objects;
    
    Node_Pairs.resize(2,-1);
    Actin_Node_Pair_list_temp.push_back(Node_Pairs);
    
    int temp_int_repeated_pair_counter_1=0;
    int temp_int_repeated_pair_counter_2=0;
    int temp_int_repeated_pair_counter_3=0;
    int temp_int_repeated_pair_counter_4=0;
    int temp_int_repeated_pair_counter_5=0;
    int temp_int_repeated_pair_counter_6=0;
    
    ifstream read;
    read.open("actin");
    
    //------------------------------------------------------------------------------------------------------
    //------------------------------------------- Improvement ----------------------------------------------
    // We Should, at some point, use a better method for reading the external files (getline for instance).
    // This will also resolve the issue with the files were we had to delete a couple of lines.
    //------------------------------------------------------------------------------------------------------
    
    read>> temp_double;//only sweep not needed elements
    for(int i=0;i<Actin_num_of_Nodes;i++)   //This loop skips (reads) the coordinates of the Nodes in the 'Actin' file.
        //We have already read the actin
    {
        read>> temp_double;//Sweeps not needed elements
        read>> temp_double;
        read>> temp_double;
        read>> temp_double;
    }
    read>> Num_of_objects;//Sweep not needed elements
    
    for(int i=0;i<Num_of_objects;i++)    // This loop reads the nettemp_double between the actin nodes
        // Remmember that we manually deleted the triangles from the 'actin' file by making an 'actin2d' file, so when we start reading the file, the triangular pyramid coordinates are written immediately after the node coordinates.
    {
        read>>temp_double;//Ignoring first five numbers
        read>>temp_double;
        read>>temp_double;
        read>>temp_double;
        read>>temp_double;

        // 4 pyramid nodes
        read>>Actin_pyramid_Node_1;
        read>>Actin_pyramid_Node_2;
        read>>Actin_pyramid_Node_3;
        read>>Actin_pyramid_Node_4;
        
        Actin_pyramid_Node_1--;   //   Lables begin from 0
        Actin_pyramid_Node_2--;
        Actin_pyramid_Node_3--;
        Actin_pyramid_Node_4--;
        
        //========================================== Very Fast, No need for modification =========================================
        for(int j=0;j<Actin_Node_Pair_list_temp.size()-1;j++)
        {
            if(  ( Actin_Node_Pair_list_temp[j][0]==Actin_pyramid_Node_1 &  Actin_Node_Pair_list_temp[j][1]==Actin_pyramid_Node_2 )  || ( Actin_Node_Pair_list_temp[j][0]==Actin_pyramid_Node_2 &  Actin_Node_Pair_list_temp[j][1]==Actin_pyramid_Node_1 )    )
            {
                temp_int_repeated_pair_counter_1=1;
            }
            
            if(  ( Actin_Node_Pair_list_temp[j][0]==Actin_pyramid_Node_2 &  Actin_Node_Pair_list_temp[j][1]==Actin_pyramid_Node_3 )  || ( Actin_Node_Pair_list_temp[j][0]==Actin_pyramid_Node_3 &  Actin_Node_Pair_list_temp[j][1]==Actin_pyramid_Node_2 )    )
            {
                temp_int_repeated_pair_counter_2=1;
            }
            
            if(  ( Actin_Node_Pair_list_temp[j][0]==Actin_pyramid_Node_1 &  Actin_Node_Pair_list_temp[j][1]==Actin_pyramid_Node_3 )  || ( Actin_Node_Pair_list_temp[j][0]==Actin_pyramid_Node_3 &  Actin_Node_Pair_list_temp[j][1]==Actin_pyramid_Node_1 )    )
            {
                temp_int_repeated_pair_counter_3=1;
            }
            
            if(  ( Actin_Node_Pair_list_temp[j][0]==Actin_pyramid_Node_1 &  Actin_Node_Pair_list_temp[j][1]==Actin_pyramid_Node_4 )  || ( Actin_Node_Pair_list_temp[j][0]==Actin_pyramid_Node_4 &  Actin_Node_Pair_list_temp[j][1]==Actin_pyramid_Node_1 )    )
            {
                temp_int_repeated_pair_counter_4=1;
            }
            
            if(  ( Actin_Node_Pair_list_temp[j][0]==Actin_pyramid_Node_2 &  Actin_Node_Pair_list_temp[j][1]==Actin_pyramid_Node_4 )  || ( Actin_Node_Pair_list_temp[j][0]==Actin_pyramid_Node_4 &  Actin_Node_Pair_list_temp[j][1]==Actin_pyramid_Node_2 )    )
            {
                temp_int_repeated_pair_counter_5=1;
            }
            
            if(  ( Actin_Node_Pair_list_temp[j][0]==Actin_pyramid_Node_3 &  Actin_Node_Pair_list_temp[j][1]==Actin_pyramid_Node_4 )  || ( Actin_Node_Pair_list_temp[j][0]==Actin_pyramid_Node_4 &  Actin_Node_Pair_list_temp[j][1]==Actin_pyramid_Node_3 )    )
            {
                temp_int_repeated_pair_counter_6=1;
            }
        }
        
        if(temp_int_repeated_pair_counter_1==0)
        {
            Actin_Node_Pair_list_temp[Actin_Node_Pair_list_temp.size()-1][0]=Actin_pyramid_Node_1;
            Actin_Node_Pair_list_temp[Actin_Node_Pair_list_temp.size()-1][1]=Actin_pyramid_Node_2;
            Actin_Node_Pair_list_temp.push_back(Node_Pairs);
            
        }
        if(temp_int_repeated_pair_counter_2==0)
        {
            Actin_Node_Pair_list_temp[Actin_Node_Pair_list_temp.size()-1][0]=Actin_pyramid_Node_2;
            Actin_Node_Pair_list_temp[Actin_Node_Pair_list_temp.size()-1][1]=Actin_pyramid_Node_3;
            Actin_Node_Pair_list_temp.push_back(Node_Pairs);
            
        }
        if(temp_int_repeated_pair_counter_3==0)
        {
            Actin_Node_Pair_list_temp[Actin_Node_Pair_list_temp.size()-1][0]=Actin_pyramid_Node_1;
            Actin_Node_Pair_list_temp[Actin_Node_Pair_list_temp.size()-1][1]=Actin_pyramid_Node_3;
            Actin_Node_Pair_list_temp.push_back(Node_Pairs);
            
        }
        
        
        if(temp_int_repeated_pair_counter_4==0)
        {
            Actin_Node_Pair_list_temp[Actin_Node_Pair_list_temp.size()-1][0]=Actin_pyramid_Node_1;
            Actin_Node_Pair_list_temp[Actin_Node_Pair_list_temp.size()-1][1]=Actin_pyramid_Node_4;
            Actin_Node_Pair_list_temp.push_back(Node_Pairs);
            
        }
        
        if(temp_int_repeated_pair_counter_5==0)
        {
            Actin_Node_Pair_list_temp[Actin_Node_Pair_list_temp.size()-1][0]=Actin_pyramid_Node_2;
            Actin_Node_Pair_list_temp[Actin_Node_Pair_list_temp.size()-1][1]=Actin_pyramid_Node_4;
            Actin_Node_Pair_list_temp.push_back(Node_Pairs);
            
        }
        
        if(temp_int_repeated_pair_counter_6==0)
        {
            Actin_Node_Pair_list_temp[Actin_Node_Pair_list_temp.size()-1][0]=Actin_pyramid_Node_3;
            Actin_Node_Pair_list_temp[Actin_Node_Pair_list_temp.size()-1][1]=Actin_pyramid_Node_4;
            Actin_Node_Pair_list_temp.push_back(Node_Pairs);
            
        }
        temp_int_repeated_pair_counter_1=0;
        temp_int_repeated_pair_counter_2=0;
        temp_int_repeated_pair_counter_3=0;
        temp_int_repeated_pair_counter_4=0;
        temp_int_repeated_pair_counter_5=0;
        temp_int_repeated_pair_counter_6=0;
        
    }
    
    cout<<"Actin_num_of_Bonds="<<int(Actin_Node_Pair_list_temp.size())-1<<endl;
    return int(Actin_Node_Pair_list_temp.size())-1;
}


void  Actin_constructor(double  Actin_Node_Position [][3],double  Actin_Node_Velocity [][3],double  Actin_Node_Force [][3],double Actin_Node_Pair_List[][3], int Actin_num_of_Bonds)
{
    int Actin_pyramid_Node_1, Actin_pyramid_Node_2, Actin_pyramid_Node_3, Actin_pyramid_Node_4;
    double temp_double;
    int Num_of_objects;
    int counter=0;
    int temp_int_repeated_pair_counter_1=0;
    int temp_int_repeated_pair_counter_2=0;
    int temp_int_repeated_pair_counter_3=0;
    int temp_int_repeated_pair_counter_4=0;
    int temp_int_repeated_pair_counter_5=0;
    int temp_int_repeated_pair_counter_6=0;
    
    for(int i=0;i<Actin_num_of_Nodes;i++)
    {
        Actin_Node_Velocity[i][0]= 0.0;
        Actin_Node_Velocity[i][1]= 0.0;
        Actin_Node_Velocity[i][2]= 0.0;
        Actin_Node_Force[i][0]= 0.0;
        Actin_Node_Force[i][1]= 0.0;
        Actin_Node_Force[i][2]= 0.0;
    }
    
    ///===========================================Find   bondslist   ( temperary files )
//    for(int i=0;i<Actin_num_of_Bonds;i++)
//    {
//        Actin_Node_Pair_List[i][0]=-1;
//        Actin_Node_Pair_List[i][1]=-1;
//    }
    ifstream read;
    read.open("actin");
    read>> temp_double;
    
    for(int i=0;i<Actin_num_of_Nodes;i++) // In this loop we will read the position coordinates of the Actin Nodes
    {
        read>> temp_double; //unnecessary numbers
        read>> Actin_Node_Position[i][0];
        read>> Actin_Node_Position[i][1];
        read>> Actin_Node_Position[i][2];
    }
    //------------------------------------------------------------------------------------------------------
    //------------------------------------------- Improvement ----------------------------------------------
    //---------------------------- This Part is exactly like what we had in the 'Actin_Node_Pair_Identifier'
    // function. It is copied down here because the 'Actin_num_of_Bonds' was needed to construct the Global
    // array'Actin_Node_Pair_List'. Since this is a very fast code, it is not worth the time to change.
    //------------------------------------------------------------------------------------------------------
    
    read>> Num_of_objects;
    for(int i=0;i<Num_of_objects;i++)    // this loop read network between nodes
    {
        read>>temp_double;//only sweep not needed elements
        read>>temp_double;//only sweep not needed elements
        read>>temp_double;//only sweep not needed elements
        read>>temp_double;//only sweep not needed elements
        read>>temp_double;//only sweep not needed elements
        
        // 4 Nodes for the triangular pyramids
        read>>Actin_pyramid_Node_1;
        read>>Actin_pyramid_Node_2;
        read>>Actin_pyramid_Node_3;
        read>>Actin_pyramid_Node_4;
        
        Actin_pyramid_Node_1--;   //Because Node labels Begin from 0
        Actin_pyramid_Node_2--;
        Actin_pyramid_Node_3--;
        Actin_pyramid_Node_4--;
        
        for(int j=0;j<counter;j++)
        {
            if(  (  Actin_Node_Pair_List[j][0]==Actin_pyramid_Node_1 &   Actin_Node_Pair_List[j][1]==Actin_pyramid_Node_2 )  || (  Actin_Node_Pair_List[j][0]==Actin_pyramid_Node_2 &   Actin_Node_Pair_List[j][1]==Actin_pyramid_Node_1 )    )
            {
                temp_int_repeated_pair_counter_1=1;
            }
            
            if(  (  Actin_Node_Pair_List[j][0]==Actin_pyramid_Node_2 &   Actin_Node_Pair_List[j][1]==Actin_pyramid_Node_3 )  || (  Actin_Node_Pair_List[j][0]==Actin_pyramid_Node_3 &   Actin_Node_Pair_List[j][1]==Actin_pyramid_Node_2 )    )
            {
                temp_int_repeated_pair_counter_2=1;
            }
            
            if(  (  Actin_Node_Pair_List[j][0]==Actin_pyramid_Node_1 &   Actin_Node_Pair_List[j][1]==Actin_pyramid_Node_3 )  || (  Actin_Node_Pair_List[j][0]==Actin_pyramid_Node_3 &   Actin_Node_Pair_List[j][1]==Actin_pyramid_Node_1 )    )
            {
                temp_int_repeated_pair_counter_3=1;
            }
            
            if(  (  Actin_Node_Pair_List[j][0]==Actin_pyramid_Node_1 &   Actin_Node_Pair_List[j][1]==Actin_pyramid_Node_4 )  || (  Actin_Node_Pair_List[j][0]==Actin_pyramid_Node_4 &   Actin_Node_Pair_List[j][1]==Actin_pyramid_Node_1 )    )
            {
                temp_int_repeated_pair_counter_4=1;
            }
            
            if(  (  Actin_Node_Pair_List[j][0]==Actin_pyramid_Node_2 &   Actin_Node_Pair_List[j][1]==Actin_pyramid_Node_4 )  || (  Actin_Node_Pair_List[j][0]==Actin_pyramid_Node_4 &   Actin_Node_Pair_List[j][1]==Actin_pyramid_Node_2 )    )
            {
                temp_int_repeated_pair_counter_5=1;
            }
            
            if(  (  Actin_Node_Pair_List[j][0]==Actin_pyramid_Node_3 &   Actin_Node_Pair_List[j][1]==Actin_pyramid_Node_4 )  || (  Actin_Node_Pair_List[j][0]==Actin_pyramid_Node_4 &   Actin_Node_Pair_List[j][1]==Actin_pyramid_Node_3 )    )
            {
                temp_int_repeated_pair_counter_6=1;
            }
        }
        
        
        if(temp_int_repeated_pair_counter_1==0)
        {
            Actin_Node_Pair_List[counter][0]=Actin_pyramid_Node_1;       //note that first node store in i and second in i+numofbonds  ---Membrane_Node_Pair_list[2*numofbonds]
            Actin_Node_Pair_List[counter][1]=Actin_pyramid_Node_2;
            counter=counter+1;
            
        }
        if(temp_int_repeated_pair_counter_2==0)
        {
            Actin_Node_Pair_List[counter][0]=Actin_pyramid_Node_2;       //note that first node store in i and second in i+numofbonds  ---Membrane_Node_Pair_list[2*numofbonds]
            Actin_Node_Pair_List[counter][1]=Actin_pyramid_Node_3;
            counter=counter+1;
            
        }
        if(temp_int_repeated_pair_counter_3==0)
        {
            Actin_Node_Pair_List[counter][0]=Actin_pyramid_Node_1;       //note that first node store in i and second in i+numofbonds  ---Membrane_Node_Pair_list[2*numofbonds]
            Actin_Node_Pair_List[counter][1]=Actin_pyramid_Node_3;
            counter=counter+1;
            
        }
        
        
        if(temp_int_repeated_pair_counter_4==0)
        {
            Actin_Node_Pair_List[counter][0]=Actin_pyramid_Node_1;       //note that first node store in i and second in i+numofbonds  ---Membrane_Node_Pair_list[2*numofbonds]
            Actin_Node_Pair_List[counter][1]=Actin_pyramid_Node_4;
            counter=counter+1;
            
        }
        
        if(temp_int_repeated_pair_counter_5==0)
        {
            Actin_Node_Pair_List[counter][0]=Actin_pyramid_Node_2;       //note that first node store in i and second in i+numofbonds  ---Membrane_Node_Pair_list[2*numofbonds]
            Actin_Node_Pair_List[counter][1]=Actin_pyramid_Node_4;
            counter=counter+1;
            
        }
        
        if(temp_int_repeated_pair_counter_6==0)
        {
            Actin_Node_Pair_List[counter][0]=Actin_pyramid_Node_3;       //note that first node store in i and second in i+numofbonds  ---Membrane_Node_Pair_list[2*numofbonds]
            Actin_Node_Pair_List[counter][1]=Actin_pyramid_Node_4;
            counter++;
            
        }
        temp_int_repeated_pair_counter_1=0;
        temp_int_repeated_pair_counter_2=0;
        temp_int_repeated_pair_counter_3=0;
        temp_int_repeated_pair_counter_4=0;
        temp_int_repeated_pair_counter_5=0;
        temp_int_repeated_pair_counter_6=0;
    }
    
    ///========================== Find xbigbeads  and  bondslist   ( temperary files )
    
    double temp_Actin_node_1_position[3],temp_Actin_node_2_position[3],temp_node_pair_distance[3];
    int node_1,node_2;
    
    
    for(int i=0;i<Actin_num_of_Bonds;i++)
    {
        node_1=Actin_Node_Pair_List[i][0];
        node_2=Actin_Node_Pair_List[i][1];
        
        temp_Actin_node_1_position[0]=Actin_Node_Position[node_1][0];
        temp_Actin_node_1_position[1]=Actin_Node_Position[node_1][1];
        temp_Actin_node_1_position[2]=Actin_Node_Position[node_1][2];
        
        temp_Actin_node_2_position[0]=Actin_Node_Position[node_2][0];
        temp_Actin_node_2_position[1]=Actin_Node_Position[node_2][1];
        temp_Actin_node_2_position[2]=Actin_Node_Position[node_2][2];
        
        temp_node_pair_distance[0]=temp_Actin_node_2_position[0]-temp_Actin_node_1_position[0];
        temp_node_pair_distance[1]=temp_Actin_node_2_position[1]-temp_Actin_node_1_position[1];
        temp_node_pair_distance[2]=temp_Actin_node_2_position[2]-temp_Actin_node_1_position[2];
        
        Actin_Node_Pair_List[i][2]=vectorlength(temp_node_pair_distance);
        
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
    //***************** It seems that the 'bond' file is created but never used *****************************
    //***************** I am going to comment it out of the code and see if I incounter any problems ********
    //*******************************************************************************************************
}











