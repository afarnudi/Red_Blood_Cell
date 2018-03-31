//-----------------------HEADERS
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <limits>
#include <cstdlib>
#include <random>
#include <string>
#include <math.h>
#include <cstdlib>

#include "General_constants.h"
#include "General_functions.hpp"
#include "General_Membrane.h"
#include "Membrane_functions.hpp"
#include "General_Actin.h"
#include "Actin_functions.hpp"
#include "Actin_membrane_shared_functions.hpp"
#include "General_Actin_Membrane_shared.h"
#include "Monte_Carlo_Bond_Flip.hpp"

using namespace std;

//-----------------------HEADERS

void  Actin_Force_calculator( double  Actin_Node_Position [][3],double  Actin_Node_Velocity [][3],double  Actin_Node_Force [][3],double Actin_Node_Pair_List[][3], int Actin_num_of_Bonds, double &Total_Potential_Energy);


void Actin_Membrane_Barrier_2(double Actin_Node_Position[][3], double Actin_Node_Velocity[][3], double Membrane_Node_Position[][3], double Membrane_Node_Velocity[][3],  vector<vector<int> > Membrane_new_triangle_list, vector<vector<int> > Membrane_Actin_shared_Node_list);


//--------------------------------ECM--------------------------------


void cellshift(double Membrane_Node_Position[][3], double Actin_Node_Position[Actin_num_of_Nodes][3], double  Membrane_Node_Velocity [][3], double Actin_Node_Velocity[Actin_num_of_Nodes][3], int Membrane_num_of_Nodes);

#define membraneshiftinXdirection 0
#define membraneshiftinZdirection 0
#define cell_downward_speed 0.0

void CellCOM( double com[3], double Membrane_Node_Position[][3], double Actin_Node_Position[][3], int Membrane_num_of_Nodes);

//-----------------------------thermostat------------------
void Thermostat_2(double Membrane_Node_Velocity[][3], double Actin_Node_Velocity[][3], int Membrane_num_of_Nodes);
//-----------------------------thermostat------------------


//---------------------------restart-----------------------

void restartsave(double Membrane_Node_Position[][3], double Membrane_Node_Velocity [][3], vector<vector<int> > Membrane_new_triangle_list, int Membrane_Triangle_Pair_Nodes[][4], int Membrane_Node_Pair_list[][2], vector<vector<int> > membrane_triangle_pair_list, double Actin_Node_Position[Actin_num_of_Nodes][3], double Actin_Node_Velocity[Actin_num_of_Nodes][3], double Actin_Node_Pair_List[][3], int Membrane_num_of_Triangle_Pairs, int Actin_num_of_Bonds, int Membrane_num_of_Nodes);


void restartread(double Membrane_Node_Position[][3], double Membrane_Node_Velocity[][3], vector<vector<int> > &Membrane_new_triangle_list, int Membrane_Triangle_Pair_Nodes[][4], int Membrane_Node_Pair_list[][2], vector<vector<int> > &membrane_triangle_pair_list, double Actin_Node_Position[Actin_num_of_Nodes][3], double Actin_Node_Velocity[Actin_num_of_Nodes][3], double Actin_Node_Pair_List[][3], int Membrane_num_of_Triangle_Pairs, int Actin_num_of_Bonds, int Membrane_num_of_Nodes);
//---------------------------restart-----------------------


//Povray------
# define export_povray_step_distance 10000  // clear!
//void povray_output_creator(int currentStep, double Membrane_Node_Position[Membrane_num_of_Nodes][3], int  Membrane_triangle_list[Membrane_num_of_Triangles][3], int Membrane_Normal_direction[Membrane_num_of_Triangles][2], int Membrane_Node_Pair_list[][2], double Actin_Node_Position[Actin_num_of_Nodes][3], double Actin_Node_Pair_List[][3], double Chromatin_Bead_Position[Chromatin_num_of_Beads][3], double  ECM_Node_Position[][3], double ECM_Node_Pair_List[][3], int ECM_surface_triangle_list[ECM_Surface_num_of_Triangles][3], int Outer_Membrane_num_of_triangles, int Membrane_num_of_Node_Pairs, int Outer_Membrane_num_of_Node_Pairs, int Actin_num_of_Bonds, int ECM_num_of_Bonds);

//Ali
void Membrane_Actin_shared_Node_Force_calculator (double Membrane_Node_Position[][3],double  Actin_Node_Position [Actin_num_of_Nodes][3], double Membrane_Node_Force [][3], double  Actin_Node_Force [Actin_num_of_Nodes][3], vector<vector<int> > Membrane_Actin_shared_Node_list, double Membrane_Node_velocity[][3], double Actin_Node_velocity[Actin_num_of_Nodes][3]);// updates forces + relavant potential energy

void generate_initial_condition_report (string initial_condition_file_name, int Membrane_num_of_Nodes);
int Membrane_Num_of_Nodes_reader (string membrane_mesh_file_name);
int Membrane_Num_of_Nodes_reader (string membrane_mesh_file_name){
    ifstream read; //This is the main ifstream that will read the Gmesh-Membrane generated file
    read.open(membrane_mesh_file_name.c_str() ); //It should be noted that the name of the file should not contain '-'. I don't know why but the memory managnet of the arrays (at the very least) in the programme will collapse when we use '-' in the file name.
    int temp_int; // This is just a temp intiger charachter that we use to read unnecessary Gmesh generated intigers. We never use these intigers in the actual programme.
    string temp_string;
    for (int i=0; i<6; i++) {
        read>> temp_string;
    }
    read>> temp_int;
    cout<<"Membrane number of Nodes is=\t"<<temp_int<<endl;
    return temp_int;
}

#define Actin_membrane_stiff_spring_coefficient 400
#define Actin_membrane_damping_coefficient 0.0
#define Membrane_ECM_Max_cos_triangle_interaction_angle 0.2
#define ECM_Membrane_Radius_of_Hard_Sphere_Interaction 1.0//0.50
#define Membrane_ECM_triangle_interaction_update_step   50 // should be 50 in contactescase but here is suspended cell
#define Membrane_ECM_triangle_unbinding_update_step   50 // should be 50 in contactescase but here is suspended cell
#define Membrane_ECM_triangle_interaction_rate 10
#define Membrane_ECM_interaction_strength 20.0
#define Membrane_ECM_Binding_prob 0.5
#define Membrane_ECM_UnBinding_prob 0.2
#define spreading_flag 0.0
#define spreading_force_magnitude -100.0
#define spreading_force_min_range 1.0
#define spreading_force_max_range 4.0
#define spreading_force_cos_triangle_interaction_angle 0.4
#define energy_calculation_flag 0.0




int main() //main**
{
    clock_t tStart = clock();//Time the programme
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    int num_of_test_particles=13;
    int test_particles[num_of_test_particles][3];
    //direction identifier
    for (int i=0; i<num_of_test_particles; i++) {
        test_particles[i][0]=12+i;
        test_particles[i][1]=-11.5+i*1;
        test_particles[i][2]=0;
        
        //        test_particles[i][0]=12-i;
        //        test_particles[i][1]=-11.5+i*1;
        //        test_particles[i][2]=i*i;
        
        //        if (i<(num_of_test_particles/2)) {
        //            test_particles[i][0]=12;
        //        }   else {
        //            test_particles[i][0]=15;
        //        }
        //        test_particles[i][1]=-11.5+i*1;
        //        test_particles[i][2]=0;
    }
    //Edges
    test_particles[num_of_test_particles-3][0]=20;
    test_particles[num_of_test_particles-3][1]=10;
    test_particles[num_of_test_particles-2][0]=-20;
    test_particles[num_of_test_particles-2][1]=10;
    
    
    test_particles[num_of_test_particles-4][2]=20;
    test_particles[num_of_test_particles-4][1]=10;
    test_particles[num_of_test_particles-4][0]=0;
    test_particles[num_of_test_particles-5][2]=-20;
    test_particles[num_of_test_particles-5][1]=10;
    test_particles[num_of_test_particles-5][0]=0;
    //Centre
    test_particles[num_of_test_particles-1][0]=0;
    test_particles[num_of_test_particles-1][1]=20;
    test_particles[num_of_test_particles-1][2]=0;
    
    //    double test_cout=123456.78910111213;
    
    char buffer [80];
    strftime (buffer,80,"%Y_%m_%d_%H:%M",now);
    //outputfiles:
    string traj_file_name, initial_condition_file_name, energy_file_name, force_file_name, com_file_name;
    
    traj_file_name="results/RBC_";
    traj_file_name +=buffer;
    initial_condition_file_name = traj_file_name;
    energy_file_name=initial_condition_file_name;
    force_file_name=initial_condition_file_name;
    com_file_name=initial_condition_file_name;
    traj_file_name +=".xyz";
    initial_condition_file_name+=".txt";
    
    ofstream trajectory;
    trajectory.open(traj_file_name.c_str() );
    trajectory << std:: fixed;
    
    
    //Energy Calculation variables-Begin
    double Membrane_total_potential_Energy=0.0, Actin_total_potential_Energy=0.0, Chromatin_total_potential_Energy=0.0, ECM_total_potential_Energy=0.0, ECM_membrane_total_potential_Energy=0.0;
    
    
    energy_file_name+="energy.txt";
    ofstream write_energy;
    write_energy.open(energy_file_name.c_str());
    write_energy<<"Time\tTotal Energy\tMembrane\tActin\tChromatin\tECM\tInteraction_4"<<endl;
    //Energy Calculation variables-End
    
    
    double Membrane_total_force[3], actin_total_force[3];
    for (int i=0; i<3; i++) {
        Membrane_total_force[i]=0;
        actin_total_force[i]=0;
    }
    
    force_file_name+="force.txt";
    ofstream write_force;
    write_force.open(force_file_name.c_str());
    write_force<<"Time\tMembrane\tActin"<<endl;
    
    bool resume=false;//If true the programme will call upon the 'restartread' function and read the positions and velocities of the nodes from an external file that is determined inside the function.
    //    bool chromatin_contant_matrix_calculation_flag=false;
    //-----------------------------------------------membrane:
    //    int counter = 0;//This counter counts the number of MD steps and wil be used to determin on which steps the preogramme is saved ('restartsave').
    string membrane_mesh_file_name="RBCMembrane";
    const int Membrane_num_of_Nodes=Membrane_Num_of_Nodes_reader(membrane_mesh_file_name);
    
    double Membrane_Node_Position[Membrane_num_of_Nodes][3], Membrane_Node_Velocity[Membrane_num_of_Nodes][3], Membrane_Node_Force[Membrane_num_of_Nodes][3];
    //    int Membrane_triangle_list[Membrane_num_of_Triangles_temp][3];
    
    generate_initial_condition_report(initial_condition_file_name, Membrane_num_of_Nodes);
    
    double Total_Kinetic_Energy, Total_Potential_Energy; // total kineti , potential ,and mechanical energy
    
    cout <<"Please double check the Membrane Radius, it should be "<< Membrane_Radius<<endl;
    cout <<"Initialising the programme ..."<<endl;
    vector<vector<int> > Membrane_new_triangle_list;
    Membrane_constructor(Membrane_Node_Position, Membrane_Node_Velocity, Membrane_Node_Force, Membrane_new_triangle_list, membrane_mesh_file_name);
    
    
    //    RBC_triangle_list_builder(Membrane_triangle_list, Membrane_triangle_list_2);
    
    
    Membrane_Normal_direction_Identifier(Membrane_Node_Position, Membrane_new_triangle_list, Membrane_num_of_Nodes);
    
    int  Membrane_num_of_Triangle_Pairs;
    Membrane_num_of_Triangle_Pairs=Membrane_triangle_pair_counter( Membrane_new_triangle_list);
    int Membrane_Triangle_Pair_Nodes[Membrane_num_of_Triangle_Pairs][4]; // pos1 pos2 pos3 and po4 of all interactions are stored here
    vector<vector<int> > membrane_triangle_pair_list;
    membrane_triangle_pair_list.resize(Membrane_new_triangle_list.size());
    
    Membrane_Triangle_Pair_Identifier(Membrane_new_triangle_list, Membrane_Triangle_Pair_Nodes, Membrane_num_of_Triangle_Pairs, membrane_triangle_pair_list); //mod**
    
    
    //    int Outer_Membrane_num_of_Node_Pairs ; // usefull in POV ray- modulate in void sortingbonds(int bondslist[][2],int tri[3][Membrane_num_of_Triangles])
    int Membrane_num_of_Node_Pairs;
    Membrane_num_of_Node_Pairs=Membrane_num_of_Node_Pair_Counter(Membrane_new_triangle_list, Membrane_num_of_Nodes);
    
    int Membrane_Node_Pair_list[Membrane_num_of_Node_Pairs][2];
    //*******************************************************************************************************
    /*BUG
     |---\   |    |  /---\
     |    |  |    |  |
     |---<   |    |  |  -\
     |    |  |    |  |   |
     |---/   \----/  \---/
     */
    //*******************************************************************************************************
    //***************** Since we use the 'Membrane_num_of_Node_Pair_Counter' just to count the number of Node pairs, we have to write a similar function to build the actual 'Membrane_Node_Pair_list', which I have assigned 'Membrane_num_of_Node_Pair_Counter' to do so. Let t be noted that not calling this function will not cause any problems if the 'resume' is on!. *****************************
    //*******************************************************************************************************
    if (resume==false) {
        Membrane_num_of_Node_Pair_Counter_2(Membrane_Node_Pair_list, Membrane_new_triangle_list, Membrane_num_of_Node_Pairs);
    }
    
    //-----------------------------------------------membrane:
    
    //*******************************************************************************************************
    /*BUG
     |---\   |    |  /---\
     |    |  |    |  |
     |---<   |    |  |  -\
     |    |  |    |  |   |
     |---/   \----/  \---/
     */
    //*******************************************************************************************************
    //***************** We spend so much energy and time to seperate the Outer  *****************************
    //***************** Membrane from the Nucleus. Why not just build 2? ************************************
    //*******************************************************************************************************
    
    
    //    Monte_carlo_bond_flip(Membrane_Node_Position, Membrane_triangle_list, membrane_triangle_pair_list, Membrane_Node_Pair_list, Membrane_num_of_Node_Pairs);
    
    double centerOfmassoftheCell[3];
    
    ofstream write_COM;
    com_file_name+="COM.txt";
    write_COM.open(com_file_name.c_str() );
    write_COM << std:: fixed;
    write_COM<<"Time\tCOM_x\tCOM_y\tCOM_z\n";
    
    
    
    //---------------------------------------------------actin:
    int Actin_num_of_Bonds=Actin_Node_Pair_Identifier(); // the value will be set automatically in actin parameter function
    
    double Actin_Node_Pair_List[Actin_num_of_Bonds][3]; //In this array 3 numbers are stored for each Actin Node pairs. The first two are of course the label (number) of the nodes that are a pair. The third element is used to store the distance between these pairs. The distance is used in the force calculations for the Maxwell spring initial length. The initial length is updated during each step so we have a diferent initial length, hence the Maxwell spring.
    
    double Actin_Node_Position[Actin_num_of_Nodes][3];
    double Actin_Node_Velocity[Actin_num_of_Nodes][3];
    double Actin_Node_Force[Actin_num_of_Nodes][3];
    //    double Actin_Node_VelocityRungKuta[Actin_num_of_Nodes][3];  //only used in  Rungekuta step
    Actin_constructor(Actin_Node_Position, Actin_Node_Velocity, Actin_Node_Force, Actin_Node_Pair_List, Actin_num_of_Bonds);
    
    
    //---------------------------------------------------actin
    
    
    //-----------------------------actin-membrane part: shared beads(suppose to have masses like mass of membrane)
    //    int Membrane_Actin_shared_Node_list[Actin_Membrane_shared_num_of_Nodes][2];// IMPORTANT: this list contains the indices of shared beads on the membrane and actin. [0] for membrane and [1] for actin.
    vector<vector<int> > Membrane_Actin_shared_Node_list;
    Membrane_Actin_shared_Node_Identifier(Membrane_Actin_shared_Node_list, Membrane_Node_Position, Actin_Node_Position, Membrane_num_of_Nodes);
    
    
    
    
    cellshift(Membrane_Node_Position, Actin_Node_Position, Membrane_Node_Velocity, Actin_Node_Velocity, Membrane_num_of_Nodes);
    if( resume==true )
    {
        restartread(Membrane_Node_Position, Membrane_Node_Velocity, Membrane_new_triangle_list, Membrane_Triangle_Pair_Nodes, Membrane_Node_Pair_list, membrane_triangle_pair_list, Actin_Node_Position, Actin_Node_Velocity, Actin_Node_Pair_List, Membrane_num_of_Triangle_Pairs, Actin_num_of_Bonds, Membrane_num_of_Nodes);
        cellshift(Membrane_Node_Position, Actin_Node_Position, Membrane_Node_Velocity, Actin_Node_Velocity, Membrane_num_of_Nodes);
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
    //***************** It would seem that the ECM has been constructed twice, once here and once the ECM section.
    // We should probably find the reason why or get rid of one.
    //*******************************************************************************************************
    
    //==========================================================================================================================
    //==========================================================================================================================
    //============================================== Beginning of the main MD loop ==============================================
    //==========================================================================================================================
    //==========================================================================================================================
    
    int counter2=0;
    cout<<"Beginning the MD loop\n";
    for(int MD_Step=0 ;MD_Step<=MD_num_of_steps ; MD_Step++)
    {
        
        for (int i=0; i<3; i++) {
            Membrane_total_force[i]=0;
            actin_total_force[i]=0;
        }
        
        
        for(int j=0 ; j<Membrane_num_of_Nodes ; j++)
        {
            Membrane_Node_Position[j][0] += Membrane_Node_Velocity[j][0]*MD_Time_Step - Membrane_Node_Force[j][0]*MD_Time_Step*MD_Time_Step/(Membrane_Node_Mass*2.0);
            Membrane_Node_Position[j][1] += Membrane_Node_Velocity[j][1]*MD_Time_Step - Membrane_Node_Force[j][1]*MD_Time_Step*MD_Time_Step/(Membrane_Node_Mass*2.0);
            Membrane_Node_Position[j][2] += Membrane_Node_Velocity[j][2]*MD_Time_Step - Membrane_Node_Force[j][2]*MD_Time_Step*MD_Time_Step/(Membrane_Node_Mass*2.0);
        }
        
        for(int j=0 ; j<Actin_num_of_Nodes ; j++)  //actin
        {
            Actin_Node_Position[j][0] += Actin_Node_Velocity[j][0]*MD_Time_Step - Actin_Node_Force[j][0]*MD_Time_Step*MD_Time_Step/(Actin_Node_Mass*2.0);
            Actin_Node_Position[j][1] += Actin_Node_Velocity[j][1]*MD_Time_Step - Actin_Node_Force[j][1]*MD_Time_Step*MD_Time_Step/(Actin_Node_Mass*2.0);
            Actin_Node_Position[j][2] += Actin_Node_Velocity[j][2]*MD_Time_Step - Actin_Node_Force[j][2]*MD_Time_Step*MD_Time_Step/(Actin_Node_Mass*2.0);
        }
        
        for(int j=0 ; j<Membrane_num_of_Nodes ; j++)  // loop to count every particle  and update its velocity
        {
            Membrane_Node_Velocity[j][0] += - Membrane_Node_Force[j][0]*MD_Time_Step/(Membrane_Node_Mass*2.0);
            Membrane_Node_Velocity[j][1] += - Membrane_Node_Force[j][1]*MD_Time_Step/(Membrane_Node_Mass*2.0);
            Membrane_Node_Velocity[j][2] += - Membrane_Node_Force[j][2]*MD_Time_Step/(Membrane_Node_Mass*2.0);
        }
        
        for(int j=0 ; j<Actin_num_of_Nodes ; j++)  //actin
        {
            
            Actin_Node_Velocity[j][0] += -Actin_Node_Force[j][0]*MD_Time_Step/(Actin_Node_Mass*2.0);
            Actin_Node_Velocity[j][1] += -Actin_Node_Force[j][1]*MD_Time_Step/(Actin_Node_Mass*2.0);
            Actin_Node_Velocity[j][2] += -Actin_Node_Force[j][2]*MD_Time_Step/(Actin_Node_Mass*2.0);
            
        }
        
        
        Total_Potential_Energy=0.0;
        Total_Kinetic_Energy=0.0;
        
        for(int j=0 ; j<Membrane_num_of_Nodes ; j++)
        {
            Membrane_Node_Force[j][0]=0.0;
            Membrane_Node_Force[j][1]=0.0;
            Membrane_Node_Force[j][2]=0.0;
        }
        for(int j=0 ; j<Actin_num_of_Nodes ; j++)
        {
            Actin_Node_Force[j][0]=0.0;
            Actin_Node_Force[j][1]=0.0;
            Actin_Node_Force[j][2]=0.0;
        }
        
        // **********calling aceel update
        
        
        Membrane_Force_Calculator(Membrane_Node_Position, Membrane_Node_Velocity, Membrane_Node_Force, Membrane_Node_Pair_list, Membrane_Triangle_Pair_Nodes, Total_Potential_Energy,  Membrane_num_of_Triangle_Pairs, Membrane_num_of_Node_Pairs, Membrane_num_of_Nodes); // **********calling aceel update
        
        
        ConstantSurfaceForceLocalTriangles( Membrane_Node_Position, Membrane_Node_Force, Membrane_new_triangle_list);
        
        if (MD_Step % mcstep == 0)// collsi
        {
            for(int s=0; s< fluidity*Membrane_num_of_Nodes ;s++)
            {
                Monte_carlo_bond_flip(Membrane_Node_Position, Membrane_new_triangle_list, membrane_triangle_pair_list, Membrane_Node_Pair_list, Membrane_num_of_Node_Pairs, Membrane_num_of_Nodes);
                //                MonteCarlo(Membrane_triangle_list, Membrane_Triangle_Pair_Nodes, Membrane_Node_Pair_list, Total_Potential_Energy, Membrane_Node_Position, Membrane_Normal_direction, Outer_Membrane_num_of_triangles, Nucleus_Membrane_num_of_triangles, Membrane_num_of_Triangle_Pairs, Membrane_num_of_Node_Pairs);
                counter2=0;
            }
            counter2++;
        }
        
        
        if (energy_calculation_flag==1.0) {
            Membrane_total_potential_Energy=Total_Potential_Energy;
            //        cout<<"Total_Potential_Energy= "<<Total_Potential_Energy<<endl;
        }
        
        
        //        for(int j=0 ; j<Actin_num_of_Nodes ; j++)  //Rung kuta update
        //        {
        //            Actin_Node_VelocityRungKuta[j][0] += -Actin_Node_Force[j][0]*MD_Time_Step/(Actin_Node_Mass);
        //            Actin_Node_VelocityRungKuta[j][1] += -Actin_Node_Force[j][1]*MD_Time_Step/(Actin_Node_Mass);
        //            Actin_Node_VelocityRungKuta[j][2] += -Actin_Node_Force[j][2]*MD_Time_Step/(Actin_Node_Mass);
        //        }
        
        
        Actin_Force_calculator(Actin_Node_Position, Actin_Node_Velocity, Actin_Node_Force, Actin_Node_Pair_List, Actin_num_of_Bonds, Total_Potential_Energy); // updates with runge kuta
        
        
        if (energy_calculation_flag==1.0) {
            Actin_total_potential_Energy=Total_Potential_Energy-Membrane_total_potential_Energy;
            //        cout<<"Total_Potential_Energy= "<<Total_Potential_Energy<<endl;
        }
        
        Membrane_Actin_shared_Node_Force_calculator(Membrane_Node_Position, Actin_Node_Position, Membrane_Node_Force,Actin_Node_Force, Membrane_Actin_shared_Node_list, Membrane_Node_Velocity, Actin_Node_Velocity);
        
        if(MD_Step%Membrane_barrier_calculation_rate==0)
        {
            Actin_Membrane_Barrier_2(Actin_Node_Position, Actin_Node_Velocity, Membrane_Node_Position, Membrane_Node_Velocity, Membrane_new_triangle_list, Membrane_Actin_shared_Node_list);
            
        }
        
        
        
        for(int j=0 ; j<Membrane_num_of_Nodes ; j++)  // loop to count every particle  and update its velocity
        {
            Membrane_Node_Velocity[j][0] += - Membrane_Node_Force[j][0]*MD_Time_Step/(Membrane_Node_Mass*2.0);
            Membrane_Node_Velocity[j][1] += - Membrane_Node_Force[j][1]*MD_Time_Step/(Membrane_Node_Mass*2.0);
            Membrane_Node_Velocity[j][2] += - Membrane_Node_Force[j][2]*MD_Time_Step/(Membrane_Node_Mass*2.0);
            
            Membrane_total_force[0]+=Membrane_Node_Force[j][0];
            Membrane_total_force[1]+=Membrane_Node_Force[j][1];
            Membrane_total_force[2]+=Membrane_Node_Force[j][2];
        }
        write_force<<MD_Step<<"\t"<<Membrane_total_force[0]<<"\t"<<Membrane_total_force[1]<<"\t"<<Membrane_total_force[2];
        
        for(int j=0 ; j<Actin_num_of_Nodes ; j++)  //actin
        {
            Actin_Node_Velocity[j][0] += -Actin_Node_Force[j][0]*MD_Time_Step/(Actin_Node_Mass*2.0);
            Actin_Node_Velocity[j][1] += -Actin_Node_Force[j][1]*MD_Time_Step/(Actin_Node_Mass*2.0);
            Actin_Node_Velocity[j][2] += -Actin_Node_Force[j][2]*MD_Time_Step/(Actin_Node_Mass*2.0);
            
            actin_total_force[0]+=Actin_Node_Force[j][0];
            actin_total_force[1]+=Actin_Node_Force[j][1];
            actin_total_force[2]+=Actin_Node_Force[j][2];
        }
        
        write_force<<"\t"<<actin_total_force[0]<<"\t"<<actin_total_force[1]<<"\t"<<actin_total_force[2]<<endl;
        
        
        if (MD_Step%10000==0 & MD_Step!=0  )
        {
            
            cout<<"Saving temp..."<<endl;
            restartsave(Membrane_Node_Position, Membrane_Node_Velocity, Membrane_new_triangle_list, Membrane_Triangle_Pair_Nodes, Membrane_Node_Pair_list, membrane_triangle_pair_list,  Actin_Node_Position, Actin_Node_Velocity, Actin_Node_Pair_List, Membrane_num_of_Triangle_Pairs, Actin_num_of_Bonds, Membrane_num_of_Nodes);
        }
        
        //--------------------------------Thermostate
        
        
        
        if(MD_Step%RunThermostatePerstep==0)
        {
            Thermostat_2(Membrane_Node_Velocity, Actin_Node_Velocity, Membrane_num_of_Nodes);
        }
        
        
        //--------------------------------------------------------------------saving and friends!
        if (MD_Step%savingstep==0  )  // SAVING AND CALCULATION SECTION
        {
            if (energy_calculation_flag==1.0) {
                write_energy<<MD_Step<<"\t"<<Total_Potential_Energy<<"\t"<<Membrane_total_potential_Energy<<"\t"<<Actin_total_potential_Energy<<"\t"<<Chromatin_total_potential_Energy<<"\t"<<ECM_total_potential_Energy<<"\t"<<ECM_membrane_total_potential_Energy<<"\n";
            }
            
            cout <<MD_Step<<endl;
            
            
            ///__________________________________________________Traj__________________________
            
            trajectory << Membrane_num_of_Nodes+Actin_num_of_Nodes+num_of_test_particles<<endl;//+Num_of_Actin_Nodes_on_Membrane_list <<endl;    // saving trajectories
            
            trajectory << " nodes  "<<endl;
            
            for(int j=0; j<Membrane_num_of_Nodes;j++) // saving trajectory
            {
                trajectory <<"membrane"  <<setprecision(5)<< setw(20)<<Membrane_Node_Position[j][0]<< setw(20)<<Membrane_Node_Position[j][1]<< setw(20)<<Membrane_Node_Position[j][2]<<endl;
            }
            
            
            for(int j=0; j< Actin_num_of_Nodes ;j++) // saving trajectory
            {
                trajectory << "actin" <<  setprecision(5)<< setw(20)<<Actin_Node_Position[j][0]<< setw(20)<<Actin_Node_Position[j][1]<< setw(20)<<Actin_Node_Position[j][2]<<endl;
                
            }
            
            
            
            for(int j=0; j<num_of_test_particles;j++) // saving trajectory
            {
                trajectory <<"test"  <<setprecision(5)<< setw(20)<<test_particles[j][0]<< setw(20)<<test_particles[j][1]<< setw(20)<<test_particles[j][2]<<endl;
            }
            
            
            
            ///__________________________________________________Traj__________________________
            
            CellCOM(centerOfmassoftheCell, Membrane_Node_Position, Actin_Node_Position, Membrane_num_of_Nodes);
            
            write_COM<<MD_Step<<"\t"<<centerOfmassoftheCell[0]<<"\t"<<centerOfmassoftheCell[1]<<"\t"<<centerOfmassoftheCell[2]<<endl;
        }
        // ______________________________________________________________SAVING AND CALCULATION SECTION
        
    }
    
    //==========================================================================================================================
    //==========================================================================================================================
    //================================================ End of the main MD loop =================================================
    //==========================================================================================================================
    //==========================================================================================================================
    
    
    
    //_______________________________restart:
    
    cout<<"Saving..."<<endl;
    restartsave(Membrane_Node_Position, Membrane_Node_Velocity, Membrane_new_triangle_list, Membrane_Triangle_Pair_Nodes, Membrane_Node_Pair_list, membrane_triangle_pair_list, Actin_Node_Position, Actin_Node_Velocity, Actin_Node_Pair_List, Membrane_num_of_Triangle_Pairs, Actin_num_of_Bonds, Membrane_num_of_Nodes);
    
    
    
    
    
    
    cout<<"done"<<endl;
    printf("Time taken: %.2f minutes\n", (double)((clock() - tStart)/CLOCKS_PER_SEC)/60.0);
}

// ==============================================================================================================
// ============================================== END OF MAIN ===================================================
// ==============================================================================================================





//______________________membrane functions

//Also the potential Energy is calculated here





void Membrane_Actin_shared_Node_Force_calculator (double Membrane_Node_Position[][3], double  Actin_Node_Position[][3], double Membrane_Node_Force[][3], double Actin_Node_Force[][3],  vector<vector<int> > Membrane_Actin_shared_Node_list, double Membrane_Node_velocity[][3], double Actin_Node_velocity[Actin_num_of_Nodes][3])
{
    double delta_x,delta_y,delta_z,temp_Node_distance,temp_force;
    
    int temp_Node_mem,temp_Node_act;
    
    for (int act_mem_temp_node=0 ; act_mem_temp_node< Membrane_Actin_shared_Node_list.size() ; act_mem_temp_node++)
    {
        temp_Node_mem = Membrane_Actin_shared_Node_list[act_mem_temp_node][0];
        temp_Node_act = Membrane_Actin_shared_Node_list[act_mem_temp_node][1];
        
        delta_x = Membrane_Node_Position[temp_Node_mem][0]-Actin_Node_Position[temp_Node_act][0];
        delta_y = Membrane_Node_Position[temp_Node_mem][1]-Actin_Node_Position[temp_Node_act][1];
        delta_z = Membrane_Node_Position[temp_Node_mem][2]-Actin_Node_Position[temp_Node_act][2];
        
        temp_Node_distance=sqrt(delta_x*delta_x + delta_y*delta_y + delta_z*delta_z);
        temp_force=0.0;
        
        if (temp_Node_distance > 0.0001 || temp_Node_distance < -0.0001) {
            temp_force = Actin_membrane_stiff_spring_coefficient*temp_Node_distance;
            
            Membrane_Node_Force[temp_Node_mem][0] += temp_force*delta_x/temp_Node_distance;
            Membrane_Node_Force[temp_Node_mem][1] += temp_force*delta_y/temp_Node_distance;
            Membrane_Node_Force[temp_Node_mem][2] += temp_force*delta_z/temp_Node_distance;
            
            Actin_Node_Force[temp_Node_act][0] += -temp_force*delta_x/temp_Node_distance;
            Actin_Node_Force[temp_Node_act][1] += -temp_force*delta_y/temp_Node_distance;
            Actin_Node_Force[temp_Node_act][2] += -temp_force*delta_z/temp_Node_distance;
        }
        Membrane_Node_Force[temp_Node_mem][0] += Actin_membrane_damping_coefficient*(Membrane_Node_velocity[temp_Node_mem][0]-Actin_Node_velocity[temp_Node_act][0]);
        Membrane_Node_Force[temp_Node_mem][1] += Actin_membrane_damping_coefficient*(Membrane_Node_velocity[temp_Node_mem][1]-Actin_Node_velocity[temp_Node_act][1]);
        Membrane_Node_Force[temp_Node_mem][2] += Actin_membrane_damping_coefficient*(Membrane_Node_velocity[temp_Node_mem][2]-Actin_Node_velocity[temp_Node_act][2]);
        Actin_Node_Force[temp_Node_act][0] += -Actin_membrane_damping_coefficient*(Membrane_Node_velocity[temp_Node_mem][0]-Actin_Node_velocity[temp_Node_act][0]);
        Actin_Node_Force[temp_Node_act][1] += -Actin_membrane_damping_coefficient*(Membrane_Node_velocity[temp_Node_mem][1]-Actin_Node_velocity[temp_Node_act][1]);
        Actin_Node_Force[temp_Node_act][2] += -Actin_membrane_damping_coefficient*(Membrane_Node_velocity[temp_Node_mem][2]-Actin_Node_velocity[temp_Node_act][2]);
        
    }
    //    exit (EXIT_FAILURE);
}











//______________________Actin functions


void  Actin_Force_calculator(double Actin_Node_Position[][3], double Actin_Node_Velocity[][3], double  Actin_Node_Force[][3], double Actin_Node_Pair_List[][3], int Actin_num_of_Bonds, double &Total_Potential_Energy)
{
    
    double deltax, deltay, deltaz, temp_distance, initial_distance;// defined below in "for loop" in detail
    int node1, node2;
    double  temp_force[3];
    
    /// Spring
    for(int i=0 ;i<Actin_num_of_Bonds ; i++ )  // all beads interaction whit the next one
    {
        node1=(int) Actin_Node_Pair_List[i][0];
        node2=(int) Actin_Node_Pair_List[i][1];
        
        initial_distance=Actin_Node_Pair_List[i][2];
        
        deltax=0.0;
        deltay=0.0;
        deltaz=0.0;
        temp_distance=0.0;
        
        deltax=Actin_Node_Position[node2][0]-Actin_Node_Position[node1][0];// delta x  betwwn i and i+1 th beads
        deltay=Actin_Node_Position[node2][1]-Actin_Node_Position[node1][1];// delta y  betwwn i and i+1 th beads
        deltaz=Actin_Node_Position[node2][2]-Actin_Node_Position[node1][2];// delta z  betwwn i and i+1 th beads
        temp_distance=sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz); // distance btween i th and i+1 th  bead
        
        if(CytoskeletonNetworkType==0)  {// if PCN is off use simple hookian network
            
            temp_force[0]=-Actin_spring_coefficient*(temp_distance-initial_distance) *  deltax/temp_distance  ;
            temp_force[1]=-Actin_spring_coefficient*(temp_distance-initial_distance) *  deltay/temp_distance  ;
            temp_force[2]=-Actin_spring_coefficient*(temp_distance-initial_distance) *  deltaz/temp_distance  ;
            
            //            cout<<temp_force[0]<<"\t"<<temp_force[1]<<"\t"<<temp_force[2]<<"\n";
            
            if (energy_calculation_flag==1.0) {
                Total_Potential_Energy += 0.5*Actin_spring_coefficient*(temp_distance-initial_distance)*(temp_distance-initial_distance);
                
            }
            
            Actin_Node_Force[node1][0] += temp_force[0]  ; // force of springs and dashes
            Actin_Node_Force[node1][1] += temp_force[1]  ;// force of springs and dashes
            Actin_Node_Force[node1][2] += temp_force[2]  ;// force of springs and dashes
            
            Actin_Node_Force[node2][0] += - temp_force[0]  ; // force of springs and dashes
            Actin_Node_Force[node2][1] += - temp_force[1]  ;// force of springs and dashes
            Actin_Node_Force[node2][2] += - temp_force[2]  ;// force of springs and dashes
            
            Actin_Node_Pair_List[i][2] +=  MD_Time_Step*   (  (temp_distance-initial_distance)/temp_distance  )   /Actin_kelvin_damping_coefficient; //We have used the Kelvin model for the Visco elasticity. Kelvin model: A spring and dashpot (in series) are in parralel with a spring.
            
        }  else if(CytoskeletonNetworkType==1  && (temp_distance-initial_distance) > 0.0)  {// passive cable network  on?
            temp_force[0]=-Actin_Passive_Cable_Network_Coefficient*(temp_distance-initial_distance)*  deltax/temp_distance  ;
            temp_force[1]=-Actin_Passive_Cable_Network_Coefficient*(temp_distance-initial_distance)*  deltay/temp_distance  ;
            temp_force[2]=-Actin_Passive_Cable_Network_Coefficient*(temp_distance-initial_distance)*  deltaz/temp_distance  ;
            if (energy_calculation_flag==1.0) {
                Total_Potential_Energy += 0.5*Actin_Passive_Cable_Network_Coefficient*(temp_distance-initial_distance)*(temp_distance-initial_distance);
            }
            
            
            Actin_Node_Force[node1][0] += temp_force[0]  ;// force of springs and dashes
            Actin_Node_Force[node1][1] += temp_force[1]  ;// force of springs and dashes
            Actin_Node_Force[node1][2] += temp_force[2]  ;// force of springs and dashes
            
            Actin_Node_Force[node2][0] += - temp_force[0]  ;// force of springs and dashes
            Actin_Node_Force[node2][1] += - temp_force[1]  ;// force of springs and dashes
            Actin_Node_Force[node2][2] += - temp_force[2]  ;// force of springs and dashes
            
        } else if(CytoskeletonNetworkType==2 )  {//active cable network
            
            if(temp_distance > initial_distance)  { //more detail in the paper: Contractile network models for adherent cells
                
                temp_force[0]=-(KActinACN_EA *(temp_distance-initial_distance) +ACN_TL0)   *  deltax/temp_distance  ;
                temp_force[1]=-(KActinACN_EA *(temp_distance-initial_distance) +ACN_TL0)   *  deltay/temp_distance  ;
                temp_force[2]=-(KActinACN_EA *(temp_distance-initial_distance) +ACN_TL0)   *  deltaz/temp_distance  ;
                if (energy_calculation_flag==1.0) {
                    Total_Potential_Energy += 0.5*KActinACN_EA*(temp_distance-initial_distance)*(temp_distance-initial_distance) + ACN_TL0*(temp_distance-initial_distance);
                }
                
            }
            else if(   (temp_distance < initial_distance)  &&  temp_distance>ACN_LC )   //more detail in the paper: Contractile network models for adherent cells
            {
                temp_force[0]=-(ACN_TL0)   *  deltax/temp_distance  ;
                temp_force[1]=-(ACN_TL0)   *  deltay/temp_distance  ;
                temp_force[2]=-(ACN_TL0)   *  deltaz/temp_distance  ;
                
                if (energy_calculation_flag==1.0) {
                    Total_Potential_Energy += ACN_TL0*(temp_distance-initial_distance);
                }
                
            }
            
            else if(temp_distance < ACN_LC)   //more detail in the paper: Contractile network models for adherent cells
                
            {
                temp_force[0]=-( ACN_TL0 *(temp_distance)/ACN_LC )   *  deltax/temp_distance  ;
                temp_force[1]=-( ACN_TL0 *(temp_distance)/ACN_LC )   *  deltay/temp_distance  ;
                temp_force[2]=-( ACN_TL0 *(temp_distance)/ACN_LC )   *  deltaz/temp_distance  ;
                if (energy_calculation_flag==1.0) {
                    Total_Potential_Energy += 0.5*ACN_TL0*(temp_distance)*(temp_distance)/ACN_LC;
                }
                
            }
            
            
            Actin_Node_Force[node1][0] +=  temp_force[0]  ;// force of springs and dashes
            Actin_Node_Force[node1][1] +=  temp_force[1]  ;// force of springs and dashes
            Actin_Node_Force[node1][2] +=  temp_force[2]  ;// force of springs and dashes
            
            Actin_Node_Force[node2][0] += - temp_force[0]  ;// force of springs and dashes
            Actin_Node_Force[node2][1] += - temp_force[1]  ;// force of springs and dashes
            Actin_Node_Force[node2][2] += - temp_force[2]  ;// force of springs and dashes
            
        }
        
        ///Damping force
        temp_force[0]= Actin_damping_Coefficient*(Actin_Node_Velocity[node1][0]-Actin_Node_Velocity[node2][0]);
        temp_force[1]= Actin_damping_Coefficient*(Actin_Node_Velocity[node1][1]-Actin_Node_Velocity[node2][1]);
        temp_force[2]= Actin_damping_Coefficient*(Actin_Node_Velocity[node1][2]-Actin_Node_Velocity[node2][2]);
        if (energy_calculation_flag==1.0) {
            Total_Potential_Energy -= 0.5*Actin_damping_Coefficient*((Actin_Node_Velocity[node1][2]-Actin_Node_Velocity[node2][2])*(Actin_Node_Velocity[node1][2]-Actin_Node_Velocity[node2][2])+(Actin_Node_Velocity[node1][1]-Actin_Node_Velocity[node2][1])*(Actin_Node_Velocity[node1][1]-Actin_Node_Velocity[node2][1])+(Actin_Node_Velocity[node1][0]-Actin_Node_Velocity[node2][0])*(Actin_Node_Velocity[node1][0]-Actin_Node_Velocity[node2][0]));
        }
        
        //        cout<<temp_force[0]<<"\t"<<temp_force[1]<<"\t"<<temp_force[2]<<"\n";
        
        
        Actin_Node_Force[node1][0] += temp_force[0]  ; // force of springs and dashes
        Actin_Node_Force[node1][1] += temp_force[1]  ;// force of springs and dashes
        Actin_Node_Force[node1][2] += temp_force[2]  ;// force of springs and dashes
        
        Actin_Node_Force[node2][0] += -temp_force[0]  ; // force of springs and dashes
        Actin_Node_Force[node2][1] += -temp_force[1]  ;// force of springs and dashes
        Actin_Node_Force[node2][2] += -temp_force[2]  ;// force of springs and dashes
        
        
        
        //
    }
    
    //    for (int i=0; i<20; i++) {
    //        cout<<Actin_Node_Position[i][0]<<"\t"<<Actin_Node_Position[i][1]<<"\t"<<Actin_Node_Position[i][2]<<"\n";
    //    }
    //        exit(EXIT_FAILURE);
    
}


void Actin_Membrane_Barrier_2( double  Actin_Node_Position [][3], double  Actin_Node_Velocity [][3], double Membrane_Node_Position[][3], double Membrane_Node_Velocity[][3], vector<vector<int> > Membrane_new_triangle_list, vector<vector<int> > Membrane_Actin_shared_Node_list)
{
    double membrane_triangle_COM_position[3]; // coordinates to the centre of mass of the triangles
    double membrane_triangle_COM_velocity[3]; // velocity of the center of mass of the triangles
    double actin_membrane_distance_amplitude;// distance between solvent particle and com of triangle
    double ABxAC[3], AB[3], AC[3], ABxAC_unit_vector[3]; // normal vector of membrane
    double actin_triangle_distance_vector[3];
    double perpendicular_distance;
    double relevant_velocity[3];
    
    for (int i =0; i < Membrane_new_triangle_list.size()  ; i++)
    {
        membrane_triangle_COM_position[0]=(Membrane_Node_Position[ Membrane_new_triangle_list[i][0]][0] + Membrane_Node_Position[Membrane_new_triangle_list[i][1]][0] + Membrane_Node_Position[ Membrane_new_triangle_list[i][2]][0])/3.0;
        membrane_triangle_COM_position[1]=(Membrane_Node_Position[ Membrane_new_triangle_list[i][0]][1] + Membrane_Node_Position[Membrane_new_triangle_list[i][1]][1] + Membrane_Node_Position[ Membrane_new_triangle_list[i][2]][1])/3.0;
        membrane_triangle_COM_position[2]=(Membrane_Node_Position[ Membrane_new_triangle_list[i][0]][2] + Membrane_Node_Position[Membrane_new_triangle_list[i][1]][2] + Membrane_Node_Position[ Membrane_new_triangle_list[i][2]][2])/3.0;
        
        membrane_triangle_COM_velocity[0]=(Membrane_Node_Velocity[ Membrane_new_triangle_list[i][0]][0] + Membrane_Node_Velocity[Membrane_new_triangle_list[i][1]][0] + Membrane_Node_Velocity[ Membrane_new_triangle_list[i][2]][0])/3.0;
        membrane_triangle_COM_velocity[1]=(Membrane_Node_Velocity[ Membrane_new_triangle_list[i][0]][1] + Membrane_Node_Velocity[Membrane_new_triangle_list[i][1]][1] + Membrane_Node_Velocity[ Membrane_new_triangle_list[i][2]][1])/3.0;
        membrane_triangle_COM_velocity[2]=(Membrane_Node_Velocity[ Membrane_new_triangle_list[i][0]][2] + Membrane_Node_Velocity[Membrane_new_triangle_list[i][1]][2] + Membrane_Node_Velocity[ Membrane_new_triangle_list[i][2]][2])/3.0;
        
        for (int actin_counter=int(Membrane_Actin_shared_Node_list.size());actin_counter<Actin_num_of_Nodes;actin_counter++)
        {
            actin_membrane_distance_amplitude = sqrt( (membrane_triangle_COM_position[0] -Actin_Node_Position[actin_counter][0] ) * (membrane_triangle_COM_position[0] - Actin_Node_Position[actin_counter][0]) + (membrane_triangle_COM_position[1] - Actin_Node_Position[actin_counter][1]) * (membrane_triangle_COM_position[1]- Actin_Node_Position[actin_counter][1]) + (membrane_triangle_COM_position[2] - Actin_Node_Position[actin_counter][2]) * (membrane_triangle_COM_position[2] - Actin_Node_Position[actin_counter][2]));
            if (  actin_membrane_distance_amplitude < sqrt(0.43*Node_radius * 0.43*Node_radius + Actin_Membrane_Radius_of_Hard_Sphere_Interaction*Actin_Membrane_Radius_of_Hard_Sphere_Interaction )  )
            {
                AB[0]=Membrane_Node_Position[ Membrane_new_triangle_list[i][1]][0]-Membrane_Node_Position[ Membrane_new_triangle_list[i][0]][0];
                AB[1]=Membrane_Node_Position[ Membrane_new_triangle_list[i][1]][1]-Membrane_Node_Position[ Membrane_new_triangle_list[i][0]][1];
                AB[2]=Membrane_Node_Position[ Membrane_new_triangle_list[i][1]][2]-Membrane_Node_Position[ Membrane_new_triangle_list[i][0]][2];
                AC[0]=Membrane_Node_Position[ Membrane_new_triangle_list[i][2]][0]-Membrane_Node_Position[ Membrane_new_triangle_list[i][0]][0];
                AC[1]=Membrane_Node_Position[ Membrane_new_triangle_list[i][2]][1]-Membrane_Node_Position[ Membrane_new_triangle_list[i][0]][1];
                AC[2]=Membrane_Node_Position[ Membrane_new_triangle_list[i][2]][2]-Membrane_Node_Position[ Membrane_new_triangle_list[i][0]][2];
                crossvector(ABxAC,AB,AC);
                //                ABxAC[0]=ABxAC[0]*Membrane_Normal_direction[i][1];
                //                ABxAC[1]=ABxAC[1]*Membrane_Normal_direction[i][1];
                //                ABxAC[2]=ABxAC[2]*Membrane_Normal_direction[i][1];
                
                ABxAC_unit_vector[0]=ABxAC[0]/vectorlength(ABxAC);
                ABxAC_unit_vector[1]=ABxAC[1]/vectorlength(ABxAC);
                ABxAC_unit_vector[2]=ABxAC[2]/vectorlength(ABxAC);
                
                actin_triangle_distance_vector[0]=Actin_Node_Position[actin_counter][0]-membrane_triangle_COM_position[0];
                actin_triangle_distance_vector[1]=Actin_Node_Position[actin_counter][1]-membrane_triangle_COM_position[1];
                actin_triangle_distance_vector[2]=Actin_Node_Position[actin_counter][2]-membrane_triangle_COM_position[2];
                
                perpendicular_distance=innerproduct(actin_triangle_distance_vector, ABxAC)/vectorlength(ABxAC);
                
                relevant_velocity[0] = Actin_Node_Velocity[actin_counter][0]-membrane_triangle_COM_velocity[0];
                relevant_velocity[1] = Actin_Node_Velocity[actin_counter][1]-membrane_triangle_COM_velocity[1];
                relevant_velocity[2] = Actin_Node_Velocity[actin_counter][2]-membrane_triangle_COM_velocity[2];
                
                if    ( (abs( perpendicular_distance )<Actin_Membrane_Radius_of_Hard_Sphere_Interaction)
                       && (innerproduct(relevant_velocity,ABxAC)>0)
                       )
                {
                    double Actin_velocity_N_new, Actin_velocity_N, Membrane_triangle_COM_velocity_N_new, Membrane_triangle_COM_velocity_N;
                    
                    Actin_velocity_N=innerproduct(Actin_Node_Velocity[actin_counter], ABxAC_unit_vector);
                    Membrane_triangle_COM_velocity_N=innerproduct(membrane_triangle_COM_velocity, ABxAC_unit_vector);
                    
                    Actin_velocity_N_new=(Actin_velocity_N*(Actin_Node_Mass-3*Membrane_Node_Mass)+2.0*3*Membrane_Node_Mass*Membrane_triangle_COM_velocity_N)/(Actin_Node_Mass+3*Membrane_Node_Mass);
                    Membrane_triangle_COM_velocity_N_new=(Membrane_triangle_COM_velocity_N*(3*Membrane_Node_Mass-Actin_Node_Mass)+2.0*Actin_Node_Mass*Actin_velocity_N)/(Actin_Node_Mass+3*Membrane_Node_Mass);
                    
                    Actin_Node_Velocity[actin_counter][0]+= (-Actin_velocity_N + Actin_velocity_N_new)*ABxAC_unit_vector[0];
                    Actin_Node_Velocity[actin_counter][1]+= (-Actin_velocity_N + Actin_velocity_N_new)*ABxAC_unit_vector[1];
                    Actin_Node_Velocity[actin_counter][2]+= (-Actin_velocity_N + Actin_velocity_N_new)*ABxAC_unit_vector[2];
                    
                    for (int i1 = 0; i1 < 3; i1++)
                    {
                        Membrane_Node_Velocity[Membrane_new_triangle_list[i][i1]][0] += (-Membrane_triangle_COM_velocity_N + Membrane_triangle_COM_velocity_N_new)*ABxAC_unit_vector[0];
                        Membrane_Node_Velocity[Membrane_new_triangle_list[i][i1]][1] += (-Membrane_triangle_COM_velocity_N + Membrane_triangle_COM_velocity_N_new)*ABxAC_unit_vector[1];
                        Membrane_Node_Velocity[Membrane_new_triangle_list[i][i1]][2] += (-Membrane_triangle_COM_velocity_N + Membrane_triangle_COM_velocity_N_new)*ABxAC_unit_vector[2];
                    }
                }//END OF:  if    ( (abs( perpendicular_distance )<Actin_Membrane_Radius_of_Hard_Sphere_Interaction) &&
            }//END OF: if (  actin_membrane_distance_amplitude < sqrt(0.43*a * 0.43*a +
            
        }//END OF: for (int actin_counter=Actin_Membrane_shared_num_of_Nodes;actin_counter<Actin_num_of_Nodes;
        
    }//END OF: for (int i =0; i < Membrane_num_of_Triangles  ; i++)
    
}



void cellshift(double Membrane_Node_Position[][3], double Actin_Node_Position[Actin_num_of_Nodes][3], double Membrane_Node_Velocity[][3], double Actin_Node_Velocity[Actin_num_of_Nodes][3], int Membrane_num_of_Nodes)
{
    for (int i=0; i<Membrane_num_of_Nodes ; i++)
    {
        Membrane_Node_Position[i][0] += membraneshiftinXdirection;
        Membrane_Node_Position[i][2] += membraneshiftinZdirection;
    }
    for (int i=0; i<Actin_num_of_Nodes ; i++)
    {
        Actin_Node_Position[i][0]= Actin_Node_Position[i][0]+membraneshiftinXdirection;
        Actin_Node_Position[i][2]= Actin_Node_Position[i][2]+membraneshiftinZdirection;
    }
    
    
    for (int i=0; i<Membrane_num_of_Nodes ; i++)
    {
        Membrane_Node_Velocity[i][1] += cell_downward_speed;
        
    }
    for (int i=0; i<Actin_num_of_Nodes ; i++)
    {
        Actin_Node_Velocity[i][1]= Actin_Node_Velocity[i][1]+cell_downward_speed;
        
    }
    
}


//___________________migration

void CellCOM( double com[3],double  Membrane_Node_Position [][3],double  Actin_Node_Position [][3], int Membrane_num_of_Nodes)
{
    com[0]=0;
    com[1]=0;
    com[2]=0;
    double sigmaM=Membrane_num_of_Nodes*Membrane_Node_Mass+Actin_num_of_Nodes*Actin_Node_Mass;
    //----------------------membrane---------------------
    for(int i=0;i<Membrane_num_of_Nodes;i++)
    {
        com [0] += Membrane_Node_Mass*Membrane_Node_Position [i][0];
        com [1] += Membrane_Node_Mass*Membrane_Node_Position [i][1];
        com [2] += Membrane_Node_Mass*Membrane_Node_Position [i][2];
        
    }
    //----------------------actin---------------------
    
    for(int i=0; i<Actin_num_of_Nodes; i++)
    {
        
        com [0] += Actin_Node_Mass* Actin_Node_Position [i] [0];
        com [1] += Actin_Node_Mass* Actin_Node_Position [i] [1];
        com [2] += Actin_Node_Mass* Actin_Node_Position [i] [2];
        //        cout<<Actin_Node_Position [i] [0]<<"\t"<<Actin_Node_Position [i] [1]<<"\t"<<Actin_Node_Position [i] [2]<<"\n";
        //        cout<<com [0]<<"\t"<<com [1]<<"\t"<<com [2]<<"\n";
    }
    //    cout<<com [0]<<"\t"<<com [1]<<"\t"<<com [2]<<"\n";
    //----------------------------------------
    com [0]=com[0]/sigmaM;
    com [1]=com[1]/sigmaM;
    com [2]=com[2]/sigmaM;
    //    cout<<com [0]<<"\t"<<com [1]<<"\t"<<com [2]<<"\n";
}









//-----------------------------thermostat------------------
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
    
    
    //    temp_V[0]=0;
    //    temp_V[1]=0;
    //    temp_V[2]=0;
    //
    //    V_com[0]=0;
    //    V_com[1]=0;
    //    V_com[2]=0;
    //    for (int i=0; i<Membrane_num_of_Nodes; i++) {
    //        V_com[0]+=Membrane_Node_Velocity[i][0];
    //        V_com[1]+=Membrane_Node_Velocity[i][1];
    //        V_com[2]+=Membrane_Node_Velocity[i][2];
    //    }
    //    V_com[0]*=Membrane_Node_Mass/Membrane_num_of_Nodes;
    //    V_com[1]*=Membrane_Node_Mass/Membrane_num_of_Nodes;
    //    V_com[2]*=Membrane_Node_Mass/Membrane_num_of_Nodes;
    //
    //    for (int i=0; i<Actin_num_of_Nodes; i++) {
    //        temp_V[0]+=Actin_Node_Velocity[i][0];
    //        temp_V[1]+=Actin_Node_Velocity[i][1];
    //        temp_V[2]+=Actin_Node_Velocity[i][2];
    //    }
    //    temp_V[0]*=Actin_Node_Mass/Actin_num_of_Nodes;
    //    temp_V[1]*=Actin_Node_Mass/Actin_num_of_Nodes;
    //    temp_V[2]*=Actin_Node_Mass/Actin_num_of_Nodes;
    //
    //    V_com[0]=(V_com[0]+temp_V[0])/(Membrane_num_of_Nodes*Membrane_Node_Mass+Actin_num_of_Nodes*Actin_Node_Mass);
    //    V_com[1]=(V_com[1]+temp_V[1])/(Membrane_num_of_Nodes*Membrane_Node_Mass+Actin_num_of_Nodes*Actin_Node_Mass);
    //    V_com[2]=(V_com[2]+temp_V[2])/(Membrane_num_of_Nodes*Membrane_Node_Mass+Actin_num_of_Nodes*Actin_Node_Mass);
    //
    //    cout<<"COM after= "<<sqrt(V_com[0]*V_com[0]+V_com[1]*V_com[1]+V_com[2]*V_com[2])<<"\nV_X="<<V_com[0]<<"\tV_Y="<<V_com[1]<<"\tV_Z="<<V_com[2]<<"\n\n";
    
}

//-----------------------------thermostat------------------






//---------------------------restart-----------------------
void restartsave(double Membrane_Node_Position[][3], double Membrane_Node_Velocity[][3], vector<vector<int> > Membrane_new_triangle_list, int Membrane_Triangle_Pair_Nodes[][4], int Membrane_Node_Pair_list[][2], vector<vector<int> > membrane_triangle_pair_list, double Actin_Node_Position[Actin_num_of_Nodes][3], double Actin_Node_Velocity[Actin_num_of_Nodes][3], double Actin_Node_Pair_List[][3], int Membrane_num_of_Triangle_Pairs, int Actin_num_of_Bonds, int Membrane_num_of_Nodes)
{
    ofstream restart;
    restart.open("restart_1.txt");
    restart << std:: fixed;
    
    for (int i=0; i<Membrane_num_of_Nodes ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart<< Membrane_Node_Position[i][j]<<endl;
        }
    }
    
    for (int i=0; i<Membrane_num_of_Nodes ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart<< Membrane_Node_Velocity[i][j]<<endl;
        }
    }
    
    
    for (int i=0; i<Membrane_new_triangle_list.size() ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart<< Membrane_new_triangle_list[i][j]<<endl;
        }
    }
    
    
    for (int i=0; i<Membrane_num_of_Triangle_Pairs ; i++)
    {
        for (int j=0; j<4 ; j++)
        {
            restart<< Membrane_Triangle_Pair_Nodes[i][j]<<endl;
        }
    }
    
    for (int i=0; i<Membrane_num_of_Triangle_Pairs ; i++)
    {
        for (int j=0; j<2 ; j++)
        {
            restart<< Membrane_Node_Pair_list[i][j]<<endl;
        }
    }
    
    for (int i=0; i<membrane_triangle_pair_list.size() ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart<< membrane_triangle_pair_list[i][j]<<endl;
        }
    }
    
    for (int i=0; i<Actin_num_of_Nodes ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart<< Actin_Node_Position[i][j]<<endl;
        }
    }
    for (int i=0; i<Actin_num_of_Nodes ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart<< Actin_Node_Velocity[i][j]<<endl;
        }
    }
    
    
    for (int i=0; i<Actin_num_of_Bonds ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart<< Actin_Node_Pair_List[i][j]<<endl;
        }
    }
    
    
}


void restartread(double Membrane_Node_Position [][3], double Membrane_Node_Velocity[][3], vector<vector<int> > &Membrane_new_triangle_list, int Membrane_Triangle_Pair_Nodes[][4], int Membrane_Node_Pair_list[][2], vector<vector<int> > &membrane_triangle_pair_list, double Actin_Node_Position[Actin_num_of_Nodes][3], double Actin_Node_Velocity[Actin_num_of_Nodes][3], double Actin_Node_Pair_List[][3], int Membrane_num_of_Triangle_Pairs, int Actin_num_of_Bonds, int Membrane_num_of_Nodes)
{
    ifstream restart;
    //    restart.open("restart-backup.txt");
    restart.open("restart_s_2.txt");
    
    for (int i=0; i<Membrane_num_of_Nodes ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart>> Membrane_Node_Position[i][j];
        }
    }
    
    for (int i=0; i<Membrane_num_of_Nodes ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart>> Membrane_Node_Velocity[i][j];
        }
    }
    
    
    for (int i=0; i<Membrane_new_triangle_list.size() ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart>> Membrane_new_triangle_list[i][j];
        }
    }
    
    for (int i=0; i<Membrane_num_of_Triangle_Pairs ; i++)
    {
        for (int j=0; j<4 ; j++)
        {
            restart>> Membrane_Triangle_Pair_Nodes[i][j];
        }
    }
    
    
    
    for (int i=0; i<Membrane_num_of_Triangle_Pairs ; i++)
    {
        for (int j=0; j<2 ; j++)
        {
            restart>> Membrane_Node_Pair_list[i][j];
        }
    }
    
    for (int i=0; i<membrane_triangle_pair_list.size() ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart>> membrane_triangle_pair_list[i][j];
        }
    }
    
    for (int i=0; i<Actin_num_of_Nodes ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart>> Actin_Node_Position[i][j];
        }
    }
    
    
    
    for (int i=0; i<Actin_num_of_Nodes ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart>> Actin_Node_Velocity[i][j];
        }
    }
    
    
    for (int i=0; i<Actin_num_of_Bonds ; i++)
    {
        for (int j=0; j<3 ; j++)
        {
            restart>> Actin_Node_Pair_List[i][j];
        }
    }
}





//---------------------------restart-----------------------



void generate_initial_condition_report (string initial_condition_file_name, int Membrane_num_of_Nodes){
    ofstream write_initial_condition;
    write_initial_condition.open(initial_condition_file_name.c_str() );
    
    write_initial_condition<<"MD_num_of_steps"<<"=\t"<<MD_num_of_steps<<endl;
    write_initial_condition<<"savingstep"<<"=\t"<<savingstep<<endl;
    write_initial_condition<<"MD_Time_Step"<<"=\t"<<MD_Time_Step<<endl;
    write_initial_condition<<"KT"<<"=\t"<<KT<<endl;
    write_initial_condition<<"RunThermostatePerstep"<<"=\t"<<RunThermostatePerstep<<endl;
    write_initial_condition<<"\n\n";
    write_initial_condition<<"Membrane_num_of_Nodes"<<"=\t"<<Membrane_num_of_Nodes<<endl;
    //    write_initial_condition<<"Membrane_num_of_Triangles"<<"=\t"<<Membrane_num_of_Triangles<<endl;
    write_initial_condition<<"Membrane_Radius"<<"=\t"<<Membrane_Radius<<endl;
    //    write_initial_condition<<"Nucleus_Membrane_radius"<<"=\t"<<Nucleus_Membrane_radius<<endl;
    write_initial_condition<<"K_surfaceConstant_local"<<"=\t"<<K_surfaceConstant_local<<endl;
    write_initial_condition<<"Membrane_spring_coefficient"<<"=\t"<<Membrane_spring_coefficient<<endl;
    write_initial_condition<<"Membrane_bending_coefficient"<<"=\t"<<Membrane_bending_coefficient<<endl;
    write_initial_condition<<"membrane_damping_coefficient"<<"=\t"<<membrane_damping_coefficient<<endl;
    write_initial_condition<<"Membrane_Node_Mass"<<"=\t"<<Membrane_Node_Mass<<endl;
    write_initial_condition<<"fluidity"<<"=\t"<<fluidity<<endl;
    write_initial_condition<<"Nucleus_Membrane_Radius_of_Hard_Sphere_Interaction"<<"=\t"<<Nucleus_Membrane_Radius_of_Hard_Sphere_Interaction<<endl;
    write_initial_condition<<"\n\n";
    write_initial_condition<<"k_actine_membrane"<<"=\t"<<k_actine_membrane<<endl;
    write_initial_condition<<"Actin_num_of_Nodes"<<"=\t"<<Actin_num_of_Nodes<<endl;
    //    write_initial_condition<<"Actin_Membrane_shared_num_of_Nodes"<<"=\t"<<Actin_Membrane_shared_num_of_Nodes<<endl;
    write_initial_condition<<"Actin_spring_coefficient"<<"=\t"<<Actin_spring_coefficient<<endl;
    write_initial_condition<<"Membrane_barrier_calculation_rate"<<"=\t"<<Membrane_barrier_calculation_rate<<endl;
    write_initial_condition<<"CytoskeletonNetworkType"<<"=\t"<<CytoskeletonNetworkType<<endl;
    write_initial_condition<<"Actin_Passive_Cable_Network_Coefficient"<<"=\t"<<Actin_Passive_Cable_Network_Coefficient<<endl;
    write_initial_condition<<"KActinACN_EA"<<"=\t"<<KActinACN_EA<<endl;
    write_initial_condition<<"ACN_TL0"<<"=\t"<<ACN_TL0<<endl;
    write_initial_condition<<"ACN_LC"<<"=\t"<<ACN_LC<<endl;
    write_initial_condition<<"Actin_kelvin_damping_coefficient"<<"=\t"<<Actin_kelvin_damping_coefficient<<endl;
    write_initial_condition<<"Actin_damping_Coefficient"<<"=\t"<<Actin_damping_Coefficient<<endl;
    write_initial_condition<<"Actin_Node_Mass"<<"=\t"<<Actin_Node_Mass<<endl;
    write_initial_condition<<"minimumlengthActin"<<"=\t"<<minimumlengthActin<<endl;
    write_initial_condition<<"maximumlengthActin"<<"=\t"<<maximumlengthActin<<endl;
    write_initial_condition<<"Actin_Membrane_Radius_of_Hard_Sphere_Interaction"<<"=\t"<<Actin_Membrane_Radius_of_Hard_Sphere_Interaction<<endl;
    write_initial_condition<<"\n\n";
    //    write_initial_condition<<"Chromatin_num_of_Beads"<<"=\t"<<Chromatin_num_of_Beads<<endl;
    //    write_initial_condition<<"Chromatin_num_of_chains"<<"=\t"<<Chromatin_num_of_chains<<endl;
    //    write_initial_condition<<"Chromatin_Bead_Mass"<<"=\t"<<Chromatin_Bead_Mass<<endl;
    //    write_initial_condition<<"Chromatin_spring_coefficient"<<"=\t"<<Chromatin_spring_coefficient<<endl;
    //    write_initial_condition<<"Chromatin_bending_coefficient"<<"=\t"<<Chromatin_bending_coefficient<<endl;
    //    write_initial_condition<<"sigmachromatin"<<"=\t"<<sigmachromatin<<endl;
    //    write_initial_condition<<"Rmaxchromatin"<<"=\t"<<Rmaxchromatin<<endl;
    //    write_initial_condition<<"Chromatin_Membrane_Radius_of_Hard_Sphere_Interaction"<<"=\t"<<Chromatin_Membrane_Radius_of_Hard_Sphere_Interaction<<endl;
    //    write_initial_condition<<"Chromatin_Scaling_Factor"<<"=\t"<<Chromatin_Scaling_Factor<<endl;
    //    write_initial_condition<<"ThermostatOnChromatin"<<"=\t"<<ThermostatOnChromatin<<endl;
    //    write_initial_condition<<"chromatin_force_cut_off"<<"=\t"<<chromatin_force_cut_off<<endl;
    //    write_initial_condition<<"chromatin_force_cut_off"<<"=\t"<<chromatin_force_cut_off<<endl;
    write_initial_condition<<"\n\n";
    //    write_initial_condition<<"ECM_num_of_Nodes"<<"=\t"<<ECM_num_of_Nodes<<endl;
    //    write_initial_condition<<"ECM_Surface_num_of_Triangles"<<"=\t"<<ECM_Surface_num_of_Triangles<<endl;
    //    write_initial_condition<<"AdhesiveIslandOffOrOn"<<"=\t"<<AdhesiveIslandOffOrOn<<endl;
    //    write_initial_condition<<"NumberOfAdhesiveIslandNodes"<<"=\t"<<NumberOfAdhesiveIslandNodes<<endl;
    //    write_initial_condition<<"ECM_Thickness"<<"=\t"<<ECM_Thickness<<endl;
    //    write_initial_condition<<"ECM_Flexibility"<<"=\t"<<ECM_Flexibility<<endl;
    //    write_initial_condition<<"sigmaECM"<<"=\t"<<sigmaECM<<endl;
    //    write_initial_condition<<"ECM_LJ_just_Repultion"<<"=\t"<<ECM_LJ_just_Repultion<<endl;
    //    write_initial_condition<<"epsilonECM"<<"=\t"<<epsilonECM<<endl;
    //    write_initial_condition<<"ECM_kelvin_damping_coefficient"<<"=\t"<<ECM_kelvin_damping_coefficient<<endl;
    //    write_initial_condition<<"miuDampECM"<<"=\t"<<miuDampECM<<endl;
    //    write_initial_condition<<"fECMcuttoff"<<"=\t"<<fECMcuttoff<<endl;
    //    write_initial_condition<<"ECM_Node_Mass"<<"=\t"<<ECM_Node_Mass<<endl;
    //    write_initial_condition<<"ECM_Min_Gradient_Coefficient"<<"=\t"<<ECM_Min_Gradient_Coefficient<<endl;
    //    write_initial_condition<<"ECM_Max_Gradient_Coefficient"<<"=\t"<<ECM_Max_Gradient_Coefficient<<endl;
    //    write_initial_condition<<"ECM_Gradient_length"<<"=\t"<<ECM_Gradient_length<<endl;
    //    write_initial_condition<<"spreading_force_magnitude"<<"=\t"<<spreading_force_magnitude<<endl;
    //    write_initial_condition<<"spreading_force_min_range"<<"=\t"<<spreading_force_min_range<<endl;
    //    write_initial_condition<<"spreading_force_max_range"<<"=\t"<<spreading_force_max_range<<endl;
    //    write_initial_condition<<"spreading_flag"<<"=\t"<<spreading_flag<<endl;
    write_initial_condition<<"spreading_force_cos_triangle_interaction_angle"<<"=\t"<<spreading_force_cos_triangle_interaction_angle<<endl;
    write_initial_condition<<"\n\n";
    write_initial_condition<<"Actin_membrane_stiff_spring_coefficient"<<"=\t"<<Actin_membrane_stiff_spring_coefficient<<endl;
    write_initial_condition<<"Actin_membrane_damping_coefficient"<<"=\t"<<Actin_membrane_damping_coefficient<<endl;
    
}






