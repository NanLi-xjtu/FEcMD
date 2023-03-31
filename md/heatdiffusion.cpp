/**************************************************************************
** This is my file with the heart of the heat diffution simulator. 
** It is where the explicit finite difference method is implemented.
**************************************************************************/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include "heatdiffusion.h"

extern "C"{
#include "main.h"
}

using namespace std;

room *temp_cell;
ofstream reading;//print temperature distribution
double dz;
double r;
double dt;
double time_elapsed = 0;
int point_nu;
int step_nu;
string file_name;

//----------- functions for block ---------------
//constructor
block::block(int x, int y, int z, double temp, double diff)
{
	x_pos = x; 
	y_pos = y;
	z_pos = z;
	
	// -1 is saying there is no block on that side
	plus_z = -1;
	minus_z = -1;
	plus_y = -1;
	minus_y = -1;
	plus_x = -1;
	minus_x = -1;
	
	number_of_connections = 0;

	heat = temp;
	new_heat = 0;

	diffution_rate = diff;
}

void block::print()
{
	cout << "x_pos = " << x_pos << ", y_pos = " << y_pos << ", z_pos = " << z_pos << ", temp = " << heat << ", new temp = " << new_heat << ",    plus_x = " << plus_x << ", minus_x = " << minus_x << ", plus_y = " << plus_y << ", minus_y = " << minus_y << ", plus_z = " << plus_z << ", minus_z = " << minus_z << ",    connections = " << number_of_connections << "\n";	
}

//----------- functions for room -----------------
// constructor
room::room(int tot_z, double T, double k)
{
	// create vector	
	int x = 0;
	int y = 0;
	int z = 0;

	for (z = 0; z < sizeTempGrid; z ++){	
		block f(x, y, z, T, k);
		room_blocks.push_back(f);
	}
	
	find_bools();
}
// prints position of data of blocks
void room::view()
{
	for(int i = 0; i < room_blocks.size(); i++)
	{
		room_blocks[i].print();
	}
}
// checks overlapping position of blocks
void room::check()
{	
	// not yet :P
}
// look up a position and return the number block in it. if none return -1
int room::look_up_pos(int x, int y, int z)
{
	for(int i = 0; i < room_blocks.size(); i++)
	{
		if((room_blocks[i].x_pos == x) && (room_blocks[i].y_pos == y) && (room_blocks[i].z_pos == z))
		{
			return i;
		}
	}
	return -1;
}
// finds the correct bool values
void room::find_bools()
{
	for(int i = 0; i < room_blocks.size(); i++)
	{
		room_blocks[i].plus_z = look_up_pos(room_blocks[i].x_pos, room_blocks[i].y_pos, room_blocks[i].z_pos + 1);
		room_blocks[i].minus_z = look_up_pos(room_blocks[i].x_pos, room_blocks[i].y_pos, room_blocks[i].z_pos - 1);
		room_blocks[i].plus_y = look_up_pos(room_blocks[i].x_pos, room_blocks[i].y_pos + 1, room_blocks[i].z_pos);
		room_blocks[i].minus_y = look_up_pos(room_blocks[i].x_pos, room_blocks[i].y_pos - 1, room_blocks[i].z_pos);
		room_blocks[i].plus_x = look_up_pos(room_blocks[i].x_pos + 1, room_blocks[i].y_pos, room_blocks[i].z_pos);
		room_blocks[i].minus_x = look_up_pos(room_blocks[i].x_pos - 1, room_blocks[i].y_pos, room_blocks[i].z_pos);
		
		if(room_blocks[i].plus_z != -1) {room_blocks[i].number_of_connections = room_blocks[i].number_of_connections + 1;}
		if(room_blocks[i].minus_z != -1) {room_blocks[i].number_of_connections = room_blocks[i].number_of_connections + 1;}
		if(room_blocks[i].plus_y != -1) {room_blocks[i].number_of_connections = room_blocks[i].number_of_connections + 1;}
		if(room_blocks[i].minus_y != -1) {room_blocks[i].number_of_connections = room_blocks[i].number_of_connections + 1;}
		if(room_blocks[i].plus_x != -1) {room_blocks[i].number_of_connections = room_blocks[i].number_of_connections + 1;}
		if(room_blocks[i].minus_x != -1) {room_blocks[i].number_of_connections = room_blocks[i].number_of_connections + 1;}
	}
}
// calc new heat for a given block
void room::new_heat(int i)
{
	r = room_blocks[i].diffution_rate*dt/(dz*dz);
	room_blocks[i].new_heat = 0;
	room_blocks[i].new_heat = room_blocks[i].heat * (1 - (room_blocks[i].number_of_connections * r));
	if(room_blocks[i].plus_z != -1) {room_blocks[i].new_heat = room_blocks[i].new_heat + (room_blocks[room_blocks[i].plus_z].heat * r);}
	if(room_blocks[i].minus_z != -1) {room_blocks[i].new_heat = room_blocks[i].new_heat + (room_blocks[room_blocks[i].minus_z].heat * r);}
	if(room_blocks[i].plus_y != -1) {room_blocks[i].new_heat = room_blocks[i].new_heat + (room_blocks[room_blocks[i].plus_y].heat * r);}
	if(room_blocks[i].minus_y != -1) {room_blocks[i].new_heat = room_blocks[i].new_heat + (room_blocks[room_blocks[i].minus_y].heat * r);}
	if(room_blocks[i].plus_x != -1) {room_blocks[i].new_heat = room_blocks[i].new_heat + (room_blocks[room_blocks[i].plus_x].heat * r);}
	if(room_blocks[i].minus_x != -1) {room_blocks[i].new_heat = room_blocks[i].new_heat + (room_blocks[room_blocks[i].minus_x].heat * r);}
	 
}
// running the new heat for the whole thing
void room::run_new_heat()
{
	for(int i = 0; i < room_blocks.size(); i++)
	{
		new_heat(i);
	}
}
// make new temp old temp and new temp 0
void room::run_once()
{
	for(int i = 0; i < room_blocks.size(); i++)
	{
		room_blocks[i].heat = room_blocks[i].new_heat;
	}
	room_blocks[room_blocks.size() - 1].heat += depositedHeat * dt / dz;
}

// run one time step
void room::one_step()
{
	run_new_heat();
	run_once();
}

//print it to a file
void room::print_open(string file_name)
{
	// create stream to file
	reading.open(file_name.c_str());
}

void room::print()
{
	// create vector	
	//reading << dz << " " << dt << "\n";
        reading << point_nu << "\n";
        reading << "SolutionReader properties=id:I:1:pos:R:3:temperature:R:1\n";
	for(int i = 0; i < point_nu; i++)
	{
		reading << i << " " <<room_blocks[i].x_pos << " " << room_blocks[i].y_pos << " " << room_blocks[i].z_pos << " " << room_blocks[i].heat * TUnit << "\n";
	}		
}

void room::print_close()
{
	reading.close();
}

void TempDistriInit()
{
	point_nu = sizeTempGrid;
	dz = region.z / (point_nu - 1);
	dt = deltaT;
	r = heatConduct * dt / ( dz * dz );
	step_nu = stepLimit;

	//build temperature cells
	room making(sizeTempGrid, temperature, heatConduct);
	temp_cell = &making;
//	temp_cell->view();

	temp_cell->print_open ("./out/tempdistri.movie");

	//tempcell
	int n, i;
	DO_TEMPCELL tempCell[i].temp = temp_cell->room_blocks[i].heat;
	CalcuTempCell ();

	/**********************Print**********************************/
	//print initial tempcells
	PrintTempCell ();
	//print initial temperature distribution
	temp_cell->print ();
	cout << "\nheat diffusion -- data\n";
	cout << "dz = " << dz * lUnit * 1.e10 << " A, dt = " << dt * tUnit * 1.e15 \
	     << " fs, r = " << r << ", total steps = " << step_nu << ", total cells = " << point_nu << "\n" << "----" << endl;
}

//onestep of calculating temperature distribution
void TempDistri()
{
	//calculate heat diffusion equation
	if (stepCount >= stepEquil) temp_cell->one_step ();

	int i;
	DO_TEMPCELL tempCell[i].temp = temp_cell->room_blocks[i].heat;

	//print temperature distribution
	temp_cell->print ();

}

void TempPrintClose ()
{
	temp_cell->print_close ();
}
