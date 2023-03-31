#ifndef HEATDIFFUSION_H
#define HEATDIFFUSION_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// what i need

#include <iostream>
#include <vector>
#include <fstream>
#include "init.h"

using namespace std;

extern double dx;
extern double dt;
extern double time_elapsed;
extern int point_nu;
extern int step_nu;

void TempDistriInit ();
void TempDistri ();
void TempPrintClose ();

class block
{

public:

// Position
int x_pos;
int y_pos;
int z_pos;

// blocks near it

int plus_z;
int minus_z;
int plus_y;
int minus_y;
int plus_x;
int minus_x;

int number_of_connections;

// HEAT

double heat;
double new_heat;

// K constant (diffition constant)
double diffution_rate;

//constructor
block(int x, int y, int z, double temp, double diff);
//print block data
void print();
// draw block
void draw_one_cube();
//

};


class room
{

public:
	
// vector of blocks
vector<block> room_blocks;
//constructor
room(int tot_z, double T, double k);
//print room 
void view();
// check to see if overlapping positions
void check();
// look up position and see if block in it. return number of block
int look_up_pos(int x, int y, int z);
// calc bools
void find_bools();
// calc new heat for one block in vector
void new_heat(int i);
// calc temp for entire vector (easyer to break it up i thought)
void run_new_heat();
// make old temp new temp and new temp 0
void run_once();
// run one time step
void one_step();
// print to file
void print_open(string file_name);
void print();
void print_close();
// testing ----------------------------
void draw();
// testing ----------------------------
void setupDrawCallback();
// double testing
//static void test();

};

#endif
