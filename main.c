//
// PP_Heat a 2D Parallel Heat Simulator
//

/*********************************************************/
/* Includes **********************************************/
/*********************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>

/*********************************************************/
/* Define Macros *****************************************/
/*********************************************************/

#define TICK 1 //length of time for a tick

/*********************************************************/
/* Global Vaiable Definitions ****************************/
/*********************************************************/

float* universe, universeNext; //X,Y,Z

int worldsize, myrank, aboveRank, belowRank;
int numTicks;
int dimX, dimY, dimZ; //Dimensions of board

/*********************************************************/
/* Function Definitions **********************************/
/*********************************************************/
void printUniverse();

float* emptyUniverse(){
    return calloc((size_t)dimX*dimY*dimZ, sizeof(float));
}

void initializeUniverse(char* filename){
    universe = (float*) calloc((size_t)dimX*dimY*dimZ, sizeof(float));

    //TODO: write the read file to read in objects into grid
}

void tick();

/*********************************************************/
/* Function Main *****************************************/
/*********************************************************/
int main(int argc, char* argv[]){

}