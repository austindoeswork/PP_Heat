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
/* Global Variable Definitions ***************************/
/*********************************************************/

typedef struct {
    float currTemp;
    float thermCond;
} object;

object* universe, universeNext; //X,Y,Z

int worldsize, myrank, aboveRank, belowRank;
int numTicks;
int dimX, dimY, dimZ; //Dimensions of board

/*********************************************************/
/* Function Definitions **********************************/
/*********************************************************/
void printUniverse();

object* emptyUniverse(){
    return calloc((size_t)dimX*dimY*dimZ, sizeof(object));
}

void initializeUniverse(char* filename){
    universe = (object*) calloc((size_t)dimX*dimY*dimZ, sizeof(object));

    //TODO: write the read file to read in objects into grid
}

void tick();

/*********************************************************/
/* Function Main *****************************************/
/*********************************************************/
int main(int argc, char* argv[]){

}