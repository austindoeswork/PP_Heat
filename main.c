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
#define cubeDim 1

/*********************************************************/
/* Global Variable Definitions ***************************/
/*********************************************************/

typedef struct {
    float currTemp;
    float thermCond;
} object;

object* universe;
object* universeNext;

MPI_Status status;

int worldsize, myrank, aboveRank, belowRank;
int numTicks, sliceSize, printOnTick;
int dimX, dimY, dimZ; //Dimensions of board

/*********************************************************/
/* Function Definitions **********************************/
/*********************************************************/
void printUniverse(){
    //TODO: Write file out to print out current universe
}

//allocates memory for next tick of universe
object* emptyUniverse(){
    return calloc((size_t)dimX*dimY*(dimZ+2), sizeof(object));
}

//initializes universe to ini file info
void initializeUniverse(char* filename){
    universe = (object*) calloc((size_t)dimX*dimY*(dimZ+2), sizeof(object));

    //TODO: write the read file to read in objects into grid
}

//completes a tick on the universe
void tick(){
    universeNext = emptyUniverse();

    MPI_Request sendGhostBelow, sendGhostAbove, receiveGhostBelow, receiveGhostAbove;
    MPI_Isend(universe+(dimX*dimY), dimX*dimY, MPI_UNSIGNED_SHORT, belowRank, 0, MPI_COMM_WORLD, &sendGhostBelow);
    MPI_Isend(universe+sliceSize, dimX*dimY, MPI_UNSIGNED_SHORT, aboveRank, 1, MPI_COMM_WORLD, &sendGhostAbove);

    MPI_Irecv(universe+(dimX*dimY+sliceSize), dimX*dimY, MPI_UNSIGNED_SHORT, aboveRank, 0, MPI_COMM_WORLD, &receiveGhostAbove);
    MPI_Irecv(universe, dimX*dimY, MPI_UNSIGNED_SHORT, belowRank, 1, MPI_COMM_WORLD, &receiveGhostBelow);

    MPI_Wait(&receiveGhostAbove, &status);
    MPI_Wait(&receiveGhostBelow, &status);

    //apply tick to universe
    for(int x = 0; x < dimX; x++){
        for(int y = 0; y < dimY; y++){
            for(int z = 1; z <= dimZ/worldsize; z++){//only as much depth as this rank handles
                //obtain pointer information
                object *target, *targetNext, *above, *below, *left, *right, *front, *back;
                target = universe + (dimX*dimY*z) + (dimX * y) + x;
                targetNext = universeNext + (dimX*dimY*z) + (dimX * y) + x;
                above = targetNext - dimX;
                below = targetNext + dimX;
                left = targetNext - 1;
                right = targetNext + 1;
                front = targetNext - (dimX*dimY);
                back = targetNext + (dimX*dimY);

                //Calculate next tick
                targetNext->thermCond = target->thermCond;
                targetNext->currTemp =
                          2/((1/target->thermCond)+(1/above->thermCond))*(target->currTemp-above->currTemp)
                        + 2/((1/target->thermCond)+(1/below->thermCond))*(target->currTemp-below->currTemp)
                        + 2/((1/target->thermCond)+(1/left->thermCond))*(target->currTemp-left->currTemp)
                        + 2/((1/target->thermCond)+(1/right->thermCond))*(target->currTemp-right->currTemp)
                        + 2/((1/target->thermCond)+(1/front->thermCond))*(target->currTemp-front->currTemp)
                        + 2/((1/target->thermCond)+(1/back->thermCond))*(target->currTemp-back->currTemp);
            }
        }
    }

    //update universe, free old memory
    free(universe);
    universe = universeNext;
}

/*********************************************************/
/* Function Main *****************************************/
/*********************************************************/
int main(int argc, char* argv[]){
    //start MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &worldsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    aboveRank = myrank - 1;
    belowRank = myrank + 1;
    if(aboveRank < 0){aboveRank = worldsize - 1;}
    if(belowRank >= worldsize){belowRank = 0;}

    //parse argv
    dimX = atoi(argv[1]);
    dimY = atoi(argv[2]);
    dimZ = atoi(argv[3]);
    char* iniFilename = argv[4];
    numTicks = atoi(argv[5]);
    printOnTick = atoi(argv[6]);

    //compute other global vars
    sliceSize = dimX*dimY*(dimZ/worldsize);

    //read in the initial universe state
    initializeUniverse(iniFilename);

    //Prep MPI_time stuff
    double start_time, total_time;
    MPI_Barrier(MPI_COMM_WORLD);
    if(myrank == 0) start_time = MPI_Wtime();

    //Run simulation
    for(int tickCount = 0; tickCount < numTicks; tickCount++){
        tick();
        if(tickCount%printOnTick == 0) printUniverse();
    }

    //Finish time and output info
    MPI_Barrier(MPI_COMM_WORLD);
    if(myrank == 0){
        total_time = MPI_Wtime() - start_time;
        printf("MPI world size: %d  Number of ticks: %d  Runtime: %lf\n", worldsize, numTicks, total_time);
    }

    //Finalize and exit
    MPI_Finalize();
    return 0;
}