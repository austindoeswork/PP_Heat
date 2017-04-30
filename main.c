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
#include<string.h>
#include <unistd.h>

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
void printUniverse(int tick){
    FILE *fp;
    fp = fopen("output.txt", "w+");
    fprintf(fp, "%d,%d,%d,%d\n", dimX,dimY,dimZ,tick);
    for (int z = 0; z < dimZ; z++) {
        for (int y = 0; y < dimY; y++) {
            for (int x = 0; x < dimX; x++) {
                if (x == dimX - 1) {fprintf(fp, "%f", (universe+x+y+z)->currTemp);}
                else {fprintf(fp, "%f,", (universe+x+y+z)->currTemp);}
            }
            fprintf(fp, "\n");
        }
    }
    fprintf(fp, "\n");
    for (int z = 0; z < dimZ; z++) {
        for (int y = 0; y < dimY; y++) {
            for (int x = 0; x < dimX; x++) {
                if (x == dimX - 1) {fprintf(fp, "%f", (universe+x+y+z)->thermCond);}
                else {fprintf(fp, "%f,", (universe+x+y+z)->thermCond);}
            }
            fprintf(fp, "\n");
        }
    }
}

//allocates memory for next tick of universe
object* emptyUniverse(){
    return calloc((size_t)dimX*dimY*(dimZ+2), sizeof(object));
}

//initializes universe to ini file info
void initializeUniverse(char* filename){
    universe = (object*) calloc((size_t)dimX*dimY*(dimZ+2), sizeof(object));
    FILE *file;
    file = fopen(filename, "r");
    if ( file != NULL )
    {
        char line[1024]; /* or other suitable maximum line size */
        while (fgets(line, sizeof line, file) != NULL) /* read a line */
        {
            if (line[0]=='/') {
                continue;
            }
            int x1, y1, z1, x2, y2, z2;
            float curTemp, thermCond;
            char* token = strtok(line, " ");
            x1 = atoi(token);
            token = strtok(NULL, " ");
            y1 = atoi(token);
            token = strtok(NULL, " ");
            z1 = atoi(token);
            token = strtok(NULL, " ");
            x2 = atoi(token);
            token = strtok(NULL, " ");
            y2 = atoi(token);
            token = strtok(NULL, " ");
            z2 = atoi(token);
            token = strtok(NULL, " ");
            curTemp = atof(token);
            token = strtok(NULL, " ");
            thermCond = atof(token);

            for (int x = x1; x < x2; x++) {
                for (int y = y1; y < y2; y++) {
                    for (int z = z1; z < z2; z++) {
                        (universe+x+y+z)->currTemp = curTemp;
                        (universe+x+y+z)->thermCond = thermCond;
                    }
                }
            }
        }
        fclose(file);
    }
}

// completes a tick on the universe
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
    // start MPI
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
    for(int i = 0; i < dimX; i++){
        for(int j = 0; j < dimY; j++){
            for(int k = 1; k <= dimZ; k++){
                printf("%lf ", universe);
            }
        }
    }
    /*universe = emptyUniverse();
    for(int i = 0; i < dimX*dimY*(dimZ/worldsize+2); i++){

    }*/

    //Prep MPI_time stuff
    double start_time, total_time;
    MPI_Barrier(MPI_COMM_WORLD);
    if(myrank == 0) start_time = MPI_Wtime();

    //Run simulation
    for(int tickCount = 0; tickCount < numTicks; tickCount++){
        tick();
        if(tickCount%printOnTick == 0) printUniverse(tickCount);
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
