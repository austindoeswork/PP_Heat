/*********************************************************/
/* Parallel Heat Transfer Simulation in Three Dimensions */
/*   by Aidan Wenzel, Austin Wilson, and Theodore Rice   */
/*********************************************************/

/*********************************************************/
/* Includes **********************************************/
/*********************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#include<string.h>
#include<unistd.h>

/*********************************************************/
/* Define Macros *****************************************/
/*********************************************************/

#define TICK 1
#define cubeDim 1
#define DIFFU 0.1

/*********************************************************/
/* Global Variable Definitions ***************************/
/*********************************************************/

// Defintion of our object which has a currTemp.
// If needed, this struct could be enhanced to
// include diffusivity, and other variables to
// make it more accurate.
typedef struct {
    double currTemp;
} object;

// The universe as it stands, and the next iteration
object* universe;
object* universeNext;

MPI_Status status;

// Global variables that show worldsize,
// ranks around, ticks, slice size, and
// dimenstions of the world
int worldsize, myrank, aboveRank, belowRank;
int numTicks, sliceSize, printOnTick;
int dimX, dimY, dimZ;

/*********************************************************/
/* Function Definitions **********************************/
/*********************************************************/
void printToConsole(int tick);

// Function that prints the current stage of the universe
// when MPI Rank 0 reaches
void printUniverse(int tick){
    if (myrank == 0) {
        printf("%d,%d,%d,%d\n", dimX,dimY,dimZ,tick);
        for (int z = 1; z < dimZ/worldsize; z++) {
            for (int y = 0; y < dimY; y++) {
                for (int x = 0; x < dimX; x++) {
                    if (x == dimX - 1) {printf("%f", (universe + (dimX*dimY*(z+1)) + (dimX * y) + x)->currTemp);}
                    else {printf("%f,", (universe + (dimX*dimY*(z+1)) + (dimX * y) + x)->currTemp);}
                }
                printf("\n");
            }
        }
    }
}

// Allocates memory for next tick of universe
object* emptyUniverse(){
    return calloc((size_t)dimX*dimY*((dimZ/worldsize)+2), sizeof(object));
}

// Opens the initialization file and set objects initial heat
void initializeUniverse(char* filename){
    universe = emptyUniverse();
    FILE *file;
    file = fopen(filename, "r");
    int localZStart = myrank * (dimZ/worldsize);
    if ( file != NULL ) {
        char line[1024];
        while (fgets(line, sizeof(line), file) != NULL) {
            if (line[0]=='/') {
                continue;
            }
            int x1, y1, z1, x2, y2, z2;
            double curTemp, thermCond;
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

            if(z1 < localZStart) z1 = 0;
            else z1 = z1 - localZStart;
            if(z2 > localZStart + dimZ/worldsize - 1) z2 = dimZ/worldsize - 1;
            else z2 = z2 - localZStart;

            if(x2 >= dimX) x2 = dimX - 1;
            if(y2 >= dimY) y2 = dimY - 1;

            for (int x = x1; x <= x2; x++) {
                for (int y = y1; y <= y2; y++) {
                    for (int z = z1; z <= z2; z++) {
                        object *target = universe + (dimX*dimY*(z+1)) + (dimX * y) + x;
                        target->currTemp = curTemp;
                    }
                }
            }
        }
        fclose(file);
    }
}

// Completes a tick on the universe. This includes receiving and sending
// information on ghost rows, and then applying the heat transfer formulas
// to objects around the target.
void tick(int tickNum){
    universeNext = emptyUniverse();

    MPI_Request sendGhostBack, sendGhostFront, receiveGhostBack, receiveGhostFront;
    MPI_Isend(universe+(dimX*dimY), dimX*dimY, MPI_DOUBLE, belowRank, 0, MPI_COMM_WORLD, &sendGhostBack);
    MPI_Isend(universe+sliceSize, dimX*dimY, MPI_DOUBLE, aboveRank, 1, MPI_COMM_WORLD, &sendGhostFront);

    MPI_Irecv(universe+(dimX*dimY+sliceSize), dimX*dimY, MPI_DOUBLE, belowRank, 1, MPI_COMM_WORLD, &receiveGhostFront);
    MPI_Irecv(universe, dimX*dimY, MPI_DOUBLE, aboveRank, 0, MPI_COMM_WORLD, &receiveGhostBack);

    MPI_Wait(&receiveGhostFront, &status);
    MPI_Wait(&receiveGhostBack, &status);

    // Apply tick to universe
    for(int z = 1; z <= dimZ/worldsize; z++){
        for(int y = 0; y < dimY; y++){
            for(int x = 0; x < dimX; x++){

                // Obtain pointer information about the target, and the objects
                // around the target.
                object *target, *targetNext, *above, *below, *left, *right, *front, *back;
                target = universe + (dimX*dimY*z) + (dimX * y) + x;
                targetNext = universeNext + (dimX*dimY*z) + (dimX * y) + x;
                above = target - dimX;
                below = target + dimX;
                left = target - 1;
                right = target + 1;
                front = target - (dimX*dimY);
                back = target + (dimX*dimY);

                targetNext->currTemp = target->currTemp;

                // Handle boundary edge cases.
                if(x != 0) targetNext->currTemp += DIFFU * (left->currTemp - target->currTemp);
                if(x != dimX-1) targetNext->currTemp += DIFFU * (right->currTemp - target->currTemp);
                if(y != 0) targetNext->currTemp += DIFFU * (above->currTemp - target->currTemp);
                if(y != dimY-1) targetNext->currTemp += DIFFU * (below->currTemp - target->currTemp);
                if(!((myrank == 0) && (z==1))) targetNext->currTemp += DIFFU * (front->currTemp - target->currTemp);
                if(!((myrank == worldsize-1)&&(z != dimZ)))
                    targetNext->currTemp += DIFFU * (back->currTemp - target->currTemp);
            }
        }
    }

    // Update the universe, and free old memory
    free(universe);
    universe = universeNext;
}

// Useful debugging function which prints every rank to console. Unused in main.
void printToConsole(int tick){
    MPI_Barrier(MPI_COMM_WORLD);
    for(int printFromRank = 0; printFromRank < worldsize; printFromRank++) {
        for (int z = 1; z <= dimZ / worldsize; z++) {
            for (int y = 0; y < dimY; y++) {
                for (int x = 0; x < dimX; x++) {
                    object *target = universe + (dimX * dimY * z) + (dimX * y) + x;
                    if(myrank == printFromRank) printf("%.2f ", target->currTemp);
                }
                if(myrank == printFromRank) printf("\n");
            }
            if(myrank == printFromRank) printf("         Printed from %d Tick:%d\n", myrank, tick);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    if(myrank == 0) printf("------------------------------------\n\n");
}

/*********************************************************/
/* Function Main *****************************************/
/*********************************************************/
int main(int argc, char* argv[]){

    // Start MPI and calculate the side ranks
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &worldsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    aboveRank = myrank - 1;
    belowRank = myrank + 1;
    if(aboveRank < 0){aboveRank = worldsize - 1;}
    if(belowRank >= worldsize){belowRank = 0;}

    // Parse arguments
    dimX = atoi(argv[1]);
    dimY = atoi(argv[2]);
    dimZ = atoi(argv[3]);
    char* iniFilename = argv[4];
    numTicks = atoi(argv[5]);
    printOnTick = atoi(argv[6]);

    // Compute slice size
    sliceSize = dimX*dimY*(dimZ/worldsize);

    // Read in the initial universe state
    initializeUniverse(iniFilename);

    // Prep MPI_time
    double start_time, total_time;
    MPI_Barrier(MPI_COMM_WORLD);
    if(myrank == 0) start_time = MPI_Wtime();

    // Run simulation
    printToConsole(-1);
    for(int tickCount = 0; tickCount < numTicks; tickCount++){
        tick(tickCount);
        if(tickCount%printOnTick == 0) printUniverse(tickCount);
    }

    // Finish time and output information
    MPI_Barrier(MPI_COMM_WORLD);
    if(myrank == 0){
        total_time = MPI_Wtime() - start_time;
        printf("MPI world size: %d  Number of ticks: %d Size of Universe: %d  Runtime: %lf\n", worldsize, numTicks,
               dimX * dimY * dimZ / worldsize, total_time);
    }

    //Finalize and exit
    MPI_Finalize();
    return 0;
}
