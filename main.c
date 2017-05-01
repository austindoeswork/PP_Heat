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
#include<unistd.h>

/*********************************************************/
/* Define Macros *****************************************/
/*********************************************************/

#define TICK 1 //length of time for a tick
#define cubeDim 1
#define DIFFU 0.1

/*********************************************************/
/* Global Variable Definitions ***************************/
/*********************************************************/

typedef struct {
    double currTemp;
    //double thermCond;
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
void printToConsole(int tick);

void printUniverse(int tick){
    MPI_File fp;
    MPI_File_open(MPI_COMM_WORLD, "output.txt", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fp);

    MPI_Status status;

    // printf("slice size %d myrank %d sizeof(double) %lu\n", sliceSize, myrank, sizeof(double));

    // printf("%lu %d\n", sizeof(double)*sliceSize*myrank, myrank);

    MPI_File_write_at_all(fp, sizeof(double)*sliceSize*myrank, universe, sizeof(double)*sliceSize, MPI_DOUBLE, &status);

    // if (myrank == 0) {
    //     FILE *fp;
    //     fp = fopen("output.txt", "w+");
    //     fprintf(fp, "%d,%d,%d,%d\n", dimX,dimY,dimZ,tick);
    //     for (int z = 0; z < dimZ; z++) {
    //         for (int y = 0; y < dimY; y++) {
    //             for (int x = 0; x < dimX; x++) {
    //                 if (x == dimX - 1) {fprintf(fp, "%f", (universe+x+y+z)->currTemp);}
    //                 else {fprintf(fp, "%f,", (universe+x+y+z)->currTemp);}
    //             }
    //             fprintf(fp, "\n");
    //         }
    //     }
    //     fprintf(fp, "\n");
    //     for (int z = 0; z < dimZ; z++) {
    //         for (int y = 0; y < dimY; y++) {
    //             for (int x = 0; x < dimX; x++) {
    //                 if (x == dimX - 1) {fprintf(fp, "%f", (universe+x+y+z)->thermCond);}
    //                 else {fprintf(fp, "%f,", (universe+x+y+z)->thermCond);}
    //             }
    //             fprintf(fp, "\n");
    //
    //         }
    //     }
    // }
}

//allocates memory for next tick of universe
object* emptyUniverse(){
    return calloc((size_t)dimX*dimY*(dimZ+2), sizeof(object));
}

void placeObjectInUniverse(int x1, int y1, int z1, int x2, int y2, int z2, double temp){
    int localZStart = myrank * (dimZ/worldsize); //Ex rank 0 starts at Z = 0

    if(z1 < localZStart) z1 = 0;
    else z1 = z1 - localZStart;
    if(z2 > localZStart + dimZ/worldsize - 1) z2 = dimZ/worldsize - 1;
    else z2 = z2 - localZStart;

    for (int x = x1; x <= x2; x++) {
        for (int y = y1; y <= y2; y++) {
            for (int z = z1; z < z2; z++) {
                object *target = universe + (dimX*dimY*(z+1)) + (dimX * y) + x;
                target->currTemp = temp;
                //                     //(universe+x+y+z)->thermCond = thermCond;
            }
        }
    }
}

//initializes universe to ini file info
void initializeUniverse(char* filename){
    universe = (object*) calloc((size_t)dimX*dimY*(dimZ+2), sizeof(object));
    FILE *file;
    file = fopen(filename, "r");
    int localZStart = myrank * (dimZ/worldsize); //Ex rank 0 starts at Z = 0
    if ( file != NULL )
    {
        char line[1024]; /* or other suitable maximum line size */
        while (fgets(line, sizeof(line), file) != NULL) /* read a line */
        {
            fprintf(stderr, "Loop starts\n");
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
    //                     //(universe+x+y+z)->thermCond = thermCond;
                    }
                }
            }
        }
        fclose(file);
    }
}

// completes a tick on the universe
void tick(int tickNum){
    universeNext = emptyUniverse();

    //printToConsole(tickNum);

    MPI_Request sendGhostBack, sendGhostFront, receiveGhostBack, receiveGhostFront;
    MPI_Isend(universe+(dimX*dimY), dimX*dimY, MPI_DOUBLE, belowRank, 0, MPI_COMM_WORLD, &sendGhostBack);
    MPI_Isend(universe+sliceSize, dimX*dimY, MPI_DOUBLE, aboveRank, 1, MPI_COMM_WORLD, &sendGhostFront);

    MPI_Irecv(universe+(dimX*dimY+sliceSize), dimX*dimY, MPI_DOUBLE, belowRank, 1, MPI_COMM_WORLD, &receiveGhostFront);
    MPI_Irecv(universe, dimX*dimY, MPI_DOUBLE, aboveRank, 0, MPI_COMM_WORLD, &receiveGhostBack);

    MPI_Wait(&receiveGhostFront, &status);
    MPI_Wait(&receiveGhostBack, &status);

    //printToConsole(tickNum);

    //apply tick to universe
    for(int z = 1; z <= dimZ/worldsize; z++){//only as much depth as this rank handles
        for(int y = 0; y < dimY; y++){
            for(int x = 0; x < dimX; x++){
                //obtain pointer information
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
                //handle edge cases
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

    //update universe, free old memory
    free(universe);
    universe = universeNext;
}

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
    //initializeUniverse(iniFilename);

    universe = emptyUniverse();
    //placeObjectInUniverse(0, 0, 0, 5000, 5000, 5000);
   // for(int i = 0; i < dimX*dimY*(dimZ/worldsize+2); i++){
   //     (universe+i)->currTemp = 0;
   //  //    (universe+i)->thermCond = 1;
   // }
   // if(myrank == 1){(universe+13)->currTemp = 50;}

    //Prep MPI_time stuff
    double start_time, total_time;
    MPI_Barrier(MPI_COMM_WORLD);
    if(myrank == 0) start_time = MPI_Wtime();

    //Run simulation
    //printToConsole(-1);
    for(int tickCount = 0; tickCount < numTicks; tickCount++){
        tick(tickCount);
        //if(tickCount%printOnTick == 0) printToConsole(tickCount);//printUniverse(tickCount);
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
