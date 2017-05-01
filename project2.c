/****************************************************************************
 * HEAT2D Example - Parallelized C Version
 * FILE: mpi_heat2D.c
 * DESCRIPTIONS: This example is based on a simplified two-dimensional heat
		* equation domain decomposition. The initial temperature is computed to be
		* high in the middle of the domain and zero at the boundaries. The
		* boundaries are held at zero throughout the simulation. During the
		* time-stepping, an array containing two domains is used; these domains
		* alternate between old data and new data.
		*
		* In this parallelized version, the grid is decomposed by the master
		* process and then distributed by rows to the worker processes. At each
		* time step, worker processes must exchange border data with neighbors,
		* because a grid point's current temperature depends upon it's previous
		* time step value plus the values of the neighboring grid points. Upon
		* completion of all time steps, the worker processes return results
		* to the master process.
		*
		* AUTHOR: Blaise Barney - adapted from D. Turner's serial version
		* CONVERTED TO MPI: George L. Gusciora (1/25/95)
		* MODIFIED BY: C. B. Connor (6/6/02)
		****************************************************************************/
#include <stdio.h>
#include <mpi.h>

#define NXPROB 50 /* x dimension of problem grid */
#define NYPROB 50 /* y dimension of problem grid */
#define TIME_STEPS 500 /* number of time steps */
#define MAXWORKER 8 /* maximum number of worker tasks */
#define MINWORKER 3 /* minimum number of worker tasks */
#define BEGIN 1 /* message type */
#define NGHBOR1 2 /* message type */
#define NGHBOR2 3 /* message type */
#define NONE 0 /* indicates no neighbor */
#define DONE 4 /* message type */
#define MASTER 0 /* taskid of first process */

struct Parms {
  float cx;
  float cy;
} diffusivity = {0.1, 0.1};

int main(argc,argv)
     int argc;
     char *argv[];
{
  void inidat(), prtdat(), update();
  float u[2][NXPROB][NYPROB]; /* array for grid */
  int taskid; /* this task's unique id */
  int   numworkers; /* number of worker processes */
  int   numtasks; /* number of tasks */
  int   min_number_rows,number_rows,offset,extra_rows;/* for sending rows of data */
  int   destination, source; /* to - from for message send-receive */
  int   worker_number, neighbor1,neighbor2; /* neighbor tasks */
  int   message_tag; /* for message types */
  int   nbytes; /* number of bytes received */
  int   rc,start,end; /* misc */
  int   i,ix,iy,iz,it; /* loop variables */
  MPI_Status status;
  
  
  /* First, find out my taskid and how many tasks are running */
  rc = MPI_Init(&argc,&argv);
  rc|= MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  rc|= MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
  if (rc != 0)
    printf ("error initializing MPI and obtaining task ID information\n");

  /* one of the processors is completely reserved for coordination */
  numworkers = numtasks-1;

  /* the MASTER node subdivides the problem by grid decomposition,
in this case by rows. The MASTER also initializes the starting values, sets the
boundary conditions on the problem , and receives results from the
nodes after TIME_STEPS */
  if (taskid == MASTER)
    {
      /************************* master code *******************************/
      /* Check if numworkers is within range - quit if not */
      if ((numworkers > MAXWORKER) || (numworkers < MINWORKER))
	{
	  printf("MP_PROCS needs to be between %d and %d for this exercise\n", MINWORKER+1,MAXWORKER+1);
	  MPI_Finalize();
	}
      
      /* Initialize grid */
      printf("Grid size: X= %d Y= %d Time steps= %d\n",NXPROB,NYPROB,TIME_STEPS);
      printf("Initializing grid and writing initial.dat file...\n");
      inidat(NXPROB, NYPROB, u);
      prtdat(NXPROB, NYPROB, u, "output/initial.dat");
      
      /* Distribute work to workers. Must first figure out how many rows to
	 send and what to do with extra rows. */
      min_number_rows = NXPROB/numworkers;
      extra_rows = NXPROB%numworkers;
      offset = 0;
      for (i=1; i<=numworkers; i++)
	{
	  
	  /* The following is a particularly compact and efficient  way of distributing the 
	     grid. It assures that the number of rows received by each node is only up to one more
	     than any other node. It can be read as:
	     if (i<=extra_rows) number_rows = min_number_rows + 1
	     else number_rows = min_number_rows			  
	  */
	  number_rows = (i <= extra_rows) ? min_number_rows+1 : min_number_rows;
	  /* Tell each worker who its neighbors are, since they must exchange
	     data with each other. */
	  if (i == 1)
	    neighbor1 = NONE;
	  else
	    neighbor1 = i - 1;
	  if (i == numworkers)
	    neighbor2 = NONE;
	  else
	    neighbor2 = i + 1;
	  /* Now send startup information to each worker Note that this 
	     information is "tagged" as the BEGIN message*/
	  destination = i;
	  worker_number = i;
	  message_tag = BEGIN;

	  /* Send the required information to each node */
	  MPI_Send(&worker_number, 1, MPI_INT, destination, message_tag, MPI_COMM_WORLD);
	  MPI_Send(&offset, 1, MPI_INT, destination, message_tag, MPI_COMM_WORLD);
	  MPI_Send(&number_rows, 1, MPI_INT, destination, message_tag, MPI_COMM_WORLD);
	  MPI_Send(&neighbor1, 1, MPI_INT, destination, message_tag, MPI_COMM_WORLD);
	  MPI_Send(&neighbor2, 1, MPI_INT, destination, message_tag, MPI_COMM_WORLD);
	  MPI_Send(&u[0][offset][0], number_rows*NYPROB, MPI_FLOAT, destination, message_tag, MPI_COMM_WORLD);

	  /*let the world know how the problem has been divided up */
	  printf("Sent to= %d offset= %d number_rows= %d neighbor1= %d neighbor2=%d\n", destination,offset,number_rows,neighbor1,neighbor2);

	  /* increment the offset by the number_rows so the next node will
	     know where its grid begins */
	  offset += number_rows;

	} /* continue doing the above for each node, i */

      /* Now wait for results from all worker tasks */
      for (i=1; i<=numworkers; i++)
	{
	  source = i;
	  message_tag = DONE;
	  MPI_Recv(&offset, 1, MPI_INT, source, message_tag, MPI_COMM_WORLD, &status);
	  MPI_Recv(&number_rows, 1, MPI_INT, source, message_tag, MPI_COMM_WORLD, &status);
	  MPI_Recv(&u[0][offset][0], number_rows*NYPROB, MPI_FLOAT, source, message_tag, MPI_COMM_WORLD, &status);
	}
      
      /* Write final output*/
      prtdat(NXPROB, NYPROB, &u[0][0][0], "output/final.dat");
      
    } /* End of master code */
  
  if (taskid != MASTER)
    {
      /************************* worker code**********************************/
      /* Initialize everything - including the borders - to zero */
      /* iz is a flag indicating one of two grids used in the analysis */
      for (iz=0; iz<2; iz++)
	for (ix=0; ix<NXPROB; ix++)
	  for (iy=0; iy<NYPROB; iy++)
	    u[iz][ix][iy] = 0.0;
      
      /* Now receive my offset, rows, neighbors and grid partition from master
       */
      source = MASTER;
      message_tag = BEGIN;
      MPI_Recv(&worker_number, 1, MPI_INT, source, message_tag, MPI_COMM_WORLD, &status);
      MPI_Recv(&offset, 1, MPI_INT, source, message_tag, MPI_COMM_WORLD, &status);
      MPI_Recv(&number_rows, 1, MPI_INT, source, message_tag, MPI_COMM_WORLD, &status);
      MPI_Recv(&neighbor1, 1, MPI_INT, source, message_tag, MPI_COMM_WORLD, &status);
      MPI_Recv(&neighbor2, 1, MPI_INT, source, message_tag, MPI_COMM_WORLD, &status);
      MPI_Recv(&u[0][offset][0], number_rows*NYPROB, MPI_FLOAT, source, message_tag, MPI_COMM_WORLD, &status);
      
      /* Determine border elements. This takes into account that the
first and last rows have fixed temperature*/

      if (offset==0)
	start=1; /* do not include row zero */
      else
	start=offset;
      if ((offset+number_rows)==NXPROB)
	end= offset + number_rows-2;  /*do not include the last row */
      else
	end = offset + number_rows-1;
      

      /* take a look at how the work is partiioned among processors */
  printf("worker number = %d offset= %d number_rows= %d start = %d end =%d\n", worker_number, offset,number_rows,start, end);

      /* Begin doing TIME_STEPS iterations. Must communicate border rows with
	 neighbors. If I have the first or last grid row, then I only need to
	 communicate with one neighbor */
      iz = 0;
      for (it = 1; it <= TIME_STEPS; it++)
	{
	   if (neighbor1 != NONE)
	    {
	      MPI_Send(&u[iz][offset][0], NYPROB, MPI_FLOAT, neighbor1, NGHBOR2, MPI_COMM_WORLD);
	      source = neighbor1;
	      message_tag = NGHBOR1;
	      MPI_Recv(&u[iz][offset-1][0], NYPROB, MPI_FLOAT, source, message_tag, MPI_COMM_WORLD, &status);
	    }
	  if (neighbor2 != NONE)
	    {
	      MPI_Send(&u[iz][offset+number_rows-1][0], NYPROB, MPI_FLOAT, neighbor2, NGHBOR1, MPI_COMM_WORLD);
	      source = neighbor2;
	      message_tag = NGHBOR2;
	      MPI_Recv(&u[iz][offset+number_rows][0], NYPROB, MPI_FLOAT, source, message_tag, MPI_COMM_WORLD, &status);
	    } 
	  /* Now call update to update the value of grid points */
	  update(start,end,NYPROB,&u[iz][0][0],&u[1-iz][0][0]);
	  iz = 1 - iz;
	}
      /* Finally, send my portion of final results back to master */
      MPI_Send(&offset, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
      MPI_Send(&number_rows, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
      MPI_Send(&u[iz][offset][0], number_rows*NYPROB, MPI_FLOAT, MASTER, DONE, MPI_COMM_WORLD);
    }
  /*gracefully exit MPI */
  MPI_Finalize();
}


/**************************************************************************
 * subroutine update
 ****************************************************************************/
void update(int start, int end, int ny, float *u1, float *u2)
{
  int ix, iy;
  for (ix = start; ix <= end; ix++)
    for (iy = 1; iy <= ny-2; iy++)
      *(u2+ix*ny+iy) = *(u1+ix*ny+iy) +
	    diffusivity.cx * (*(u1+(ix+1)*ny+iy) +
		    *(u1+(ix-1)*ny+iy) -
		    2.0 * *(u1+ix*ny+iy)) +
	    diffusivity.cy * (*(u1+ix*ny+iy+1) +
		    *(u1+ix*ny+iy-1) -
		    2.0 * *(u1+ix*ny+iy));
}

/*****************************************************************************
 * subroutine inidat
 *****************************************************************************/
void inidat(int nx, int ny, float *u) {
  int ix, iy;
  
  for (ix = 0; ix <= nx-1; ix++)
    for (iy = 0; iy <= ny-1; iy++)
      *(u+ix*ny+iy) = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1))/(4*nx*ny);
}

/**************************************************************************
 * subroutine prtdat
 **************************************************************************/
void prtdat(int nx, int ny, float *u1, char *fnam) {
  int ix, iy;
  FILE *fp;
  
  fp = fopen(fnam, "w");
  for (iy = ny-1; iy >= 0; iy--) {
    for (ix = 0; ix <= nx-1; ix++) {
      fprintf(fp, "%d %d %6.1f\n", iy, ix, *(u1+ix*ny+iy));
    }
  }
  fclose(fp);
}
