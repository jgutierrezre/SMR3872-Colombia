/* connect mpi processes to a ring and send a message around */

#include <stdio.h>
#include <mpi.h>

int main(int argc, char **argv)
{
  /* Variables:
   ! me   : rank of current process
   ! prev : rank of process that we are receiving data from
   ! next : rank of process that we are sending data to
   ! mesg : storage for message data
   ! tmp  : holding buffer for incoming data
   */
  int me, ncpu, prev, next, mesg, i, tmp;

  /* EDIT: declare a variable to store the non-blocking MPI request handle. */
  MPI_Request req;

  /* set up MPI environment */
  i = MPI_Init(&argc,&argv);
  if (i != MPI_SUCCESS) {
    puts("problem initializing MPI");
    return 1;
  }

  MPI_Comm_size(MPI_COMM_WORLD,&ncpu);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);

  /* EDIT: determne previous and next rank to the current rank */
  prev = me - 1;
  if (prev < 0) prev = ncpu - 1;

  next = me + 1;
  if (next == ncpu) next = 0;
  
  /* EDIT: at the beginning assign the message (a number != 0) to rank 0
   * all other ranks have nothing (= 0) */
  if (me == 0) {
    mesg = 1;
  } else {
    mesg = 0;
  }

  /* run loop to communicate the message buffer from one rank to the next.
   to avoid a deadlock we use a non-blocking receive with a blocking send. */

  for (i=0; i < 20; ++i) {
    /* EDIT: first post the non-blocking receive from the previous rank.
     *       then use a blocking send to the next rank.
     *       finally wait until the receive is completed. */
    MPI_Irecv(&tmp,1,MPI_INT,prev,0,MPI_COMM_WORLD,&req);
    MPI_Send(&mesg,1,MPI_INT,next,0,MPI_COMM_WORLD);
    MPI_Wait(&req,MPI_STATUS_IGNORE);

    /* communication complete for this step. copy incoming buffer to message */
    mesg = tmp;
    if (mesg != 0) printf("step %d: rank %d has message\n",i,me);
  }

  /* shut down MPI environment */
  MPI_Finalize();
  return 0;
}
