#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <omp.h>

int main(int argc, char **argv) {
  int rank, size;
  long long numParts;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  numParts = atoll(argv[1]);

  double compute(int numberP, int start, int bl){ 
	double result = 0;
	int end = start+bl;
	double temp;
	double temp2;
  #pragma omp parallel for reduction(+:result)
	for(int i=start; i<end; i++){
		temp = (double)i/numberP;
		temp = (double)1 - pow(temp, 2);
		temp2 = (double)(1)/numberP;
		result = result + (temp2*sqrt(temp));
	}
	return result;
  }

  int block;
  int local_start;
  double local_result;
  double global_result;

  block = numParts/size;
  local_start = rank*block;
  local_result = compute(numParts, local_start, block);

  MPI_Reduce(&local_result, &global_result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  if (rank == 0){
        printf("%.12lf", global_result*4);
  }

  MPI_Finalize();
  
  return 0;
}
