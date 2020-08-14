#include <sched.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>

int** e;
int* b;

int main(int argc, char** argv) {
    MPI_Init(&argc,&argv);
    /* detect how many CPUs are available */
    cpu_set_t cpu_set;
    sched_getaffinity(0, sizeof(cpu_set), &cpu_set);
    int num_thread = CPU_COUNT(&cpu_set);
    if(num_thread>8) num_thread=8;
    
    MPI_File f;
	
    int* v_e = (int*)malloc(2*sizeof(int));
    MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &f);
    MPI_File_read_at(f, 0, v_e, 2, MPI_INT, MPI_STATUS_IGNORE);

    int v_n = v_e[0];
    int e_n = v_e[1];

    e = (int**)malloc(v_n*sizeof(int*));
    for(int i=0; i<v_n; i++)
        e[i] = (int*)malloc(v_n*sizeof(int));

    b = (int*)malloc((e_n*3)*sizeof(int));
    MPI_File_read_at(f, sizeof(int)*2, b, e_n*3, MPI_INT, MPI_STATUS_IGNORE);


    // initialize the 2-D array to store the edge weight
    for(int i=0; i<v_n ;i++){
        for(int j=0; j<v_n ;j++)
            e[i][j] = 1073741823;
        e[i][i] = 0;
    }

    // give the original edge weight 
    int edge = e_n*3;
    for(int i=0; i<edge; i+=3){
        int src = b[i];
        int dst = b[i+1];
        e[src][dst] = b[i+2];
    }

    // Floyd Warshall
    for(int k=0; k<v_n; k++){
        #pragma omp parallel for
        for(int i=0; i<v_n; i++){
            for(int j=0; j<v_n; j++){
                if(e[i][j] > (e[i][k] + e[k][j]))
                    e[i][j] = e[i][k] + e[k][j];
            }
        }
    }

    MPI_File_close(&f);

    MPI_File_open(MPI_COMM_SELF, argv[2], MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &f);
    #pragma omp parallel for num_threads(num_thread)
    for(int i=0; i<v_n ;i++){
        MPI_File_write_at(f, i*v_n*sizeof(int), &e[i][0], v_n*sizeof(int), MPI_BYTE, MPI_STATUS_IGNORE);
    }

    MPI_File_close(&f);
	MPI_Finalize();
}