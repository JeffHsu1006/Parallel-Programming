#include <cstdio>
#include <omp.h>
#include <mpi.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <stdexcept>
int main(int argc, char** argv){
    // assert(argc == 3);
    // char *input = argv[1];
    // char *output = argv[2];
    // FILE *fp = fopen(input, "r");
    MPI_Init(&argc,&argv);
    #if 1
    std::ifstream f(argv[1]);
    
    int V, E;

    // f.read((char*)&V, sizeof V);
    // f.read((char*)&E, sizeof E);
    MPI_File ff;
    MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &ff);
    MPI_File_read_at(ff, 0, &V, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read_at(ff, 1, &E, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_close(&ff);
    int* dist = (int*)malloc(V * V * sizeof(int));
    int i, j, k;
    // int *dist[V];
    // for (i = 0; i < V; i++) {
    //     dist[i] = (int*)malloc(V * sizeof(int));
    // }

    for (i = 0; i < V; i++) {
        for (j = 0; j < V; j++) {
            if (i == j) {
                dist[i * V + j] = 0;
            } else {
                dist[i * V + j] = 1073741823;
            }
        }
    }
    float read_file_start = omp_get_wtime();
    int *e = (int*)malloc(sizeof(int) * 3 * E);
    f.read((char*)e, sizeof(int) * 3 * E);
    for (i = 0; i < 3*E; i+=3) {
        dist[e[i] * V + e[i+1]] = e[i+2];
        // printf("%4d -> %4d, dist=%4d\n", e[0], e[1], e[2]);
    }
    f.close();
    float read_file_end = omp_get_wtime();
    printf("read file time is: %f\n", read_file_end - read_file_start);
    double start = MPI_Wtime();
    for(k=0; k<V; k++){
        #pragma omp parallel for
        for(i=0; i<V; i++){
            for(j=0; j<V; j++){
                if(dist[i*V+j] > (dist[i*V+k] + dist[k*V+j]))
                    dist[i*V+j] = dist[i*V+k] + dist[k*V+j];
            }
        }
    }
    float end_cal = omp_get_wtime();
    printf("write file time: %f\n", MPI_Wtime()-start);
    // std::ofstream ggoutoutder(argv[2], std::ios::out | std::ios::binary);
    // for (i = 0; i < V; i++) {
    //     ggoutoutder.write(reinterpret_cast<const char*>(dist[i]),V * sizeof(int));
    // }
    std::ofstream ggoutoutder(argv[2], std::ios::out | std::ios::binary);
    
    ggoutoutder.write(reinterpret_cast<const char*>(dist),V * V *sizeof(int));
    
    ggoutoutder.close();
    MPI_Finalize();
    #else
    std::ifstream f(argv[1]);
    int i;
    for (i = 0; i <100 ; i++) {
        int e;
        f.read((char*)&e, sizeof e);
        printf("%d ", e);
    }
    printf("\n");
    #endif
    return 0;
}