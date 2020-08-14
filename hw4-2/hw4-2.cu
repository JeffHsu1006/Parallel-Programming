#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <cuda_profiler_api.h>

const int INF = ((1 << 30) - 1);
void input(char* inFileName);
void output(char* outFileName);

int ceil(int a, int b);
__global__ void phase1(int round, int n, int V, int* Dist, int B);
__global__ void phase2(int round, int n, int V, int* Dist, int B);
__global__ void phase3(int round, int n, int V, int* Dist, int B, int R, int rank);
extern __shared__ int S[];

int n, m, num_devices;
int *d_Dist_0, *d_Dist_1, *d_n;
int *Dist;


int main(int argc, char* argv[]) {

    input(argv[1]);
    int B = 32;
    int round = ceil(n, B);

    cudaGetDeviceCount(&num_devices);
    int *d_Dist[num_devices];
    
    #pragma omp parallel num_threads(2)
    {
        int thread_id = omp_get_thread_num();
        cudaSetDevice(thread_id);

        cudaMalloc((void **)&d_Dist[thread_id], n * n * sizeof(int));
        cudaMemcpy(d_Dist[thread_id], Dist, n * n * sizeof(int), cudaMemcpyHostToDevice);
        
        cudaMalloc((void **)&d_n, sizeof(int));
        cudaMemcpy(d_n, &n, sizeof(int), cudaMemcpyHostToDevice);

        dim3 grid1(1, 1);
        dim3 grid2(round, 2);
        dim3 grid3((round/2)+1, round);
        dim3 blk(B, B);

        for (int r = 0; r < round; ++r) {
            #pragma omp barrier
            if(n > B && r < (round/2) && thread_id == 1){
                cudaMemcpyPeer((void*) &d_Dist[1][r * B * n], 1, (void*) &d_Dist[0][r * B * n], 0, B * n * sizeof(int));
  
            }else if(n > B && r >= (round/2) && thread_id == 0){
                if(r == (round-1))
                    cudaMemcpyPeer((void*) &d_Dist[0][r * B * n], 0, (void*) &d_Dist[1][r * B * n], 1, (n - r * B) * n * sizeof(int));
                else
                    cudaMemcpyPeer((void*) &d_Dist[0][r * B * n], 0, (void*) &d_Dist[1][r * B * n], 1, B * n * sizeof(int));
            }
            #pragma omp barrier
            phase1<<<grid1, blk, B*B*sizeof(int)>>>(r, n, n, d_Dist[thread_id], B);
            phase2<<<grid2, blk, 2*B*B*sizeof(int)>>>(r, n, n, d_Dist[thread_id], B);
            phase3<<<grid3, blk, 2*B*B*sizeof(int)>>>(r, n, n, d_Dist[thread_id], B, round, thread_id);
        }

        if(thread_id == 0)
            cudaMemcpy(Dist, d_Dist[0], (round/2) * B * n * sizeof(int), cudaMemcpyDeviceToHost);
        else if(n > B && thread_id == 1)
            cudaMemcpy(&Dist[(round/2) * B * n], &d_Dist[1][(round/2) * B * n], (n - (round/2) * B) * n * sizeof(int), cudaMemcpyDeviceToHost);

    }

    output(argv[2]);

    return 0;
}

void input(char* infile) {
    FILE* file = fopen(infile, "rb");
    fread(&n, sizeof(int), 1, file);
    fread(&m, sizeof(int), 1, file);
    Dist = (int*)malloc(n*n*sizeof(int));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                Dist[i * n + j] = 0;
            } else {
                Dist[i * n + j] = INF;
            }
        }
    }

    int pair[3];
    for (int i = 0; i < m; ++i) {
        fread(pair, sizeof(int), 3, file);
        Dist[pair[0] * n + pair[1]] = pair[2];
    }
    fclose(file);
}

void output(char* outFileName) {
    FILE* outfile = fopen(outFileName, "w");
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (Dist[i * n + j] >= INF) Dist[i * n + j] = INF;
        }
    }
    fwrite(Dist, sizeof(int), n*n, outfile);
    fclose(outfile);
}

int ceil(int a, int b) { return (a + b - 1) / b; }

__global__ void phase1(int round, int n, int V, int* Dist, int B){
    int s_i = threadIdx.y;
    int s_j = threadIdx.x;
    int i = round * B + s_i;
    int j = round * B + s_j;
    
    if((i < n && j < n))
        S[s_i * B + s_j] = Dist[i * V + j];
    __syncthreads();

    int tt = round * B;
    int ss = s_i * B;
    #pragma unroll
    for (int k = 0; k < B && tt + k < n; ++k) {
        if (S[ss + k] + S[k * B + s_j] < S[ss + s_j])
            S[ss + s_j] = S[ss + k] + S[k * B + s_j];
        
        __syncthreads();
    }
    if (i < n && j < n) Dist[i * V + j] = S[ss + s_j];
    __syncthreads();

}

__global__ void phase2(int round, int n, int V, int* Dist, int B){

    if (blockIdx.x == round) return;

    int* S_pivot = &S[0];
    int* S_dist = &S[B * B];

    int s_i = threadIdx.y;
    int s_j = threadIdx.x;
    int i = round * B + s_i;
    int j = round * B + s_j;
    
    int ss = s_i * B;

    if((i < n && j < n))
        S_pivot[ss + s_j] = Dist[i * V + j];
    __syncthreads();

    if (blockIdx.y == 0)
        j = blockIdx.x * B + s_j;
    else
        i = blockIdx.x * B + s_i;

    if (i >= n || j >= n) return;

    if((i < n && j < n))
        S_dist[ss + s_j] = Dist[i * V + j];
    __syncthreads();

    int tt = round * B;
    if(blockIdx.y == 1){
        #pragma unroll
        for (int k = 0; k < B && tt + k < n; ++k) {
            if (S_dist[ss + k] + S_pivot[k * B + s_j] < S_dist[ss + s_j])
                S_dist[ss + s_j] = S_dist[ss + k] + S_pivot[k * B + s_j];
        }
    }else{
        #pragma unroll
        for (int k = 0; k < B && tt + k < n; ++k) {
            if (S_pivot[ss + k] + S_dist[k * B + s_j] < S_dist[ss + s_j])
                S_dist[ss + s_j] = S_pivot[ss + k] + S_dist[k * B + s_j];
        }
    }
    
    if (i < n && j < n) Dist[i * V + j] = S_dist[ss + s_j];
    __syncthreads();
}

__global__ void phase3(int round, int n, int V, int* Dist, int B, int R, int rank){

    int block_i = blockIdx.x;
    int block_j = blockIdx.y;

    if(rank == 1)
        block_i += (R/2);

    if (block_i == round || block_j == round) return;

    int* S_pivot_row = &S[0];
    int* S_pivot_col= &S[B * B];

    int s_i = threadIdx.y;
    int s_j = threadIdx.x;
    int i = block_i * B + s_i;
    int j = block_j * B + s_j;
    int b_i = round * B + s_i;
    int b_j = round * B + s_j;

    int ss = s_i * B;
    
    if(i < n && b_j < n) S_pivot_row[ss + s_j] = Dist[i * V + b_j];
    if(j < n && b_i < n) S_pivot_col[ss + s_j] = Dist[b_i * V + j];
    __syncthreads();

    if (i >= n || j >= n) return;

    int dst = Dist[i * V + j];

    int tt = round * B;
    #pragma unroll
    for (int k = 0; k < B && tt + k < n; ++k) {
        int tmp_result = S_pivot_row[ss + k] + S_pivot_col[k * B + s_j];
        if (tmp_result < dst) dst = tmp_result;
    }
    
    Dist[i * V + j] = dst;
}
