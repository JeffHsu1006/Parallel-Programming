#include <cstdlib>
#include <cassert>
#include <cmath>
#include <chrono>
#include <iostream>
#include "cuda_fp16.h"
#include "utils/utils.h"

#define Z 2
#define Y 5
#define X 5
#define xBound X / 2
#define yBound Y / 2
#define SCALE 8
#define num_blocks 2048
#define num_threads 256

// const int tx = 128;
// const int ty = 5;

__constant__ int filter[Z][Y][X] = { { { -1, -4, -6, -4, -1 },
                                        { -2, -8, -12, -8, -2 },
                                        { 0, 0, 0, 0, 0 },
                                        { 2, 8, 12, 8, 2 },
                                        { 1, 4, 6, 4, 1 } },
                                      { { -1, -2, 0, 2, 1 },
                                        { -4, -8, 0, 8, 4 },
                                        { -6, -12, 0, 12, 6 },
                                        { -4, -8, 0, 8, 4 },
                                        { -1, -2, 0, 2, 1 } } };


inline __device__ int bound_check(int val, int lower, int upper) {
    if (val >= lower && val < upper)
        return 1;
    else
        return 0;
}


__global__ void sobel(unsigned char *s, unsigned char *t, unsigned height,
                      unsigned width, unsigned channels) {

    // int tid = blockIdx.x * blockDim.x + threadIdx.x;
    float val[Z][3];
    // if (tid >= height) return;

     const int size = 48*1024;
     __shared__ unsigned char ss[size];
    // __shard__ unsigned char sp[ty + 4][(tx + 4) * channels];
    
     int end = blockIdx.x*size+size;
     int all_size = height*width*channels;
     if(end > all_size) end = all_size;
     for(int i=blockIdx.x*(num_threads); i<end; i++)
        ss[i] = s[i];  
    
     __syncthreads();
    // int y = tid; 
    for(int y = blockIdx.x; y < height; y+=gridDim.x) {
    // for (int x = 0; x < width; ++x) {
    for(int x = threadIdx.x; x < width; x+=blockDim.x){

	// for (int sy = threadIdx.y; sy < ty + 4; sy += blockDim.y)
	// int gy; if (bound_check)
	// for (int sx = threadIdx.x; sx < tx + 4; sx += blockDim.x)
	// int gx; if (bound
	// sp[sy][sx * channels + 0] = ss[...]
	// sp[sy][sx * channels + 1] = ss[...]
	// sp[sy][sx * channels + 2] = ss[...]
	


        /* Z axis of filter */
        for (int i = 0; i < Z; ++i) {

            val[i][2] = 0.;
            val[i][1] = 0.;
            val[i][0] = 0.;

            /* Y and X axis of filter */
            for (int v = -yBound; v <= yBound; ++v) {
                for (int u = -xBound; u <= xBound; ++u) {
                    if (bound_check(x + u, 0, width) &&
                        bound_check(y + v, 0, height)) {
                        const unsigned char R =
                            s[channels * (width * (y + v) + (x + u)) + 2];
                        const unsigned char G =
                            s[channels * (width * (y + v) + (x + u)) + 1];
                        const unsigned char B =
                            s[channels * (width * (y + v) + (x + u)) + 0];
                        val[i][2] += R * filter[i][u + xBound][v + yBound];
                        val[i][1] += G * filter[i][u + xBound][v + yBound];
                        val[i][0] += B * filter[i][u + xBound][v + yBound];
                    }
                }
            }
        }
        float totalR = 0.;
        float totalG = 0.;
        float totalB = 0.;
        for (int i = 0; i < Z; ++i) {
            totalR += val[i][2] * val[i][2];
            totalG += val[i][1] * val[i][1];
            totalB += val[i][0] * val[i][0];
        }
        totalR = sqrt(totalR) / SCALE;
        totalG = sqrt(totalG) / SCALE;
        totalB = sqrt(totalB) / SCALE;
        const unsigned char cR = (totalR > 255.) ? 255 : totalR;
        const unsigned char cG = (totalG > 255.) ? 255 : totalG;
        const unsigned char cB = (totalB > 255.) ? 255 : totalB;
        t[channels * (width * y + x) + 2] = cR;
        t[channels * (width * y + x) + 1] = cG;
        t[channels * (width * y + x) + 0] = cB;
    }
    }
}

int main(int argc, char **argv) {
    assert(argc == 3);
    unsigned height, width, channels;
    unsigned char *src = NULL, *dst;
    unsigned char *dsrc, *ddst;

    /* read the image to src, and get height, width, channels */
    if (read_png(argv[1], &src, &height, &width, &channels)) {
        std::cerr << "Error in read png" << std::endl;
        return -1;
    }
    
    auto Tstart = std::chrono::high_resolution_clock::now();

    dst = (unsigned char *)malloc(height * width * channels *
                                  sizeof(unsigned char));
    cudaHostRegister(src, height * width * channels * sizeof(unsigned char), cudaHostRegisterDefault);

    // cudaMalloc(...) for device src and device dst
    cudaMalloc(&dsrc, height * width * channels * sizeof(unsigned char));
    cudaMalloc(&ddst, height * width * channels * sizeof(unsigned char));

    // cudaMemcpy(...) copy source image to device (filter matrix if necessary)
    cudaMemcpy(dsrc, src, height * width * channels * sizeof(unsigned char),
               cudaMemcpyHostToDevice);

    // decide to use how many blocks and threads
    // const int num_threads = 256;
    // const int num_blocks = 2048;

    // launch cuda kernel
    sobel << <num_blocks, num_threads>>> (dsrc, ddst, height, width, channels);

    // cudaMemcpy(...) copy result image to host
    cudaMemcpy(dst, ddst, height * width * channels * sizeof(unsigned char),
               cudaMemcpyDeviceToHost);

    auto Tend = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = Tend - Tstart;
    std::cout << "Measured Compute Time: " << duration.count() << "s" << std::endl;

    write_png(argv[2], dst, height, width, channels);
    free(src);
    free(dst);
    cudaFree(dsrc);
    cudaFree(ddst);
    return 0;
}
