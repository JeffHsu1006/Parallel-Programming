#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#define PNG_NO_SETJMP
#include <sched.h>
#include <assert.h>
#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <mpi.h>
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

int num_parts, last_id, last_part, iters, width, height;
double upper, lower, right ,left;
int* image;
int current_height = 0;
int num_threads;
double computing_time = 0;

void *cal_man_set (void *input);
static double getDoubleTime();

void write_png(const char* filename, int iters, int width, int height, const int* buffer) {
    FILE* fp = fopen(filename, "wb");
    assert(fp);
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    assert(png_ptr);
    png_infop info_ptr = png_create_info_struct(png_ptr);
    assert(info_ptr);
    png_init_io(png_ptr, fp);
    png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_set_filter(png_ptr, 0, PNG_NO_FILTERS);
    png_write_info(png_ptr, info_ptr);
    png_set_compression_level(png_ptr, 1);
    size_t row_size = 3 * width * sizeof(png_byte);
    png_bytep row = (png_bytep)malloc(row_size);
    for (int y = 0; y < height; ++y) {
        memset(row, 0, row_size);
        for (int x = 0; x < width; ++x) {
            int p = buffer[(height - 1 - y) * width + x];
            png_bytep color = row + x * 3;
            if (p != iters) {
                if (p & 16) {
                    color[0] = 240;
                    color[1] = color[2] = p % 16 * 16;
                } else {
                    color[0] = p % 16 * 16;
                }
            }
        }
        png_write_row(png_ptr, row);
    }
    free(row);
    png_write_end(png_ptr, NULL);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
}

int main(int argc, char** argv) {
    int rank, size, local_start;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* detect how many CPUs are available */
    cpu_set_t cpu_set;
    sched_getaffinity(0, sizeof(cpu_set), &cpu_set);
    num_threads = CPU_COUNT(&cpu_set);

    /* argument parsing */
    assert(argc == 9);
    const char* filename = argv[1];
    iters = strtol(argv[2], 0, 10);
    left = strtod(argv[3], 0);
    right = strtod(argv[4], 0);
    lower = strtod(argv[5], 0);
    upper = strtod(argv[6], 0);
    width = strtol(argv[7], 0, 10);
    height = strtol(argv[8], 0, 10);

    /* allocate memory for image */
    image = (int*)calloc(width * height,  sizeof(int));
    assert(image);

    /* creat thread */
    pthread_t threads[num_threads];

    int t;
    double start = MPI_Wtime();
    for(t = 0; t < num_threads; t++){
        pthread_create(&threads[t], NULL, cal_man_set, (void *)&t);
    }
    
    
    for(int i=0; i<num_threads; i++)
        pthread_join(threads[i], NULL);
    computing_time += MPI_Wtime() - start;
    printf("computing time: %f\n", computing_time);

    /* draw and cleanup */
    write_png(filename, iters, width, height, image);
    free(image);
    MPI_Finalize();
}


void* cal_man_set(void* input){
    double y0;
    pthread_mutex_lock(&mutex);
    int j = current_height;
    current_height+=1;
    pthread_mutex_unlock(&mutex);

    /* mandelbrot set */
    while(current_height <= height){
        
        
        y0 = j * ((upper - lower) / height) + lower;
        for (int i = 0; i < width; ++i) {
            double x0 = i * ((right - left) / width) + left;

            int repeats = 0;
            double x = 0;
            double y = 0;
            double length_squared = 0;
            while (repeats < iters && length_squared < 4) {
                double temp = x * x - y * y + x0;
                y = 2 * x * y + y0;
                x = temp;
                length_squared = x * x + y * y;
                ++repeats;
            }

            image[j * width + i] = repeats;

        }

        pthread_mutex_lock(&mutex);
        j = current_height;
        current_height+=1;
        pthread_mutex_unlock(&mutex);
    }

    pthread_exit(NULL);
}
