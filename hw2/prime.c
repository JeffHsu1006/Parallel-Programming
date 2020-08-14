#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <omp.h>
#define CHUNKSIZE 1000


int isPrime(int i) {
    int sqrt_i = (int)sqrt((double)i);
    int j;
    for (j = 2; j <= sqrt_i; ++ j) {
        if (i % j == 0) return 0;
    }
    return 1;
}
int main(int argc, char** argv) {
    assert(argc == 2);
    int N = atoi(argv[1]);
    int i;
    int count = 0;
    int num_threads = omp_get_num_threads();
    int thread_id;

    printf("%d",num_threads);
    # pragma omp parallel shared(count)
    {   
	int local_count = 0;
	# pragma omp for schedule(dynamic, CHUNKSIZE) nowait
    	for (i = 2; i <= N; i ++) {
          local_count += isPrime(i);
        }

	count += local_count;
    }
    printf("There are %d prime numbers <= %d\n", count, N);
    return 0;
}
