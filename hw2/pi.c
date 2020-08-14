#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

void *cal_pi (void *a);

struct args {
  double num_part;
  double id;
  double pi;
};


int main(int argc, char *argv[]){
  long long numParts = atoll(argv[2]);
  int num_threads = atoi(argv[1]);
  pthread_t threads[num_threads];
  struct args *a = (struct args *)malloc(sizeof (struct args));
  double num_part = numParts/num_threads;
  a->num_part = num_part;

  double ID[num_threads];
  

  for(t = 0; t < num_threads; t++){
     a->id = ID[t];
     pthread_create(&threads[t], NULL, cal_pi, (void *)a);
  }

  for(int i=0; i<num_threads; i++)
     pthread_join(threads[i], NULL);

    

  printf("%.12lf\n", a->pi * 4 / numParts);

  return 0;
}

void* cal_pi(void* input){
   double pii = 0;
   double num_parts = ((struct args*)input)->num_part;
   double tid = ((struct args*)input)->id;
   int start=tid*num_parts;
   int end=tid*num_parts+num_parts;
   for (int i = start; i < end; i++) {
       pii += sqrt(1 - ((double)i / num_parts) * ((double)i / num_parts));
   }
   pthread_mutex_lock(&mutex);
   ((struct args*)input)->pi += pii;
   pthread_mutex_unlock(&mutex);
   
   pthread_exit(NULL);
}
