#include <cstdio>
#include <cstdlib>
#include <mpi.h>
#include <algorithm>
using namespace std;
#define SWAP(x,y) {float t; t = x; x = y; y = t;}

int partition(float[], int, int); 
void quickSort(float[], int, int);
void mergeSort(float[], int, int);

void merge(float[] ,int ,int ,int ,int);
double cpu_time = 0;
double com_time = 0;
double io_time = 0;
double total_time = 0;
double total_cpu_time, total_com_time, total_io_time, start_total, end_total;

int main(int argc, char** argv) {
	MPI_Init(&argc,&argv);
	start_total = MPI_Wtime();
	int rank, size, o_size, nsub_arr, n_sub_arr, n;
	double start, end;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &o_size);

	n = atoi(argv[1]);
	size=o_size;
	if((n%size) == 0){
		nsub_arr = n/size;
		n_sub_arr = n/size;
	}else{
		nsub_arr = n/size+1;
		n_sub_arr = n/size+1;
	}

	if((nsub_arr*size-n) >= nsub_arr){
		if(n/size==1)
			size=n/nsub_arr+1;
		else
			size=n/nsub_arr;
	}
		
	if(rank==(size-1)){
		if(size!=1)
			nsub_arr=n-(size-1)*nsub_arr;
	}
	if((rank+1)>size) nsub_arr=1;
	float* data = (float*) malloc (nsub_arr*sizeof(float));
	
	MPI_File f;
	start = MPI_Wtime();
	MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_RDONLY, MPI_INFO_NULL, &f);
	end = MPI_Wtime();
	io_time+=(end-start);
	if((rank+1)<=size){
		if(rank==(size-1)){
			start = MPI_Wtime();
			MPI_File_read_at(f, sizeof(float) * rank * n_sub_arr, data, nsub_arr, MPI_FLOAT, MPI_STATUS_IGNORE);
			end = MPI_Wtime();
			io_time+=(end-start);
		}else{
			start = MPI_Wtime();
			MPI_File_read_at(f, sizeof(float) * rank * nsub_arr, data, nsub_arr, MPI_FLOAT, MPI_STATUS_IGNORE);
			end = MPI_Wtime();
			io_time+=(end-start);
		}
		start = MPI_Wtime();
		sort(data, data+(nsub_arr));
		end = MPI_Wtime();
		cpu_time+=(end-start);

		float* data_trans = (float*) malloc (nsub_arr*sizeof(float));
		float* temp = (float*) malloc (nsub_arr*sizeof(float));
		float arr_biggest;
		long long i=0;
		int r=0, round=1, count_loc, count_trans, count;
		
		while(round!=(size+2)){
			
			if((((rank+r)%2)==0) && (rank!=(size-1))){
				start = MPI_Wtime();
				MPI_Send(&data[nsub_arr-1], 1, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD);
				MPI_Recv(&i, 1, MPI_INT, rank+1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				MPI_Send(&data[nsub_arr-i-1], i+1, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD);
				MPI_Recv(data_trans, i+1, MPI_FLOAT, rank+1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				end = MPI_Wtime();
				com_time+=(end-start);

				count_loc=0, count_trans=0, count=0;

				start = MPI_Wtime();
				for(int j=0;j<nsub_arr;j++){
					temp[j] = data[j];
				}

				while(count<=(nsub_arr-1)){
					if(data_trans[count_trans] < temp[count_loc]){			
						data[count] = data_trans[count_trans];
						count_trans++;
					}else{
						data[count] = temp[count_loc];
						count_loc++;
					}
					if(count_trans>i && (count+1)<nsub_arr){
						for(;count<nsub_arr;count_loc++){
							count++;
							data[count] = temp[count_loc];
						}	
						break;
					}
					count++;
				}
				end = MPI_Wtime();
				cpu_time+=(end-start);

			}else if(rank!=0){
				if((rank==(size-1)) && (((rank+r)%2)==(size%2)));
				else{
					start = MPI_Wtime();
					MPI_Recv(&arr_biggest, 1, MPI_FLOAT, rank-1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					end = MPI_Wtime();
					com_time+=(end-start);
					i = 0;
					
					int q = (n%size);
					start = MPI_Wtime();
					while(true){
						if(i==(nsub_arr-1)){
							break;
						}
						if(arr_biggest < data[i]){
							i--;
							if(i<0) i=0;
							break;
						}
						i++;
					}
					end = MPI_Wtime();
					cpu_time+=(end-start);
					
					start = MPI_Wtime();
					MPI_Send(&i, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD);
					end = MPI_Wtime();
					com_time+=(end-start);
					
					if(i>=0){
						start = MPI_Wtime();
						MPI_Recv(data_trans, i+1, MPI_FLOAT, rank-1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						MPI_Send(data, i+1, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD);
						end = MPI_Wtime();
						com_time+=(end-start);
						
						count_loc=i, count_trans=i, count=i;;
						
						start = MPI_Wtime();
						for(int j=0;j<(i+1);j++){
							temp[j] = data[j];
						}
		
						while(count>=0){	
							if(data_trans[count_trans] > temp[count_loc]){
								
								data[count] = data_trans[count_trans];
								count_trans--;
							}else{
								data[count] = temp[count_loc];
								count_loc--;
							}
							if(count_trans<0 && (count-1)>=0){
								for(;count>=0;count_loc--){
									count--;
									if(count<0) break;
									data[count] = temp[count_loc];
								}	
								break;
							}
							
							count--;
						}
						end = MPI_Wtime();
						cpu_time+=(end-start);
					}
				}
			}
			r++;
			round++;
		}
		start = MPI_Wtime();
		MPI_File_close(&f);
		MPI_File_open(MPI_COMM_SELF, argv[3], MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &f);
		if(rank==(size-1)){
			MPI_File_write_at(f, sizeof(float) * rank * n_sub_arr, data, nsub_arr,  MPI_FLOAT, MPI_STATUS_IGNORE);	
		}else{
			MPI_File_write_at(f, sizeof(float) * rank * n_sub_arr, data, nsub_arr,  MPI_FLOAT, MPI_STATUS_IGNORE);
		}
		end = MPI_Wtime();
		io_time+=(end-start);
		
	}

	MPI_File_close(&f);
	end_total = MPI_Wtime();
	total_time = end_total - start_total;
	if(rank == 0){
		printf("total_cpu_time: %f\n", cpu_time);
		printf("total_com_time: %f\n", com_time);
		printf("total_io_time: %f\n", io_time);
		printf("total_time: %f\n", total_time);
	}

	MPI_Finalize();
}
