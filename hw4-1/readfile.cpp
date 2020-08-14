#include <stdio.h>
#include <stdlib.h>

void input(char* inFileName);
const int INF = ((1 << 30) - 1);
const int V = 40000;
int Dist[V][V];
int n, m;

int main(int argc, char* argv[]) {
    FILE* file = fopen(argv[1], "rb");
    fread(&n, sizeof(int), 1, file);
    fread(&m, sizeof(int), 1, file);
    printf("vector: %d\n", n);
    printf("edge: %d\n", m);
    // input(argv[1]);

    // for (int i = 0; i < n; ++i) {
    //     for (int j = 0; j < n; ++j) 
    //         printf("%d\n", Dist[i][j]);
    // }
    return 0;
}


void input(char* infile) {
    FILE* file = fopen(infile, "rb");
    fread(&n, sizeof(int), 1, file);
    fread(&m, sizeof(int), 1, file);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                Dist[i][j] = 0;
            } else {
                Dist[i][j] = INF;
            }
        }
    }

    int pair[3];
    for (int i = 0; i < m; ++i) {
        fread(pair, sizeof(int), 3, file);
        Dist[pair[0]][pair[1]] = pair[2];
    }
    fclose(file);
}