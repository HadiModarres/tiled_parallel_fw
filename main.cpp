#include <cstdio>
#include <omp.h>

using namespace std;

#define MAX 2048
#define NOT_CONNECTED -1
#define BLOCK_SIZE 256

int distance[MAX][MAX];

int nodesCount;

void Initialize() {
    for (int i = 0; i < MAX; ++i) {
        for (int j = 0; j < MAX; ++j) {
            distance[i][j] = NOT_CONNECTED;

        }
        distance[i][i] = 0;
    }
}


void FW(int i, int j, int BlockSize, int k) {

    // printf("FW for i,j: %d,%d blocksize: %d k:%d\n",i,j,BlockSize,k);

    for (int k1 = k; k1 <= k + BlockSize; ++k1) {
        //  printf("k1: %d",k1);
        //  printf ("now computing: %d %d block:%d, k:%d",i,j,BlockSize,k);
        for (int i1 = i; i1 <= i + BlockSize - 1; ++i1) {
            //  printf("i1: %d",i1);
            if (distance[i1][k1] != NOT_CONNECTED) {
                for (int j1 = j; j1 <= j + BlockSize - 1; ++j1) {
                    //        printf("j1: %d",j1);

                    if (distance[k1][j1] != NOT_CONNECTED &&
                        (distance[i1][j1] == NOT_CONNECTED || distance[i1][k1] + distance[k1][j1] < distance[i1][j1])) {
                        distance[i1][j1] = distance[i1][k1] + distance[k1][j1];
                    }
                }
            }
        }

    }
    //   printf ("computing finished: %d %d block:%d, k:%d\n",i,j,BlockSize,k);
    //   printf("fw done\n");

}


int main(int argc, char **argv) {

    int i = 2;
    for (; i <= 5; ++i) {
        printf("i");
    }

    if (argc != 2) {
        printf("The path to the input file is not specified as a parameter.\n");
        return -1;
    }
    FILE *in_file = fopen(argv[1], "r");
    if (in_file == NULL) {
        printf("Can't open file for reading.\n");
        return -1;
    }

    double start = omp_get_wtime();
    Initialize();

    fscanf(in_file, "%d", &nodesCount);

    printf("nodes count read: %d", nodesCount);
    int a, b, c;
    while (fscanf(in_file, "%d %d %d", &a, &b, &c) != EOF) {
        if (a > nodesCount || b > nodesCount) {
            printf("Vertex index out of boundary.");
            return -1;
        }
        distance[a][b] = c;
    }

    //   printf("done reading");
    for (int k = 0; k < MAX / BLOCK_SIZE; k++) {
        printf("big loop\n");
//#pragma omp parallel
//        {

//#pragma omp single


        FW(BLOCK_SIZE * k, BLOCK_SIZE * k, BLOCK_SIZE, k * BLOCK_SIZE);
#pragma omp parallel for schedule(static,4)
        for (int i = 0; i < BLOCK_SIZE * k; i = i + BLOCK_SIZE) {
            //   printf("here1");


            FW(i, k * BLOCK_SIZE, BLOCK_SIZE, k * BLOCK_SIZE);


        }
#pragma omp parallel for schedule(static,4)

        for (int i = BLOCK_SIZE * (k + 1); i < MAX; i = i + BLOCK_SIZE) {


            //        printf("here2");

            FW(i, k * BLOCK_SIZE, BLOCK_SIZE, k * BLOCK_SIZE);


        }
#pragma omp parallel for schedule(guided)

        for (int j = 0; j < BLOCK_SIZE * k; j = j + BLOCK_SIZE) {
            //         printf("here3");


            FW(k * BLOCK_SIZE, j, BLOCK_SIZE, k * BLOCK_SIZE);


        }
#pragma omp parallel for schedule(static,4)
        for (int j = BLOCK_SIZE * (k + 1); j < MAX; j = j + BLOCK_SIZE) {
            //       printf("here4");

            FW(k * BLOCK_SIZE, j, BLOCK_SIZE, k * BLOCK_SIZE);


        }
#pragma omp barrier
#pragma omp parallel for collapse(2) schedule(static,4)
        for (int i = 0; i < MAX; i = i + BLOCK_SIZE) {
            for (int j = 0; j < MAX; j = j + BLOCK_SIZE) {
                if ((i != k * BLOCK_SIZE) && (j != k * BLOCK_SIZE)) {
                    //             printf("here5");


                    FW(i, j, BLOCK_SIZE, k * BLOCK_SIZE);


                }
            }
        }
#pragma omp barrier
//    }

    }
    //Floyd-Warshall
//    for (int k=1;k<=nodesCount;++k){
//        for (int i=1;i<=nodesCount;++i){
//            if (distance[i][k]!=NOT_CONNECTED){
//                for (int j=1;j<=nodesCount;++j){
//                    if (distance[k][j]!=NOT_CONNECTED && (distance[i][j]==NOT_CONNECTED || distance[i][k]+distance[k][j]<distance[i][j])){
//                        distance[i][j]=distance[i][k]+distance[k][j];
//                    }
//                }
//            }
//        }
//    }

    int diameter = -1;

    //look for the most distant pair
    for (int i = 1; i <= nodesCount; ++i) {
        for (int j = 1; j <= nodesCount; ++j) {
            if (diameter < distance[i][j]) {
                diameter = distance[i][j];
                //      printf("%d-%d-%d\n", i, diameter, j);
            }
        }
    }

    printf("%d\n", diameter);

    double end = omp_get_wtime();
    printf("total running time: %f", (end - start));
    return 0;

}

