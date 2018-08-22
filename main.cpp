#include <cstdio>
#include <omp.h>
#include <cstdlib>

using namespace std;

#define MAX 10000
#define NOT_CONNECTED -1
int BLOCK_SIZE=29;

int distance[MAX][MAX];

int nodesCount;

void Initialize() {
    printf("nodes count: %d",nodesCount);
    for (int i = 0; i < MAX; ++i) {
        for (int j = 0; j < MAX; ++j) {
            distance[i][j] = NOT_CONNECTED;

        }
        distance[i][i] = 0;
    }
}




void FW_parallel(int i, int j, int BlockSize, int k) {
    if ((i>nodesCount)||(j>nodesCount))
        return;
    // printf("FW for i,j: %d,%d blocksize: %d k:%d\n",i,j,BlockSize,k);

    int step = BLOCK_SIZE/omp_get_num_procs();

    for (int k1 = k; k1 <= k + BlockSize; ++k1) {
      //    printf("step: %d",step);
        //  printf ("now computing: %d %d block:%d, k:%d",i,j,BlockSize,k);
//#pragma omp parallel for schedule(static,step)
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

void FW(int i, int j, int BlockSize, int k) {
    if ((i>nodesCount)||(j>nodesCount))
        return;
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

void parallel_diameter(){
    int stat_step = nodesCount/(BLOCK_SIZE*omp_get_num_procs());
    printf("stat step: %d\n",stat_step);

    double timeBegin, timeEnd;
    timeBegin = omp_get_wtime();

    for (int k = 0; k < MAX / BLOCK_SIZE; k++) {
        if (k*BLOCK_SIZE>nodesCount){
            break;
        }
        //     printf("big loop\n");
//#pragma omp parallel
//        {

//#pragma omp single


        FW(BLOCK_SIZE * k, BLOCK_SIZE * k, BLOCK_SIZE, k * BLOCK_SIZE);
#pragma omp parallel for schedule(dynamic,stat_step)
        for (int i = 0; i < BLOCK_SIZE * k; i = i + BLOCK_SIZE) {
            //   printf("here1");

//            if ((i>nodesCount)||(k*BLOCK_SIZE>nodesCount)){
//                continue;
//            }

            FW(i, k * BLOCK_SIZE, BLOCK_SIZE, k * BLOCK_SIZE);


        }
#pragma omp parallel for schedule(dynamic,stat_step)
        for (int i = BLOCK_SIZE * (k + 1); i < MAX; i = i + BLOCK_SIZE) {


            //        printf("here2");
            if ((i>nodesCount)||(k*BLOCK_SIZE>nodesCount)){
                continue;
            }
            FW(i, k * BLOCK_SIZE, BLOCK_SIZE, k * BLOCK_SIZE);


        }
#pragma omp parallel for schedule(dynamic,stat_step)
        for (int j = 0; j < BLOCK_SIZE * k; j = j + BLOCK_SIZE) {
            //         printf("here3");


            FW(k * BLOCK_SIZE, j, BLOCK_SIZE, k * BLOCK_SIZE);


        }
#pragma omp parallel for schedule(dynamic,stat_step)
        for (int j = BLOCK_SIZE * (k + 1); j < MAX; j = j + BLOCK_SIZE) {
            //       printf("here4");
            if ((j>nodesCount)||(k*BLOCK_SIZE>nodesCount)){
                continue;
            }
            FW(k * BLOCK_SIZE, j, BLOCK_SIZE, k * BLOCK_SIZE);


        }
//#pragma omp barrier
#pragma omp parallel for schedule(dynamic,stat_step)
        for (int i = 0; i < MAX; i = i + BLOCK_SIZE) {
            for (int j = 0; j < MAX; j = j + BLOCK_SIZE) {
                if ((i != k * BLOCK_SIZE) && (j != k * BLOCK_SIZE)) {
                    //             printf("here5");
                    if ((i>nodesCount)||(j>nodesCount)){
                        continue;
                    }

                    FW(i, j, BLOCK_SIZE, k * BLOCK_SIZE);


                }
            }
        }

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


    timeEnd = omp_get_wtime();
    printf("%d %.16g\n", diameter,(timeEnd-timeBegin));
}

int optimizeBlockSize(){
    double minimumTime=30000;
    int optimumBlockSize=-1;
    nodesCount=2000;
    for (int i=0;i<2000;i++){
        for (int j=0;j<2000;j++){
            distance[i][j]= rand()* 100;
        }
    }
    for (int i=50;i>15;i--){
        printf("now checking tile size: %d\n",i);
        BLOCK_SIZE=i;
        double avgTime=0 ;
        for (int j=0;j<6;j++){
            double start = omp_get_wtime();
            parallel_diameter();
            double end = omp_get_wtime();
            avgTime+=(end-start);
        }

        if (avgTime<minimumTime){
            minimumTime=avgTime ;
            optimumBlockSize=i;
            printf("new optimum block size: %d\n",optimumBlockSize);
        }

    }
    printf("optimization done -> blocksize:%d\n",optimumBlockSize);


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

//    double start = omp_get_wtime();
    fscanf(in_file, "%d", &nodesCount);

    Initialize();
    optimizeBlockSize();

    return 0;


    printf("nodes count read: %d", nodesCount);
    int a, b, c;
  //  int count=0;
    while (fscanf(in_file, "%d %d %d", &a, &b, &c) != EOF) {
        if (a > nodesCount || b > nodesCount) {
            printf("Vertex index out of boundary.");
            return -1;
        }
        distance[a][b] = c;
       // count++;
       /// printf("read %d\n",count);
    }

       printf("done reading");

        parallel_diameter();



    return 0;

}

