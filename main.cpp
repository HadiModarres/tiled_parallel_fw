#include <cstdio>
#include <omp.h>
#include <malloc.h>

using namespace std;

#define MAX 32
#define NOT_CONNECTED -1
#define BLOCK_SIZE 16

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

typedef struct
{
//    int dist[BLOCK_SIZE][BLOCK_SIZE];
    int *dist = (int *)malloc(BLOCK_SIZE * BLOCK_SIZE * sizeof(int));


}Tile;


int tile_new(Tile *out,int source[][MAX], int i,int j,int tile_size){
    for (int i1=0;i1<tile_size;i1++){
        for (int j1=0;j1<tile_size;j1++){
            out->dist[i1+j1*BLOCK_SIZE]= source[i1+i][j1+j];

        }
    }
}

void makeTiles(Tile tiles[MAX/BLOCK_SIZE][MAX/BLOCK_SIZE],int source[MAX][MAX],int sourceSize){
    for (int i=0;i<sourceSize;i+=BLOCK_SIZE){
        for (int j=0;j<sourceSize;j+=BLOCK_SIZE){
            printf("new tile i:%d j:%d",i,j);

            tile_new(&tiles[i][j],source,i,j,MAX);
        }
    }
}

void printTile(Tile tile){
    printf("print tile start\n");
    for (int i=0;i<BLOCK_SIZE;i++){
        for (int j=0;j<BLOCK_SIZE;j++){
            printf("%d ",tile.dist[i+j*BLOCK_SIZE]);
        }
        printf("\n");
    }
    printf("print tile end\n");
}


//void FW_tile(Tile *tile,int k){
//
//    for (int k1 = k; k1 <= k + BLOCK_SIZE; ++k1) {
//        //  printf("k1: %d",k1);
//        //  printf ("now computing: %d %d block:%d, k:%d",i,j,BlockSize,k);
//        for (int i1 = 0; i1 <= BLOCK_SIZE - 1; ++i1) {
//            //  printf("i1: %d",i1);
//            if (distance[i1][k1] != NOT_CONNECTED) {
//                for (int j1 = 0; j1 <= BLOCK_SIZE - 1; ++j1) {
//                    //        printf("j1: %d",j1);
//
//                    if (distance[k1][j1] != NOT_CONNECTED &&
//                        (distance[i1][j1] == NOT_CONNECTED || distance[i1][k1] + distance[k1][j1] < distance[i1][j1])) {
//                        tile->dist[i1+j1*BLOCK_SIZE] = distance[i1][k1] + distance[k1][j1];
//                    }
//                }
//            }
//        }
//
//    }
//
//}
void FW_tile(Tile *tile,int i, int j, int BlockSize, int k) {

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
//                        tile->dist[i1+j1*BlockSize] = distance[i1][k1] + distance[k1][j1];
                        tile->dist[(i1-i)+(j1-j)*BlockSize] = distance[i1][k1] + distance[k1][j1];
                    }
                }
            }
        }

    }
    //   printf ("computing finished: %d %d block:%d, k:%d\n",i,j,BlockSize,k);
    //   printf("fw done\n");

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


void printDistanceMatrix(int distance[MAX][MAX]){
    printf("print distance start\n");
    for (int i=0;i<MAX;i++){
        for (int j=0;j<MAX;j++){
            printf("%d ",distance[i][j]);
        }
        printf("\n");
    }
    printf("print distance end\n");
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

    printf("nodes count read: %d\n", nodesCount);
    int a, b, c;
    while (fscanf(in_file, "%d %d %d", &a, &b, &c) != EOF) {
        if (a > nodesCount || b > nodesCount) {
            printf("Vertex index out of boundary.");
            return -1;
        }
        distance[a][b] = c;
    }

    Tile tiles[MAX/BLOCK_SIZE][MAX/BLOCK_SIZE];
    for (int i=0;i<MAX/BLOCK_SIZE;i++){
        for (int j=0;j<MAX/BLOCK_SIZE;j++){
            for (int a=0;a<BLOCK_SIZE;a++){
                for (int b=0;b<BLOCK_SIZE;b++){
                    tiles[i][j].dist[a+b*BLOCK_SIZE] = distance[a+i*BLOCK_SIZE][b+j*BLOCK_SIZE];
                }
            }
        }
    }

//    for (int i=0;i<MAX/BLOCK_SIZE;i++) {
//        for (int j = 0; j < MAX / BLOCK_SIZE; j++) {
//            printTile(tiles[i][j]);
//        }
//    }
//
//    printDistanceMatrix(distance);

//            printf("i1:%d j1:%d\n",distance[2][17],tiles[0][1].dist[2][1]);

    printf("tiles created");
    //   printf("done reading");
    for (int k = 0; k < MAX / BLOCK_SIZE; k++) {
     //   printf("big loop\n");
//#pragma omp parallel
//        {

//#pragma omp single

//        FW_tile(&(tiles[k][k]),k*BLOCK_SIZE);
        printf("here");
        FW_tile(&(tiles[k][k]),BLOCK_SIZE * k, BLOCK_SIZE * k, BLOCK_SIZE, k * BLOCK_SIZE);
//#pragma omp parallel for schedule(guided)
        for (int i = 0; i < k; i = i + 1) {
//               printf("here1\n");

//            FW_tile(&(tiles[i][k]),k*BLOCK_SIZE);
            FW_tile(&(tiles[i][k]),i*BLOCK_SIZE, k * BLOCK_SIZE, BLOCK_SIZE, k * BLOCK_SIZE);


        }
//#pragma omp parallel for schedule(guided)

        for (int i = (k + 1); i < MAX/BLOCK_SIZE; i = i + 1) {


//                    printf("here2\n");

            FW_tile(&(tiles[i][k]),i*BLOCK_SIZE, k * BLOCK_SIZE, BLOCK_SIZE, k * BLOCK_SIZE);
//            FW_tile(&(tiles[i][k]),k*BLOCK_SIZE);

        }
//#pragma omp parallel for schedule(guided)

        for (int j = 0; j <  k; j = j + 1) {
//                     printf("here3\n");


            FW_tile(&(tiles[k][j]),k * BLOCK_SIZE, j*BLOCK_SIZE, BLOCK_SIZE, k * BLOCK_SIZE);
//            FW_tile(&(tiles[k][j]),k*BLOCK_SIZE);

        }
//#pragma omp parallel for schedule(guided)
        for (int j =  (k + 1); j < MAX/BLOCK_SIZE; j = j + 1) {
//                   printf("here4\n");

            FW_tile(&(tiles[k][j]),k * BLOCK_SIZE, j*BLOCK_SIZE, BLOCK_SIZE, k * BLOCK_SIZE);
//            FW_tile(&(tiles[k][j]),k*BLOCK_SIZE);

        }
//#pragma omp barrier
//#pragma omp parallel for schedule(guided)
        for (int i = 0; i < MAX/BLOCK_SIZE; i = i + 1) {
            for (int j = 0; j < MAX/BLOCK_SIZE; j = j + 1) {
                if ((i != k ) && (j != k )) {
//                                 printf("here5\n");

//                    FW_tile(&(tiles[i][j]),k*BLOCK_SIZE);
                    FW_tile(&(tiles[i][j]),i*BLOCK_SIZE, j*BLOCK_SIZE, BLOCK_SIZE, k * BLOCK_SIZE);


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


    printf("finished fw\n");
    for (int i=0;i<MAX/BLOCK_SIZE;i++){
        for (int j=0;j<MAX/BLOCK_SIZE;j++){
            for (int a=0;a<BLOCK_SIZE;a++){
                for (int b=0;b<BLOCK_SIZE;b++){
                    distance[a+i*BLOCK_SIZE][b+j*BLOCK_SIZE]= tiles[i][j].dist[a+b*BLOCK_SIZE] ;
                }
            }
        }
    }

    double endfw = omp_get_wtime();

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
    printf("fw time: %f\n",(endfw-start));
    printf("total running time: %f\n", (end - start));

//    for (int i=0;i<MAX/BLOCK_SIZE;i++) {
//        for (int j = 0; j < MAX / BLOCK_SIZE; j++) {
//            printTile(tiles[i][j]);
//        }
//    }

//    for (int i=0;i<MAX/BLOCK_SIZE;i++) {
//        for (int j = 0; j < MAX / BLOCK_SIZE; j++) {
//            printTile(tiles[i][j]);
//        }
//    }
    return 0;

}

