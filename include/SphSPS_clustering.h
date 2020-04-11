#include<vector>
#include<algorithm>
#include<string.h>
#include <pthread.h>
#include <semaphore.h>
#include <sys/sysinfo.h>

#include"shortest_path.h"

#ifndef pi
#define pi 3.141593
#endif

#ifndef EPSILON
#define EPSILON 0.0000001
#endif

#ifndef MAX
#define MAX(a, b) ((a)>(b)?(a):(b))
#endif

#ifndef MIN
#define MIN(a, b) ((a)<(b)?(a):(b))
#endif


void xy_to_xyz_sphere(float x, float y, int h, int w, float* x3, float* y3, float* z3) {
    
    //Conversion to 3D sphere
    *x3 = sin(y*pi/h)*cos(x*2*pi/w);
    *y3 = sin(y*pi/h)*sin(x*2*pi/w);
    *z3 = cos(y*pi/h);
    
}


typedef struct{
    
    float * dist;
    float *L1;
    float *L2;
    float *a1;
    float *a2;
    float *b1;
    float *b2;
    float *x1;
    float *x2;
    float *y1;
    float *y2;
    float *centerL1;
    float *centerL2;
    float *centera1;
    float *centera2;
    float *centerb1;
    float *centerb2;
    float * Lab_2;
    unsigned short int* label;
    float * seedsxyz;
    int seedNum;
    int nRows;
    int nCols;
    float lambda;
    float compactness;
    float *contour;
    int path_color;
    int path_contour;
    sem_t* mutex;
    unsigned short int * line_area;
    float S;
    unsigned char * N_path_vect;
    float * xyz;
    float* WSum;
    int* clusterSize;
    int ini;
    int fin;
    int iter;
    int iter_max;
}sphsps_core_arg;



void* SphSPS_core(void*  pArguments)
{
    
    sphsps_core_arg inputs;
    inputs = *(sphsps_core_arg*) pArguments;
    
    float* L1 = inputs.L1;
    float* L2 = inputs.L2;
    float* a1 = inputs.a1;
    float* a2 = inputs.a2;
    float* b1 = inputs.b1;
    float* b2 = inputs.b2;
    float* Lab_2 = inputs.Lab_2;
    float* centerL1 = inputs.centerL1;
    float* centerL2 = inputs.centerL2;
    float* centera1 = inputs.centera1;
    float* centera2 = inputs.centera2;
    float* centerb1 = inputs.centerb1;
    float* centerb2 = inputs.centerb2;
    unsigned short int* label = inputs.label;
    float* seedsxyz = inputs.seedsxyz;
    int seedNum = inputs.seedNum;
    int nRows = inputs.nRows;
    int nCols = inputs.nCols;
    int path_color = inputs.path_color;
    int path_contour = inputs.path_contour;
    float * contour = inputs.contour;
    float compactness = inputs.compactness;
    float lambda = inputs.lambda;
    float *dist = inputs.dist;
    int ini = inputs.ini;
    int fin = inputs.fin;
    sem_t*mutex = inputs.mutex;
    float S = inputs.S;
    unsigned char * N_path_vect = inputs.N_path_vect;
    unsigned short int* line_area = inputs.line_area;
    float* xyz = inputs.xyz;
    float* WSum = inputs.WSum;
    int* clusterSize = inputs.clusterSize;
    int iter = inputs.iter;
    int iter_max = inputs.iter_max;
    
    int nRows_nCols = nRows*nCols;
    float nRows_PI = nRows/pi;
    float nCols_2PI = nCols/(2*pi);
    
    //To store information of already crossed pixels
    float * stack_cont = (float *)malloc(nRows_nCols*sizeof(float));
    float * stack_col = (float *)malloc(nRows_nCols*sizeof(float));
    float * stack_col_total = (float *)malloc(nRows_nCols*sizeof(float));
    unsigned char* visited = (unsigned char*)malloc(nRows_nCols*sizeof(unsigned char));
    
    int minX,minY,maxX,maxY,minX1,maxX1,minX2,maxX2, y, x, pos_i;
    float D, D_s, D_contour, D_color, x3c, y3c, z3c;
    
    
    if (iter > 0) {
        for(int i=ini;i<fin;i++)
        {
            WSum[i]=(WSum[i]==0)?1:WSum[i];
            clusterSize[i]=(clusterSize[i]==0)?1:clusterSize[i];
        }
        for(int i=ini;i<fin;i++)
        {
            centerL1[i]/=WSum[i];
            centerL2[i]/=WSum[i];
            centera1[i]/=WSum[i];
            centera2[i]/=WSum[i];
            centerb1[i]/=WSum[i];
            centerb2[i]/=WSum[i];
            seedsxyz[i]/=clusterSize[i];
            seedsxyz[i+seedNum]/=clusterSize[i];
            seedsxyz[i+seedNum*2]/=clusterSize[i];
        }
    }
    
    
    for(unsigned short int i=ini;i<fin;i++) {
        
        
        float Lab_cc = centerL1[i]*centerL1[i] + centerL2[i]*centerL2[i] + centera1[i]*centera1[i] +
                centera2[i]*centera2[i] + centerb1[i]*centerb1[i] + centerb2[i]*centerb2[i];
        
        
        //////// REGION ////////////
        //sphere coordinates
        x3c = seedsxyz[i];
        y3c = seedsxyz[i+seedNum];
        z3c = seedsxyz[i+seedNum*2];
        
        float phi = acos(z3c);
        float theta = atan2(y3c,x3c);
        y = (int) ((phi*nRows/(pi)));
        x = (int) ((theta)*nCols/(2*pi));
        if (x<0)
            x = nCols+x;
        
        minY = MIN(MAX(round(y - S),0), nRows-1);
        maxY = MIN(MAX(round(y + S),0), nRows-1);
        
        memset(visited,0,nRows*nCols);
        
        for (unsigned short int n = minY; n <=maxY; n++) {
            
            //Deformation for search area
            minX = x - line_area[n];
            maxX = x + line_area[n];
            
            //Handle mirror effect
            minX1 = 0;
            maxX1 = -1;
            minX2 = minX;
            maxX2 = maxX;
            if (minX<0) {
                minX1 = 0;
                maxX1 = maxX;
                minX2 = nCols+minX;
                maxX2 = nCols-1;
            }
            if (maxX>=nCols) {
                minX1 = 0;
                maxX1 = maxX-nCols;
                minX2 = minX;
                maxX2 = nCols-1;
            }
            
            for (int mirror=0; mirror<=1; mirror++) {
                
                if (mirror == 0) {
                    minX = minX1;
                    maxX = maxX1;
                }
                else {
                    minX = minX2;
                    maxX = maxX2;
                }
                
                for(unsigned short int m=minX;m<=maxX;m++) {
                    
                    pos_i = m*nRows+n;
                    
                    D_contour = 1;
                    if (path_color || path_contour) {  //Computation of features on the spherical shortest path
                        
                        if (visited[pos_i]==0) {
                            
                            D_color = geodesic_spherical_path_rec(xyz[pos_i], xyz[pos_i+nRows_nCols], xyz[pos_i+nRows_nCols*2], x3c, y3c, z3c,
                                    N_path_vect[n], nRows, nCols, nRows_PI, nCols_2PI,
                                    L1, L2, a1, a2, b1, b2, Lab_2,
                                    centerL1[i], centerL2[i], centera1[i], centera2[i], centerb1[i], centerb2[i], Lab_cc,
                                    path_contour, contour, &D_contour, stack_col_total, stack_col, stack_cont, visited, &D_s);
                            
                            stack_col_total[pos_i] = D_color;
                            stack_cont[pos_i] = D_contour;
                            
                            //Color distance
                            stack_col[pos_i] =  Lab_cc + Lab_2[pos_i] -
                                    2*(centerL1[i]*L1[pos_i] + centerL2[i]*L2[pos_i] +
                                    centera1[i]*a1[pos_i] + centera2[i]*a2[pos_i] +
                                    centerb1[i]*b1[pos_i] + centerb2[i]*b2[pos_i]);
                            
                        }
                        else {
                            D_color = stack_col_total[pos_i];
                            D_contour = stack_cont[pos_i];
                            //Cosine dissimilarity spatial distance
                            D_s = 1 - (x3c*xyz[pos_i] + y3c*xyz[pos_i+nRows_nCols] + z3c*xyz[pos_i+nRows_nCols*2]);
                            
                        }
                        D = stack_col[pos_i]*lambda;
                        
                    }
                    else {
                        D = Lab_cc + Lab_2[pos_i] -
                                2*(centerL1[i]*L1[pos_i] + centerL2[i]*L2[pos_i] +
                                centera1[i]*a1[pos_i] + centera2[i]*a2[pos_i] +
                                centerb1[i]*b1[pos_i] + centerb2[i]*b2[pos_i]);
                        
                        //Cosine dissimilarity spatial distance
                        D_s = 1 - (x3c*xyz[pos_i] + y3c*xyz[pos_i+nRows_nCols] + z3c*xyz[pos_i+nRows_nCols*2]);
                        
                    }
                    
                    D += D_s*compactness;
                    
                    //Color distance on path
                    if (path_color) {
                        D += D_color*(1-lambda);
                    }
                    
                    //Contour distance on path
                    if (path_contour) {
                        D *= D_contour;
                    }
                    
                    sem_wait(&(mutex[pos_i]));
                    if(D<dist[pos_i])
                    {
                        label[pos_i]=i;
                        dist[pos_i]=D;
                    }
                    sem_post(&(mutex[pos_i]));
                    
                }
                
            }
            
            
            
            
        }
        
        
    }
    
    if (iter < iter_max) {
        for(unsigned short int i=ini;i<fin;i++)
        {
            centerL1[i]=0;
            centerL2[i]=0;
            centera1[i]=0;
            centera2[i]=0;
            centerb1[i]=0;
            centerb2[i]=0;
            WSum[i]=0;
            clusterSize[i]=0;
            seedsxyz[i] = 0;
            seedsxyz[i+seedNum] = 0;
            seedsxyz[i+seedNum*2] = 0;
        }
    }
    
    free(visited);
    free(stack_cont);
    free(stack_col);
    free(stack_col_total);
    
    
    return NULL;
}


void SphSPS_clustering(
        float* L1,
        float* L2,
        float* a1,
        float* a2,
        float* b1,
        float* b2,
        float* x1,
        float* x2,
        float* y1,
        float* y2,
        float* Lab_2,
        float* xy_2,
        float* centerL1,
        float* centerL2,
        float* centera1,
        float* centera2,
        float* centerb1,
        float* centerb2,
        float* centerx1,
        float* centerx2,
        float* centery1,
        float* centery2,
        float* W,
        unsigned short int* label,
        float* seedsx,
        float* seedsy,
        int seedNum,
        int nRows,
        int nCols,
        int StepX,
        int StepY,
        int iterationNum,
        int path_color,
        int path_contour,
        float * contour,
        float compactness,
        float lambda,
        unsigned char* R,
        unsigned char* G,
        unsigned char* B,
        int N_path
        )
{
    
    
    int nRows_nCols = nRows*nCols;
    float* dist = (float *)calloc(nRows_nCols,sizeof(float));
    float* WSum = (float *)calloc(seedNum,sizeof(float));
    int* clusterSize = (int *)calloc(seedNum,sizeof(int));
    
    
    //Pre conversion to spherical coordinates
    float *xyz, *seedsxyz;
    seedsxyz = (float *)calloc(seedNum*3,sizeof(float));
    xyz  = (float *) malloc(nRows_nCols*3*sizeof(float));
    for (int j=0; j<nCols; j++) {
        for (int i=0; i<nRows; i++) {
            xy_to_xyz_sphere(j,i,nRows,nCols,&xyz[i+j*nRows],&xyz[i+j*nRows+nRows_nCols],&xyz[i+j*nRows+nRows_nCols*2]);
        }
    }
    for (int i=0; i<seedNum; i++) {
        xy_to_xyz_sphere(seedsx[i],seedsy[i],nRows,nCols,&seedsxyz[i],&seedsxyz[i+seedNum],&seedsxyz[i+seedNum*2]);
    }
    
    
    //Mean area superpixel
    float S = nCols/(sqrt(seedNum*pi));
    
    //Deform line for region search
    unsigned short int* line_area = (unsigned short int*) malloc(nRows*sizeof(unsigned short int));
    unsigned char* N_path_vect = (unsigned char*) malloc(nRows*sizeof(unsigned char));
    for (int i=0; i<nRows; i++) {
        line_area[i] = (unsigned short int) round(MIN(S/sin(pi*i/nRows),nCols/2));
        N_path_vect[i] = N_path;
    }
    
    // Thread stuff //
    int thread_nbr = get_nprocs();
    
    int step = floor(seedNum/thread_nbr);
    int * slice_vect = (int *) calloc(thread_nbr+1, sizeof(int));
    
    slice_vect[0] = 0;
    for (int i=1; i<thread_nbr; i++)
        slice_vect[i] = i*step;
    slice_vect[thread_nbr] = seedNum;
    
    //Mutex
    sem_t* mutex= (sem_t*) malloc(nRows_nCols*sizeof(sem_t));
    for (int i=0;i<nCols*nRows;i++)    {
        sem_init(&(mutex[i]), 0, 1);
    }
    
    //THREAD BUILDING
    pthread_t *thread_list_th = (pthread_t *) calloc(thread_nbr, sizeof(pthread_t));
    sphsps_core_arg * ThreadArgs = (sphsps_core_arg*)calloc(thread_nbr, sizeof(sphsps_core_arg));
    
    
    //Iterative clustering process
    for(int iteration=0;iteration<=iterationNum;iteration++) {
        
        for(int i=0;i<nCols*nRows;i++)
            dist[i] = FLT_MAX;
        
        /*Launching of the THREADS*/
        for (int i=0; i < thread_nbr; i++) {
            
            /*Thread arguments*/
            ThreadArgs[i].L1 = L1;
            ThreadArgs[i].L2 = L2;
            ThreadArgs[i].a1 = a1;
            ThreadArgs[i].a2 = a2;
            ThreadArgs[i].b1 = b1;
            ThreadArgs[i].b2 = b2;
            ThreadArgs[i].Lab_2 = Lab_2;
            ThreadArgs[i].centerL1 = centerL1;
            ThreadArgs[i].centerL2 = centerL2;
            ThreadArgs[i].centera1 = centera1;
            ThreadArgs[i].centera2 = centera2;
            ThreadArgs[i].centerb1 = centerb1;
            ThreadArgs[i].centerb2 = centerb2;
            ThreadArgs[i].label = label;
            ThreadArgs[i].seedsxyz = seedsxyz;
            ThreadArgs[i].seedNum = seedNum;
            ThreadArgs[i].xyz = xyz;
            
            ThreadArgs[i].mutex = mutex;
            ThreadArgs[i].S = S;
            ThreadArgs[i].N_path_vect = N_path_vect;
            ThreadArgs[i].line_area = line_area;
            
            ThreadArgs[i].path_color = path_color;
            ThreadArgs[i].path_contour = path_contour;
            ThreadArgs[i].contour = contour;
            ThreadArgs[i].compactness = compactness;
            ThreadArgs[i].lambda = lambda;
            ThreadArgs[i].iter = iteration;
            ThreadArgs[i].iter_max = iterationNum;
            
            ThreadArgs[i].ini = slice_vect[i];
            ThreadArgs[i].fin = slice_vect[i+1];
            ThreadArgs[i].dist = dist;
            ThreadArgs[i].nCols = nCols;
            ThreadArgs[i].nRows = nRows;
            ThreadArgs[i].clusterSize = clusterSize;
            ThreadArgs[i].WSum = WSum;
            
            if (pthread_create(&thread_list_th[i], NULL,
                    SphSPS_core, &ThreadArgs[i]))
                printf("Error creating a thread!\n");
        }
        
        /*Wait for all threads to end*/
        for (int j=0; j<thread_nbr; j++) {
            pthread_join(thread_list_th[j],NULL);
        }
        
        //Update clusters
        if (iteration<iterationNum) {
            
            for(unsigned short int j=0;j<nCols;j++)
            {
                for(unsigned short int i=0;i<nRows;i++)
                {
                    int pos = i+nRows*j;
                    int L=label[pos];
                    float Weight=W[pos];
                    
                    centerL1[L]+=Weight*L1[pos];
                    centerL2[L]+=Weight*L2[pos];
                    centera1[L]+=Weight*a1[pos];
                    centera2[L]+=Weight*a2[pos];
                    centerb1[L]+=Weight*b1[pos];
                    centerb2[L]+=Weight*b2[pos];
                    
                    clusterSize[L]++;
                    WSum[L]+=Weight;
                    seedsxyz[L] += xyz[pos];
                    seedsxyz[L+seedNum] += xyz[pos+nRows_nCols];
                    seedsxyz[L+seedNum*2] += xyz[pos+nRows_nCols*2];
                    
                }
            }
            
        }
        
        
        
    }
    
    
    //Free
    free(slice_vect);
    free(mutex);
    free(N_path_vect);
    free(xyz);
    free(seedsxyz);
    free(line_area);
    free(WSum);
    free(clusterSize);
    free(dist);
    
    return;
}

