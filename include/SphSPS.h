#include "myrgb2lab.h"
#include "hammersley_sampling.h"
#include "initialize.h"
#include "bilateral_lab_filtering.h"
#include "SphSPS_clustering.h"
#include "enforce_connectivity.h"


void SphSPS_setup(unsigned char* R, unsigned char* G, unsigned char* B, int nRows, int nCols, int seedNum, float compactness, unsigned short* label,
        int path_color, int path_contour, float*contour)
{
    //Setting Parameters
    float colorCoefficient=20;
    float distCoefficient=colorCoefficient*compactness;
    int thresholdCoef=4;
    int iteration = 5;
    int N_path = 15;
    
    unsigned char *L, *a, *b;
    L=new unsigned char[nRows*nCols];
    a=new unsigned char[nRows*nCols];
    b=new unsigned char[nRows*nCols];
    
    //Conversion to LAB
    myrgb2lab1(R,G,B,L,a,b,nRows,nCols);
    
    /////////////// Produce Seeds ///////////// 
    int newSeedNum = seedNum;
    float *seedsx = (float *) calloc(seedNum,sizeof(float));
    float *seedsy = (float *) calloc(seedNum,sizeof(float));
    int ColNum,RowNum,StepY,StepX;
    ColNum=sqrt(float(seedNum*nCols/nRows));
    RowNum=seedNum/ColNum;
    StepY=nRows/RowNum;
    StepX=nCols/ColNum;
    //Hammersely spherical sampling
    seeds_sp_sampling_hammersley(seedNum, nRows, nCols, seedsx, seedsy);

        
    compactness = 1.0/((StepY/compactness)*(StepX/compactness));
    float lambda = 1;
    if (path_color == 1)
        lambda = 0.5;
    
    
    //Initialization 6-Lab features
    float *L1,*L2,*a1,*a2,*b1,*b2,*x1,*x2,*y1,*y2;
    float *W;
    L1=new float[nRows*nCols];
    L2=new float[nRows*nCols];
    a1=new float[nRows*nCols];
    a2=new float[nRows*nCols];
    b1=new float[nRows*nCols];
    b2=new float[nRows*nCols];
    x1=new float[nRows*nCols];
    x2=new float[nRows*nCols];
    y1=new float[nRows*nCols];
    y2=new float[nRows*nCols];
    W =new float[nRows*nCols];
   
    Initialize(L,a,b,L1,L2,a1,a2,b1,b2,x1,x2,y1,y2,W,nRows,nCols,StepX,StepY,colorCoefficient,distCoefficient);
    delete [] L;
    delete [] a;
    delete [] b;
    
    
    //Filtered and average features - set pw = 0 for no action 
    float* Lab_2=new float[nCols*nRows]();
    float* xy_2=new float[nCols*nRows]();
    int pw = 2;
    bilateral_lab_features(L1, L2, a1, a2, b1, b2, x1, x2, y1, y2, Lab_2, xy_2, nRows, nCols, pw, R, G, B);
    
    
    //Superpixel seeds initialization
    float* centerL1=new float[newSeedNum]();
    float* centerL2=new float[newSeedNum]();
    float* centera1=new float[newSeedNum]();
    float* centera2=new float[newSeedNum]();
    float* centerb1=new float[newSeedNum]();
    float* centerb2=new float[newSeedNum]();
    float* centerx1=new float[newSeedNum]();
    float* centerx2=new float[newSeedNum]();
    float* centery1=new float[newSeedNum]();
    float* centery2=new float[newSeedNum]();
    int step_init = 2;
    for(int i=0;i<newSeedNum;i++)
    {
        int x = seedsx[i];
        int y = seedsy[i];
        int minX=(x-step_init<=0)?0:x-step_init;
        int minY=(y-step_init<=0)?0:y-step_init;
        int maxX=(x+step_init>=nCols-1)?nCols-1:x+step_init;
        int maxY=(y+step_init>=nRows-1)?nRows-1:y+step_init;
        int Count=0;
        
        //Centered on the seeds
        for(int j=minX;j<=maxX;j++)
            for(int k=minY;k<=maxY;k++)
            {
                Count++;
                
                int pos = k+j*nRows;
                //Remplacer par features filtrés
                centerL1[i]+=L1[pos];
                centerL2[i]+=L2[pos];
                centera1[i]+=a1[pos];
                centera2[i]+=a2[pos];
                centerb1[i]+=b1[pos];
                centerb2[i]+=b2[pos];
                centerx1[i]+=x1[pos];
                centerx2[i]+=x2[pos];
                centery1[i]+=y1[pos];
                centery2[i]+=y2[pos];
                
            }
        centerL1[i]/=Count;
        centerL2[i]/=Count;
        centera1[i]/=Count;
        centera2[i]/=Count;
        centerb1[i]/=Count;
        centerb2[i]/=Count;
        centerx1[i]/=Count;
        centerx2[i]/=Count;
        centery1[i]/=Count;
        centery2[i]/=Count;
    }
    
    
    
    //Superpixel Clustering
    SphSPS_clustering(L1, L2, a1, a2, b1, b2, x1, x2, y1, y2, Lab_2, xy_2,
                     centerL1, centerL2, centera1, centera2, centerb1, centerb2, centerx1, centerx2, centery1, centery2, 
                     W, label, seedsx, seedsy, newSeedNum, nRows, nCols, StepX, StepY, iteration,
                    path_color, path_contour, contour, compactness, lambda, R, G, B, N_path);
 
    
    //Enforce Connectivity
    int threshold= (nCols*nRows)/(newSeedNum*thresholdCoef);
    preEnforceConnectivity(label,nRows,nCols);
    EnforceConnectivity(L1,L2,a1,a2,b1,b2,x1,x2,y1,y2,W,label,threshold,nRows,nCols);
 
    
    //Clear Memory
    delete []L1;
    delete []L2;
    delete []a1;
    delete []a2;
    delete []b1;
    delete []b2;
    delete []x1;
    delete []x2;
    delete []y1;
    delete []y2;
    delete []W;
    free(seedsx);
    free(seedsy);
        
    //Clear Memory
    delete []centerL1;
    delete []centerL2;
    delete []centera1;
    delete []centera2;
    delete []centerb1;
    delete []centerb2;
    delete []centerx1;
    delete []centerx2;
    delete []centery1;
    delete []centery2;
    
    delete []Lab_2;
    delete []xy_2;
    
    
}

