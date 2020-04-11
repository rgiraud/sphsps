#ifndef INITIALIZE
#define INITIALIZE

#include<cmath>
using namespace std;

//map pixels into ten dimensional feature space

const float PI=3.1415926;

void Initialize(
        unsigned char* L,
        unsigned char* a,
        unsigned char* b,
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
        float* W,
        int nRows,
        int nCols,
        int StepX,
        int StepY,
        float Color,
        float Distance
        )
{
    
    
    float thetaL,thetaa,thetab,thetax,thetay;
    for(int j=0;j<nCols;j++)
        for(int i=0;i<nRows;i++)
        {
            int pos = i+j*nRows;
            thetaL=((float)L[pos]/(float)255)*PI/2;
            thetaa=((float)a[pos]/(float)255)*PI/2;
            thetab=((float)b[pos]/(float)255)*PI/2;
            thetax=((float)j/(float)StepX)*PI/2;
            thetay=((float)i/(float)StepY)*PI/2;
            L1[pos]=Color*cos(thetaL);
            L2[pos]=Color*sin(thetaL);
            a1[pos]=Color*cos(thetaa)*2.55;
            a2[pos]=Color*sin(thetaa)*2.55;
            b1[pos]=Color*cos(thetab)*2.55;
            b2[pos]=Color*sin(thetab)*2.55;
            x1[pos]=Distance*cos(thetax);
            x2[pos]=Distance*sin(thetax);
            y1[pos]=Distance*cos(thetay);
            y2[pos]=Distance*sin(thetay);
        }
    float sigmaL1=0,sigmaL2=0,sigmaa1=0,sigmaa2=0,sigmab1=0,sigmab2=0,sigmax1=0,sigmax2=0,sigmay1=0,sigmay2=0;
    float size=nRows*nCols;
    for(int j=0;j<nCols;j++)
        for(int i=0;i<nRows;i++)
        {
            int pos = i+j*nRows;
            sigmaL1+=L1[pos];
            sigmaL2+=L2[pos];
            sigmaa1+=a1[pos];
            sigmaa2+=a2[pos];
            sigmab1+=b1[pos];
            sigmab2+=b2[pos];
            sigmax1+=x1[pos];
            sigmax2+=x2[pos];
            sigmay1+=y1[pos];
            sigmay2+=y2[pos];
        }
    sigmaL1/=size;
    sigmaL2/=size;
    sigmaa1/=size;
    sigmaa2/=size;
    sigmab1/=size;
    sigmab2/=size;
    sigmax1/=size;
    sigmax2/=size;
    sigmay1/=size;
    sigmay2/=size;
    
    for(int j=0;j<nCols;j++)
        for(int i=0;i<nRows;i++){
            int pos = i+j*nRows;
            W[pos]=L1[pos]*sigmaL1+
                    L2[pos]*sigmaL2+
                    a1[pos]*sigmaa1+
                    a2[pos]*sigmaa2+
                    b1[pos]*sigmab1+
                    b2[pos]*sigmab2+
                    x1[pos]*sigmax1+
                    x2[pos]*sigmax2+
                    y1[pos]*sigmay1+
                    y2[pos]*sigmay2;
            L1[pos]/=W[pos];
            L2[pos]/=W[pos];
            a1[pos]/=W[pos];
            a2[pos]/=W[pos];
            b1[pos]/=W[pos];
            b2[pos]/=W[pos];
            x1[pos]/=W[pos];
            x2[pos]/=W[pos];
            y1[pos]/=W[pos];
            y2[pos]/=W[pos];
        }
    return;
}

#endif