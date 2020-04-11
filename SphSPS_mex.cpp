#include<mex.h>
#include<matrix.h>
#include"./include/SphSPS.h"


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    unsigned char* R,* G,* B;
    int SuperpixelNum = 0;
    int nRows = 0; int nCols = 0;
    unsigned short* label;
    unsigned short* output;
    unsigned char* img;
    
    //Default parameters
    float compactness = 0.12;
    int path_color = 1;
    int path_contour = 0;
    
    if(nrhs>=2)
    {
        if(mxGetNumberOfDimensions(prhs[0])!=3)
            mexErrMsgTxt("The input image must be in CIERGB form");
        if(mxGetClassID(prhs[0])!=mxUINT8_CLASS)
            mexErrMsgTxt("The input image must be in CIERGB form");
        nRows = (int) mxGetM(prhs[0]);
        nCols = (int) mxGetN(prhs[0])/3;
        
        img=(unsigned char*)mxGetPr(prhs[0]);
        SuperpixelNum= (int) mxGetScalar(prhs[1]);
    }
    else if(nrhs>4)
        mexErrMsgTxt("Too many inputs!");
    
    
    float * contour = (float*)malloc(nRows*nCols*sizeof(float));
    
    if(nrhs>=3)
        compactness= (float) mxGetScalar(prhs[2]);
    
    if (nrhs==4) {
        path_contour = 1;
        float * contour_init = (float*) mxGetPr(prhs[3]);
        for (int i=0; i<nRows*nCols; i++) {
            contour[i] = 1 + 10*contour_init[i];
        }
    }
    else {
        for (int i=0; i<nRows*nCols; i++)
            contour[i] = 1;
    }
    
    
    
    int pixel=nRows*nCols;
    R=new unsigned char[pixel];
    G=new unsigned char[pixel];
    B=new unsigned char[pixel];
    label = (unsigned short*) calloc(pixel,sizeof(unsigned short));
    
    
    for(int i=0;i<pixel;i++)
    {
        R[i]=img[i];
        G[i]=img[i+pixel];
        B[i]=img[i+pixel+pixel];
    }
    
    
    //SphSPS
    SphSPS_setup(R,G,B,nRows,nCols,SuperpixelNum,compactness,label,path_color,path_contour,contour);
    
    
    //Output label map
    plhs[0]=mxCreateNumericMatrix(nRows,nCols,mxUINT16_CLASS,mxREAL);
    output=(unsigned short*)mxGetPr(plhs[0]);
    for(int i=0;i<pixel;i++)
        output[i]=label[i];
    
    delete [] R;
    delete [] G;
    delete [] B;
    free(label);
    free(contour);
}
