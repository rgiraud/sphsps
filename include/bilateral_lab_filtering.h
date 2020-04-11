
#ifndef MAX
#define MAX(a, b) ((a)>(b)?(a):(b))
#endif


void lookup_table2(float * g_s, float * g_r, int max_hw, float sigma_s, float sigma_c) {
    
    float  s2  = sigma_s*sigma_s;
    for(int i=0; i< max_hw; i++){
        float v = exp(-0.5*(i*i)/(s2));
        if(v<0.1){
            break;
        }
        g_s[i]=v;
    }
    //pre-compute range gaussian
    float r2 = sigma_c*sigma_c;
    for(int i=0; i<256; i++){
        g_r[i] = exp(-0.5*(i*i)/r2);
    }
    
}


//BILATERAL
void bilateral_lab_features(
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
        int nRows,
        int nCols,
        int pw,
        unsigned char* R,
        unsigned char* G,
        unsigned char* B
        )
{
    
    float* L1_1=new float[nCols*nRows]();
    float* L2_1=new float[nCols*nRows]();
    float* a1_1=new float[nCols*nRows]();
    float* b1_1=new float[nCols*nRows]();
    float* a2_1=new float[nCols*nRows]();
    float* b2_1=new float[nCols*nRows]();
    
    float sigma_s = FLT_MAX;
    float sigma_c = 40;
    float *g_s = (float*) calloc(MAX(nRows,nCols),sizeof(float));
    float *g_r = (float*) calloc(256,sizeof(float));
    lookup_table2(g_s,g_r,MAX(nRows,nCols),sigma_s,sigma_c);
    
    for(int i=0;i<nCols;i++) {
        for(int j=0;j<nRows;j++) {
            float count = 0;
            
            int pos_i = j+i*nRows;
            
            for (int dx=-pw; dx<=pw; dx++) {
                for (int dy=-pw; dy<=pw; dy++) {
                    if ((i+dx<nCols)&&(i+dx>=0)&(j+dy>=0)&(j+dy<nRows)){
                        
                        int pos_d = (i+dx)*nRows+(j+dy);
                        float d_s = g_s[abs(dx)+abs(dy)];
                        
                        //New RGB
                        float d_rr = g_r[abs(R[pos_i]-R[pos_d])];
                        float d_rg = g_r[abs(G[pos_i]-G[pos_d])];
                        float d_rb = g_r[abs(B[pos_i]-B[pos_d])];
                        d_rr = (d_rr+d_rg+d_rb)/3;
                        
                        float kk = d_rr*d_s;
                        
                        L1_1[pos_i] += (float) (L1[pos_d]*kk);
                        L2_1[pos_i] += (float) (L2[pos_d]*kk);
                        a1_1[pos_i] += (float) (a1[pos_d]*kk);
                        b1_1[pos_i] += (float) (b1[pos_d]*kk);
                        a2_1[pos_i] += (float) (a2[pos_d]*kk);
                        b2_1[pos_i] += (float) (b2[pos_d]*kk);
                        Lab_2[pos_i] += (float) (L1[pos_d]*L1[pos_d] +
                                L2[pos_d]*L2[pos_d] +
                                a1[pos_d]*a1[pos_d] +
                                a2[pos_d]*a2[pos_d] +
                                b1[pos_d]*b1[pos_d] +
                                b2[pos_d]*b2[pos_d])*kk;
                        
                        count += (float) kk;
                    }
                }
            }
            L1_1[pos_i]  /= count;
            L2_1[pos_i]  /= count;
            a1_1[pos_i]  /= count;
            b1_1[pos_i]  /= count;
            a2_1[pos_i]  /= count;
            b2_1[pos_i]  /= count;
            
            Lab_2[pos_i] /= count;
            
        }
    }
    
    //Output recopy
    for (int i=0; i<nCols*nRows; i++) {
        L1[i] = L1_1[i];
        L2[i] = L2_1[i];
        a1[i] = a1_1[i];
        a2[i] = a2_1[i];
        b1[i] = b1_1[i];
        b2[i] = b2_1[i];
    }
    
    delete []L1_1;
    delete []L2_1;
    delete []a1_1;
    delete []a2_1;
    delete []b1_1;
    delete []b2_1;
    
}
