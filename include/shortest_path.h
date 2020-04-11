#ifndef EPSILON
#define EPSILON 0.00000001
#endif

#ifndef MAX
#define MAX(a, b) ((a)>(b)?(a):(b))
#endif

#ifndef MIN
#define MIN(a, b) ((a)<(b)?(a):(b))
#endif


void orthogonalization_xyz(float x1, float y1, float z1, float *x2, float *y2, float *z2, float proj)
{
    
    *x2 = *x2 - proj*x1;
    *y2 = *y2 - proj*y1;
    *z2 = *z2 - proj*z1;

    float norm = sqrt(1 - proj*proj);
    
    if (norm > EPSILON) {
        *x2 /= norm;
        *y2 /= norm;
        *z2 /= norm;
    }
    
}


float geodesic_spherical_path_rec_core(
        int nRows, int nCols, float nRows_PI, float nCols_2PI,
        float* L1, float* L2, float* a1, float* a2, float* b1, float* b2, float* Lab_2,
        float centerL1, float centerL2, float centera1, float centera2, float centerb1, float centerb2, float Lab_cc,
        float xs1, float ys1, float zs1, float xs2, float ys2, float zs2,
        float * contour, float * D_contour, float * stack_col_total, float *stack_col,  float * stack_cont, unsigned char *visited,
        int step, int pos_1, float angle_, unsigned short int N_path) {
    
    int pos = 0;
    
    //New intermediate angle between Xs1 and Xs2
    float phi = ((float)step)*angle_;
    float cphi = cos(phi);
    float sphi = sin(phi);
    
    //Phi - longitude
    float   ccos = cphi*zs1 + sphi*zs2;
    if (fabs(ccos)>0.99999){
        if (ccos > 0)
            ccos = 0.99999;
        else
            ccos = -0.99999;
    }
    float phib = (float) acos(ccos);
    
    //Theta - latitude
    float thetab = (float) atan2(cphi*ys1 + sphi*ys2, cphi*xs1 + sphi*xs2);
    int xx = (int) (thetab*nCols_2PI);
    if (xx<0)
        pos = ((int) (phib*nRows_PI)) + (xx+nCols)*nRows;
    else
        pos = ((int) (phib*nRows_PI)) + xx*nRows;
    
    
    float Dc = contour[pos];
    
    if (pos != pos_1) {
        
        if (visited[pos]==0) {
            stack_col[pos] = Lab_cc + Lab_2[pos] - 2*(centerL1*L1[pos] + centerL2*L2[pos] +
                    centera1*a1[pos] + centera2*a2[pos] +
                    centerb1*b1[pos] + centerb2*b2[pos]);
            
            
            if (step+1<N_path) {
                stack_col_total[pos] = stack_col[pos] +
                        geodesic_spherical_path_rec_core(nRows, nCols, nRows_PI, nCols_2PI,
                        L1, L2, a1, a2, b1, b2, Lab_2,
                        centerL1, centerL2, centera1, centera2, centerb1, centerb2, Lab_cc,
                        xs1, ys1, zs1, xs2, ys2, zs2,
                        contour, &Dc, stack_col_total, stack_col, stack_cont, visited, step+1, pos, angle_, N_path);
                
                stack_col_total[pos] /= ((float)N_path-step);
                stack_cont[pos] = Dc;
                
            }
            else {
                stack_col_total[pos] = stack_col[pos];
                stack_cont[pos] = contour[pos];
            }
            
            visited[pos] = 1;
            
        }
    }
    else {
        if (step+1<N_path) {
            stack_col_total[pos] = stack_col[pos] + geodesic_spherical_path_rec_core(nRows, nCols, nRows_PI, nCols_2PI,
                    L1, L2, a1, a2, b1, b2, Lab_2,
                    centerL1, centerL2, centera1, centera2, centerb1, centerb2, Lab_cc,
                    xs1, ys1, zs1, xs2, ys2, zs2,
                    contour, &Dc, stack_col_total, stack_col, stack_cont, visited, step+1, pos, angle_, N_path);
            
            stack_col_total[pos] /= ((float)N_path-step);
            stack_cont[pos] = Dc;
            
        }
        else {
            return 0;
        }
    }
    
    *D_contour = MAX(stack_cont[pos],*D_contour);
    
    return stack_col_total[pos]*((float)N_path-step);
    
    
}



float geodesic_spherical_path_rec(float xs1, float ys1, float zs1, float xs2, float ys2, float zs2,
        unsigned char N_path, int nRows, int nCols, float nRows_PI, float nCols_2PI,
        float* L1, float* L2, float* a1, float* a2, float* b1, float* b2, float* Lab_2,
        float centerL1, float centerL2, float centera1, float centera2, float centerb1, float centerb2, float Lab_cc,
        int bres_contour, float * contour, float * D_contour, float * stack_col_total, float *stack_col,  float * stack_cont, unsigned char *visited,
        float*D_s) {
    
    float proj = (xs1*xs2 + ys1*ys2 + zs1*zs2) - 0.000001;  //To avoid error in acos(proj)
    float angle_ = (float) acos(proj)/((float)N_path-1);
    
    //Orthogonalization process to build the orthogonal coordinate system
    orthogonalization_xyz(xs1, ys1, zs1, &xs2, &ys2, &zs2, proj);
    
    //Recursive shortest path route
    float D = geodesic_spherical_path_rec_core(nRows, nCols, nRows_PI, nCols_2PI,
            L1, L2, a1, a2, b1, b2, Lab_2,
            centerL1, centerL2, centera1, centera2, centerb1, centerb2, Lab_cc,
            xs1, ys1, zs1, xs2, ys2, zs2,
            contour, D_contour, stack_col_total, stack_col, stack_cont, visited, 
            0, -1, angle_, (int) N_path);
    
    *D_s = 1 - proj; //Spatial distance ds
    
    return D/((float)N_path);
    
}







