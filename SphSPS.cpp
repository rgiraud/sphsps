#include "./include/CImg.h"
#include <iostream>
#include <string>

#include"./include/SphSPS.h"

using namespace cimg_library;

int main(int argc, char* argv[])
{
    if (argc < 3) {
        printf("Not enough inputs!\nUsage : ./SphSPS -i img_name [-outm output_map_name] [-outb output_border_name]  [-k superpixel_nbr] [-m compactness] [-c contour]\n");
        return -1;
    }

    //Inputs
    string img_name = string(cimg_option("-i","","Input image file"));
    int SuperpixelNum = cimg_option("-k",1000,"Number of desired superpixels");
    float compactness = cimg_option("-m", 0.12, "Compactness value");
    string contour_input = string(cimg_option("-c","","Input contour image file"));

    //Outputs
    string output_map_name = string(cimg_option("-outm","labelmap.png","Output Labeled image file"));
    string output_border_name = string(cimg_option("-outb","borders.png","Output borders of superPixels with original image as background"));

    unsigned char* R,* G,* B;
    unsigned short* label;

    //Default parameters
    int path_color = 1;
    int path_contour = 0;

    //Image loading
    cout << img_name.c_str() << "\n";
    CImg<unsigned char> img_in(img_name.c_str());
    CImg<unsigned char> img = img_in;
    int nCols = img.width();
    int nRows = img.height();

    float * contour = (float*)malloc(nRows*nCols*sizeof(float));

    if(contour_input!="") {

        path_contour = 1;
        CImg<int> contour_img(contour_input.c_str());
        if ( (contour_img.width() != nCols) || (contour_img.height() != nRows) ) {
            printf("Contour dimensions do not match image size!\n");
            return -1;
        }
        for (int j=0; j<nCols; j++) {
            for (int i=0; i<nRows; i++) {
                contour[i+j*nRows] = 1 + 10*contour_img(j,i)/255;
            }
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


    for (int j=0; j<nCols; j++) {
        for (int i=0; i<nRows; i++)  {
            R[i+j*nRows] = img(j,i,0,0);
            G[i+j*nRows] = img(j,i,0,1);
            B[i+j*nRows] = img(j,i,0,2);
        }
    }

    //SphSPS
    SphSPS_setup(R,G,B,nRows,nCols,SuperpixelNum,compactness,label,path_color,path_contour,contour);

    // SAVE OUTPUTS

    // Output Label image (main output for many applications)
    //int max_sp = 0;
    CImg<int> output(nCols,nRows);
    for (int i=0; i<nCols; i++){
        for (int j=0; j<nRows; j++) {
            output(i,j) = (int) label[j+i*nRows];
//             if (output(i,j) > max_sp)
//                 max_sp = output(i,j);
        }
    }

//     printf("max %d\n", max_sp);

    char str[100];
    strcpy(str, "res/");
    strcat(str, output_map_name.c_str());
    output.save(str);


    // Output borders of SuperPixels with original image as background
    // Draw borders in a slice by slice fashion for 3d images
    CImg<> output_border = img;
    int v4x[]={-1,0,1,0};
    int v4y[]={0,-1,0,1};

    cimg_forZ(img,z) {
        cimg_forXY(output_border,x,y) {
            int lab1=output(x,y,z);
            for(int k=0;k<4;k++)
                if(output_border.containsXYZC(x+v4x[k],y+v4y[k],z))
                    if(lab1 != output(x+v4x[k],y+v4y[k],z))
                        cimg_forC(output_border,c)
                        output_border(x,y,z,c)=0;
        }
    }

    char str2[100];
    strcpy(str2, "res/");
    strcat(str2, output_border_name.c_str());
    output_border.save(str2);


    delete [] R;
    delete [] G;
    delete [] B;
    free(label);
    free(contour);

    return 0;
}
