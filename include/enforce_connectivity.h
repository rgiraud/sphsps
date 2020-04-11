#include<queue>
#include<vector>
#include<float.h>

#ifndef MAX
#define MAX(a, b) ((a)>(b)?(a):(b))
#endif

#ifndef MIN
#define MIN(a, b) ((a)<(b)?(a):(b))
#endif

class Superpixel
{
public:
    int Label;
    int Size;
    vector<int> Neighbor;
    Superpixel(int L=0,int S=0):Label(L),Size(S){}
    vector<int> xLoc;
    vector<int> yLoc;
    friend bool operator==(Superpixel& S,int L);
    friend bool operator==(int L,Superpixel& S);
};

bool operator==(Superpixel& S,int L)
{
    return S.Label==L;
}

bool operator==(int L,Superpixel& S)
{
    return S.Label==L;
}




//Enforce Connectivity by merging very small superpixels with their neighbors

void preEnforceConnectivity(unsigned short int* label, int nRows, int nCols)
{
    const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
    const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};
    int adj=0;
    unsigned int Bond=20;
    bool *mask=new bool[nRows*nCols];
    for(int i=0;i<nRows*nCols;i++)
        mask[i]=0;
    
    vector<unsigned short> xLoc;
    vector<unsigned short> yLoc;
    for(int i=0;i<nRows;i++)
        for(int j=0;j<nCols;j++)
        {
            if(mask[i+j*nRows]==0)
            {
                int L=label[i+j*nRows];
                for(int k=0;k<8;k++)
                {
                    int x=j+dx8[k];
                    int y=i+dy8[k];
                    if(x>=0&&x<=nCols-1&&y>=0&&y<=nRows-1)
                    {
                        if(mask[y+x*nRows]==1&&label[x*nRows+y]!=L) {
                            adj=label[x*nRows+y];break;
                        }
                        
                        
                    }
                    
                }
                mask[i+j*nRows]=1;
                xLoc.insert(xLoc.end(),j);
                yLoc.insert(yLoc.end(),i);
                unsigned int indexMarker=0;
                while(indexMarker<xLoc.size())
                {
                    int x=xLoc[indexMarker];int y=yLoc[indexMarker];
                    indexMarker++;
                    int minX=x-1;
                    int maxX=x+1;
                    int minY=(y-1<=0)?0:y-1;
                    int maxY=(y+1>=nRows-1)?nRows-1:y+1;
                    for(int m=MAX(minX,0);m<=MIN(maxX,nCols-1);m++){
                        for(int n=minY;n<=maxY;n++)
                        {
                            if(mask[n+m*nRows]==0&&label[n+m*nRows]==L)
                            {
                                mask[n+m*nRows]=1;
                                xLoc.insert(xLoc.end(),m);
                                yLoc.insert(yLoc.end(),n);
                            }
                        }
                    }
                    if (minX<0) {
                        int m = nCols-1;
                        for(int n=minY;n<=maxY;n++)
                        {
                            if(mask[n+m*nRows]==0&&label[n+m*nRows]==L)
                            {
                                mask[n+m*nRows]=1;
                                xLoc.insert(xLoc.end(),m);
                                yLoc.insert(yLoc.end(),n);
                            }
                        }
                    }
                    if (maxX>=nCols) {
                        int m = 0;
                        for(int n=minY;n<=maxY;n++)
                        {
                            if(mask[n+m*nRows]==0&&label[n+m*nRows]==L)
                            {
                                mask[n+m*nRows]=1;
                                xLoc.insert(xLoc.end(),m);
                                yLoc.insert(yLoc.end(),n);
                            }
                        }
                    }
                    
                    
                }
                if(indexMarker<Bond)
                {
                    for(int k=0;k<(int) xLoc.size();k++)
                    {
                        int x=xLoc[k];int y=yLoc[k];
                        label[x*nRows+y]=adj;
                    }
                }
                xLoc.clear();
                yLoc.clear();
            }
        }
    
    delete [] mask;
}



void EnforceConnectivity(
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
        unsigned short int* label,
        int threshold,
        int nRows,
        int nCols
        )
{
    unsigned char* mask=new unsigned char[nRows*nCols];
    for(int i=0;i<nRows*nCols;i++)
        mask[i]=0;
    
    vector<unsigned short>strayX;
    vector<unsigned short>strayY;
    vector<unsigned short>Size;
    queue<unsigned short> xLoc;
    queue<unsigned short> yLoc;
    vector<float>centerL1;
    vector<float>centerL2;
    vector<float>centera1;
    vector<float>centera2;
    vector<float>centerb1;
    vector<float>centerb2;
    vector<float>centerx1;
    vector<float>centerx2;
    vector<float>centery1;
    vector<float>centery2;
    vector<float>centerW;
    int sLabel=-1;
    int pos = 0;
    int L;
    for(int i=0;i<nRows;i++)
        for(int j=0;j<nCols;j++)
        {
            pos = i+j*nRows;
            if(mask[pos]==0)
            {
                sLabel++;
                int Count=1;
                centerL1.insert(centerL1.end(),0);
                centerL2.insert(centerL2.end(),0);
                centera1.insert(centera1.end(),0);
                centera2.insert(centera2.end(),0);
                centerb1.insert(centerb1.end(),0);
                centerb2.insert(centerb2.end(),0);
                centerx1.insert(centerx1.end(),0);
                centerx2.insert(centerx2.end(),0);
                centery1.insert(centery1.end(),0);
                centery2.insert(centery2.end(),0);
                centerW.insert(centerW.end(),0);
                strayX.insert(strayX.end(),j);
                strayY.insert(strayY.end(),i);
                float Weight=W[pos];
                centerL1[sLabel]+=L1[pos]*Weight;
                centerL2[sLabel]+=L2[pos]*Weight;
                centera1[sLabel]+=a1[pos]*Weight;
                centera2[sLabel]+=a2[pos]*Weight;
                centerb1[sLabel]+=b1[pos]*Weight;
                centerb2[sLabel]+=b2[pos]*Weight;
                centerx1[sLabel]+=x1[pos]*Weight;
                centerx2[sLabel]+=x2[pos]*Weight;
                centery1[sLabel]+=y1[pos]*Weight;
                centery2[sLabel]+=y2[pos]*Weight;
                centerW[sLabel]+=W[pos];
                L=label[pos];
                label[pos]=sLabel;
                mask[pos]=1;
                xLoc.push(j);yLoc.push(i);
                while(!xLoc.empty())
                {
                    int x=xLoc.front();xLoc.pop();
                    int y=yLoc.front();yLoc.pop();
                    int minX=x-1;
                    int maxX=x+1;
                    int minY=(y-1<=0)?0:y-1;
                    int maxY=(y+1>=nRows-1)?nRows-1:y+1;
                    for(int m=MAX(minX,0);m<=MIN(maxX,nCols-1);m++) {
                        for(int n=minY;n<=maxY;n++)
                        {
                            pos = m*nRows+n;
                            if(mask[pos]==0&&label[pos]==L)
                            {
                                Count++;
                                xLoc.push(m);
                                yLoc.push(n);
                                mask[pos]=1;
                                label[pos]=sLabel;
                                Weight=W[pos];
                                centerL1[sLabel]+=L1[pos]*Weight;
                                centerL2[sLabel]+=L2[pos]*Weight;
                                centera1[sLabel]+=a1[pos]*Weight;
                                centera2[sLabel]+=a2[pos]*Weight;
                                centerb1[sLabel]+=b1[pos]*Weight;
                                centerb2[sLabel]+=b2[pos]*Weight;
                                centerx1[sLabel]+=x1[pos]*Weight;
                                centerx2[sLabel]+=x2[pos]*Weight;
                                centery1[sLabel]+=y1[pos]*Weight;
                                centery2[sLabel]+=y2[pos]*Weight;
                                centerW[sLabel]+=W[pos];
                            }
                        }
                    }
                    if (minX<0) {
                        int m = nCols-1;
                        for(int n=minY;n<=maxY;n++)
                        {
                            pos = m*nRows+n;
                            if(mask[pos]==0&&label[pos]==L)
                            {
                                Count++;
                                xLoc.push(m);
                                yLoc.push(n);
                                mask[pos]=1;
                                label[pos]=sLabel;
                                Weight=W[pos];
                                centerL1[sLabel]+=L1[pos]*Weight;
                                centerL2[sLabel]+=L2[pos]*Weight;
                                centera1[sLabel]+=a1[pos]*Weight;
                                centera2[sLabel]+=a2[pos]*Weight;
                                centerb1[sLabel]+=b1[pos]*Weight;
                                centerb2[sLabel]+=b2[pos]*Weight;
                                centerx1[sLabel]+=x1[pos]*Weight;
                                centerx2[sLabel]+=x2[pos]*Weight;
                                centery1[sLabel]+=y1[pos]*Weight;
                                centery2[sLabel]+=y2[pos]*Weight;
                                centerW[sLabel]+=W[pos];
                            }
                        }
                    }
                    if (maxX>=nCols) {
                        int m = 0;
                        for(int n=minY;n<=maxY;n++)
                        {
                            pos = m*nRows+n;
                            if(mask[pos]==0&&label[pos]==L)
                            {
                                Count++;
                                xLoc.push(m);
                                yLoc.push(n);
                                mask[pos]=1;
                                label[pos]=sLabel;
                                Weight=W[pos];
                                centerL1[sLabel]+=L1[pos]*Weight;
                                centerL2[sLabel]+=L2[pos]*Weight;
                                centera1[sLabel]+=a1[pos]*Weight;
                                centera2[sLabel]+=a2[pos]*Weight;
                                centerb1[sLabel]+=b1[pos]*Weight;
                                centerb2[sLabel]+=b2[pos]*Weight;
                                centerx1[sLabel]+=x1[pos]*Weight;
                                centerx2[sLabel]+=x2[pos]*Weight;
                                centery1[sLabel]+=y1[pos]*Weight;
                                centery2[sLabel]+=y2[pos]*Weight;
                                centerW[sLabel]+=W[pos];
                            }
                        }
                    }
                    
                    
                }
                Size.insert(Size.end(),Count);
                centerL1[sLabel]/=centerW[sLabel];
                centerL2[sLabel]/=centerW[sLabel];
                centera1[sLabel]/=centerW[sLabel];
                centera2[sLabel]/=centerW[sLabel];
                centerb1[sLabel]/=centerW[sLabel];
                centerb2[sLabel]/=centerW[sLabel];
                centerx1[sLabel]/=centerW[sLabel];
                centerx2[sLabel]/=centerW[sLabel];
                centery1[sLabel]/=centerW[sLabel];
                centery2[sLabel]/=centerW[sLabel];
            }
        }
    
    sLabel=sLabel+1;
    
    vector<int>::iterator Pointer;
    vector<Superpixel> Sarray;
    for(int i=0;i<sLabel;i++)
    {
        if(Size[i]<threshold)
        {
            int x=strayX[i];int y=strayY[i];
            L=label[x*nRows+y];
            mask[y+x*nRows]=0;
            int indexMark=0;
            Superpixel S(L,Size[i]);
            S.xLoc.insert(S.xLoc.end(),x);
            S.yLoc.insert(S.yLoc.end(),y);
            while(indexMark<(int) S.xLoc.size())
            {
                x=S.xLoc[indexMark]; y=S.yLoc[indexMark];
                indexMark++;
                int minX=x-1;
                int maxX=x+1;
                int minY=(y-1<=0)?0:y-1;
                int maxY=(y+1>=nRows-1)?nRows-1:y+1;
                for(int m=MAX(minX,0);m<=MIN(maxX,nCols-1);m++) {
                    for(int n=minY;n<=maxY;n++)
                    {
                        if(mask[n+m*nRows]==1&&label[m*nRows+n]==L)
                        {
                            mask[n+m*nRows]=0;
                            S.xLoc.insert(S.xLoc.end(),m);
                            S.yLoc.insert(S.yLoc.end(),n);
                        }
                        else if(label[m*nRows+n]!=L)
                        {
                            int NewLabel=label[m*nRows+n];
                            Pointer=find(S.Neighbor.begin(),S.Neighbor.end(),NewLabel);
                            if(Pointer==S.Neighbor.end())
                            {
                                S.Neighbor.insert(S.Neighbor.begin(),NewLabel);
                            }
                        }
                    }
                }
                if (minX<0) {
                    int m = nCols-1;
                    for(int n=minY;n<=maxY;n++)   {
                        if(mask[n+m*nRows]==1&&label[m*nRows+n]==L)
                        {
                            mask[n+m*nRows]=0;
                            S.xLoc.insert(S.xLoc.end(),m);
                            S.yLoc.insert(S.yLoc.end(),n);
                        }
                        else if(label[m*nRows+n]!=L)
                        {
                            int NewLabel=label[m*nRows+n];
                            Pointer=find(S.Neighbor.begin(),S.Neighbor.end(),NewLabel);
                            if(Pointer==S.Neighbor.end())
                            {
                                S.Neighbor.insert(S.Neighbor.begin(),NewLabel);
                            }
                        }
                    }
                }
                if (maxX>=nCols) {
                    int m = 0;
                    for(int n=minY;n<=maxY;n++)   {
                        if(mask[n+m*nRows]==1&&label[m*nRows+n]==L)
                        {
                            mask[n+m*nRows]=0;
                            S.xLoc.insert(S.xLoc.end(),m);
                            S.yLoc.insert(S.yLoc.end(),n);
                        }
                        else if(label[m*nRows+n]!=L)
                        {
                            int NewLabel=label[m*nRows+n];
                            Pointer=find(S.Neighbor.begin(),S.Neighbor.end(),NewLabel);
                            if(Pointer==S.Neighbor.end())
                            {
                                S.Neighbor.insert(S.Neighbor.begin(),NewLabel);
                            }
                        }
                    }
                }
                
            }
            Sarray.insert(Sarray.end(),S);
        }
    }
    
    if (1) {
        vector<Superpixel>::iterator S;
        vector<int>::iterator I;
        vector<int>::iterator I2;
        S=Sarray.begin();
        while(S!=Sarray.end())
        {
            float MinDist=DBL_MAX;
            int Label1=(*S).Label;
            int Label2=-1;
            for(I=(*S).Neighbor.begin();I!=(*S).Neighbor.end();I++)
            {
                float D=(centerL1[Label1]-centerL1[*I])*(centerL1[Label1]-centerL1[*I])+
                        (centerL2[Label1]-centerL2[*I])*(centerL2[Label1]-centerL2[*I])+
                        (centera1[Label1]-centera1[*I])*(centera1[Label1]-centera1[*I])+
                        (centera2[Label1]-centera2[*I])*(centera2[Label1]-centera2[*I])+
                        (centerb1[Label1]-centerb1[*I])*(centerb1[Label1]-centerb1[*I])+
                        (centerb2[Label1]-centerb2[*I])*(centerb2[Label1]-centerb2[*I])+
                        (centerx1[Label1]-centerx1[*I])*(centerx1[Label1]-centerx1[*I])+
                        (centerx2[Label1]-centerx2[*I])*(centerx2[Label1]-centerx2[*I])+
                        (centery1[Label1]-centery1[*I])*(centery1[Label1]-centery1[*I])+
                        (centery2[Label1]-centery2[*I])*(centery2[Label1]-centery2[*I]);
                if(D<MinDist)
                {
                    MinDist=D;
                    Label2=(*I);
                }
            }
            float W1=centerW[Label1];
            float W2=centerW[Label2];
            float W=W1+W2;
            centerL1[Label2]=(W2*centerL1[Label2]+W1*centerL1[Label1])/W;
            centerL2[Label2]=(W2*centerL2[Label2]+W1*centerL2[Label1])/W;
            centera1[Label2]=(W2*centera1[Label2]+W1*centera1[Label1])/W;
            centera2[Label2]=(W2*centera2[Label2]+W1*centera2[Label1])/W;
            centerb1[Label2]=(W2*centerb1[Label2]+W1*centerb1[Label1])/W;
            centerb2[Label2]=(W2*centerb2[Label2]+W1*centerb2[Label1])/W;
            centerx1[Label2]=(W2*centerx1[Label2]+W1*centerx1[Label1])/W;
            centerx2[Label2]=(W2*centerx2[Label2]+W1*centerx2[Label1])/W;
            centery1[Label2]=(W2*centery1[Label2]+W1*centery1[Label1])/W;
            centery2[Label2]=(W2*centery2[Label2]+W1*centery2[Label1])/W;
            centerW[Label2]=W;
            for(int i=0;i<(int)(*S).xLoc.size();i++)
            {
                int x=(*S).xLoc[i];int y=(*S).yLoc[i];
                label[x*nRows+y]=Label2;
            }
            vector<Superpixel>::iterator Stmp;
            Stmp=find(Sarray.begin(),Sarray.end(),Label2);
            if(Stmp!=Sarray.end())
            {
                Size[Label2]=Size[Label1]+Size[Label2];
                if(Size[Label2]>=threshold)
                {
                    Sarray.erase(Stmp);
                    Sarray.erase(S);
                }
                else
                {
                    (*Stmp).xLoc.insert((*Stmp).xLoc.end(),(*S).xLoc.begin(),(*S).xLoc.end());
                    (*Stmp).yLoc.insert((*Stmp).yLoc.end(),(*S).yLoc.begin(),(*S).yLoc.end());
                    (*Stmp).Neighbor.insert((*Stmp).Neighbor.end(),(*S).Neighbor.begin(),(*S).Neighbor.end());
                    sort((*Stmp).Neighbor.begin(),(*Stmp).Neighbor.end());
                    I=unique((*Stmp).Neighbor.begin(),(*Stmp).Neighbor.end());
                    (*Stmp).Neighbor.erase(I,(*Stmp).Neighbor.end());
                    I=find((*Stmp).Neighbor.begin(),(*Stmp).Neighbor.end(),Label1);
                    (*Stmp).Neighbor.erase(I);
                    I=find((*Stmp).Neighbor.begin(),(*Stmp).Neighbor.end(),Label2);
                    (*Stmp).Neighbor.erase(I);
                    Sarray.erase(S);
                }
            }
            else
            {
                Sarray.erase(S);
            }
            for(int i=0;i<(int)Sarray.size();i++)
            {
                I=find(Sarray[i].Neighbor.begin(),Sarray[i].Neighbor.end(),Label1);
                I2=find(Sarray[i].Neighbor.begin(),Sarray[i].Neighbor.end(),Label2);
                if(I!=Sarray[i].Neighbor.end()&&I2!=Sarray[i].Neighbor.end())
                    Sarray[i].Neighbor.erase(I);
                else if(I!=Sarray[i].Neighbor.end()&&I2==Sarray[i].Neighbor.end())
                    (*I)=Label2;
            }
            S=Sarray.begin();
        }
    }
    
    delete [] mask;
    return;
}



