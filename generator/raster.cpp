// Minimal z-buffer renderer. Input: int32 nf, float[nf*9] tris, uint8[nf*4] rgba,
// float[16] view (row-major camera: right,up,fwd,eye + ortho scale), int W,H. Output PPM.
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cstdint>
#include <algorithm>
int main(int argc,char**argv){
  FILE*fp=fopen(argv[1],"rb"); int nf; if(fread(&nf,4,1,fp)!=1) return 1;
  std::vector<float> T(nf*9); if(fread(T.data(),4,nf*9,fp)!=(size_t)nf*9) return 1;
  std::vector<uint8_t> C(nf*4); if(fread(C.data(),1,nf*4,fp)!=(size_t)nf*4) return 1;
  float cam[16]; if(fread(cam,4,16,fp)!=16) return 1; int WH[2]; if(fread(WH,4,2,fp)!=2) return 1; fclose(fp);
  int W=WH[0],H=WH[1];
  float rx=cam[0],ry=cam[1],rz=cam[2], ux=cam[3],uy=cam[4],uz=cam[5], fx=cam[6],fy=cam[7],fz=cam[8], ex=cam[9],ey=cam[10],ez=cam[11], sc=cam[12];
  float bg[3]={cam[13],cam[14],cam[15]};
  std::vector<float> zb(W*H,1e30f); std::vector<float> img(W*H*3); for(int i=0;i<W*H;i++){img[3*i]=bg[0];img[3*i+1]=bg[1];img[3*i+2]=bg[2];}
  std::vector<float> nb(W*H*3,0.f);
  float L1[3]={-0.45f,0.55f,0.70f}, L2[3]={0.6f,-0.3f,0.3f};
  auto nrm=[](float*v){float l=std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); v[0]/=l;v[1]/=l;v[2]/=l;}; nrm(L1); nrm(L2);
  for(int f=0;f<nf;f++){
    float*t=&T[9*f]; float sx[3],sy[3],sz[3];
    for(int k=0;k<3;k++){ float dx=t[3*k]-ex,dy=t[3*k+1]-ey,dz=t[3*k+2]-ez;
      sx[k]=(dx*rx+dy*ry+dz*rz)*sc+W/2.f; sy[k]=H/2.f-(dx*ux+dy*uy+dz*uz)*sc; sz[k]=dx*fx+dy*fy+dz*fz; }
    float e1[3]={t[3]-t[0],t[4]-t[1],t[5]-t[2]}, e2[3]={t[6]-t[0],t[7]-t[1],t[8]-t[2]};
    float n[3]={e1[1]*e2[2]-e1[2]*e2[1],e1[2]*e2[0]-e1[0]*e2[2],e1[0]*e2[1]-e1[1]*e2[0]};
    float ln=std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]); if(ln<1e-12f) continue; n[0]/=ln;n[1]/=ln;n[2]/=ln;
    float facing = n[0]*fx+n[1]*fy+n[2]*fz; // <0 faces camera
    float wn[3]={n[0],n[1],n[2]}; if(facing>0){wn[0]=-wn[0];wn[1]=-wn[1];wn[2]=-wn[2];}
    float d1=std::max(0.f,wn[0]*L1[0]+wn[1]*L1[1]+wn[2]*L1[2]); float d2=std::max(0.f,wn[0]*L2[0]+wn[1]*L2[1]+wn[2]*L2[2]);
    float vd=std::fabs(facing);
    float shade=0.42f+0.50f*d1+0.18f*d2+0.12f*vd; if(facing>0) shade*=0.6f;
    float a=C[4*f+3]/255.f;
    float col[3]={C[4*f]/255.f*shade,C[4*f+1]/255.f*shade,C[4*f+2]/255.f*shade};
    int x0=std::max(0,(int)std::floor(std::min({sx[0],sx[1],sx[2]}))), x1=std::min(W-1,(int)std::ceil(std::max({sx[0],sx[1],sx[2]})));
    int y0=std::max(0,(int)std::floor(std::min({sy[0],sy[1],sy[2]}))), y1=std::min(H-1,(int)std::ceil(std::max({sy[0],sy[1],sy[2]})));
    float area=(sx[1]-sx[0])*(sy[2]-sy[0])-(sx[2]-sx[0])*(sy[1]-sy[0]); if(std::fabs(area)<1e-9f) continue;
    for(int y=y0;y<=y1;y++) for(int x=x0;x<=x1;x++){
      float px=x+0.5f,py=y+0.5f;
      float w0=((sx[1]-px)*(sy[2]-py)-(sx[2]-px)*(sy[1]-py))/area;
      float w1=((sx[2]-px)*(sy[0]-py)-(sx[0]-px)*(sy[2]-py))/area;
      float w2=1-w0-w1; if(w0<0||w1<0||w2<0) continue;
      float z=w0*sz[0]+w1*sz[1]+w2*sz[2]; int i=y*W+x;
      if(a<0.99f){ // transparent: blend, no depth write
        if(z<zb[i]){ for(int c=0;c<3;c++) img[3*i+c]=img[3*i+c]*(1-a)+col[c]*a; } continue; }
      if(z<zb[i]){ zb[i]=z; img[3*i]=col[0];img[3*i+1]=col[1];img[3*i+2]=col[2]; nb[3*i]=wn[0];nb[3*i+1]=wn[1];nb[3*i+2]=wn[2]; }
    }
  }
  // silhouette / crease lines
  std::vector<float> out=img;
  for(int y=1;y<H-1;y++) for(int x=1;x<W-1;x++){ int i=y*W+x; float z=zb[i]; if(z>1e29f) continue;
    float m=0; int nbrs[4]={i-1,i+1,i-W,i+W};
    for(int k=0;k<4;k++){ int j=nbrs[k]; if(zb[j]>1e29f){m=1;break;} float dz=std::fabs(zb[j]-z); if(dz>1.5f/sc*4) m=std::max(m,0.8f);
      float dn=nb[3*i]*nb[3*j]+nb[3*i+1]*nb[3*j+1]+nb[3*i+2]*nb[3*j+2]; if(dn<0.5f) m=std::max(m,0.5f);}
    for(int c=0;c<3;c++) out[3*i+c]=out[3*i+c]*(1-0.75f*m); }
  FILE*o=fopen(argv[2],"wb"); fprintf(o,"P6\n%d %d\n255\n",W,H);
  for(int i=0;i<W*H*3;i++){ int v=(int)std::lround(std::min(1.f,std::max(0.f,out[i]))*255); fputc(v,o);} fclose(o); return 0; }
