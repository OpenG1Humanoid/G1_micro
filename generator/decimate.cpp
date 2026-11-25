// Quadric-error-metric mesh decimation (Garland & Heckbert) for closed manifold meshes.
// Binary I/O: int32 nv, int32 nf, double[nv*3], int32[nf*3]
// usage: decimate in.bin out.bin max_error target_faces
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <queue>
#include <algorithm>
#include <cstdint>

struct Q { double a[10]; Q(){ for(int i=0;i<10;i++) a[i]=0; }
  void add(const Q&o){ for(int i=0;i<10;i++) a[i]+=o.a[i]; }
  void plane(double nx,double ny,double nz,double d,double w){
    a[0]+=w*nx*nx; a[1]+=w*nx*ny; a[2]+=w*nx*nz; a[3]+=w*nx*d;
    a[4]+=w*ny*ny; a[5]+=w*ny*nz; a[6]+=w*ny*d;
    a[7]+=w*nz*nz; a[8]+=w*nz*d; a[9]+=w*d*d; }
  double err(double x,double y,double z) const {
    return a[0]*x*x+2*a[1]*x*y+2*a[2]*x*z+2*a[3]*x+a[4]*y*y+2*a[5]*y*z+2*a[6]*y+a[7]*z*z+2*a[8]*z+a[9]; }
};

struct V3 { double x,y,z; };
static inline V3 sub(V3 a,V3 b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline V3 cross(V3 a,V3 b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double dot(V3 a,V3 b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline double len(V3 a){return std::sqrt(dot(a,a));}

struct Edge { double cost; int a,b; unsigned va,vb; V3 p; };
struct Cmp { bool operator()(const Edge&x,const Edge&y)const{return x.cost>y.cost;} };

std::vector<V3> P; std::vector<int> F; std::vector<char> fdead, vdead;
std::vector<std::vector<int>> VF; std::vector<Q> QV; std::vector<unsigned> ver;

static bool solve(const Q&q, V3&out){
  double m[3][4]={{q.a[0],q.a[1],q.a[2],-q.a[3]},{q.a[1],q.a[4],q.a[5],-q.a[6]},{q.a[2],q.a[5],q.a[7],-q.a[8]}};
  double det = m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1]) - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0]) + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
  if (std::fabs(det) < 1e-10) return false;
  auto D3=[&](int c){ double t[3][3]; for(int i=0;i<3;i++)for(int j=0;j<3;j++) t[i][j]= (j==c)? m[i][3]: m[i][j];
    return t[0][0]*(t[1][1]*t[2][2]-t[1][2]*t[2][1]) - t[0][1]*(t[1][0]*t[2][2]-t[1][2]*t[2][0]) + t[0][2]*(t[1][0]*t[2][1]-t[1][1]*t[2][0]); };
  out = {D3(0)/det, D3(1)/det, D3(2)/det};
  return true;
}

static Edge make_edge(int a,int b){
  Q q=QV[a]; q.add(QV[b]);
  V3 best; double bc;
  V3 cand[4]; int nc=0;
  V3 s; if (solve(q,s)) {
    // keep optimum only if close to the edge (avoid spikes)
    V3 mid={(P[a].x+P[b].x)/2,(P[a].y+P[b].y)/2,(P[a].z+P[b].z)/2};
    if (len(sub(s,mid)) < 2.0*len(sub(P[a],P[b]))+1e-9) cand[nc++]=s; }
  cand[nc++]=P[a]; cand[nc++]=P[b];
  cand[nc++]={(P[a].x+P[b].x)/2,(P[a].y+P[b].y)/2,(P[a].z+P[b].z)/2};
  bc=1e300; best=cand[0];
  for(int i=0;i<nc;i++){ double e=q.err(cand[i].x,cand[i].y,cand[i].z); if(e<bc){bc=e;best=cand[i];} }
  Edge E; E.cost=std::max(bc,0.0); E.a=a; E.b=b; E.va=ver[a]; E.vb=ver[b]; E.p=best; return E;
}

static void neighbours(int v, std::vector<int>&out){
  out.clear();
  for(int f: VF[v]) if(!fdead[f]) for(int k=0;k<3;k++){ int w=F[3*f+k]; if(w!=v) out.push_back(w);}
  std::sort(out.begin(),out.end()); out.erase(std::unique(out.begin(),out.end()),out.end());
}

int main(int argc,char**argv){
  if(argc<5){fprintf(stderr,"usage\n");return 1;}
  double maxerr=atof(argv[3]); long target=atol(argv[4]);
  FILE*fp=fopen(argv[1],"rb"); int nv,nf; fread(&nv,4,1,fp); fread(&nf,4,1,fp);
  P.resize(nv); std::vector<double> buf(nv*3); fread(buf.data(),8,nv*3,fp);
  for(int i=0;i<nv;i++) P[i]={buf[3*i],buf[3*i+1],buf[3*i+2]};
  F.resize(nf*3); fread(F.data(),4,nf*3,fp); fclose(fp);
  fdead.assign(nf,0); vdead.assign(nv,0); VF.assign(nv,{}); QV.assign(nv,Q()); ver.assign(nv,0);
  for(int f=0;f<nf;f++){
    int a=F[3*f],b=F[3*f+1],c=F[3*f+2];
    VF[a].push_back(f); VF[b].push_back(f); VF[c].push_back(f);
    V3 n=cross(sub(P[b],P[a]),sub(P[c],P[a])); double l=len(n); if(l<1e-14) continue;
    double area=0.5*l; n={n.x/l,n.y/l,n.z/l}; double d=-dot(n,P[a]);
    Q q; q.plane(n.x,n.y,n.z,d,area); QV[a].add(q); QV[b].add(q); QV[c].add(q);
  }
  std::priority_queue<Edge,std::vector<Edge>,Cmp> heap;
  for(int f=0;f<nf;f++) for(int k=0;k<3;k++){ int a=F[3*f+k], b=F[3*f+(k+1)%3]; if(a<b) heap.push(make_edge(a,b)); }
  long alive=nf; std::vector<int> na,nb,common;
  while(!heap.empty() && alive>target){
    Edge e=heap.top(); heap.pop();
    if(vdead[e.a]||vdead[e.b]||ver[e.a]!=e.va||ver[e.b]!=e.vb) continue;
    if(e.cost>maxerr) break;
    int a=e.a,b=e.b;
    // link condition
    neighbours(a,na); neighbours(b,nb); common.clear();
    std::set_intersection(na.begin(),na.end(),nb.begin(),nb.end(),std::back_inserter(common));
    int shared=0; for(int f: VF[a]) if(!fdead[f]){ int*t=&F[3*f]; if(t[0]==b||t[1]==b||t[2]==b) shared++; }
    if(shared!=2 || common.size()!=2) continue;
    // flip / sliver check
    bool bad=false;
    for(int pass=0;pass<2 && !bad;pass++){ int v= pass? b:a;
      for(int f: VF[v]){ if(fdead[f]) continue; int*t=&F[3*f]; if(t[0]==(pass?a:b)||t[1]==(pass?a:b)||t[2]==(pass?a:b)) continue;
        V3 p0=P[t[0]],p1=P[t[1]],p2=P[t[2]]; V3 n0=cross(sub(p1,p0),sub(p2,p0));
        V3 q0=p0,q1=p1,q2=p2; if(t[0]==v) q0=e.p; if(t[1]==v) q1=e.p; if(t[2]==v) q2=e.p;
        V3 n1=cross(sub(q1,q0),sub(q2,q0)); double l0=len(n0), l1=len(n1);
        if(l1<1e-12){bad=true;break;}
        if(dot(n0,n1) < 0.3*l0*l1){bad=true;break;}
        // reject very skinny results
        double e0=len(sub(q1,q0)),e1=len(sub(q2,q1)),e2=len(sub(q0,q2)); double mx=std::max(e0,std::max(e1,e2));
        if(l1/(mx*mx) < 0.02){bad=true;break;}
      } }
    if(bad) continue;
    // collapse b into a
    P[a]=e.p; QV[a].add(QV[b]); vdead[b]=1;
    for(int f: VF[b]){ if(fdead[f]) continue; int*t=&F[3*f];
      if(t[0]==a||t[1]==a||t[2]==a){ fdead[f]=1; alive--; continue; }
      for(int k=0;k<3;k++) if(t[k]==b) t[k]=a; VF[a].push_back(f); }
    std::vector<int> keep; for(int f: VF[a]) if(!fdead[f]) keep.push_back(f);
    std::sort(keep.begin(),keep.end()); keep.erase(std::unique(keep.begin(),keep.end()),keep.end()); VF[a]=keep;
    ver[a]++;
    neighbours(a,na); for(int w: na){ ver[w]++; }
    for(int w: na){ heap.push(make_edge(std::min(a,w),std::max(a,w)));
      // refresh edges around neighbours (their version changed)
      std::vector<int> nw; neighbours(w,nw); for(int u: nw) if(u!=a) heap.push(make_edge(std::min(w,u),std::max(w,u))); }
  }
  std::vector<int> remap(nv,-1); int cnt=0;
  for(int f=0;f<nf;f++) if(!fdead[f]) for(int k=0;k<3;k++){ int v=F[3*f+k]; if(remap[v]<0) remap[v]=cnt++; }
  std::vector<double> ov(cnt*3); for(int v=0;v<nv;v++) if(remap[v]>=0){ ov[3*remap[v]]=P[v].x; ov[3*remap[v]+1]=P[v].y; ov[3*remap[v]+2]=P[v].z; }
  std::vector<int> of; for(int f=0;f<nf;f++) if(!fdead[f]) for(int k=0;k<3;k++) of.push_back(remap[F[3*f+k]]);
  int onf=of.size()/3; fp=fopen(argv[2],"wb"); fwrite(&cnt,4,1,fp); fwrite(&onf,4,1,fp); fwrite(ov.data(),8,cnt*3,fp); fwrite(of.data(),4,onf*3,fp); fclose(fp);
  return 0;
}
