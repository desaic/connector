#include "connector.hpp"
#include <fstream>
#include <iostream>
#include <set>
#include <algorithm>
#include <cmath>
using ClipperLib::long64;
using ClipperLib::Clipper;
using ClipperLib::Polygons;
using ClipperLib::Polygon;
using ClipperLib::IntPoint;
#define PI 3.1415926
real_t teethRatio=1;
real_t min_depth_ratio=0.2;

/**@brief
offset a little so that tv0 is not exactly on the line
*/
real_t normalOffset=10;
real_t unitRatio=2;
real_t slotUnit=1;
real_t startRatio=2;
real_t testExtend=1;
real_t reserveRatio=0.07;
real_t extraRatio=0.06;
#define CLIP_VAL(x,y) ((x)<(y)?(y):(x))
void make_rect(Vec3 &start, const Vec3 & lineDir, const Vec3 &normal,
               real_t t, real_t len,std::vector<Vert> & rect);
void vert2poly(const std::vector<Vert> & rect, Polygon & p);
void poly2vert(const Polygon & polySeg, std::vector<Vert>&lineseg);
void fix_dir(Vec3 & dir, const Vec3 & lineDir, const Vec3 lineNormal);
void PolyMesh::chopAlongEdge(const Edge&e,const EdgeVal & ev, int ii,real_t depth, Polygon & rect)
{
  int planeIdx = ev.p[ii];
  VertIdx vi0 = vertp[e.id[0]][planeIdx];
  VertIdx vi1 = vertp[e.id[1]][planeIdx];
  if(ii==1){
    vi0= vertp[e.id[1]][planeIdx];
    vi1= vertp[e.id[0]][planeIdx];
  }
  std::vector<Vert> & lineseg =  planes[vi0.i][vi0.j];
  Vec3 v0  = lineseg[vi0.k].v;
  Vec3 v1  = lineseg[vi1.k].v;
  //vertex before v0
  Vec3 v_1 = vi0.k==0 ?
             lineseg[lineseg.size()-1].v : lineseg[vi0.k-1].v;
  //vertex after v1
  Vec3 v2  = lineseg[(vi1.k+1)%lineseg.size()].v;

  Vec3 lineDir=v1-v0;
  real_t len = lineDir.norm();
  lineDir/=len;
  //(a,b)-->(b,-a) is perpendicular on the right of the edge
  Vec3 lineNormal = Vec3(lineDir[1],-lineDir[0],0);

  Vec3 dir0 = v_1-v0;
  dir0/=dir0.norm();
  fix_dir(dir0,lineDir,lineNormal);
  real_t depth0=depth/(-dir0.dot(lineNormal));

  Vec3 dir1 = v2-v1;
  dir1/=dir1.norm();
  fix_dir(dir1,-lineDir, lineNormal);
  real_t depth1=depth/(-dir1.dot(lineNormal));

  real_t alpha = -0.001;
  Vec3 tv0=(1-alpha)*v0+alpha*v1;
  Vec3 tv01=tv0+depth0*dir0;
  //offset a little so that tv0 is not exactly on the line
  tv0 += lineNormal*normalOffset;
  rect.push_back(IntPoint((long64)tv0[0],(long64)tv0[1]));
  rect.push_back(IntPoint((long64)tv01[0],(long64)tv01[1]));

  alpha = 1.001;
  Vec3 tv1=(1-alpha)*v0+alpha*v1;
  Vec3 tv11=tv1+depth1*dir1;
  tv1+=lineNormal*normalOffset;
  rect.push_back(IntPoint((long64)tv11[0],(long64)tv11[1]));
  rect.push_back(IntPoint((long64)tv1[0],(long64)tv1[1]));
}

/**@brief 1.chop away both sides
2. add teeth
*/
void PolyMesh::teeth()
{
  std::map<Edge,EdgeVal>::iterator it;
  real_t teethWidth = teethRatio*t;
  for(it=eset.begin(); it!=eset.end(); it++) {
    if(it->second.hasConn){
      continue;
    }
    int pid[2]= {it->second.p[0],it->second.p[1]};
    Vec3 n1=planes[pid[0]].n, n2=planes[pid[1]].n;
    real_t cosine=n1.dot(n2);
    real_t chop=0;
    if(cosine>0) {
      real_t halftan=half_tan(cosine);
      chop=t*halftan*(cosine+0.5);
    } else {
      real_t halftan=half_tan(-cosine);
      halftan=CLIP_VAL(halftan, min_depth_ratio);
      chop=(t/2)/halftan;
    }

    for(int ii=0; ii<2; ii++) {
      Polygon rect;
      chopAlongEdge(it->first, it->second, ii,chop,rect);
      chopPoly(rect, pid[ii]);
    }
  }

  for(it=eset.begin(); it!=eset.end(); it++) {
    if(it->second.hasConn){
      continue;
    }
    int pid[2]= {it->second.p[0],it->second.p[1]};
    Vec3 n1=planes[pid[0]].n, n2=planes[pid[1]].n;
    real_t cosine=n1.dot(n2);
    real_t chop=0,extend=0;
    if(cosine>0) {
      real_t halftan=half_tan(cosine);
      chop=t*halftan*(cosine+0.5);
      extend=t*halftan*(1+cosine);
    } else {
      real_t halftan=half_tan(-cosine);
      halftan=CLIP_VAL(halftan, min_depth_ratio);
      chop=(t/2)/halftan;
      extend=chop;
    }

    for(int ii=0; ii<2; ii++) {
      Vec3 v0,v1;
      edgeVertPos(it->first, pid[ii],v0,v1);
      Vec3 lineDir = v1-v0;
      real_t len = lineDir.norm();
      lineDir/=len;
      Vec3 lineNormal=Vec3(lineDir[1],-lineDir[0],0);
      if(ii!=0){
        lineNormal=-lineNormal;
      }
      real_t start=0;
      int nteeth=0;
      if(ii!=0){
        nteeth=1;
      }
      std::vector<Vert>lineseg;
      poly2vert(poly[pid[ii]][0],lineseg);
      while(start+teethWidth<len){
        start = nteeth*teethWidth;
        Vec3 r0=v0+start*lineDir-chop*lineNormal;
        r0+=normalOffset*lineNormal;
        std::vector<Vert>rect;
        make_rect(r0,lineDir,-lineNormal,teethWidth,extend,rect);

        bool valid=true;
        size_t jj0=lineseg.size()-1;
        for(size_t jj=0;jj<lineseg.size();jj++){
          if(  lineIntersect(lineseg[jj0].v,lineseg[jj].v,rect[0].v,rect[1].v)
             ||lineIntersect(lineseg[jj0].v,lineseg[jj].v,rect[2].v,rect[3].v)){
            valid=false;
            goto ENDTEETHCHECK;
          }
          jj0=jj;
        }
        r0 -= 2*normalOffset*lineNormal;
        rect.clear();
        make_rect(r0,lineDir,-lineNormal,teethWidth,extend,rect);
        if(   pnpoly(lineseg,rect[2])
           || pnpoly(lineseg,rect[1])
           ||!pnpoly(lineseg,rect[0])
           ||!pnpoly(lineseg,rect[3])){
            valid=false;
        }
      ENDTEETHCHECK:
        if(valid){
          Polygon p;
          vert2poly(rect,p);
          chopPoly(p,pid[ii],ClipperLib::ctUnion);
        }
        nteeth+=2;
      }
    }
  }
}

int spots=3;
real_t polyArea(std::vector<Vert> & l)
{
  int ii0=l.size()-1;
  real_t A=0;
  for(size_t ii=0; ii<l.size(); ii++) {
    A+= l[ii0].v[0]*l[ii].v[1]-l[ii].v[0]*l[ii0].v[1];
    ii0=ii;
  }
  return A;
}

void PolyMesh::save_result(const char * filename)
{
  std::ofstream out;
  out.open(filename);
  if(!out.good()) {
    std::cout<<"cannot open file"<<filename<<"\n";
    return;
  }
  out<<planes.size()+c.size()<<"\n";
  for(size_t ii=0; ii<poly.size(); ii++) {
    out<<poly[ii].size()<<"\n";
    if(poly[ii].size()<1) {
      out<<"\n";
      continue;
    }
    for(size_t jj=0; jj<poly[ii].size(); jj++) {
      out<<poly[ii][jj].size()<<"\n";
      for(size_t kk=0; kk<poly[ii][jj].size(); kk++) {
        IntPoint point = poly[ii][jj][kk];
        real_t x = (real_t)point.X;
        real_t y = (real_t)point.Y;
        x=(x/intscale)*obj_scale;
        y=(y/intscale)*obj_scale;

        out<<x<<","<<y<<" ";
      }
      out<<"\n";
    }
    out<<"\n";
  }

  for(size_t ii=0; ii<c.size(); ii++) {
    out<<"1 "<<c[ii].pid[0]<<" "<<c[ii].pid[1]<<"\n";
    out<<c[ii].l.size()<<"\n";
    for(size_t jj=0; jj<c[ii].l.size(); jj++) {
      real_t x = c[ii].l[jj][0];
      real_t y = c[ii].l[jj][1];
      x=(x/intscale)*obj_scale;
      y=(y/intscale)*obj_scale;
      out<<x<<","<<y<<" ";
    }
    out<<"\n\n";
  }

  out.close();
}


PolyMesh::PolyMesh(const char * filename):intscale(1),obj_scale(1),
  t(0) {
  std::ifstream in;
  in.open(filename);
  if(!in.good()) {
    std::cout<<"cannot open"<<filename<<"\n";
    return;
  }
  size_t nplane=0;
  in>>nplane;


  int nvert=0;
  int pcnt=0;
  for(size_t ii=0; ii<nplane; ii++) {
    int nseg = 0 ;
    in>>nseg;
    if(nseg <= 0) {
      continue;
    }

    planes.push_back(Plane());
    Plane &p=planes[pcnt];
    pcnt++;
    in>> p.n[0];
    in>> p.n[1];
    in>> p.n[2];
    in>> p.ax[0];
    in>> p.ax[1];
    in>> p.ax[2];
    in>> p.ay[0];
    in>> p.ay[1];
    in>> p.ay[2];
    p.l.resize(nseg);
    for(size_t jj=0; jj<p.l.size(); jj++) {
      //number of points in a segment
      int npt = 0;
      in >> npt;
      p[jj].resize(npt);
      for(size_t kk=0; kk<p[jj].size(); kk++) {
        in >> p[jj][kk].v[0];
        char c;
        //comma;
        in >>c;
        in >> p[jj][kk].v[1];
      }
      for(size_t kk=0; kk<p[jj].size(); kk++) {
        int id;
        in >> id;
        if(vid.find(id)==vid.end()) {
          vid[id]=nvert;
          nvert++;
        }
        id=vid[id];
        p[jj][kk].id=id;
      }
    }
  }

  //wind polygons in counter-clock-wise order
  //inner line segments have clock-wise order
  //area is always on left of line segment
  //assume planes have no self-intersections
  //also it cannot have anything else inside a hole
  //first find a line segment that's the contour of the plane
  //then deal with holes

  //for line segments i and i-1
  //if its point is inside the other polygon then that polygon is the most outside one
  //otherwise if there are an odd number of polygons the last one is the outmost
  for(size_t ii=0; ii<planes.size(); ii++) {
    if(planes[ii].size()<2) {
      continue;
    }
    //to deal with odd number of line segments
    int outmost = planes[ii].size()-1;
    for(size_t jj=1; jj<planes[ii].size(); jj+=2) {
      //jj is inside jj-1
      if(pnpoly(planes[ii][jj-1],planes[ii][jj][0])) {
        outmost=jj-1;
        break;
      }//jj-1 is inside jj
      else if(pnpoly(planes[ii][jj],planes[ii][jj-1][0])) {
        outmost=jj;
        break;
      }
    }
    if(outmost!=0) {
      //swap the outmost line segment to the first segment
      std::vector<Vert>  tmp=planes[ii][0];
      planes[ii][0]=planes[ii][outmost];
      planes[ii][outmost]=tmp;
    }
  }

  //now fix winding orders
  for(size_t ii=0; ii<planes.size(); ii++) {

    real_t A=polyArea(planes[ii][0]);
    if(A<0) {
      std::reverse(planes[ii][0].begin(),planes[ii][0].end());
    }
    for(size_t jj=1; jj<planes[ii].size(); jj++) {
      A=polyArea(planes[ii][jj]);
      if(A>0) {
        std::reverse(planes[ii][jj].begin(),planes[ii][jj].end());
      }
    }
  }
  in.close();
};

void PolyMesh::scale(real_t s)
{
  intscale=s;
  for(size_t ii=0; ii<planes.size(); ii++) {
    for(size_t jj=0; jj<planes[ii].size(); jj++) {
      for(size_t kk=0; kk<planes[ii][jj].size(); kk++) {
        planes[ii][jj][kk].v*=s;
      }
    }
  }
}

void PolyMesh::buildEdge() {
  vertp.resize(vid.size());
  for(size_t ii=0; ii<planes.size(); ii++) {
    for(size_t jj=0; jj<planes[ii].size(); jj++) {
      size_t kk1=planes[ii][jj].size()-1;
      for(size_t kk=0; kk<planes[ii][jj].size(); kk++) {
        int v0=planes[ii][jj][kk].id;
        int v1=planes[ii][jj][kk1].id;
        vertp[v0][ii]=VertIdx(ii,jj,kk);
        Edge e(v0,v1);
        std::map<Edge,EdgeVal>::iterator it=eset.find(e);
        if(it!=eset.end()) {
          eset[e].p[1]=ii;
        }
        else {
          eset[e].p[0]=ii;
          eset[e].p[1]=-1;
        }
        kk1=kk;
      }
    }
  }

  std::map<Edge,EdgeVal>::iterator it=eset.begin();
  std::map<Edge,EdgeVal>::iterator nextit;

  for(; it!=eset.end(); it=nextit) {
    nextit=it;
    nextit++;
    if(it->second.p[1]<0 || it->second.p[0]==it->second.p[1]) {
      eset.erase(it);
    }
    int pid[2]= {it->second.p[0],it->second.p[1]};
    int v0pIdx=vertp[it->first.id[0]][pid[0]].k;
    int v1pIdx=vertp[it->first.id[1]][pid[0]].k;
    int npt = planes[pid[0]][vertp[it->first.id[0]][pid[0]].j].size();
    if(v1pIdx!=(v0pIdx+1)%npt) {
      it->second.p[0]=pid[1];
      it->second.p[1]=pid[0];
    }
  }
}

void make_rect(Vec3 &start, const Vec3 & lineDir, const Vec3 &normal,
               real_t t, real_t len,std::vector<Vert> & rect)
{
  Vec3 v = start;
  //id is not used here
  rect.push_back(Vert(v));

  v -= normal*len;
  rect.push_back(Vert(v));

  v += lineDir*t;
  rect.push_back(Vert(v));

  v += normal*len;
  rect.push_back(Vert(v));
}

void PolyMesh::edgeVertPos(const Edge &e, int planeIdx, Vec3 &v0,Vec3 & v1)
{
  VertIdx & vi0 = vertp[e.id[0]][planeIdx];
  VertIdx & vi1 = vertp[e.id[1]][planeIdx];
  std::vector<Vert>&lineseg =  planes[vi0.i][vi0.j];
  v0  = lineseg[vi0.k].v;
  v1  = lineseg[vi1.k].v;
}

void PolyMesh::slot(real_t frac)
{
  std::map<Edge,EdgeVal>::iterator it;
  real_t unitlen = unitRatio*t*frac;
  real_t start = t*startRatio;
  real_t slot_len=unitlen*slotUnit;
  real_t testLen=testExtend*unitlen;
  //  real_t end = start+slot_len;
  for(it=eset.begin(); it!=eset.end(); it++) {
    if(it->second.hasConn) {
      continue;
    }
    int pid[2]= {it->second.p[0],it->second.p[1]};
    //if concave shift the slots outwards by (t/2) * or / tan(theta/2);
    //if convex shift inwards
    //positive value for shifting inwards
    real_t shift =0;
    Vec3 &n1=planes[pid[0]].n;
    Vec3 &n2 = planes[pid[1]].n;

    real_t cosine = n1.dot(n2);
    real_t halftangent = half_tan(std::abs(cosine));
    if(cosine>0) {
      //obtuse angle
      shift = (t/2)*halftangent;
    } else {
      //acute angle
      halftangent =CLIP_VAL(halftangent , min_depth_ratio);
      shift=(t/2)/halftangent;
    }

    if(!isConvex(it->first, it->second)) {
      //concave
      shift=-shift;
    }
    bool possible = true;
    real_t alpha=1.0/(spots+1);
    for(int spot=1; spot<=spots; spot++) {
      possible=true;
      alpha=spot/(spots+1.0);
      for(size_t ii=0; ii<2; ii++) {
        int planeIdx = pid[ii];
        Vec3 v0, v1;
        edgeVertPos(it->first, planeIdx, v0,v1);
        Vec3 lineDir=v1-v0;
        real_t len = lineDir.norm();
        lineDir/=len;
        Vec3 lineNormal = Vec3(lineDir[1],-lineDir[0],0);
        if(ii!=0){
          lineNormal=-lineNormal;
        }
        Vec3 mid=alpha*v0+(1-alpha)*v1;
        mid-=lineNormal*(shift+start);
        mid-=lineDir*(testLen+t)/2;
        std::vector<Vert>rect;

        make_rect(mid,lineDir,lineNormal,t+testLen,
                  slot_len+testLen, rect);
        size_t jj0=rect.size()-1;
        for(size_t jj=0; jj<rect.size(); jj++) {
          bool expected=true;
          for(size_t seg=0; seg<poly[planeIdx].size(); seg++) {
            Polygon polySeg = poly[planeIdx][seg];
            std::vector<Vert>lineseg;
            poly2vert(polySeg, lineseg);
            bool ret = pnpoly(lineseg, rect[jj]);
            if(ret!=expected) {
              possible = false;
              goto ENDEDGELOOP;
            }
            size_t kk0=lineseg.size()-1;
            for(size_t kk=0; kk<lineseg.size(); kk++) {
              ret = pnpoly(rect, lineseg[kk]);
              if(ret) {
                possible = false;
                goto ENDEDGELOOP;
              }
              bool intersect = lineIntersect(rect[jj].v,rect[jj0].v,
                                             lineseg[kk0].v,lineseg[kk].v);
              if(intersect) {
                possible=false;
                goto ENDEDGELOOP;
              }
              kk0=kk;
            }
            //slot should not cut into inner holes
            expected = false;
          }
          jj0=jj;
        }
      }//for two pieces
      if(possible) {
        break;
      }
ENDEDGELOOP:
      ;
    }//for each spots


    if(possible) {
      it->second.hasConn=true;
      it->second.connSize=slot_len;
      for(size_t ii=0; ii<2; ii++) {
        int planeIdx = pid[ii];
        Vec3 v0, v1;
        edgeVertPos(it->first, planeIdx, v0,v1);
        Vec3 mid=alpha*v0+(1-alpha)*v1;

        Vec3 lineDir=v1-v0;
        real_t len = lineDir.norm();
        lineDir/=len;
        Vec3 lineNormal = Vec3(lineDir[1],-lineDir[0],0);
        if(ii!=0){
          lineNormal=-lineNormal;
        }
        real_t reserveLen=t*reserveRatio;
        mid+=(-lineNormal*(start+reserveLen+shift));
        mid-=lineDir*(t)/2;
        std::vector<Vert> rect;
        make_rect(mid, lineDir, lineNormal,t,slot_len-2*reserveLen,rect);
        Polygon p ;
        vert2poly(rect,p);
        poly[planeIdx].push_back(p);
        planes[planeIdx].l.push_back(rect);
      }
    }
    std::cout<<"possible"<<possible<<"\n";
  }
}

Vec3 rotate(const Vec3 & v,real_t cosine)
{
  real_t sine= std::sqrt(1-cosine*cosine);
  Vec3 ret(cosine*v.get(0)+sine*v.get(1), -sine*v.get(0)+cosine*v.get(1));
  return ret;
}

void PolyMesh::connector()
{
  //for laser cutter
  //leave a bit extra material for friction fit
  real_t unitlen = unitRatio*t;
  real_t u = unitlen;
  std::map<Edge,EdgeVal>::iterator it;
  real_t extra=t*extraRatio;
  for(it = eset.begin(); it!=eset.end(); it++) {
    if( !it->second.hasConn) {
      continue;
    }
    std::vector<Vec3> baseShape(6);
    real_t slot_len=it->second.connSize;
    baseShape[0]=Vec3(0,0,0);
    baseShape[1]=Vec3(-u+extra/2,0,0);
    baseShape[2]=Vec3(-u,-t,0);
    baseShape[3]=Vec3(-u-slot_len,-t,0);
    baseShape[4]=Vec3(-u-extra/2-slot_len,0,0);
    baseShape[5]=Vec3(-u-extra/2-slot_len,u,0);

    int pid[2]= {it->second.p[0],it->second.p[1]};

    Connector conn;
    Vec3 &n1 = planes[pid[0]].n;
    Vec3 &n2 = planes[pid[1]].n;

    conn.l.insert(conn.l.end(),baseShape.begin(),baseShape.end());
    conn.pid[0]=pid[0];
    conn.pid[1]=pid[1];

    real_t cosine=n1.dot(n2);
    if(isConvex(it->first,it->second)) {
      //convex
      cosine = -cosine;
      std::vector<Vec3>rot(6);
      for(size_t ii=0; ii<rot.size(); ii++) {
        size_t ri=rot.size()-1-ii;
        rot[ri]=Vec3(baseShape[ii][0],-baseShape[ii][1],0);
        Vec3 tmp = rotate(rot[ri],cosine);
        rot[ri]=tmp;
      }
      conn.l.insert(conn.l.end(),rot.begin(),rot.end());
    } else {
      //concave
      std::vector<Vec3>rot(6);
      for(size_t ii=0; ii<rot.size(); ii++) {
        size_t ri=rot.size()-1-ii;
        rot[ri]=Vec3(-baseShape[ii][0],baseShape[ii][1],0);
        Vec3 tmp = rotate(rot[ri],cosine);
        rot[ri]=tmp;
      }
      if(cosine>0) {
        real_t halftan=half_tan(cosine);
        real_t extend=u*halftan;
        Vec3 exv = Vec3(extend,u,0);
        conn.l.push_back(exv);
        conn.l.insert(conn.l.end(),rot.begin(),rot.end());
      }
      else {
        real_t extend=u;
        Vec3 exv = Vec3(extend,u,0);
        conn.l.push_back(exv);
        conn.l.insert(conn.l.end(),rot.begin(),rot.end());
      }
    }
    c.push_back(conn);
  }
}

void fix_dir(Vec3 & dir, const Vec3 & lineDir, const Vec3 lineNormal)
{
  if(dir.dot(lineNormal)>=0) {
    dir=-lineNormal;
  }
  if(dir.dot(lineNormal)>-0.3) {
    dir=-lineNormal;
  }
  if(dir.dot(lineDir)>0) {
    dir=-lineNormal;
  }
}

void PolyMesh::zz(real_t _t)
{
  buildEdge();
  t=(_t/obj_scale)*intscale;
  std::map<Edge,EdgeVal>::iterator it;
  poly.resize(planes.size());

  for(size_t ii=0; ii<planes.size(); ii++) {
    poly[ii].resize(planes[ii].size());
    for(size_t jj=0; jj<planes[ii].size(); jj++) {
      vert2poly(planes[ii][jj],poly[ii][jj]);
    }
  }
  for(it=eset.begin(); it!=eset.end(); it++) {
    real_t depth=0;
    int pid[2]= {it->second.p[0],it->second.p[1]};
    Vec3 n1=planes[pid[0]].n, n2=planes[pid[1]].n;
    real_t cosine=n1.dot(n2);
    if(cosine>=0) {
      //obtuse angle
      depth = std::sqrt(1-cosine*cosine)*t/2;
    } else {
      //acute angle
      real_t halftangent = half_tan(std::abs(cosine));
      halftangent =CLIP_VAL(halftangent , min_depth_ratio);
      depth = (t/2)/halftangent;
    }

    if(depth<1) {
      continue;
    }
    for(size_t ii=0; ii<1; ii++) {
      Polygon rect;
      chopAlongEdge(it->first,it->second, ii,depth, rect);
      chopPoly(rect,pid[ii]);
    }
  }
}

bool PolyMesh::isConvex(const Edge & e, const EdgeVal&ev)
{
  int pid[2]= {ev.p[0],ev.p[1]};
  int planeIdx = pid[0];
  VertIdx& vi0 = vertp[e.id[0]][planeIdx];
  VertIdx& vi1 = vertp[e.id[1]][planeIdx];
  std::vector<Vert> & lineseg =  planes[vi0.i][vi0.j];
  Vec3 v0  = lineseg[vi0.k].v;
  Vec3 v1  = lineseg[vi1.k].v;
  Vec3 lineDir=v1-v0;
  Vec3 &n1=planes[pid[0]].n;
  Vec3 &ax = planes[pid[0]].ax;
  Vec3 &ay = planes[pid[0]].ay;
  Vec3 n2 = planes[pid[1]].n;
  n2 = Vec3(n2.dot(ax),n2.dot(ay),n2.dot(n1));
  Vec3 zaxis = n2.cross(Vec3(0,0,1));
  return zaxis.dot(lineDir)<=0;
}

bool lineIntersect(Vec3 la0,Vec3 la1,Vec3 lb0, Vec3 lb1) {
  Vec3 n1=la0-la1;
  n1=Vec3(n1[1],-n1[0])/n1.norm();
  real_t prod0=(lb0-la0).dot(n1);
  real_t prod1=(lb1-la0).dot(n1);
  if( (prod0>0 && prod1>0) || (prod0<0 && prod1<0) ) {
    return false;
  }

  Vec3 n2=lb0-lb1;
  n2=Vec3(n2[1],-n2[0])/n2.norm();
  prod0=(la0-lb0).dot(n2);
  prod1=(la1-lb0).dot(n2);
  if( (prod0>0 && prod1>0) || (prod0<0 && prod1<0) ) {
    return false;
  }
  return true;
}

void PolyMesh::chopPoly(const Polygon & rect, int pid,ClipperLib::ClipType ct)
{
  ClipperLib::Clipper c;
  c.AddPolygon(poly[pid][0],ClipperLib::ptSubject);
  c.AddPolygon(rect,ClipperLib::ptClip);
  Polygons solution;
  bool ret = c.Execute(ct,solution);
  if(solution.size()<1) {
    std::cout<<"ret "<<ret<<" "<<pid<<" is completely clipped\n";
  } else {
    std::cout<<"plane "<<pid<<"\n";
    std::cout<<"soln size "<<solution.size()<<"\n";
    poly[pid][0]=solution[0];
  }
}

void vert2poly(const std::vector<Vert> & rect, Polygon & p)
{
  for(size_t ii=0;ii<rect.size();ii++){
    p.push_back(IntPoint((long64)rect[ii].v.get(0),
                        (long64)rect[ii].v.get(1)));
  }
}
void poly2vert(const Polygon & polySeg, std::vector<Vert>&lineseg){
  for(size_t kk=0; kk<polySeg.size(); kk++) {
              lineseg.push_back(Vec3((real_t)polySeg[kk].X,
                                     (real_t)polySeg[kk].Y,0));
            }
}
