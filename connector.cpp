#include "connector.hpp"
#include <fstream>
#include <iostream>
#include <set>
#include <algorithm>

using ClipperLib::long64;
using ClipperLib::Clipper;
using ClipperLib::Polygons;
using ClipperLib::Polygon;
using ClipperLib::IntPoint;
real_t teethRatio=5;
real_t min_depth_ratio=2/teethRatio;

/**@brief
offset a little so that tv0 is not exactly on the line
*/
real_t normalOffset=3;
real_t unitRatio=2;
real_t slotUnit=1;
real_t testSlotUnit=2.0;
real_t testSlotStart=0.5;
real_t polyArea(std::vector<Vert> & l)
{
	int ii0=l.size()-1;
	real_t A=0;
	for(size_t ii=0;ii<l.size();ii++){
		A+= l[ii0].v[0]*l[ii].v[1]-l[ii].v[0]*l[ii0].v[1];
		ii0=ii;
	}
	return A;
}

void PolyMesh::save_result(const char * filename)
{
	std::ofstream out;
	out.open(filename);
	if(!out.good()){
		std::cout<<"cannot open file"<<filename<<"\n";
		return;
	}
	out<<planes.size()<<"\n";
  for(size_t ii=0; ii<poly.size(); ii++) {
    out<<poly[ii].size()<<"\n";
    if(poly[ii].size()<1){
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
  out.close();
}


PolyMesh::PolyMesh(const char * filename):intscale(1),obj_scale(1),
                                          t(0){
	std::ifstream in;
	in.open(filename);
	if(!in.good()){
		std::cout<<"cannot open"<<filename<<"\n";
		return;
	}
	size_t nplane=0;
	in>>nplane;

	
	int nvert=0;
	int pcnt=0;
	for(size_t ii=0;ii<nplane;ii++){
		int nseg = 0 ;
		in>>nseg;
		if(nseg <= 0){
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
		for(size_t jj=0;jj<p.l.size();jj++){
			//number of points in a segment
			int npt = 0;
			in >> npt;
			p[jj].resize(npt);
			for(size_t kk=0;kk<p[jj].size();kk++){
				in >> p[jj][kk].v[0];
				char c;
				//comma;
				in >>c;
				in >> p[jj][kk].v[1];
			}
			for(size_t kk=0;kk<p[jj].size();kk++){
				int id;
				in >> id;
				if(vid.find(id)==vid.end()){
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
	for(size_t ii=0;ii<planes.size();ii++){
		if(planes[ii].size()<2){
			continue;
		}
		//to deal with odd number of line segments
		int outmost = planes[ii].size()-1;
		for(size_t jj=1;jj<planes[ii].size();jj+=2){
			//jj is inside jj-1
			if(pnpoly(planes[ii][jj-1],planes[ii][jj][0])){
				outmost=jj-1;
				break;
			}//jj-1 is inside jj
			else if(pnpoly(planes[ii][jj],planes[ii][jj-1][0])){
				outmost=jj;
				break;
			}
		}
		if(outmost!=0){
			//swap the outmost line segment to the first segment
			std::vector<Vert>  tmp=planes[ii][0];
			planes[ii][0]=planes[ii][outmost];
			planes[ii][outmost]=tmp;
		}
	}

	//now fix winding orders
	for(size_t ii=0;ii<planes.size();ii++){
	
		real_t A=polyArea(planes[ii][0]);
		if(A<0){
			std::reverse(planes[ii][0].begin(),planes[ii][0].end());
		}
		for(size_t jj=1;jj<planes[ii].size();jj++){
			A=polyArea(planes[ii][jj]);
			if(A>0){
				std::reverse(planes[ii][jj].begin(),planes[ii][jj].end());
			}
		}
	}
	in.close();
};

void PolyMesh::scale(real_t s)
{
	intscale=s;
	for(size_t ii=0;ii<planes.size();ii++){
		for(size_t jj=0;jj<planes[ii].size();jj++){
			for(size_t kk=0;kk<planes[ii][jj].size();kk++){
				planes[ii][jj][kk].v*=s;
			}
		}
	}
}

void PolyMesh::buildEdge(){
	vertp.resize(vid.size());
	for(size_t ii=0;ii<planes.size();ii++){
		for(size_t jj=0;jj<planes[ii].size();jj++){
			size_t kk1=planes[ii][jj].size()-1;
			for(size_t kk=0;kk<planes[ii][jj].size();kk++){
				int v0=planes[ii][jj][kk].id;
				int v1=planes[ii][jj][kk1].id;
				vertp[v0][ii]=VertIdx(ii,jj,kk);
				Edge e(v0,v1);
				std::map<Edge,EdgeVal>::iterator it=eset.find(e);
				if(it!=eset.end()){
					eset[e].p[1]=ii;
				}
				else{
					eset[e].p[0]=ii;
					eset[e].p[1]=-1;
				}
				kk1=kk;
			}
		}
	}
	
}

void make_rect(Vec3 &start, Vec3 & lineDir, Vec3 &normal,
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

void PolyMesh::slot()
{
  std::map<Edge,EdgeVal>::iterator it;
  real_t unitlen = unitRatio*t;
  real_t start = unitlen;
  real_t slot_len=unitlen*slotUnit;
  real_t test_start= testSlotStart*unitlen;
  real_t test_slot_len=unitlen*testSlotUnit;
  //  real_t end = start+slot_len;
	for(it=eset.begin();it!=eset.end();it++){
    int pid[2]={it->second.p[0],it->second.p[1]};
		int v0pIdx=vertp[it->first.id[0]][pid[0]].k;
		int v1pIdx=vertp[it->first.id[1]][pid[0]].k;
    int npt = planes[pid[0]][vertp[it->first.id[0]][pid[0]].j].size();
		if(v1pIdx!=(v0pIdx+1)%npt){
			int tmpPIdx = pid[0];
			pid[0]=pid[1];
			pid[1]=tmpPIdx;
		}
    bool possible = true;
		for(size_t ii=0;ii<2;ii++){
			int planeIdx = pid[ii];
			VertIdx& vi0 = vertp[it->first.id[0]][planeIdx];
			VertIdx& vi1 = vertp[it->first.id[1]][planeIdx];
      std::vector<Vert> & lineseg =  planes[vi0.i][vi0.j];
			Vec3 v0  = lineseg[vi0.k].v;
			Vec3 v1  = lineseg[vi1.k].v;

			Vec3 lineDir=v1-v0;
			real_t len = lineDir.norm();
			lineDir/=len;
			Vec3 lineNormal = Vec3(lineDir[1],-lineDir[0],0);
      if(ii>0){
        lineNormal = -lineNormal;
      }
      Vec3 mid=(v0+v1)/2;
      mid-=lineNormal*test_start;
      mid-=lineDir*t/2;
      std::vector<Vert>rect;

      make_rect(mid,lineDir,lineNormal,t,test_slot_len, rect);
      size_t jj0=rect.size()-1;
      for(size_t jj=0;jj<rect.size();jj++){
        //TODO:
        //reject if any vertex of the rectangle is not inside the polygon
        //or any edge of the rectangle intersects with the polygon
        bool expected=true;
        for(size_t seg=0;seg<poly[planeIdx].size();seg++){
          Polygon polySeg = poly[planeIdx][seg];
          std::vector<Vert>lineseg;
          for(size_t kk=0;kk<polySeg.size();kk++){
            lineseg.push_back(Vec3((real_t)polySeg[kk].X,(real_t)polySeg[kk].Y,0));
          }
          bool ret = pnpoly(lineseg, rect[jj]);
          if(ret!=expected){
            //std::cout<<"outside "<<planeIdx<<"\n";
            possible = false;
            goto ENDEDGELOOP;
          }
          size_t kk0=lineseg.size()-1;
          for(size_t kk=0;kk<lineseg.size();kk++){
            bool intersect = lineIntersect(rect[jj].v,rect[jj0].v,
                                           lineseg[kk0].v,lineseg[kk].v);
            if(intersect){
              //std::cout<<"intersect"<<planeIdx<<"\n";
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
    }
    if(possible){
      for(size_t ii=0;ii<2;ii++){
        int planeIdx = pid[ii];
        VertIdx& vi0 = vertp[it->first.id[0]][planeIdx];
        VertIdx& vi1 = vertp[it->first.id[1]][planeIdx];
        std::vector<Vert> & lineseg =  planes[vi0.i][vi0.j];
        Vec3 v0  = lineseg[vi0.k].v;
        Vec3 v1  = lineseg[vi1.k].v;
        Vec3 mid=(v0+v1)/2;

        Vec3 lineDir=v1-v0;
        real_t len = lineDir.norm();
        lineDir/=len;
        Vec3 lineNormal = Vec3(lineDir[1],-lineDir[0],0);
        if(ii>0){
          lineNormal = -lineNormal;
        }
        mid=(v0+v1)/2;
        mid-=lineNormal*start;
        mid-=lineDir*t/2;
        std::vector<Vert> rect;
        make_rect(mid, lineDir, lineNormal,t,slot_len,rect);
        Polygon p ;
        for(size_t jj=0;jj<rect.size();jj++){
          p.push_back(IntPoint((long64)rect[jj][0],(long64)rect[jj][1]));
        }
        poly[planeIdx].push_back(p);
        planes[planeIdx].l.push_back(rect);
      }
    }
  ENDEDGELOOP:
    ;
    std::cout<<"possible"<<possible<<"\n";
  }
}

void PolyMesh::zz(real_t _t)
{
	buildEdge();
  t=(_t/obj_scale)*intscale;
	real_t width =t*teethRatio;
	//	std::map<Edge, bool> processed;
	std::map<Edge,EdgeVal>::iterator it;
	poly.resize(planes.size());

	for(size_t ii=0;ii<planes.size();ii++){
		poly[ii].resize(planes[ii].size());
		for(size_t jj=0;jj<planes[ii].size();jj++){
			for(size_t kk=0;kk<planes[ii][jj].size();kk++){
				poly[ii][jj].push_back(ClipperLib::IntPoint(
					(long64)planes[ii][jj][kk].v[0],
					(long64)planes[ii][jj][kk].v[1]));
			}
		}
	}
	for(it=eset.begin();it!=eset.end();it++){
		real_t depth=0;
		int pid[2]={it->second.p[0],it->second.p[1]};
		Vec3 n1=planes[pid[0]].n, n2=planes[pid[1]].n;
		real_t cosine=n1.dot(n2);
		if(cosine>=0){
			//obtuse angle
			depth = std::sqrt(1-cosine*cosine)*t/2;
		}else{
			//acute angle
			real_t sine = std::sqrt(1-cosine*cosine);
			real_t tangent = sine/cosine;
			//tangent half angle formula
			real_t halftangent = tangent/(1+std::sqrt(1+tangent*tangent));
			if(halftangent<min_depth_ratio){
				halftangent=min_depth_ratio;
			}
			depth = (t/2)/halftangent;
		}

    //Only truncate concave connections
    if(depth<1){
			continue;
		}
		//check which plane is on the left and swap if p[0] is on the right
		//use winding order
		//index of the first edge vertex in the plane
		int v0pIdx=vertp[it->first.id[0]][pid[0]].k;
		int v1pIdx=vertp[it->first.id[1]][pid[0]].k;
    int npt = planes[pid[0]][vertp[it->first.id[0]][pid[0]].j].size();
		//plane ids
		if(v1pIdx!=(v0pIdx+1)%npt){
			int tmpPIdx = pid[0];
			pid[0]=pid[1];
			pid[1]=tmpPIdx;
		}
		for(size_t ii=0;ii<1;ii++){
			real_t s=0;
			int planeIdx = pid[ii];
			VertIdx& vi0 = vertp[it->first.id[0]][planeIdx];
			VertIdx& vi1 = vertp[it->first.id[1]][planeIdx];
      std::vector<Vert> & lineseg =  planes[vi0.i][vi0.j];
			Vec3 v0  = lineseg[vi0.k].v;
			Vec3 v1  = lineseg[vi1.k].v;
      //vertex before v0
      Vec3 v_1 = vi0.k==0 ? lineseg[lineseg.size()-1].v : lineseg[vi0.k-1].v;
      //vertex after v1
      Vec3 v2  = lineseg[(vi1.k+1)%lineseg.size()].v;

			Vec3 lineDir=v1-v0;
			real_t len = lineDir.norm();
			lineDir/=len;
			//(a,b)-->(b,-a) is perpendicular on the right of the edge
			//(-b,a) is on the other side
			Vec3 lineNormal = Vec3(-lineDir[1],lineDir[0],0);
			//second time teeth grow the other way
			if(ii==1){
				lineNormal=-lineNormal;
			}

      Vec3 dir0 = v_1-v0;
      dir0/=dir0.norm();
      if(dir0.dot(lineNormal)>=0){
        dir0=-dir0;
      }
      if(dir0.dot(lineNormal)>-0.3){
        dir0=-lineNormal;
      }
      real_t depth0=depth/(-dir0.dot(lineNormal));

      Vec3 dir1 = v2-v1;
      dir1/=dir1.norm();
      if(dir1.dot(lineNormal)>=0){
        dir1=-dir1;
      }
      if(dir1.dot(lineNormal)>-0.3){
        dir1=-lineNormal;
      }
      real_t depth1=depth/(-dir1.dot(lineNormal));


			//first start at 0
			//second time start at 1
			int nTeeth=ii;
			//teeth vertices
			// <-- lineNormal
			//          line direction
			// tv01--tv0   |
			//    |  |     |
			//    |  |    \|/
			//tv11|__|tv1

      width=len;
			s=nTeeth*width;
      //      while(s<len){
				ClipperLib::Clipper c;
				c.AddPolygons(poly[pid[ii]],ClipperLib::ptSubject);
				
				Polygons rect(1);
				
				real_t alpha = -0.0001;//s/len;
				Vec3 tv0=(1-alpha)*v0+alpha*v1;
        //				Vec3 tv01=tv0+depth*lineNormal;
        Vec3 tv01=tv0+depth0*dir0;
				//offset a little so that tv0 is not exactly on the line
				tv0 -= lineNormal*normalOffset;
				rect[0].push_back(IntPoint((long64)tv0[0],(long64)tv0[1]));
				rect[0].push_back(IntPoint((long64)tv01[0],(long64)tv01[1]));

				alpha = 1.00001;//(nTeeth+1)*width/len;
        //  if(alpha>1){alpha=1.0000001;}
				Vec3 tv1=(1-alpha)*v0+alpha*v1;
				//Vec3 tv11=tv1+lineNormal*depth;
        Vec3 tv11=tv1+depth1*dir1;
				tv1-=lineNormal*normalOffset;
				rect[0].push_back(IntPoint((long64)tv11[0],(long64)tv11[1]));
				rect[0].push_back(IntPoint((long64)tv1[0],(long64)tv1[1]));

				c.AddPolygons(rect,ClipperLib::ptClip);

				Polygons solution;
				bool ret = c.Execute(ClipperLib::ctDifference,solution);
        if(solution.size()<1){
           std::cout<<"ret "<<ret<<" "<<pid[ii]<<" is completely clipped\n";
        }else{
          std::cout<<"plane "<<pid[ii]<<"\n";
          std::cout<<"soln size "<<solution.size()<<"\n";
          poly[pid[ii]]=solution;
          
        }
				nTeeth+=2;
				s=nTeeth*width;
        //	}
    }
	}
}

bool lineIntersect(Vec3 la0,Vec3 la1,Vec3 lb0, Vec3 lb1){
  Vec3 n1=la0-la1;
  n1=Vec3(n1[1],-n1[0]);
  if( (lb0-la0).dot(n1) * (lb1-la0).dot(n1)>0){
    return false;
  }

  Vec3 n2=lb0-lb1;
  n2=Vec3(n2[1],-n2[0]);
  if( (la0-lb0).dot(n2) * (la1-lb0).dot(n2)>0){
    return false;
  }
  return true;
}
