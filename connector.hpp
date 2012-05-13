#ifndef CONNECTOR_HPP
#define CONNECTOR_HPP
#include <vector>
#include <set>
#include "math.hpp"
#include <string>
#include <map>
#include <iostream>
#include "clipper.hpp"
//teeth width with respect to thickness
extern real_t teethRatio;
using ClipperLib::Polygon;
struct Vert{
	int id;
	Vec3 v;
	Vert(const Vec3 & _v=Vec3(0,0,0) ,int _id=0):id(_id),v(_v){}
  real_t operator[](size_t i){return v[i];}
};

struct VertIdx{
	//polygon index
	//line segment index
	//point index
	int i,j,k;
	VertIdx(int _i=0,int _j=0,int _k=0):i(_i),j(_j),k(_k){
	}

};
struct Edge{
	//id of its vertices
	int id[2];

	Edge(int id0=0,int id1=0){
		if(id0>id1){
			id[0]=id1;
			id[1]=id0;
		}else{
			id[0]=id0;
			id[1]=id1;
		}
	}
	bool operator<(const Edge e)const{
		if(id[0]<e.id[0]){
			return true;
		}
		if(id[0]>e.id[0]){
			return false;
		}
		if(id[1]<e.id[1]){
			return true;
		}
		return false;
	}
};
/**@brief used in map*/
struct EdgeVal{
	//planes at its two sides
	//two planes may share more than one edge if there is a hole between the two planes
	int p[2];
  bool hasConn;
  Vec3 v0, norm,dir;
  real_t connSize;
  EdgeVal():hasConn(false){
  }
};

struct Connector{
  std::vector<Vec3>l;
  std::vector<Vec3>world;
  int pid[2];
  /**@brief point attached to pid[0]*/
  Vec3 v0, norm,dir;
  Vec3 local2plane(const Vec3 & v);
  size_t size()const {
    return l.size();
  }
  Vec3 & operator [] (int idx){
    return l[idx];
  }
  Vec3 operator [] (int idx)const {
    return l[idx];
  }
};

//a plane may have more one line segment if it has holes inside
class Plane{
public:
	//normal  two axis  position of the first vertex
	Vec3 n,ax,ay,v0;
	Vec3 local2world(const Vec3& v);
	//line segments
	std::vector<std::vector<Vert> >l;
	//planes in world coordinates
	std::vector<std::vector<Vec3> >world_p;
	std::vector<Vert> &operator[](size_t index){
		return l[index];
	};
	size_t size()const{
		return l.size();
	}

};

using ClipperLib::Polygons;
class PolyMesh{
public:
	std::vector<Plane>planes;

	std::vector<Polygons>poly;
	PolyMesh(const char * filename);
	std::map<Edge,EdgeVal> eset;
	real_t intscale;
	/**@brief length of the longest side of the bounding box of the object*/
	real_t obj_scale;
  real_t t;
	/**@brief what planes are at each vertex
	maps from vertex indices to indices into planes
	*/
	std::vector<std::map<int, VertIdx> > vertp;

	void scale (real_t s);
	/**@param t thickness*/
	void zz(real_t t);
	bool isConvex(const Edge & e, const EdgeVal&ev);
	void slot(real_t frac);
  void connector(const Edge & e,const  EdgeVal & ev,Connector & conn);
  bool intersect(const Connector & conn);
  bool linePlaneInter(Plane & p, const Vec3& v0, const Vec3& v1);
  void teeth();
	void buildEdge();
	void draw();
	void edgeVertPos(const Edge &e, int planeIdx, Vec3 &v0,Vec3 & v1);
	//map from input vertex id to [1..n]
	std::map<int,int>vid;
	void save_result(const char * filename);
  std::vector<Connector> conns;
  void chopAlongEdge(const Edge&e, const EdgeVal & ev, int ii,real_t depth,  Polygon & rect);
  void chopPoly(const Polygon & poly, int pid, ClipperLib::ClipType ct=ClipperLib::ctDifference);
  /**@param len length to chop away on each piece

  */
  void chopLen(const Edge&e, const EdgeVal & ev, real_t & len1, real_t & len2);
};

template<class T>
bool pnpoly(std::vector<T>& l,T & test)
{
  size_t i, j;
  bool c= false;
  for (i = 0, j = l.size()-1; i < l.size(); j = i++) {
	  if ( ((l[i][1]>test[1]) != (l[j][1]>test[1])) &&
		(test[0] < (l[j][0]-l[i][0]) * (test[1]-l[i][1])
		/ (l[j][1]-l[i][1]) + l[i][0]) ){
		c = !c;
	  }
  }
  return c;
}
/**@brief two line segments intersect if
   for both segments,
   the end points are at different sides of the other segment

*/
bool lineIntersect(Vec3 la0,Vec3 la1,Vec3 lb0, Vec3 lb1);

#endif
