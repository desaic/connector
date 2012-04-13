#ifndef CONNECTOR_HPP
#define CONNECTOR_HPP
#include <vector>
#include <set>
#include "math.hpp"
#include <string>
#include <map>
#include <iostream>
//teeth width with respect to thickness
extern real_t teethRatio;
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
};

typedef std::vector<Vec3> Connector;
//a plane may have more one line segment if it has holes inside
class Plane{
public:
	//normal  two axis
	Vec3 n,ax,ay;
	//line segments
	std::vector<std::vector<Vert> >l;
	std::vector<Vert> &operator[](size_t index){
		return l[index];
	};
	size_t size(){
		return l.size();
	}
	
};
#include "clipper.hpp"
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
	void slot();
	void buildEdge();
	//map from input vertex id to [1..n]
	std::map<int,int>vid;
	void save_result(const char * filename);
  std::vector<Connector> conns;
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
