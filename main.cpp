#include <iostream>
#include "connector.hpp"
#define INT_SCALE 1000000
real_t size = 8; //inches
real_t t=0.125;//thickness
int main(int argc, char * argv[])
{
  const char * filename = "plane.txt";
  if(argc>1){
    filename=argv[1];
  }
	PolyMesh pm(filename);
	pm.scale(INT_SCALE);
	pm.obj_scale=size;
	pm.zz(t);
  pm.slot();
	pm.save_result("zz.txt");
	return 0;
}
