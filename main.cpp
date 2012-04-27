#include <iostream>
#include "connector.hpp"
#define INT_SCALE 1000000
real_t size = 12; //inches
real_t t=0.12;//thickness
int main(int argc, char * argv[])
{
  const char * filename = "bunny_30.txt";
  size=8;
  if(argc>1){
    filename=argv[1];
  }
  if(argc>2){
    size=atoi(argv[2]);
  }
	PolyMesh pm(filename);
	pm.scale(INT_SCALE);
	pm.obj_scale=size;
	pm.zz(t);
  pm.slot(1);
  pm.slot(0.5);
  pm.connector();
  pm.teeth();
	pm.save_result("zz.txt");
	return 0;
}
