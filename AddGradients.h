#ifndef __ADFRcc__AddGradients__
#define __ADFRcc__AddGradients__

#include "Common.h"
#include "GridMap.h"
#include "DensityClustering.h"
#include <iostream>
#include <list>
#include <vector>

namespace ADFR {

class AddGradients
{
 public:
  AddGradients(float gradL=5.0);
  ~AddGradients(){};
  float gradLow;
  
  //std::vector<std::vector<std::vector <int> > > clusters;
  std::map<std::string, std::vector<std::vector <int> > > clusters;

  void processMaps(std::vector<GridMap *>& maps, float spacing[3], int nneighb, int cutOffValue, float errorCut, const char *filename);
  void testGridOrder(Float *gridMap, int npoints[3]);
  std::vector<std::vector <int> > getMapClusters(const std::string& mapType);
  private:
  void onepass(Float* gridMap, int npoints[3], std::vector<Vec>& indices, float generation, std::map<Vec, bool>& seen);
};

}  // namespace ADFR
#endif
