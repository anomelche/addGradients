//DensityClustering algorithm
#ifndef __ADFRcc__DensityClustering__
#define __ADFRcc__DensityClustering__

#include <iostream>
#include <list>
#include <vector>

typedef std::vector<int> Vec;
//class to represent a grid point 
class Point
{
  friend std::ostream &operator<<(std::ostream &, const Point &);
  
  public:

      // grid point indices
      int x;
      int y;
      int z;
      // index of the point in a list of indices 
      int index;
      // indicates if the point is in cluster 
      bool incluster;
      Point();
      Point(int i, int j, int k); 
      Point(const Point &);
      ~Point(){};
      Point &operator=(const Point &rhs);
      int operator==(const Point &rhs) const;
      int operator<(const Point &rhs) const;
};


class DensityClustering
{
 public:

  int neighborpts;
  float spacing[3];
  std::vector<std::vector <int> > clusters;
  DensityClustering();
  DensityClustering(float spacing[3], int neighborpts=14);
  ~DensityClustering(){};
  void findClusters(std::vector<Vec> &gridIndices, int cVolcut=1);
  };


#endif
