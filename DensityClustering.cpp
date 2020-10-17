//

#include "DensityClustering.h"
#include <map>
#include <list>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <functional>
#include <string>


#if defined(_MSC_VER) && _MSC_VER < 1900
#include <stdarg.h>  
#define snprintf c99_snprintf
#define vsnprintf c99_vsnprintf

__inline int c99_vsnprintf(char *outBuf, size_t size, const char *format, va_list ap)
{
    int count = -1;

    if (size != 0)
        count = _vsnprintf_s(outBuf, size, _TRUNCATE, format, ap);
    if (count == -1)
        count = _vscprintf(format, ap);

    return count;
}

__inline int c99_snprintf(char *outBuf, size_t size, const char *format, ...)
{
    int count;
    va_list ap;

    va_start(ap, format);
    count = c99_vsnprintf(outBuf, size, format, ap);
    va_end(ap);

    return count;
}

#endif



struct greater_than
{
  inline bool operator() (const std::vector<int>& v1, const std::vector<int>& v2)
  {
    return (v1.size() > v2.size());
  }
};

Point::Point()   // Default Constructor
{
   x = 0;
   y = 0;
   z = 0;
   incluster = false;
}

Point::Point(int i, int j, int k) //Constructor
{
   x = i;
   y = j;
   z = k;
   incluster = false;
}

Point::Point(const Point &copyin)   // Copy constructor to handle pass by value.
{                             
   x = copyin.x;
   y = copyin.y;
   z = copyin.z;
   incluster = copyin.incluster;  
}

std::ostream &operator<<(std::ostream &output, const Point &point)
{
  output << point.x << ' ' << point.y << ' ' << point.z  << std::endl;
   return output;
}

Point& Point::operator=(const Point &rhs)
{
   this->x = rhs.x;
   this->y = rhs.y;
   this->z = rhs.z;
   return *this;
}

int Point::operator==(const Point &rhs) const
{
   if( this->x != rhs.x) return 0;
   if( this->y != rhs.y) return 0;
   if( this->z != rhs.z) return 0;
   return 1;
}

// This function is required for built-in STL list functions like sort
int Point::operator<(const Point &rhs) const
{
   if( this->x == rhs.x && this->y == rhs.y && this->z < rhs.z) return 1;
   if( this->x == rhs.x && this->y < rhs.y) return 1;
   if( this->x < rhs.x ) return 1;
   return 0;
}

class PointLessThan
{
public:
    bool operator( )(const Point& p1, const Point& p2) const {
      if(p1.x == p2.x && p1.y == p2.y && p1.z < p2.z) return(true);
      if( p1.x == p2.x && p1.y < p2.y) return(true);
      if( p1.x < p2.x ) return(true);
      return(false);
   }
};

DensityClustering::DensityClustering()
{
  neighborpts=14;
  //spacing
}

DensityClustering::DensityClustering(float sp[3], int npts)
{
  for (int i=0; i<3; i++)  spacing[i]=sp[i];
  neighborpts = npts;
} 


void DensityClustering::findClusters(std::vector<Vec> &gridIndices, int cVolcut)
{
  // gridIndices is a vector of (i, j, k) grid indices of negative grid points:

  int si, sj, sk;
  si = sj = sk = 1;
  std::vector<Vec>::iterator indIter;
  // std::map ptCache containing all points in form {"i j k": Point(i,j,k, ii) }
  // where ii is in range(0, gridIndices.size()-1)
  
  //std::map<std::string, Point> ptCache;
  //std::map<std::string, Point>::iterator cacheIter;
  std::map<std::string, int> ptCache;
  std::map<std::string, int>::iterator cacheIter;
  
  // std::map noClusterPoints contains points that have not been clustered yet
  //std::map<int, Point> noClusterPts;
  std::map<int, bool> noClusterPts;
  std::map<int, bool>::iterator ptrIter;
  std::vector<Point> pointVec;
  // single cluster (vector of  ii indices) 
  std::vector<int> cluster;
  // a list (vector) of clusters 
  std::vector<std::vector<int> > _clusters;

  // std::map inCluster contains indices of clustered points in form {ii : "i j k"}
  std::map<int, std::string> inCluster;
  // format to create std::string object containing i, j ,k indices
  const char *format = "%d %d %d";
  std::string indStr;
  
  int i,j,k, ii, size;
  ii = 0;
  //for (ii=0; ii<ninds; ii++)
  for (indIter = gridIndices.begin(); indIter != gridIndices.end(); ++indIter)
    {
      i = (*indIter)[0]; j = (*indIter)[1]; k = (*indIter)[2];
      Point p(i, j, k);
      // create a string for i, j, k. This will be a key in ptCache map.
       
      size = snprintf(NULL, 0, format, i,j,k);
      //std::cout << "create string for " << i<< " "<<j<<" "<<k << " " << format << " " << size << std::endl;
      indStr = std::string(size + 1, '\0');
      #if defined(_MSC_VER) && _MSC_VER < 1900
      sprintf_s(&indStr[0], size+1, format, i,j,k);
      #else
      sprintf(&indStr[0], format, i,j,k);
      #endif
      //std::cout << "kword " << indStr << std::endl;
      //if (i % si == 0 and j % sj == 0 and k % sk == 0)
      pointVec.push_back(p);
      noClusterPts[ii] = true;
      ptCache[indStr] = ii;
      ii += 1;
    }
  int ci = 0; //index of point in current cluster
  Point point;
  
  // add first point to cluster
  ptrIter = noClusterPts.begin();
  // first point in the list of points that have not been clustered yet
  ii  = ptrIter->first;
  point = pointVec[ii];
  cluster.push_back(ii);
  pointVec[ii].incluster = true;
  noClusterPts.erase(ptrIter); // remove this point from map containing points that have not been included an any clusters
  int ncl = 0;
  while(true)
    {
      // search for neighbors of the cluster[ci] point  
      for (i=point.x-si; i<=point.x+si; i+=si)
	for (j=point.y-sj; j<=point.y+sj; j+=sj)
	  for (k=point.z-sk; k<=point.z+sk; k+=sk)
	    {
	     
	     // create an "i j k" key from neighbor grid indices 
	     size = snprintf(NULL, 0, format, i,j,k);
	     indStr = std::string(size + 1, '\0');
	     #if defined(_MSC_VER) && _MSC_VER < 1900
	     sprintf_s(&indStr[0], size+1, format, i,j,k);
             #else
	     sprintf(&indStr[0], format, i,j,k);
             #endif
	     //sprintf(&indStr[0], format, i,j,k);
	     // look up this key in PtCache map
	     //std::cout << "Looking for neighbor " << i<< " "<<j<<" "<<k << " " << indStr << std::endl; 
	     cacheIter = ptCache.find(indStr);
	     if (cacheIter != ptCache.end())
		{
		  int nind = cacheIter->second;
		  if (i == point.x && j == point.y && k == point.z){
		    // this is the same point , skip
		    pointVec[nind].incluster = true; // this is not necessary. Incluster should be True without this call. Well , just in case.
		  }
		  else if (pointVec[nind].incluster == false)
		    {
		      // we have a neighbor point that is not in any cluster 
		      //int nind = cacheIter->second.index;
		      cluster.push_back(nind);
		      //inCluster[nind] = indStr;
		      //cacheIter->second.incluster = true;
		      pointVec[nind].incluster = true;
		      ptrIter = noClusterPts.find(nind);
		      if (ptrIter != noClusterPts.end())
			{
			  noClusterPts.erase(ptrIter);
			}
		    }
		}
	    } // closes for (k= ...) loop
      ci += 1;
      if (ci < cluster.size()) // more points in this cluster can be checked for neighbors
	{
	  ii = cluster[ci];
	  point = pointVec[ii];
	}
      else //all points in cluster have been checked for neighbors
	{
	  //std::cout << "cluster: ";
	  //for (int cc=0; cc<cluster.size(); cc++)
	  //  std::cout << cluster[cc] << " ";
	  //std::cout << std::endl;
	  _clusters.push_back(cluster);
	  if  (noClusterPts.size() == 0)
	    {
	      ptCache.clear();
	      pointVec.clear();
	      break; //Done clustering
	    }
	  cluster.clear();
	  ci = 0;
	  ncl = ncl +1;
	  ptrIter = noClusterPts.begin();
	  //point = ptrIter->second;
	  ii = ptrIter->first; 
	  point = pointVec[ii];
	  cluster.push_back(ii);
	  noClusterPts.erase(ptrIter);
	}
    } //closes while(true)
  std::sort(_clusters.begin(), _clusters.end(), greater_than());

  for (i=0; i<_clusters.size(); i++)
     {
       int s1 = _clusters[i].size();
       if (s1 >= cVolcut)
	 {
	   this->clusters.push_back(_clusters[i]);
	 }
     }
}


