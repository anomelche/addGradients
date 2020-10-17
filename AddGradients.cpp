#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <map>
#include "AddGradients.h"
//# include <omp.h>
namespace ADFR {
AddGradients::AddGradients(float gradL)
{
  gradLow = gradL;
}

  void AddGradients::processMaps(std::vector<GridMap*>& maps, float spacing[3], int nneighb, int cutOffValue, float errorCut, const char *filename)
  {
    /**cutOffValue can be either -1 to keep only the largest cluster of points as 
       the "outside",  or a positive integer to keep clusters with more than that many points 
       as the "outside" when computing the gradient for the receptor "inside" **/
    
    //for (std::vector<GridMap>::iterator nm=maps.begin(); nm != maps.end(); ++nm)
    //std::cout << "Num of procs: " << omp_get_num_procs() << std::endl;
    //std::cout << "Max num of threads: " << omp_get_max_threads() << std::endl;

    /** Open file (if filename is not NULL) for writing. 
    The purpose of this temporary file is to write the
    map type name of each processed map on a separate line.
    The calling python function reads this file (in a separate thread) and 
    uses it to display the progress of computing the gradient **/
    std::ofstream outfile;
    if (filename)
      {
	outfile.open(filename);
      }
    #pragma omp parallel for
    for (int nm = 0; nm < maps.size(); nm++) 
      {
	//std::cout << "IN OMP FOR LOOP: Num of threads: " << omp_get_num_threads() << std::endl;
	//std::cout << "Thread number: " << omp_get_thread_num() << std::endl;
	// get the data array
	Float* mapData = maps[nm]->gridData;
	int* npoints = maps[nm]->numGridPoints;
	//std::cout << maps[nm]->getMapType() << std::endl;
        std::string mapType = maps[nm]->getMapType();
	int nx, ny, nz, i, j, k, numPoints;
	nx = npoints[0];
	ny = npoints[1];
	nz = npoints[2];
	numPoints = nx*ny*nz;
	Vec vv(3,0);
	std::vector<Vec> gridInds;  // this vector obj will contain grid indices (i, j, k) of grid points with negative values
	float generation = this->gradLow;
	float nzy = nz*ny;
	for (int nn=0; nn < numPoints; nn++)
	  {
	    if (mapData[nn] <= generation)
	      {
		// compute  3D index into 1D grid data array C-order
		i = nn/nzy;
		j = (nn-nzy*i)/nz;
		k = nn-nz*j-nzy*i;
		vv[0] = i; vv[1] = j; vv[2] = k;
		gridInds.push_back(vv);
	      }
	  }
	//compute clusters
	DensityClustering dc(spacing, nneighb);
	dc.findClusters(gridInds);
	//this->clusters.push_back(dc.clusters);
	this->clusters[mapType] = dc.clusters;
	// dc.clusters is a vector of integer vectors [ [ii_0 , ii_1 ...], [....]]
	// where ii_# is an int in range (0, len(gridInds)-1)
	std::map<Vec, bool> seen;
	for (i=0; i<nx; i++)
	  for (j=0; j<ny; j++)
	    for (k=0; k<nz; k++)
	      {
		vv[0] = i; vv[1] = j; vv[2] = k;
		seen[vv] = false;
	      }
	Vec cluster0;
	std::vector<Vec> outsideIndices;
	for (i=0; i<dc.clusters.size(); i++)
	  {
	    //std::cout << "##### cluster and size " << i << " " <<  dc.clusters[i].size() << " "<< cutOffValue << "\n";
	    if (cutOffValue < 0)
	      {
		// keep only the largest cluster of points as the "outside" 
		if (i==0)
		  cluster0 = dc.clusters[0];
		else
		  break;
	      }
	    else // cutOffValue is a positive number
	      {
		//keep clusters with more than cutOffValue points as the "outside"
		if (dc.clusters[i].size() > cutOffValue)
		  cluster0 = dc.clusters[i];
		else
		  continue;
	      }
	    //int ct = 0;
	    for (Vec::iterator it = cluster0.begin(); it != cluster0.end(); ++it)
	      {
		vv = gridInds[*it];
		outsideIndices.push_back(vv);
		seen[vv] = true;
		//ct += 1;
	      }
	    //std::cout << "##### Points less equal gradlow" << ct << "\n";
	  }
	while(outsideIndices.size() > 0)
	  {
	    generation += 1.0 ;
	    onepass(mapData, npoints, outsideIndices, generation, seen);
	  }
	if (filename)
	  {
            #pragma omp critical
	    outfile << maps[nm]->getMapType() << std::endl;
	  }
      } // end of for loop over maps
    if (filename)
      outfile.close();
  } // end of processMaps


  void AddGradients::onepass(Float *gridMap, int *npoints, std::vector<Vec>& indices, float generation, std::map<Vec, bool>& seen)
  {
    std::vector<Vec> edge;
    int a,b,c, a1, b1, c1, nx, ny, nz;
    nx = npoints[0];
    ny = npoints[1];
    nz = npoints[2];
    Vec vv(3, 0);
    std::map<Vec, bool>::iterator iter;
    for (std::vector<Vec>::iterator it = indices.begin(); it != indices.end(); ++it)
      {
	a = (*it)[0]; b = (*it)[1]; c = (*it)[2];
	//int nbrinds[18] = {a-1,b,c,  a+1,b,c,  a,b-1,c,  a,b+1,c,  a,b,c-1,  a,b,c+1};
	/**
        int nbrinds[78] = {
	 a-1, b-1, c-1,   a-1, b-1, c,   a-1, b-1, c+1,   a-1, b, c-1,   a-1, b, c,   a-1, b, c+1,
	 a-1, b+1, c-1,   a-1, b+1, c,   a-1, b+1, c+1,   a, b-1, c-1,   a, b-1, c,   a, b-1, c+1,
	 a, b, c-1,       a, b, c+1,   a, b+1, c-1,   a, b+1, c,   a, b+1, c+1,   a+1, b-1, c-1,
	 a+1, b-1, c,   a+1, b-1, c+1,   a+1, b, c-1,   a+1, b, c,   a+1, b, c+1,   a+1, b+1, c-1,
	 a+1, b+1, c,   a+1, b+1, c+1 };
        **/
	int nbrinds[54] = {
	  a-1, b-1, c,   a-1, b, c-1,   a-1, b, c,   a-1, b, c+1,
	  a-1, b+1, c,   a, b-1, c-1,   a, b-1, c,   a, b-1, c+1,
	  a, b, c-1,     a, b, c+1,     a, b+1, c-1,   a, b+1, c,   a, b+1, c+1,
	  a+1, b-1, c,   a+1, b, c-1,   a+1, b, c,   a+1, b, c+1,   a+1, b+1, c };
	
	for (int i=0; i<18; i++)
	  {
	    a1 = nbrinds[3*i]; b1 = nbrinds[3*i+1]; c1 = nbrinds[3*i+2];
	    vv[0] = a1; vv[1] = b1; vv[2] = c1;
	    iter = seen.find(vv);
	    if (iter != seen.end()) // found this point in seen
	      {
		if (iter->second == false)
		  {
		    edge.push_back(vv);
		    iter->second = true;
		    gridMap[a1 * ny * nz + b1 * nz + c1] = generation;
		    
		  }
	      }
	  }
      }
    indices = edge;
  }

  
  //for debugging
  void AddGradients::testGridOrder(Float *gridMap, int npoints[3])
  {
    int i, j, k , nx, ny, nz;
    nx = npoints[0];
    ny = npoints[1];
    nz = npoints[2];
    std::cout << "nx, ny, nz == " << nx << " " << ny << " " << nz << std::endl;
    for (i=0; i<nx; i++)
      for (j=0; j<ny; j++)
	for (k=0; k<nz; k++)
	  {
	    int ind = i * ny * nz + j * nz + k;
	    std::cout << i << " " << j << " " << k << " (" << ind << ") is " << gridMap[ind] << std::endl; 
	  }
  } 

  
  std::vector<std::vector <int> > AddGradients::getMapClusters(const std::string& mapType)
  {
    return this->clusters[mapType];
  }
}  // namespace ADFR
