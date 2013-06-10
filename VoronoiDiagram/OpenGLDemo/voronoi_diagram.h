#ifndef VORONOI_DIAGRAM_H_
#define VORONOI_DIAGRAM_H_

#include <list>
using std::list;
#include <vector>
using std::vector;

#include "basic_types.h"
#include "GeometryTool.hpp"

class VoronoiDiagram
{

public:
    VoronoiDiagram( void );
    ~VoronoiDiagram( void );
    
    bool ToLeft( Site* s1, Site* s2, Site* new_site );
    bool IsIntersectionPointOnHalfedge( double x, double y, Halfedge* he );
    Face* LocateInFace( Point* p );
    bool UpdateConvexHull( Site* new_site );
    Site* FindClosestSite( Point* p );
    void CalcIntersectionVertices( Halfedge* he, Site* site, bool is_in_convex_hull, Vertex* vbegin, Vertex* vend );
    void GenerateNewHalfedges( Site* new_site, Site* site, Halfedge* he1, Halfedge* he2 );
    VoronoiDiagram* IncrementalConstruction( vector< Point* > points );
    VoronoiDiagram* DevideConquerConstruction( vector< Point > &points );

 

public:
    vector< Site* > sites;
    vector< Site* > convex_hull;
    vector< Vertex* > vertices;
    vector< Halfedge* > halfedges;
    vector< Face* > faces;
};


#endif