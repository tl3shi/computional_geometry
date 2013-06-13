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
    ~VoronoiDiagram( void )
    {
        for (unsigned int i = 0; i < sites.size(); i++)
        {
            delete sites[i];
        }
        
        for (unsigned int i = 0; i < halfedges.size(); i++)
        {
            delete halfedges[i];
        }
        for (unsigned int i = 0; i < convex_hull.size(); i++)
        {
            delete convex_hull[i];
        }
        for (unsigned int i = 0; i < vertices.size(); i++)
        {
            delete vertices[i];
        }
        for (unsigned int i = 0; i < faces.size(); i++)
        {
            delete faces[i];
        }
    }
    
    bool ToLeft( Site* s1, Site* s2, Site* new_site );
    bool IsIntersectionPointOnHalfedge( double x, double y, Halfedge* he );
    Face* LocateInFace( VoronoiDiagram* vd, Point* p );
    bool UpdateConvexHull( Site* new_site );
    Site* FindClosestSite( VoronoiDiagram* vd, Point* p );
    vector<Halfedge*> CalcIntersectionVertices( Halfedge* he, Site* site, bool is_in_convex_hull, Vertex* vbegin, Vertex* vend, Halfedge* he_intersection1, Halfedge* he_intersection_1 );
    void CalcIntersectionVertex( Halfedge* he, Site* site, bool is_in_convex_hull, Vertex* vbegin, Vertex* vend, Halfedge* he_intersection1, Halfedge* he_intersection_1 );
    void GenerateNewHalfedges( Site* new_site, Site* site, Halfedge* he1, Halfedge* he2 );
    VoronoiDiagram* IncrementalConstruction( vector< Point* > &points );
    void PrintVoronoiDiagram( VoronoiDiagram* vd );
    void DeleteUselessHalfedges( VoronoiDiagram* vd );
    VoronoiDiagram* DevideConquerConstruction( vector< Point > &points );
public:
    vector< Site* > sites;
    vector< Site* > convex_hull;
    vector< Vertex* > vertices;
    vector< Halfedge* > halfedges;
    vector< Face* > faces;

    int method; // 0 for IncrementalConstruction, 1 for DevideConquerConstruction
};

#endif