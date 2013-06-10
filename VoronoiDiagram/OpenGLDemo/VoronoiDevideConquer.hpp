#include "stdafx.h"
#include "voronoi_diagram.h"
#include "GeometryTool.hpp"
#include <algorithm>
#include <assert.h>
#include <limits>
#include <iostream>

using namespace std;
#define INFINITE_LENGTH 1000

#define TOLERANCE 0.001

#ifndef Vector
    typedef Point Vector ;
#endif

    Vertex* infiniteVertex = new Vertex( 
        DBL_MAX, 
        DBL_MAX, 
        NULL );
    Point infinitePoint = Point( DBL_MAX, DBL_MAX);

    VoronoiDiagram*   operator + (const VoronoiDiagram &u, const VoronoiDiagram  &v)
    {
        VoronoiDiagram* result = new VoronoiDiagram();
        for (unsigned int i = 0; i < u.faces.size(); i++)
        {
            result->faces.push_back(u.faces[i]);
        }
        for (unsigned int i = 0; i < v.faces.size(); i++)
        {
            result->faces.push_back(v.faces[i]);
        }
        for (unsigned int i = 0; i < u.sites.size(); i++)
        {
            result->sites.push_back(u.sites[i]);
        }
        for (unsigned int i = 0; i < v.sites.size(); i++)
        {
            result->sites.push_back(v.sites[i]);
        }
        for (unsigned int i = 0; i < u.halfedges.size(); i++)
        {
            result->halfedges.push_back(u.halfedges[i]);
        }
        for (unsigned int i = 0; i < v.halfedges.size(); i++)
        {
            result->halfedges.push_back(v.halfedges[i]);
        }
        for (unsigned int i = 0; i < u.faces.size(); i++)
        {
            result->faces.push_back(u.faces[i]);
        }
        for (unsigned int i = 0; i < v.faces.size(); i++)
        {
            if( v.faces[i]->site() == NULL && v.faces[i]->incEdge() == NULL)
                continue;//the whole site has added
            result->faces.push_back(v.faces[i]);
        }

        return result;

    }

VoronoiDiagram* smallVD(vector<Point> &origin, int left, int right)
{
    VoronoiDiagram * result = new VoronoiDiagram();
    if(right - left == 2)
    {
       //three points
       Point p0 = origin.at(left);
       Point p1 = origin.at(left+1);
       Point p2 = origin.at(right);
       
       Site *s0 , *s1, *s2;// ccw
       if(GeometryTool::to_left(p0, p1, p2))
       {
            s0 = new Site(p0);
            s1 = new Site(p1);
            s2 = new Site(p2);
       }else
       {
           s0 = new Site(p0);
           s1 = new Site(p2);
           s2 = new Site(p1);
       }
       Point* mid_01 = new Point((s0->p + s1->p) / 2.0);
       Point* mid_12 = new Point((s1->p + s2->p) / 2.0);
       Point* mid_20 = new Point((s0->p + s2->p) / 2.0);

       //direction_perpendicular
       Vector dir_per_01 = Vector(s1->p - s0->p).getPerpendicularVector();
       Vector dir_per_12 = Vector(s2->p - s1->p).getPerpendicularVector();
       Vector dir_per_20 = Vector(s0->p - s2->p).getPerpendicularVector();
       Point * center = GeometryTool::intersectPointVector(*mid_01, dir_per_01, *mid_12, dir_per_12);
       if(center == NULL) 
       {
        //parallel
           Halfedge *e1, *e2, *e3, *e4;

           e1 = new Halfedge();
           e2 = new Halfedge();
           e1->SetDirection(dir_per_01);
           e1->SetOriVertex(infiniteVertex);

           e2->SetDirection(-dir_per_01);
           e2->SetOriVertex(infiniteVertex);
           e1->SetTwinEdge(e2);
           e2->SetTwinEdge(e1);

           e1->SetMidPoint(mid_01);
           e2->SetMidPoint(mid_01);

           e3 = new Halfedge();
           e4 = new Halfedge();
           e3->SetDirection(dir_per_12);
           e3->SetOriVertex(infiniteVertex);

           e4->SetDirection(-dir_per_12);
           e4->SetOriVertex(infiniteVertex);
           e3->SetTwinEdge(e4);
           e4->SetTwinEdge(e3);
           e3->SetMidPoint(mid_12);
           e4->SetMidPoint(mid_12);

           Face *f0, *f1, *f2;
           f0 = new Face(s0, e1);
           f1 = new Face(s1, e3);
           f2 = new Face(s2, e4);

           s0->SetFace(f0);
           s1->SetFace(f1);
           s2->SetFace(f2);

           result->sites.push_back(s0);
           result->sites.push_back(s1);
           result->sites.push_back(s2);
           result->halfedges.push_back(e1);
           result->halfedges.push_back(e2);
           result->halfedges.push_back(e3);
           result->halfedges.push_back(e4);
           result->faces.push_back(f0);
           result->faces.push_back(f1);
           result->faces.push_back(f2);
           return result;
       }
       Vertex* v_venter = new Vertex(center);

       Halfedge *e1, *e2, *e3, *e4, *e5, *e6;

       e1 = new Halfedge();
       e2 = new Halfedge();
       e1->SetDirection(dir_per_01);
       e1->SetOriVertex(infiniteVertex);

       e2->SetDirection(-dir_per_01);
       e2->SetOriVertex(v_venter);
       e1->SetTwinEdge(e2);
       e2->SetTwinEdge(e1);
       e1->SetMidPoint(mid_01);
       e2->SetMidPoint(mid_01);


       e3 = new Halfedge();
       e4 = new Halfedge();
       e3->SetDirection(dir_per_12);
       e3->SetOriVertex(infiniteVertex);
       e3->SetTwinEdge(e4);
       e4->SetDirection(-dir_per_12);
       e4->SetOriVertex(v_venter);
       e4->SetTwinEdge(e3);
       e3->SetNextEdge(e2);
       e2->SetPrevEdge(e3);

       e3->SetMidPoint(mid_12);
       e4->SetMidPoint(mid_12);


       e5 = new Halfedge();
       e6 = new Halfedge();
       e5->SetDirection(dir_per_20);
       e5->SetOriVertex(infiniteVertex);
       e5->SetTwinEdge(e6);
       e6->SetDirection(-dir_per_20);
       e6->SetOriVertex(v_venter);
       e6->SetTwinEdge(e5);
       e5->SetNextEdge(e4);
       e4->SetPrevEdge(e5);
       
       e6->SetMidPoint(mid_20);
       e5->SetMidPoint(mid_20);

       e1->SetNextEdge(e6);
       e6->SetPrevEdge(e1);


       Face *f0, *f1, *f2;
       f0 = new Face(s0, e1);
       f1 = new Face(s1, e3);
       f2 = new Face(s2, e5);
       
       e1->SetIncFace(f0);
       e2->SetIncFace(f1);
       e3->SetIncFace(f1);
       e4->SetIncFace(f2);
       e5->SetIncFace(f2);
       e6->SetIncFace(f0);

       s0->SetFace(f0);
       s1->SetFace(f1);
       s2->SetFace(f2);

       result->sites.push_back(s0);
       result->sites.push_back(s1);
       result->sites.push_back(s2);
       result->vertices.push_back(v_venter);
       result->halfedges.push_back(e1);
       result->halfedges.push_back(e2);
       result->halfedges.push_back(e3);
       result->halfedges.push_back(e4);
       result->halfedges.push_back(e5);
       result->halfedges.push_back(e6);
       result->faces.push_back(f0);
       result->faces.push_back(f1);
       result->faces.push_back(f2);

       result->convex_hull.push_back(s0);
       result->convex_hull.push_back(s1);
       result->convex_hull.push_back(s2);
    }else if(right - left == 1)
    {
        //two points
        Site* s1 = new Site(origin.at(left));
        Site* s2 = new Site(origin.at(right));
        result->sites.push_back(s1);
        result->sites.push_back(s2);

        Point* mid = new Point((origin.at(left) + origin.at(right)) / 2);
        Vector* direction = new Vector(origin.at(right) - origin.at(left));
        direction->mf_normalize();
        //assert(direction != 0);// direction == 0, TODO
         
        Vector* norDirection = new Vector(direction->getPerpendicularVector());


        Halfedge* e1 = new Halfedge();
        Halfedge* e2 = new Halfedge();

        e1->SetDirection(norDirection);
        e1->SetOriVertex(infiniteVertex);
        e1->SetMidPoint(mid);
        e1->SetTwinEdge(e2);

        e2->SetDirection((-*norDirection));
        e2->SetOriVertex(infiniteVertex);
        e2->SetMidPoint(mid);
        e2->SetTwinEdge(e1);

        result->halfedges.push_back(e1);
        result->halfedges.push_back(e2);
        Face* f1 = new Face();
        Face* f2 = new Face();
        f1->SetSite(s1);
        f2->SetSite(s2);
        f1->SetIncEdge(e1);
        f2->SetIncEdge(e2);

        s1->SetFace(f1);
        s2->SetFace(f2);

        e1->SetIncFace(f1);
        e2->SetIncFace(f2);

        result->faces.push_back(f1);
        result->faces.push_back(f2);

        result->convex_hull.push_back(s1);
        result->convex_hull.push_back(s2);

    }
    return result;
}

 
void tangentLine(vector<Site*> left, vector<Site*> right, Site &leftMax, Site  &leftMin, Site &rightMax, Site  &rightMin)
{
    /*Point leftMax = Point(numeric_limits<double>::min(),numeric_limits<double>::min());
    Point leftMin = Point(numeric_limits<double>::max(),numeric_limits<double>::max());
    Point rightMax = Point(numeric_limits<double>::min(),numeric_limits<double>::min());
    Point rightMin = Point(numeric_limits<double>::max(),numeric_limits<double>::max());*/
    
    double min =  DBL_MAX;
    double max = DBL_MIN;

    for (unsigned int i = 0; i < left.size(); i++)
    {
        Point t = left[i]->p;
        if(t.y() < min)
        {
            leftMin = *left[i];
            min = t.y();
        }
        if(t.y() > max)
        {
            leftMax = *left[i];
            max = t.y();
        }
    }
    
    min =  DBL_MAX;
    max = DBL_MIN;

    for (unsigned int i = 0; i < right.size(); i++)
    {
        Point t = right[i]->p;
        if(t.y() < min)
        {
            rightMin = *right[i];
            min = t.y();
        }
        if(t.y() > max)
        {
            rightMax = *right[i];
            max = t.y();
        }
    }
}
vector<Site*> processConvexHull(vector<Site*> &convexhull)
{
    vector<Site*> result;
    Point p = convexhull.at(0)->p;
    int lowestandleftest = 0;
    for (unsigned int i = 0; i < convexhull.size(); i++)
    {
        Point t = convexhull.at(i)->p;
        if((t.y() < p.y()) ||
            (t.y() == p.y() && t.x() < p.x()))
        {
            lowestandleftest = i;
        }
    }
    for (int k = lowestandleftest; k >= 0; k--)
    {
        result.push_back(convexhull.at(k));
    }
    for (int k = convexhull.size() - 1; k > lowestandleftest; k--)
    {
        result.push_back(convexhull.at(k));
    }

    return result;
}


//set the edge to whose previous edge
void connectLeftEdge(Halfedge *edge, Face* face)
{
    Vertex * end = edge->twinEdge()->oriVertex();
    if(end == infiniteVertex)// the edge's next NULL
        return;


    Halfedge * initialLeft = face->incEdge();
    //set next
    while(initialLeft!= NULL)
    {
        if(initialLeft->oriVertex()->p == end->p)
        {
            edge->SetNextEdge(initialLeft);
            initialLeft->SetPrevEdge(edge);
            return;
        }
        initialLeft = initialLeft->prevEdge();
    }

    initialLeft = face->incEdge();
    while (initialLeft != NULL)
    {
        if(initialLeft->twinEdge()->oriVertex()->p == edge->oriVertex()->p)
        {
            edge->SetNextEdge(initialLeft);
            initialLeft->SetPrevEdge(edge);
            return;
        }
        initialLeft = initialLeft->nextEdge();
    }
}

//set the edge to whose next edge
void connectRightEdge(Halfedge *edge, Face *face)
{
    if(edge->oriVertex() == infiniteVertex)// the edge's prev NULL
        return;

    Halfedge * initialRight = face->incEdge();

    while (initialRight !=NULL )
    {
        if (initialRight->twinEdge()->oriVertex()->p == edge->oriVertex()->p)
        {
            (*initialRight).SetNextEdge(edge);
            edge->SetPrevEdge(initialRight);
            return;
        }
        initialRight = initialRight->nextEdge();
    }

    initialRight = face->incEdge();
    while (initialRight !=NULL )
    {
        if (initialRight->oriVertex()->p == edge->twinEdge()->oriVertex()->p)
        {
            (*initialRight).SetNextEdge(edge);
            edge->SetPrevEdge(initialRight);
            return;
        }
        initialRight = initialRight->prevEdge();
    }
     
}

class DevideChain
{
public:
    bool left;
    Vertex * intersectionV;
    Halfedge * intersectionEdge;
    Halfedge * downEdge;//the downEdge  of the chain, store midpoint,directioin
    DevideChain(bool paraLeft, Vertex * v, Halfedge* l, Halfedge * d):left(paraLeft),intersectionV(v),intersectionEdge(l),downEdge(d){};

    ~DevideChain()
    {
    }
};

void deleteNextEdge(vector<Halfedge*> &edges,  Halfedge * edge)
{
    Halfedge * nextEdge = edge->nextEdge();
    if (nextEdge != NULL)
    {
        vector<Halfedge*>::iterator it;
        it = edges.begin(); 
        while (it!= edges.end())
        {
            if((*nextEdge) == *(*it) || (*nextEdge->twinEdge()) == *(*it))
            {
             it = edges.erase(it);    
            }else
            {
                it++;
            }
        }
        //edge->nextEdge()->SetPrevEdge(NULL);No need, has delete, none can access it
        edge->SetNextEdge(NULL);
        delete nextEdge;
    }

}

void deletePrevEdge(vector<Halfedge*> edges,  Halfedge * edge)
{
    Halfedge * prevEdge = edge->prevEdge();
    if (prevEdge != NULL)
    {
        vector<Halfedge*>::iterator it;
        it = edges.begin(); 
        while (it!= edges.end())
        {
            if((*prevEdge) == *(*it) || (*prevEdge->twinEdge()) == *(*it))
            {
                it = edges.erase(it);    
            }else
            {
                it++;
            }
            //edge->prevEdge()->SetNextEdge(NULL);//no need, none can access it.
            edge->SetPrevEdge(NULL);
            //delete prevEdge;
        }
    }
}

void connectWithChainOld(vector<DevideChain> &devideChain, VoronoiDiagram* left, VoronoiDiagram* right, Site* leftMin, Site * rightMin)
{
    Halfedge * lastE1 = new Halfedge();
    Halfedge * lastE2 = new Halfedge();
    bool lastLeft = false;//flag the last intersection with right or left
    //clip
    for (unsigned int i = 0; i < devideChain.size(); i++)
    {
        DevideChain chain = devideChain.at(i);
        Halfedge * e1 = new Halfedge();
        Halfedge * e2 = new Halfedge();

        Vector * downDirection = chain.downEdge->direction();
        Point * mid = chain.downEdge->midPoint();
        e2->SetDirection(-(*downDirection));
        e1->SetDirection(downDirection);
        e1->SetTwinEdge(e2);
        e2->SetTwinEdge(e1);
        e2->SetMidPoint(mid);
        e1->SetMidPoint(mid);

        e2->SetOriVertex(chain.intersectionV);

        if(i == 0)//the bisector vertical line from infinite
            e1->SetOriVertex(infiniteVertex);
        else
            e1->SetOriVertex(devideChain.at(i-1).intersectionV);



        //the last bisector to add begin
        if(i == devideChain.size()-1)//the last bisector
        {
            {
                if (lastLeft)//both this intersection and the last intersection are with the left sub voronoi diagram
                {
                    e2->SetNextEdge(devideChain.at(i-1).intersectionEdge->twinEdge());
                    devideChain.at(i-1).intersectionEdge->twinEdge()->SetPrevEdge(e2);
                    //lastE2 has set the face,but not set the next edge
                    e1->SetPrevEdge(lastE1);
                    lastE1->SetNextEdge(e1);
                }else
                {
                    e2->SetNextEdge(lastE2);
                    lastE2->SetPrevEdge(e2);
                    lastE2->SetIncFace(leftMin->incFace());

                    e1->SetPrevEdge(lastE1->nextEdge()->twinEdge());
                    lastE1->nextEdge()->twinEdge()->SetNextEdge(e1);
                }
            }

            left->halfedges.push_back(e2);
            right->halfedges.push_back(e1);

            break;
        }
        //the last bisector to end

        Halfedge * edge = chain.intersectionEdge;
        Face *  face = edge->incFace();

        if (chain.left)
        {
            e2->SetIncFace(edge->incFace());
            //before set,should delete whether the edge's next edge 
            //should have the condition: next's orgin is on the left of e1 --- to left is not ok....
            //if(GeometryTool::to_left(Point(edge->nextEdge()->oriVertex()->p.x() - TOLERANCE, edge->nextEdge()->oriVertex()->p.x())
            //                         ,e1->oriVertex()->p, e2->oriVertex()->p))
            //if edge.next.orgin != last last intersection
            if(i>0 && edge->nextEdge() != NULL && edge->nextEdge()->oriVertex()->p != devideChain.at(i-1).intersectionV->p)
                deleteNextEdge(left->halfedges, edge);

            edge->SetNextEdge(e2);//edge.twin.setorgin = intersecionV
            edge->twinEdge()->SetOriVertex(chain.intersectionV);
            e2->SetPrevEdge(edge);

            if(i != 0)
            {
                if (lastLeft)//both this intersection and the last intersection are with the left sub voronoi diagram
                {
                    
                    e2->SetNextEdge(devideChain.at(i-1).intersectionEdge->twinEdge());
                    devideChain.at(i-1).intersectionEdge->twinEdge()->SetPrevEdge(e2);
                    //lastE2 has set the face but not set next edge
                    e1->SetPrevEdge(lastE1);
                    lastE1->SetNextEdge(e1);
                }else
                {
                    e2->SetNextEdge(lastE2);
                    lastE2->SetPrevEdge(e2);
                    lastE2->SetIncFace(e2->incFace());// the same site
                    
                    e1->SetPrevEdge(lastE1->nextEdge()->twinEdge());
                    lastE1->nextEdge()->twinEdge()->SetNextEdge(e1);
                }
            }
            lastLeft = true;
            left->vertices.push_back(chain.intersectionV);
        }else//right 
        {
            e1->SetIncFace(edge->incFace());
            //before e1.set next, should delete the edge.prev
            if(i>0 && edge->prevEdge() != NULL && edge->prevEdge()->twinEdge()->oriVertex()->p != devideChain.at(i-1).intersectionV->p)
                deletePrevEdge(right->halfedges, edge);

            edge->SetPrevEdge(e1);
            e1->SetNextEdge(edge);//edge .set orivertex to e1.ending
            edge->SetOriVertex(chain.intersectionV);

            if(i != 0)
            {
                if(!lastLeft)//both this intersection and last intersection are with the right sub voronoi diagram
                {
                    e1->SetPrevEdge(devideChain.at(i-1).intersectionEdge->twinEdge());
                    devideChain.at(i-1).intersectionEdge->twinEdge()->SetNextEdge(e1);
                    //lastE1->SetIncFace(edge->twinEdge()->incFace());//should be the opposite site, infact has been set face above
                    e2->SetNextEdge(lastE2);
                    lastE2->SetPrevEdge(e2);
                }else
                {
                    e1->SetPrevEdge(lastE1);
                    lastE1->SetNextEdge(e1);
                    lastE1->SetIncFace(e1->incFace());// the same site

                    e2->SetNextEdge(lastE2->prevEdge()->twinEdge());
                    lastE2->prevEdge()->twinEdge()->SetPrevEdge(e2);
                }
            }
            lastLeft = false;
            right->vertices.push_back(chain.intersectionV);
        }
 
        left->halfedges.push_back(e2);
        right->halfedges.push_back(e1);
        lastE1 = e1;
        lastE2 = e2;
    }
}
bool checkHalfedge(Halfedge * edge)
{
    return false;
};
bool checkDCEL(vector<Halfedge*> &edges)
{
    for (unsigned int i = 0; i < edges.size(); i++)
    {
        Halfedge *edge = edges.at(i);
        if(edge->incFace() == NULL)
            return false;
        if(edge->twinEdge() == NULL)
            return false;
        if(edge->twinEdge()->incFace() == NULL)
            return false;
        if(edge->nextEdge() != NULL)
            if (!(edge->nextEdge()->oriVertex()->p == edge->twinEdge()->oriVertex()->p)) return false;
        if(edge->prevEdge() != NULL)
            if (!(edge->oriVertex()->p == edge->prevEdge()->twinEdge()->oriVertex()->p)) return false;

    }
    return true;
}

void connectWithChain(vector<DevideChain> &devideChain, VoronoiDiagram* left, VoronoiDiagram* right, Site* leftMin, Site * rightMin)
{
    connectWithChainOld(devideChain, left, right, leftMin, rightMin);
   
   
}



VoronoiDiagram* mergeVD(VoronoiDiagram* left, VoronoiDiagram* right)
{
     //check every edge has face....for further use in change leftmax and rightmax
    assert(checkDCEL(left->halfedges) == true);
    assert(checkDCEL(right->halfedges) == true);

    VoronoiDiagram * result = new VoronoiDiagram();
    vector<DevideChain> devideChain;

    //the convex hull is cw
    //left->convex_hull = GeometryTool::getConvexHullUseGrahamScan(left->sites);
    //right->convex_hull = GeometryTool::getConvexHullUseGrahamScan(right->sites);
    //process the convex hull to ccw, and the first one is lowest then leftest
    
    Site *leftMax = new Site();
    Site *leftMin = new Site();
    Site *rightMax = new Site();
    Site *rightMin = new Site(); 
    tangentLine(left->sites, right->sites, *leftMax,  *leftMin,  *rightMax, *rightMin);
   
    bool tobottom = false;
    //not zig zag, but the site opposite
    do
    {
        tobottom = ((*leftMax).p == (*leftMin).p && (*rightMax).p == (*rightMin).p);

        Point leftp = leftMax->p;
        Point rightp = rightMax->p;

        Point* mid = new Point(( leftp + rightp ) / 2.0);
        Vector d = Vector(rightp - leftp);
        d.mf_normalize();//direction:point to the up
        d = -d.getPerpendicularVector();//direction:point to the down 
        
        Vertex* lastV;
        if(devideChain.size() > 0)
             lastV = devideChain.at(devideChain.size()-1).intersectionV;
        else
            lastV = infiniteVertex;
        //the intersection should below the mid point

        Point* leftIntersectP = new Point(DBL_MIN, DBL_MIN);
 
        Halfedge * leftIntersectionEdge;
        Halfedge * rightIntersectionEdge;
         
        Halfedge* initalEdge = leftMax->incFace()->incEdge();
        Halfedge*  edge;
        edge = initalEdge;
        bool hasDetectAll = false;
        while(edge != NULL)
        {
            //to get the intersection halfege with ccw 
            //Point edge_p = left->halfedges[i]->oriVertex()->p;//the point may be infinite

            Point* edge_p = edge->midPoint();//*(left->halfedges[i]->midPoint());
            Vector* edge_d = edge->direction();//*(left->halfedges[i]->direction());
            
            Point * t = GeometryTool::intersectPointVector(*mid, d, *edge_p, *edge_d);
            //to get a higher t
            if (t!= NULL && t->y() > leftIntersectP->y() + TOLERANCE 
                && t->y() + TOLERANCE < lastV->y())//ensure the intersection point is the NEW highest 
            {
                leftIntersectP = t;
                leftIntersectionEdge = edge;
            }
            edge = edge->nextEdge();
            if(edge != NULL && *edge == *initalEdge)
            {
                hasDetectAll = true;
                break;
            }
        }
         
        //if there is..but no new intersection s.t condition ok..the edge is changed...
        edge = initalEdge->prevEdge();
        while(!hasDetectAll && edge != NULL)//prev 
        {
            Point* edge_p = edge->midPoint();  
            Vector* edge_d = edge->direction();  

            Point * t = GeometryTool::intersectPointVector(*mid, d, *edge_p, *edge_d);
            if (t!= NULL && t->y() > leftIntersectP->y() + TOLERANCE 
                 && t->y() + TOLERANCE < lastV->y())
            {
                leftIntersectP = t;
                leftIntersectionEdge = edge;
            }
            edge = edge->prevEdge();
        }
         

        Point* rightIntersectP = new Point(DBL_MIN, DBL_MIN);

        initalEdge = rightMax->incFace()->incEdge();
        edge = initalEdge;
        hasDetectAll = false;
        while(edge != NULL)
        {
            //no need all , twin no need, only need the site's face's edges.
            //Point edge_p = right->halfedges[i]->oriVertex()->p;//the p may be infinte
            Point* edge_p = edge->midPoint(); //*(right->halfedges[i]->midPoint());
            Vector* edge_d = edge->direction(); //*(right->halfedges[i]->direction());

            Point * t = GeometryTool::intersectPointVector(*mid, d, *edge_p, *edge_d);
            if (t!= NULL && t->y() > rightIntersectP->y() + TOLERANCE 
                &&t->y() + TOLERANCE < lastV->y())
            {
                rightIntersectP = t;
                rightIntersectionEdge = edge;
            }
            edge = edge->nextEdge();
            if(edge != NULL && (*edge) == *(initalEdge))
            {
                break;
                hasDetectAll = true;
            }
        }
        
        edge = initalEdge->prevEdge();
        while(!hasDetectAll && edge != NULL)//prev 
        {
            Point* edge_p = edge->midPoint();  
            Vector* edge_d = edge->direction();  

            Point * t = GeometryTool::intersectPointVector(*mid, d, *edge_p, *edge_d);
            if (t!= NULL && t->y() > rightIntersectP->y() + TOLERANCE 
                 && t->y() + TOLERANCE < lastV->y())
            {
                rightIntersectP = t;
                rightIntersectionEdge = initalEdge;
            }
            edge = edge->prevEdge();
        }
        

        Halfedge * downEdge = new Halfedge();
        downEdge->SetDirection(d);
        downEdge->SetMidPoint(mid);

        //intersect with left vorinoi first
        if(leftIntersectP->y() > rightIntersectP->y() + TOLERANCE)
        {

            Vertex * v = new Vertex(leftIntersectP);
            v->SetIncEdge(leftIntersectionEdge);

            if(!tobottom)
                leftMax = leftIntersectionEdge->twinEdge()->incFace()->site();

            devideChain.push_back(DevideChain(true, v, leftIntersectionEdge, downEdge));

        }else if (leftIntersectP->y() + TOLERANCE < rightIntersectP->y() )//intersect with right voronoi first
        {
            Vertex * v = new Vertex(rightIntersectP);
            v->SetIncEdge(rightIntersectionEdge);

            if(!tobottom)
                rightMax = rightIntersectionEdge->twinEdge()->incFace()->site();

            devideChain.push_back(DevideChain(false, v, rightIntersectionEdge, downEdge));
            //clicp update
        }else// no intersection below the last intersection point, 
        {
            tobottom = true;

            devideChain.push_back(DevideChain(false, infiniteVertex, NULL, downEdge));
 
        }
        //chain to the lower bound
    }while(!tobottom); 

    connectWithChain(devideChain, left, right, leftMin, rightMin);
    //left right to result
    result = (*right) + (*left);
    return result;
}

VoronoiDiagram* partition(vector<Point> &origin, int left, int right)
{
    if(right - left < 3)
        return smallVD(origin, left, right);
    else
        return mergeVD(partition(origin, left, (left+right)/2), partition(origin, (left+right)/2 + 1, right));
}   

VoronoiDiagram* VoronoiDiagram::DevideConquerConstruction( vector< Point > &points )
{
    //remove TOLERANCE same points should.
    sort(points.begin(), points.end(), &GeometryTool::compareByX);

    return partition(points, 0, points.size()-1);

}