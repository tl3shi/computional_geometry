#include "stdafx.cpp"
#include "basic_types.h"
#include <stack>
#include <vector>
#ifndef MAX_INT
    #define  MAX_INT 65536
#endif // !MAX_INT
#define  TOLERANCE 0.000001

#include <algorithm>
#include <cmath>

using namespace std;

class AngleComparer
{
    Point leftest_and_lowest;
    //static Site  leftest_and_lowest;
public:
    AngleComparer(Point _param)
    {
        leftest_and_lowest = _param;
    }
    //leftmost, theta in [0, pi], increasing 
    bool operator()(const Point& p1, const Point& p2) const
    {
        if(p1 == leftest_and_lowest) return 1;//leftmost is the min means p1 < p2
        if(p2 == leftest_and_lowest) return 0;//leftmost is the min means p1 > p2

        double arccos_angle1 = (Point(MAX_INT, leftest_and_lowest.y()) - leftest_and_lowest) * (p1 - leftest_and_lowest) / (Point(MAX_INT, leftest_and_lowest.y()) - leftest_and_lowest).mf_getLength() / (p1 - leftest_and_lowest).mf_getLength();
        double arccos_angle2 = (Point(MAX_INT, leftest_and_lowest.y()) - leftest_and_lowest) * (p2 - leftest_and_lowest) / (Point(MAX_INT, leftest_and_lowest.y()) - leftest_and_lowest).mf_getLength() / (p2 - leftest_and_lowest).mf_getLength();
        //cos  decreasing
        return arccos_angle1 > arccos_angle2 ;
    }
 
};

class SiteAngleComparer
{
    Point leftest_and_lowest;
    //static Site  leftest_and_lowest;
public:
    SiteAngleComparer(Point _param)
    {
        leftest_and_lowest = _param;
    }
 
    //leftmost, theta in [0, pi], increasing 
    bool operator()(const Site * p1, const Site * p2) const
    {
        if(p1->p == leftest_and_lowest) return 1;//leftmost is the min means p1 < p2
        if(p2->p == leftest_and_lowest) return 0;//leftmost is the min means p1 > p2

        double arccos_angle1 = (Point(MAX_INT, leftest_and_lowest.y()) - leftest_and_lowest) * (p1->p - leftest_and_lowest) / (Point(MAX_INT, leftest_and_lowest.y()) - leftest_and_lowest).mf_getLength() / (p1->p - leftest_and_lowest).mf_getLength();
        double arccos_angle2 = (Point(MAX_INT, leftest_and_lowest.y()) - leftest_and_lowest) * (p2->p - leftest_and_lowest) / (Point(MAX_INT, leftest_and_lowest.y()) - leftest_and_lowest).mf_getLength() / (p2->p - leftest_and_lowest).mf_getLength();
        //cos  decreasing
        return arccos_angle1 > arccos_angle2 ;
    }

};

class GeometryTool
{
public:

    //test k is in triangle(p,q,r)
    bool static in_triangle(const Point &k, const Point &p, const  Point &q, const  Point &r)
    {
        //this case, p q r must be ccw
        //return to_left(k, p, q) && to_left(k, q, r) && to_left(k, r, p);
        //not mattter ccw or cw
        bool left1 = to_left(k, p, q);
        bool left2 = to_left(k, q, r);
        bool left3 = to_left(k, r, p);
        return (left1 == left2) && (left2 == left3);
    }

    static int compareByY(Point &p0, Point &p1)
    {
        if (p0.y() == p1.y())
            return p0.x() < p1.x();
        return p0.y() < p1.y();
    }

    static int comparePointByX(Point &p0, Point &p1)
    {
        if (p0.x() == p1.x())
            return p0.y() < p1.y();
        return p0.x() < p1.x();
    }

    static int compareSiteByX(Site* &p0, Site* &p1)
    {
        if (p0->x() == p1->x())
            return p0->y() < p1->y();
        return p0->x() < p1->x();
    }

    //test if po is on the left of p2p3
    bool static to_left(const Point &p0, const Point &p1, const Point &p2)
    {
        //p0p1 * p2p1 * sin(theta)
        return ((p2 - p1) ^ (p0 - p1)) > 0;
    }

    bool static to_left(Site* p0,  Site *p1,  Site *p2)
    {
        //p0p1 * p2p1 * sin(theta)
        return ((p2->p - p1->p) ^ (p0->p - p1->p)) > 0;
    }

    bool static isParallel(Vector &v1, Vector &v2)
    {
        v1.mf_normalize();
        v2.mf_normalize();
        return abs(v1 * v2) < 1 + TOLERANCE && abs(v1 * v2) > 1 - TOLERANCE; 
    }

    //return line(p0,p1) intersect with line(p2, p3)
   static Point* intersect(const Point &p0, const Point &p1, const Point &p2, const Point &p3, bool isSegment=false)
    {
        Vector d1 = Vector(p0 - p1);
        Vector d2 = Vector(p2 - p3);
        d1.mf_normalize();
        d2.mf_normalize();
        if (isParallel(d1, d2))
            return NULL;
        Vector p = p2 - p0;
        double delta = d1 ^ d2;
        double parameter_on_l1 = (p ^ d2) / delta;
        double parameter_on_l2 = (p ^ d1) / delta;
        Point * retval = new Point((d1 * parameter_on_l1 + p1 + d2 * parameter_on_l2 + p2) / 2.0);
        return retval;
    }

    //intersection line p1 ,direction d1 
    static Point*  intersectPointVector(const Point &p1,  Vector &d1, const Point &p2, Vector &d2, bool isSegment=false)
    {
        d1.mf_normalize();
        d2.mf_normalize();
        if (isParallel(d1, d2))
            return NULL;
        Vector p = p2 - p1;
        double delta = d1 ^ d2;
        double parameter_on_l1 = (p ^ d2) / delta;
        double parameter_on_l2 = (p ^ d1) / delta;
        Point * retval = new Point((d1 * parameter_on_l1 + p1 + d2 * parameter_on_l2 + p2) / 2.0);
        return retval;
    }
    /*
    //O(nlgn)
    vector<Point> static getConvexHullUseGrahamScan(vector<Point> points)
    {
        int ltl = 0;//find the lowest-and-leftmost point
        for (unsigned int i = 1; i < points.size(); i++)
        {
            if(points[i].y() < points[ltl].y() || 
                (points[i].y() == points[ltl].y() && points[i].x() < points[ltl].x()))
                ltl = i;
        }

        leftest_and_lowest = points[ltl];
        //compareByAngle should use the ltl point
        sort(points.begin(), points.end(), GeometryTool::compareByAngle);

        vector<Point> convexPoints;
        stack<Point> S, T;
        S.push(points[0]);S.push(points[1]);
        for (int i = points.size()-1; i >= 2;i--)
        {
            T.push(points[i]);
        }
        Point s0, s1;
        while (!T.empty())
        {
            s0 = S.top();
            S.pop();//s0 is gone
            if(S.empty())
            {
                S.push(T.top());
                break;
            }
            s1 = S.top();//get s1
            S.push(s0);
            while (!to_left(T.top(), s1, s0))
            { 
                S.pop(); 
                s0 = S.top();
                S.pop();//s0 is gone
                s1 = S.top();//get s1
                S.push(s0);
            }
            S.push(T.top());
            T.pop();
        }
        int size = S.size();
        for(int i = 0; i < size; i++)
        {
            convexPoints.push_back(S.top());
            S.pop();//size is change, cannot use i < S.size()
        }
        return convexPoints;
    }
    */


    //O(nlgn)
    vector<Site*> static getConvexHullUseGrahamScan(vector<Site*> sites)
    {
        int ltl = 0;//find the lowest-and-leftmost point
        for (unsigned int i = 1; i < sites.size(); i++)
        {
            if(sites[i]->y() < sites[ltl]->y() || 
                (sites[i]->y() == sites[ltl]->y() && sites[i]->x() < sites[ltl]->x()))
                ltl = i;
        }

        Point leftest_and_lowest = (*(sites[ltl])).p;

        SiteAngleComparer comparer(leftest_and_lowest);
        //compareByAngle should use the ltl point
        sort(sites.begin(), sites.end(), comparer);

        vector<Site*> convexPoints;
        stack<Site*> S, T;
        S.push(sites[0]);S.push(sites[1]);
        for (int i = sites.size()-1; i >= 2;i--)
        {
            T.push(sites[i]);
        }
        Site* s0 = new Site(), *s1 = new Site();
        while (!T.empty())
        {
            s0 = S.top();
            S.pop();//s0 is gone
            if(S.empty())
            {
                S.push(T.top());
                break;
            }
            s1 = S.top();//get s1
            S.push(s0);
            while (!to_left(T.top(), s1, s0))
            { 
                S.pop(); 
                s0 = S.top();
                S.pop();//s0 is gone
                s1 = S.top();//get s1
                S.push(s0);
            }
            S.push(T.top());
            T.pop();
        }
        int size = S.size();
        for(int i = 0; i < size; i++)
        {
            convexPoints.push_back(S.top());
            S.pop();//size is change, cannot use i < S.size()
        }
        return convexPoints;
    }
};
