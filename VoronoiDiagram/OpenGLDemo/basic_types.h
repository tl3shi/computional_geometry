#ifndef BASIC_TYPES_H_
#define BASIC_TYPES_H_

#include <cmath>
#include <limits>

class Point;
class Site;
class Vertex;
class Face;
class Halfedge;

typedef Point Vector;

class Point
{
public:
    Point( void ) { p_[0] = 0.0f; p_[1] = 0.0f; };
    Point( double x, double y ) 
    { p_[0] = x; p_[1] = y; 
    };
    ~Point( void ) {};

    inline double x( void ) const { return p_[0]; };
    inline double y( void ) const { return p_[1]; };
    inline void SetX( double x ) { p_[0] = x; };
    inline void SetY( double y ) { p_[1] = y; };
    inline double* data() { return p_; }

    void operator += ( const Point& p ) { p_[0] += p.x(); p_[1] += p.y(); };
    void operator += ( double c ) { p_[0] += c; p_[1] += c; };
    void operator -= ( const Point& p ) { p_[0] -= p.x(); p_[1] -= p.y(); };
    void operator -= ( double c ) { p_[0] -= c; p_[1] -= c; };
    void operator *= ( double c ) { p_[0] *= c; p_[1] *= c; };
    void operator /= ( double c ) { p_[0] /= c; p_[1] /= c; };
   
    CString toString()
    {
        CString str;
        str.Format(L"(%.2f, %.2f) ", p_[0], p_[1]);
        return str;
    }

    //the direction perpendicular to this, ccw 90 degree
    Vector getPerpendicularVector()
    {
        Vector v =  Vector(-p_[1], p_[0]);
        v.mf_normalize();
        return v;
    }
    double mf_getLength()
    {
        return sqrt(p_[0]*p_[0] + p_[1]*p_[1]); 
    }
    
    void mf_normalize( )
    {
        double a = mf_getLength( );
        (*this) /= a; // 注意: 这里没有处理除数为0的情况
    } //成员函数mf_normalize结束


    inline double& operator []( int index ) { return p_[index]; };

    const Point& operator = ( const Point& p ) { p_[0] = p.x(); p_[1] = p.y(); return p; };

    Point operator - () 
    {
        return Point(-p_[0], -p_[1]);
    }

    //
    friend inline bool operator == ( const Point& p1, const Point& p2 ){

        return ( ( p1.p_[0] == p2.p_[0] ) && ( p1.p_[1] == p2.p_[1] ) );
    };
    //

    //
    friend inline bool operator != ( const Point& p1, const Point& p2 ){

        return ( ( p1.p_[0] != p2.p_[0] ) || ( p1.p_[1] != p2.p_[1] ) );
    };
    //

    //
    friend inline  Point operator + ( const Point& p1, const Point& p2 ){

        return Point( p1.x() + p2.x(), p1.y() + p2.y() );
    };
    //

    // 
    friend inline  Point operator + ( const Point& p, double c ){

        return Point( p.x() + c, p.y() + c );
    };
    //

    //
    friend inline  Point operator - ( const Point& p1, const Point& p2 ){

        return Point( p1.x() - p2.x(), p1.y() - p2.y() );
    };
    //

    // 
    friend inline  Point operator - ( const Point& p, double c ){

        return Point( p.x() - c, p.y() - c );
    };
    //

    // 
    friend inline  Point operator * ( const Point& p, double c ){

        return Point( p.x() * c, p.y() * c );
    };
    //

    //
    friend inline  Point operator / ( const Point& p, double c ){
    
        return Point( p.x() / c, p.y() / c );
    };
    //

private:
    double p_[2];
};

// 点积
static double  operator * (const Point& u, const Point& v)
{
    return u.x() * v.x() + u.y() * v.y();
} //函数operator *结束

static double operator ^ (const Point& u, const Point& v)
{
    return(u.x()* v.y() - u.y()* v.x() );
} //函数operator ^结束


// ! Point

class Site
{

public:
    Site( void ) { x_ = 0.0f; y_ = 0.0f; incFace_ = NULL; };
    Site( double x, double y ) { x_ = x; y_ = y; p.SetX(x); p.SetY(y);incFace_ = NULL; };
    Site( Point* p ) { x_ = p->x(); y_ = p->y(); this->p = *p; incFace_ = NULL; };
    Site( Point p ) { x_ = p.x(); y_ = p.y(); this->p = p; incFace_ = NULL; };
    Site( double x, double y , Face* incFace ) { x_ = x; y_ = y; p.SetX(x); p.SetY(y); incFace_ = incFace; };
    ~Site( void ) {};

    inline double x( void ) const { return x_; };
    inline double y( void ) const { return y_; };
    inline Face* incFace( void ) const { return incFace_; };
    inline void SetX( double x ) { x_ = x;  p.SetX(x);};
    inline void SetY( double y ) { y_ = y;  p.SetY(y);};
    inline void SetFace( Face* incFace ) { incFace_ = incFace; };

public:
    Point p;
private:
    double x_, y_;  //coordinates of this site
    Face* incFace_; //the face that this site belongs to
    
};

class Vertex
{

public:
    Vertex( void ) { x_ = 0.0f; y_ = 0.0f; incEdge_ = NULL; };
    Vertex( double x, double y ) { x_ = x; y_ = y; incEdge_ = NULL; };
    Vertex( Point  p) { x_ = p.x(); y_ = p.y(); this->p = p;incEdge_ = NULL; };
    Vertex( Point*  p) { x_ = p->x(); y_ = p->y(); this->p = *p;incEdge_ = NULL; };
    Vertex( Point*  p, Halfedge* incEdge ) { x_ = p->x(); y_ = p->y(); this->p = *p; incEdge_ = incEdge; };
    Vertex( double x, double y , Halfedge* incEdge ) { x_ = x; y_ = y; p = Point(x, y); incEdge_ = incEdge; };
    ~Vertex( void ) {};

    inline double x( void ) const { return x_; };
    inline double y( void ) const { return y_; };
    inline Halfedge* incEdge( void ) const { return incEdge_; };
    inline void SetX( double x ) { x_ = x; };
    inline void SetY( double y ) { y_ = y; };
    inline void SetIncEdge( Halfedge* incEdge ) { incEdge_ = incEdge; };
public:
    Point p;
private:
    double x_, y_;      //coordinates of v
    Halfedge* incEdge_; //pionter to any outgoing incident halfedge
    
};

class Face
{

public:
    Face( void ) { site_ = NULL; incEdge_ = NULL; };
    Face( Site* site ) { site_ = site; incEdge_ = NULL; };
    Face( Halfedge* incEdge ) { site_ = NULL; incEdge_ = incEdge; };
    Face( Site* site, Halfedge* incEdge ) { site_ = site; incEdge_ = incEdge; };
    ~Face( void ) {};

    inline Site* site( void ) const { return site_; };
    inline Halfedge* incEdge( void ) const { return incEdge_; };
    inline void SetIncEdge( Halfedge* incEdge ) { incEdge_ = incEdge; };
    inline void SetSite( Site* site ) { site_ = site; };

private:
    Site* site_;       //the site of this face
    Halfedge* incEdge_; //any incident halfedge
    
};

class Halfedge
{

public:
    Halfedge( void ) { twinEdge_ = NULL; oriVertex_ = NULL; incFace_ = NULL; prevEdge_ = NULL; nextEdge_ = NULL; 
                       midPoint_ = NULL; direction_ = NULL; endVertex = NULL;};
    Halfedge( Halfedge* twinEdge, Vertex* oriVertex, Face* incFace, Halfedge* prevEdge, Halfedge* nextEdge, 
              Point* midPoint, Vector* direction ) 
            { twinEdge_ = twinEdge; oriVertex_ = oriVertex; incFace_ = incFace; prevEdge_ = prevEdge; nextEdge_ = nextEdge; 
              midPoint_ = midPoint; direction_ = direction; };
    ~Halfedge( void ) {};

    inline Halfedge* twinEdge( void ) const { return twinEdge_; };
    inline Vertex* oriVertex( void ) const { return oriVertex_; };
    inline Face* incFace( void ) const { return incFace_; };
    inline Halfedge* prevEdge( void ) const { return prevEdge_; };
    inline Halfedge* nextEdge( void ) const { return nextEdge_; };
    inline Point* midPoint( void ) const { return midPoint_; };
    inline Vector* direction( void ) const { direction_->mf_normalize(); return direction_; };
    inline void SetTwinEdge( Halfedge* twinEdge ) { twinEdge_ = twinEdge; };
    inline void SetOriVertex( Vertex* oriVertex ) { oriVertex_ = oriVertex; };
    inline void SetIncFace( Face* incFace ) { incFace_ = incFace; };
    inline void SetPrevEdge( Halfedge* prevEdge ) { prevEdge_ = prevEdge; };
    inline void SetNextEdge( Halfedge* nextEdge ) { nextEdge_ = nextEdge; };
    inline void SetMidPoint( Point* midPoint ) { midPoint_ = midPoint; };
    inline void SetDirection( Vector* direction ) { direction_ = direction; };
    inline void SetDirection( Vector direction ) { direction_ = new Vector(direction.x(), direction.y());};
    inline void SetEndVertex(Vertex* vertex){endVertex = vertex;};
    //
    friend inline bool operator == ( const Halfedge& p1, const Halfedge& p2 )
    {
        //cannot use nextedge == nextedge, may be in infinite loop
        return (p1.direction() == p2.direction() && p1.oriVertex()->p == p2.oriVertex()->p
            &&(*p1.midPoint()) == *(p2.midPoint()) );


    };
    //
public:
    Vertex*   endVertex;

private:
    Halfedge* twinEdge_;  //pointer to twin halfedge
    Vertex*   oriVertex_; //pointer to origin vertex
    Face*     incFace_;   //pointer to left incident face
    Halfedge* prevEdge_;  //pointer to CCW previous halfedge
    Halfedge* nextEdge_;  //pointer to CCW next halfedge
    Point*    midPoint_;  //the midpoint of the two sites of this halfedge
    Vector*   direction_; //the direction of this halfedge
   
};


class Line
{
public:
    Point begin;
    Point end;
    Line(Point b, Point e):begin(b),end(e){}
    ~Line( void ) {};
};

#endif
