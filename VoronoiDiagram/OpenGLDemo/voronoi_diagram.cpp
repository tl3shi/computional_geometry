#include "stdafx.h"

#include <limits>

#include "basic_types.h"
#include "voronoi_diagram.h"

VoronoiDiagram::VoronoiDiagram( void ) {
    sites = vector< Site* >();
    vertices = vector< Vertex* >();
    halfedges = vector< Halfedge* >();
    faces = vector< Face* >();

    vertices.push_back( new Vertex( 
        DBL_MAX , 
         DBL_MAX , 
        NULL ) );
    faces.push_back( new Face( NULL, NULL ) );
}

VoronoiDiagram::~VoronoiDiagram( void ) {

}

bool VoronoiDiagram::ToLeft( Site* s1, Site* s2, Site* new_site ) {
    double temp = ( s2->x() - s1->x() ) * ( new_site->y() - s1->y() ) - 
                  ( new_site->x() - s1->x() ) * ( s2->y() - s1->y() );
    if ( temp >= 0.0 ) { // if new_site lies on the line ( s1, s2 ), I say it's on the left side.
        return true;
    } else {
        return false;
    }
}

bool VoronoiDiagram::IsIntersectionPointOnHalfedge( double x, double y, Halfedge* he ) {
    if ( he->oriVertex()->x() == DBL_MAX &&
         he->nextEdge()->oriVertex()->x() == DBL_MAX ) {
        // he is a straight line
        return true;
    } else if ( he->oriVertex()->x() == DBL_MAX ) {
        // he is a ray, he's start vertex is an infinite point
        if ( ( x - he->nextEdge()->oriVertex()->x() ) * ( -he->direction()->x() ) + 
            ( y - he->nextEdge()->oriVertex()->y() ) * ( -he->direction()->y() ) > 0 ) {
            return true;
        } else {
            return false;
        }
    } else if ( he->nextEdge()->oriVertex()->x() == DBL_MAX ) {
        // he is a ray, he's end vertex is an infinite point
        if ( ( x - he->oriVertex()->x() ) * he->direction()->x() + 
            ( y - he->oriVertex()->y() ) * he->direction()->y() > 0 ) {
            return true;
        } else {
            return false;
        }
    } else {
        // he is a segment
        if ( ( x - he->nextEdge()->oriVertex()->x() ) * ( -he->direction()->x() ) + 
            ( y - he->nextEdge()->oriVertex()->y() ) * ( -he->direction()->y() ) > 0 && 
            ( x - he->oriVertex()->x() ) * he->direction()->x() + 
            ( y - he->oriVertex()->y() ) * he->direction()->y() > 0 ) {
                return true;
        } else {
            return false;
        }
    }
}

Face* VoronoiDiagram::LocateInFace( VoronoiDiagram* vd, Point* p ) {
    Site* closest_site = FindClosestSite( vd, p );
    return closest_site->incFace();
}

bool VoronoiDiagram::UpdateConvexHull( Site* new_site ) {
    if ( this->sites.size() <= 1 ) {
        this->convex_hull.push_back( new_site );
        return true;
    } else {
        vector< Site* >::iterator iter = this->convex_hull.begin();
        int i;
        bool is_inserted = false;
        for ( i = 0; i < this->convex_hull.size(); i++ ) {
            if ( !ToLeft( *( iter + i ), *( iter + (i==this->convex_hull.size()-1?0:i) ), new_site ) ) {
                // indicate where this new site should be inserted in the convex hull site vector
                // however, after inserting the new site into this positon, 
                // the convex hull list site should be checeked because
                // some sites may become non-convex hull sites after the insertion.
                is_inserted = true;
                iter = this->convex_hull.insert( iter + i + 1, new_site );
            }
        }
        if ( is_inserted ) {
            while ( true ) {
                vector< Site* >::iterator s1 = iter;
                vector< Site* >::iterator s2;
                vector< Site* >::iterator n_s;
                if ( iter + 1 == this->convex_hull.end() ) {
                    s2 = this->convex_hull.begin() + 1;
                    n_s = this->convex_hull.begin();
                } else if ( iter + 2 == this->convex_hull.end() ) {
                    s2 = this->convex_hull.begin();
                    n_s = iter + 1;
                } else {
                    s2 = iter + 2 ;
                    n_s = iter + 1 ;
                }
                if ( ToLeft( *s1, *s2, *n_s ) ) {
                    this->convex_hull.erase( n_s );
                } else {
                    break;
                }
            }
            while ( true ) {
                vector< Site* >::iterator s1 = iter;
                vector< Site* >::iterator s2;
                vector< Site* >::iterator n_s;
                if ( iter == this->convex_hull.begin() ) {
                    s2 = this->convex_hull.end() - 2;
                    n_s = this->convex_hull.end() - 1;
                } else if ( iter - 1 == this->convex_hull.begin() ) {
                    s2 = this->convex_hull.end() - 1;
                    n_s = iter - 1;
                } else {
                    s2 = iter - 2 ;
                    n_s = iter - 1 ;
                }
                if ( !ToLeft( *s1, *s2, *n_s ) ) {
                    this->convex_hull.erase( n_s );
                } else {
                    break;
                }
            }
        }
        return is_inserted;
    }
}

Site* VoronoiDiagram::FindClosestSite( VoronoiDiagram* vd, Point* p ) {
    if ( vd->sites.size() == 0 ) {
        return NULL;
    } else {
        vector< Site* >::iterator iter = vd->sites.begin();
        double min_distance = (*iter)->x() * p->x() + (*iter)->y() * p->y();
        Site* closest_site = *iter;
        iter++;
        for ( ; iter != vd->sites.end(); iter++ ) {
            double temp_distance = (*iter)->x() * p->x() + (*iter)->y() * p->y();
            if ( temp_distance < min_distance ) {
                closest_site = *iter;
                min_distance = temp_distance;
            }
        }
        return closest_site;
    }
}

vector<Halfedge*> VoronoiDiagram::CalcIntersectionVertices( Halfedge* new_he, Site* site, bool is_in_convex_hull, Vertex* vbegin, Vertex* vend, Halfedge* he_intersection1, Halfedge* he_intersection_1 ) {
    vector<Halfedge*> result;
    Halfedge* hebegin = site->incFace()->incEdge();
    double a0, b0, c0, a1, b1, c1, d, x, y;
    a1 = -new_he->direction()->y();
    b1 = new_he->direction()->x();
    c1 = new_he->midPoint()->x() * ( new_he->midPoint()->y() + new_he->direction()->y() ) - 
        new_he->midPoint()->y() * ( new_he->midPoint()->x() + new_he->direction()->x() );
    if ( hebegin->nextEdge() == hebegin ) { // this face has only one halfedge
        a0 = -hebegin->direction()->y();
        b0 = hebegin->direction()->x();
        c0 = hebegin->midPoint()->x() * ( hebegin->midPoint()->y() + hebegin->direction()->y() ) - 
                    hebegin->midPoint()->y() * ( hebegin->midPoint()->x() + hebegin->direction()->x() );
        d = a0 * b1 - a1 * b0;
        if ( d > -0.0000001 && d < 0.0000001 ) { // new halfedge and the old one are parallel
            vbegin->SetX( DBL_MAX );
            vbegin->SetY( DBL_MAX );
            vend->SetX( DBL_MAX );
            vend->SetY( DBL_MAX );
            result.push_back(he_intersection1);
            result.push_back(he_intersection_1);
            return result;
        } else { // new halfedge and the old one are not parallel
            x = ( b0 * c1 - b1 * c0 ) / d;
            y = ( a1 * c0 - a0 * c1 ) / d;
            if ( ( x - new_he->midPoint()->x() ) * new_he->direction()->x() + 
                ( y - new_he->midPoint()->y() ) * new_he->direction()->y() > 0 ) {
                    // ( x, y ) lies on the direction of new_he, it's the end vertex
                    vbegin->SetX( DBL_MAX );
                    vbegin->SetY( DBL_MAX );
                    vend->SetX( x );
                    vend->SetY( y );
                    he_intersection1 = hebegin;
                    result.push_back(he_intersection1);
                    result.push_back(he_intersection_1);
                    return result;
            } else {
                // it's the begin vertex
                vbegin->SetX( x );
                vbegin->SetY( y );
                vend->SetX( DBL_MAX );
                vend->SetY( DBL_MAX );
                he_intersection_1 = hebegin;
                result.push_back(he_intersection1);
                result.push_back(he_intersection_1);
                return result;
            }
        }
    } else { // this face has more than one halfedge
        if ( is_in_convex_hull ) { // new halfedge must have two intersection points
            Halfedge* temp_he = hebegin;
            bool first_vertex = true;
            while ( true ) {
                a0 = -temp_he->direction()->y();
                b0 = temp_he->direction()->x();
                c0 = temp_he->midPoint()->x() * ( temp_he->midPoint()->y() + temp_he->direction()->y() ) - 
                    temp_he->midPoint()->y() * ( temp_he->midPoint()->x() + temp_he->direction()->x() );
                d = a0 * b1 - a1 * b0;
                if ( !( d > -0.0000001 && d < 0.0000001 ) ) { // not parallel
                    x = ( b0 * c1 - b1 * c0 ) / d;
                    y = ( a1 * c0 - a0 * c1 ) / d;
                    if ( !IsIntersectionPointOnHalfedge( x, y, temp_he ) ) {
                        temp_he = temp_he->nextEdge();
                        continue;
                    }
                    if ( ( x - new_he->midPoint()->x() ) * new_he->direction()->x() + 
                        ( y - new_he->midPoint()->y() ) * new_he->direction()->y() > 0 ) {
                            vend->SetX( x );
                            vend->SetY( y );
                            he_intersection_1 = temp_he;
                            if ( !first_vertex ) {
                                result.push_back(he_intersection1);
                                result.push_back(he_intersection_1);
                                return result;
                            } else {
                                first_vertex = false;
                            }
                    } else {
                        vbegin->SetX( x );
                        vbegin->SetY( y );
                        he_intersection1 = temp_he;
                        if ( !first_vertex ) {
                            result.push_back(he_intersection1);
                            result.push_back(he_intersection_1);
                            return result;
                        } else {
                            first_vertex = false;
                        }
                    }
                }
                temp_he = temp_he->nextEdge();
                if ( temp_he == hebegin ) {
                    break;
                }
            }
        } else { // new halfedge may have one intersection point and one infinite point
            Halfedge* temp_he = hebegin;
            bool first_vertex = true;
            bool find_begin_vertex = false;
            bool find_end_vertex = false;
            while ( true ) {
                a0 = -temp_he->direction()->y();
                b0 = temp_he->direction()->x();
                c0 = temp_he->midPoint()->x() * ( temp_he->midPoint()->y() + temp_he->direction()->y() ) - 
                    temp_he->midPoint()->y() * ( temp_he->midPoint()->x() + temp_he->direction()->x() );
                d = a0 * b1 - a1 * b0;
                if ( !( d > -0.0000001 && d < 0.0000001 ) ) { // not parallel
                    x = ( b0 * c1 - b1 * c0 ) / d;
                    y = ( a1 * c0 - a0 * c1 ) / d;
                    if ( !IsIntersectionPointOnHalfedge( x, y, temp_he ) ) {
                        temp_he = temp_he->nextEdge();
                        continue;
                    }
                    if ( ( x - new_he->midPoint()->x() ) * new_he->direction()->x() + 
                        ( y - new_he->midPoint()->y() ) * new_he->direction()->y() > 0 ) {
                            vend->SetX( x );
                            vend->SetY( y );
                            he_intersection1 = temp_he;
                            find_end_vertex = true;
                            if ( !first_vertex ) {
                                result.push_back(he_intersection1);
                                result.push_back(he_intersection_1);
                                return result;
                            } else {
                                first_vertex = false;
                            }
                    } else {
                        vbegin->SetX( x );
                        vbegin->SetY( y );
                        he_intersection_1 = temp_he;
                        find_begin_vertex = true;
                        if ( !first_vertex ) {
                            result.push_back(he_intersection1);
                            result.push_back(he_intersection_1);
                            return result;
                        } else {
                            first_vertex = false;
                        }
                    }
                }
                temp_he = temp_he->nextEdge();
                if ( temp_he == hebegin ) {
                    break;
                }
            }
            if ( !find_begin_vertex ) {
                vbegin->SetX( DBL_MAX );
                vbegin->SetY( DBL_MAX );
                result.push_back(he_intersection1);
                result.push_back(he_intersection_1);
                return result;
            }
            if ( !find_end_vertex ) {
                vend->SetX( DBL_MAX );
                vend->SetY( DBL_MAX );
                result.push_back(he_intersection1);
                result.push_back(he_intersection_1);
                return result;
            }
        }
    }
}

VoronoiDiagram* VoronoiDiagram::IncrementalConstruction( vector< Point* > &points ) {
    VoronoiDiagram* vd = new VoronoiDiagram();
    vector< Point* >::iterator iter;
    for( iter = points.begin(); iter != points.end(); iter++) {
        if ( vd->sites.size() == 0 ) {
            // adding the first site, it takes the whole plane
            Site* s = new Site( (*iter)->x(), (*iter)->y(), *(vd->faces.begin()) );
            vd->sites.push_back( s );
            (*(vd->faces.begin()))->SetSite( s );
        } else if ( vd->sites.size() == 1 ) {
            // adding the second site, generate two halfedge (line)
            Halfedge* he1 = new Halfedge();
            Halfedge* he2 = new Halfedge();
            Face* f = new Face;
            Site* s = new Site( (*iter)->x(), (*iter)->y(), f );

            he1->SetOriVertex( *(vd->vertices.begin()) );
            (*(vd->vertices.begin()))->SetIncEdge( he1 );
            he1->SetIncFace( *(vd->faces.begin()) );
            (*(vd->faces.begin()))->SetIncEdge( he1 );
            he1->SetPrevEdge( he1 );
            he1->SetNextEdge( he1 );
            he1->SetTwinEdge( he2 );

            he2->SetOriVertex( *(vd->vertices.begin()) );
            he2->SetIncFace( f );
            f->SetIncEdge( he2 );
            he2->SetPrevEdge( he2 );
            he2->SetNextEdge( he2 );
            he2->SetTwinEdge( he1 );

            f->SetSite( s );

            Point* midp = new Point( ( he1->incFace()->site()->x() + he2->incFace()->site()->x() ) / 2.0, 
                                     ( he1->incFace()->site()->y() + he2->incFace()->site()->y() ) / 2.0 );
            Vector* v1 = new Vector( he1->incFace()->site()->y() - he2->incFace()->site()->y(), 
                                     -he1->incFace()->site()->x() + he2->incFace()->site()->x() );
            Vector* v2 = new Vector( -he1->incFace()->site()->y() + he2->incFace()->site()->y(), 
                                     he1->incFace()->site()->x() - he2->incFace()->site()->x() );
            he1->SetMidPoint( midp );
            he2->SetMidPoint( midp );
            he1->SetDirection( v1 );
            he2->SetDirection( v2 );

            vd->sites.push_back( s );
            vd->faces.push_back( f );
            vd->halfedges.push_back( he1 );
            vd->halfedges.push_back( he2 );
        } else {
            // adding the third and after site
            Site* new_site = new Site( *iter );
            Face* new_face = new Face( new_site );
            new_site->SetFace( new_face );
            vd->sites.push_back( new_site );
            vd->faces.push_back( new_face );

            Site* closest_site = FindClosestSite( vd, *iter );

            Halfedge* he1 = new Halfedge();
            Halfedge* he_1 = new Halfedge();
            he1->SetIncFace( new_face );
            he_1->SetIncFace( closest_site->incFace() );

            Point* midp1 = new Point( ( he1->incFace()->site()->x() + he_1->incFace()->site()->x() ) / 2.0, 
                ( he1->incFace()->site()->y() + he_1->incFace()->site()->y() ) / 2.0 );
            Vector* dir1 = new Vector( he1->incFace()->site()->y() - he_1->incFace()->site()->y(), 
                -he1->incFace()->site()->x() + he_1->incFace()->site()->x() );
            Vector* dir_1 = new Vector( -he1->incFace()->site()->y() + he_1->incFace()->site()->y(), 
                he1->incFace()->site()->x() - he_1->incFace()->site()->x() );
            he1->SetMidPoint( midp1 );
            he_1->SetMidPoint( midp1 );
            he1->SetDirection( dir1 );
            he_1->SetDirection( dir_1 );

            Vertex* v1 = new Vertex();
            Vertex* v_1 = new Vertex();

            Halfedge* he_intersection1 = NULL;
            Halfedge* he_intersection_1 = NULL;

            bool is_on_convex_hull = UpdateConvexHull( new_site );
            vector<Halfedge*> result = CalcIntersectionVertices( he1, closest_site, !is_on_convex_hull, v1, v_1, 
                he_intersection1, he_intersection_1 );
            he_intersection1 = result[0];
            he_intersection_1 = result[1];

            he1->SetOriVertex( v1 );
            he_1->SetOriVertex( v_1 );
            he1->SetTwinEdge( he_1 );
            he_1->SetTwinEdge( he1 );
            v1->SetIncEdge( he1 );
            v_1->SetIncEdge( he_1 );

            if ( he_intersection1 == NULL && he_intersection_1 == NULL ) {

            } else if ( he_intersection1 != NULL && he_intersection_1 == NULL ) {
                Halfedge* tmp = he_intersection1->nextEdge();
                while ( tmp != he_intersection_1 ) {
                    tmp->SetPrevEdge( NULL );
                    tmp->SetNextEdge( NULL );
                    tmp = tmp->nextEdge();
                }
                he_intersection1->nextEdge()->SetPrevEdge( he_1 );
                he_1->SetNextEdge( he_intersection1->nextEdge() );
                he_intersection1->SetNextEdge( he_1 );
                he_1->SetPrevEdge( he_intersection1 );
            } else if ( he_intersection1 == NULL && he_intersection_1 != NULL ) {
                Halfedge* tmp = he_intersection_1;
                while ( tmp->oriVertex()->x() != DBL_MAX ) {
                }
                he_intersection_1->prevEdge()->SetNextEdge( he_1 );
                he_1->SetPrevEdge( he_intersection_1->prevEdge() );
                he_intersection_1->SetPrevEdge( he_1 );
                he_1->SetNextEdge( he_intersection_1 );
            } else {
                Halfedge* tmp = he_intersection1->nextEdge();
                while ( tmp != he_intersection_1 ) {
                    tmp->SetPrevEdge( NULL );
                    tmp->SetNextEdge( NULL );
                    tmp = tmp->nextEdge();
                }
                he_intersection1->SetNextEdge( he_1 );
                he_1->SetPrevEdge( he_intersection1 );
                he_intersection_1->SetPrevEdge( he_1 );
                he_1->SetNextEdge( he_intersection_1 );
            }

            vd->halfedges.push_back( he1 );
            vd->halfedges.push_back( he_1 );
            vd->vertices.push_back( v1 );
            vd->vertices.push_back( v_1 );

            Halfedge* first_he = he1;
            Halfedge* first_he_ = he_1;

            if ( !is_on_convex_hull ) { // in convex hull
                Site* next_site;
                Halfedge* he2;
                Halfedge* he_2;
                Halfedge* he_intersection2;
                Halfedge* he_intersection_2;

                while ( true ) {
                    next_site = he_intersection1->twinEdge()->incFace()->site();

                    he2 = new Halfedge();
                    he_2 = new Halfedge();
                    he2->SetIncFace( new_face );
                    he_2->SetIncFace( next_site->incFace() );

                    Point* midp2 = new Point( ( he2->incFace()->site()->x() + he_2->incFace()->site()->x() ) / 2.0, 
                        ( he2->incFace()->site()->y() + he_2->incFace()->site()->y() ) / 2.0 );
                    Vector* dir2 = new Vector( he2->incFace()->site()->y() - he_2->incFace()->site()->y(), 
                        -he2->incFace()->site()->x() + he_2->incFace()->site()->x() );
                    Vector* dir_2 = new Vector( -he2->incFace()->site()->y() + he_2->incFace()->site()->y(), 
                        he2->incFace()->site()->x() - he_2->incFace()->site()->x() );
                    he2->SetMidPoint( midp2 );
                    he_2->SetMidPoint( midp2 );
                    he2->SetDirection( dir2 );
                    he_2->SetDirection( dir_2 );

                    Vertex* v2 = v_1;
                    Vertex* v_2 = new Vertex();

                    he_intersection2 = NULL;
                    he_intersection_2 = he_intersection1->twinEdge();

                    he1->SetNextEdge( he2 );
                    he2->SetPrevEdge( he1 );
                    he2->SetTwinEdge( he_2 );
                    he_2->SetTwinEdge( he2 );
                    v_2->SetIncEdge( he_2 );


                    //bool is_on_convex_hull = UpdateConvexHull( new_site );
                    //CalcIntersectionVertex( he2, next_site, !is_on_convex_hull, v2, v_2, 
                    //    he_intersection2, he_intersection_2 );
                    vector<Halfedge*> result = CalcIntersectionVertices( he2, next_site, !is_on_convex_hull, v2, v_2, 
                        he_intersection2, he_intersection_2 );
                    he_intersection2 = result[0];
                    he_intersection_2 = result[1];

                    if ( he_intersection2 != NULL ) {
                        Halfedge* tmp = he_intersection2->nextEdge();
                        while ( tmp != he_intersection_2 ) {
                            tmp->SetPrevEdge( NULL );
                            tmp->SetNextEdge( NULL );
                            tmp = tmp->nextEdge();
                        }
                        he_intersection2->SetNextEdge( he_2 );
                    }

                    if ( he_intersection_2 != NULL ) {
                        he_intersection_2->SetPrevEdge( he_2 );
                    }

                    he2->SetOriVertex( v2 );
                    he_2->SetOriVertex( v_2 );

                    vd->halfedges.push_back( he2 );
                    vd->halfedges.push_back( he_2 );
                    vd->vertices.push_back( v_2 );

                    if ( he_intersection2->twinEdge()->incFace()->site() == closest_site ) {
                        he2->SetNextEdge( first_he );
                        first_he->SetPrevEdge( he2 );
                        he_2->SetPrevEdge( first_he_ );
                        first_he_->SetNextEdge( he_2 );
                        break;
                    } else {
                        he1 = he2;
                        he_1 = he_2;
                        he_intersection1 = he_intersection_2;
                        he_intersection_1 = he_intersection_1;
                        continue;
                    }
                }
            } else { // not in convex hull
                if ( he_intersection_1 == NULL ) {
                    Site* next_site;
                    Halfedge* he2;
                    Halfedge* he_2;
                    Halfedge* he_intersection2;
                    Halfedge* he_intersection_2;

                    while ( true ) {
                        next_site = he_intersection1->twinEdge()->incFace()->site();

                        he2 = new Halfedge();
                        he_2 = new Halfedge();
                        he2->SetIncFace( new_face );
                        he_2->SetIncFace( next_site->incFace() );

                        Point* midp2 = new Point( ( he2->incFace()->site()->x() + he_2->incFace()->site()->x() ) / 2.0, 
                            ( he2->incFace()->site()->y() + he_2->incFace()->site()->y() ) / 2.0 );
                        Vector* dir2 = new Vector( he2->incFace()->site()->y() - he_2->incFace()->site()->y(), 
                            -he2->incFace()->site()->x() + he_2->incFace()->site()->x() );
                        Vector* dir_2 = new Vector( -he2->incFace()->site()->y() + he_2->incFace()->site()->y(), 
                            he2->incFace()->site()->x() - he_2->incFace()->site()->x() );
                        he2->SetMidPoint( midp2 );
                        he_2->SetMidPoint( midp2 );
                        he2->SetDirection( dir2 );
                        he_2->SetDirection( dir_2 );

                        Vertex* v2 = v_1;
                        Vertex* v_2 = new Vertex();

                        he_intersection2 = NULL;
                        he_intersection_2 = he_intersection1->twinEdge();

                        he1->SetNextEdge( he2 );
                        he2->SetPrevEdge( he1 );
                        he2->SetTwinEdge( he_2 );
                        he_2->SetTwinEdge( he2 );
                        v_2->SetIncEdge( he_2 );

                        //bool is_on_convex_hull = UpdateConvexHull( new_site );
                        //CalcIntersectionVertex( he2, next_site, !is_on_convex_hull, v2, v_2, 
                        //    he_intersection2, he_intersection_2 );
                        vector<Halfedge*> result = CalcIntersectionVertices( he2, next_site, !is_on_convex_hull, v2, v_2, 
                            he_intersection2, he_intersection_2 );
                        he_intersection2 = result[0];
                        he_intersection_2 = result[1];

                        if ( he_intersection2 != NULL ) {
                            Halfedge* tmp = he_intersection2->nextEdge();
                            while ( tmp != he_intersection_2 ) {
                                tmp->SetPrevEdge( NULL );
                                tmp->SetNextEdge( NULL );
                                tmp = tmp->nextEdge();
                            }
                            he_intersection2->SetNextEdge( he_2 );
                        }
                        if ( he_intersection_2 != NULL ) {
                            he_intersection_2->SetPrevEdge( he_2 );
                        }

                        he2->SetOriVertex( v2 );
                        he_2->SetOriVertex( v_2 );

                        vd->halfedges.push_back( he2 );
                        vd->halfedges.push_back( he_2 );
                        vd->vertices.push_back( v_2 );

                        if ( he_intersection2->twinEdge()->incFace()->site() == closest_site ) {
                            he2->SetNextEdge( first_he );
                            first_he->SetPrevEdge( he2 );
                            he_2->SetPrevEdge( first_he_ );
                            first_he_->SetNextEdge( he_2 );
                            break;
                        } else {
                            he1 = he2;
                            he_1 = he_2;
                            he_intersection1 = he_intersection_2;
                            he_intersection_1 = he_intersection_1;
                            continue;
                        }
                    }
                } else {
                    Site* next_site;
                    Halfedge* he2;
                    Halfedge* he_2;
                    Halfedge* he_intersection2;
                    Halfedge* he_intersection_2;

                    while ( true ) {
                        next_site = he_intersection_1->twinEdge()->incFace()->site();

                        he2 = new Halfedge();
                        he_2 = new Halfedge();
                        he2->SetIncFace( new_face );
                        he_2->SetIncFace( next_site->incFace() );

                        Point* midp2 = new Point( ( he2->incFace()->site()->x() + he_2->incFace()->site()->x() ) / 2.0, 
                            ( he2->incFace()->site()->y() + he_2->incFace()->site()->y() ) / 2.0 );
                        Vector* dir2 = new Vector( he2->incFace()->site()->y() - he_2->incFace()->site()->y(), 
                            -he2->incFace()->site()->x() + he_2->incFace()->site()->x() );
                        Vector* dir_2 = new Vector( -he2->incFace()->site()->y() + he_2->incFace()->site()->y(), 
                            he2->incFace()->site()->x() - he_2->incFace()->site()->x() );
                        he2->SetMidPoint( midp2 );
                        he_2->SetMidPoint( midp2 );
                        he2->SetDirection( dir2 );
                        he_2->SetDirection( dir_2 );

                        Vertex* v2 = new Vertex();
                        Vertex* v_2 = v1;

                        he_intersection2 = he_intersection_1->twinEdge();
                        he_intersection_2 = NULL;

                        he1->SetPrevEdge( he2 );
                        he2->SetNextEdge( he1 );
                        he2->SetTwinEdge( he_2 );
                        he_2->SetTwinEdge( he2 );
                        v2->SetIncEdge( he2 );

                        //bool is_on_convex_hull = UpdateConvexHull( new_site );
                        //CalcIntersectionVertex( he2, next_site, !is_on_convex_hull, v2, v_2, 
                        //    he_intersection2, he_intersection_2 );
                        vector<Halfedge*> result = CalcIntersectionVertices( he2, next_site, !is_on_convex_hull, v2, v_2, 
                            he_intersection2, he_intersection_2 );
                        he_intersection2 = result[0];
                        he_intersection_2 = result[1];

                        if ( he_intersection2 != NULL ) {
                            Halfedge* tmp = he_intersection2->nextEdge();
                            while ( tmp != he_intersection_2 ) {
                                tmp->SetPrevEdge( NULL );
                                tmp->SetNextEdge( NULL );
                                tmp = tmp->nextEdge();
                            }
                            he_intersection2->SetNextEdge( he_2 );
                        }
                        if ( he_intersection2 != NULL ) {
                            he_intersection_2->SetPrevEdge( he_2 );
                        }

                        he2->SetOriVertex( v2 );
                        he_2->SetOriVertex( v_2 );

                        vd->halfedges.push_back( he2 );
                        vd->halfedges.push_back( he_2 );
                        vd->vertices.push_back( v2 );

                        if ( he_intersection2->twinEdge()->incFace()->site() == closest_site ) {
                            he2->SetPrevEdge( first_he );
                            first_he->SetNextEdge( he2 );
                            he_2->SetNextEdge( first_he_ );
                            first_he_->SetPrevEdge( he_2 );
                            break;
                        } else {
                            he1 = he2;
                            he_1 = he_2;
                            he_intersection1 = he_intersection_2;
                            he_intersection_1 = he_intersection_1;
                            continue;
                        }
                    }
                }
            }
        }
    }
    
    DeleteUselessHalfedges();

    return vd;
}

void VoronoiDiagram::DeleteUselessHalfedges() {
    vector< Halfedge* >::iterator iter = this->halfedges.begin();
    while ( iter != this->halfedges.end() ) {
        if ( (*iter)->prevEdge() == NULL ) {
            iter = this->halfedges.erase( iter );
        } else {
            iter++;
        }
    }
}
void VoronoiDiagram::PrintVoronoiDiagram( VoronoiDiagram* vd ) {

}