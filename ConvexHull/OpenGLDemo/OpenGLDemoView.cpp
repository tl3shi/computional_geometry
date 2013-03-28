// OpenGLDemoView.cpp : COpenGLDemoView 类的实现
//

#include "stdafx.h"
// SHARED_HANDLERS 可以在实现预览、缩略图和搜索筛选器句柄的
// ATL 项目中进行定义，并允许与该项目共享文档代码。
#ifndef SHARED_HANDLERS
#include "OpenGLDemo.h"
#endif


#include "OpenGLDemoDoc.h"
#include "OpenGLDemoView.h"
#include <math.h>
#ifdef _DEBUG
#define new DEBUG_NEW
#endif

#include <algorithm>
#include <vector>
#include <cmath>
#include <stack>

#include "GL/GLU.H" // 已经包含GL/GL.h
#include  "CP_PointVector.h"

using namespace std;



// COpenGLDemoView

IMPLEMENT_DYNCREATE(COpenGLDemoView, CView)

    BEGIN_MESSAGE_MAP(COpenGLDemoView, CView)
        // 标准打印命令
        ON_COMMAND(ID_FILE_PRINT, &CView::OnFilePrint)
        ON_COMMAND(ID_FILE_PRINT_DIRECT, &CView::OnFilePrint)
        ON_COMMAND(ID_FILE_PRINT_PREVIEW, &COpenGLDemoView::OnFilePrintPreview)
        ON_WM_CONTEXTMENU()
        ON_WM_RBUTTONUP()
        ON_WM_LBUTTONUP()
        ON_WM_CREATE()
        ON_WM_DESTROY()
        ON_WM_SIZE()
        ON_WM_MOUSEMOVE()
        ON_WM_ERASEBKGND()
        ON_WM_LBUTTONDBLCLK()
        ON_COMMAND(ID_Nsqr4, &COpenGLDemoView::OnNsqr4)
        ON_COMMAND(ID_Nsqr3, &COpenGLDemoView::OnNsqr3)
        ON_COMMAND(ID_GiftWrapping, &COpenGLDemoView::OnGiftwrapping)
        ON_COMMAND(ID_GrahamScan, &COpenGLDemoView::OnGrahamscan)
    END_MESSAGE_MAP()

    const int MAX_INT = 65536;
    vector<CP_Vector2D> points;
    vector<CP_Vector2D> convexHullResult;
    CP_Vector2D leftest_and_lowest;
    enum Algorithm{ ExtreamPoint_Method, ExtremeEdge_Method} method;

    class MyPoint
    {
    public:
        CP_Vector2D point;
        bool isExtreme;
        MyPoint(CP_Vector2D p, bool extreme):point(p),isExtreme(extreme)
        {}
    };

    class ExtremeEdge
    {
    public:
        CP_Vector2D start;
        CP_Vector2D end;
        ExtremeEdge(CP_Vector2D p1, CP_Vector2D p2):start(p1),end(p2)
        {}
    };


    //test if po is on the left of p2p3
    bool to_left(const CP_Vector2D &p0, const CP_Vector2D &p1, const CP_Vector2D &p2)
    {
        //p0p1 * p2p1 * sin(theta)
        return ((p2 - p1) ^ (p0 - p1)) > 0;
    }
    //test k is in triangle(p,q,r)
    bool in_triangle(const CP_Vector2D &k, const CP_Vector2D &p, const  CP_Vector2D &q, const  CP_Vector2D &r)
    {
        //this case, p q r must be ccw
        //return to_left(k, p, q) && to_left(k, q, r) && to_left(k, r, p);
        //not mattter ccw or cw
        bool left1 = to_left(k, p, q);
        bool left2 = to_left(k, q, r);
        bool left3 = to_left(k, r, p);
        return (left1 == left2) && (left2 == left3);
    }
    int compare(CP_Vector2D &p0, CP_Vector2D &p1)
    {
        if (p0.m_y == p1.m_y)
            return p0.m_x < p1.m_x;
        return p0.m_y < p1.m_y;
    }
    
    //leftmost, theta in [0, pi], increasing 
    int compareByAngle(CP_Vector2D &p1, CP_Vector2D &p2)
    {
        if(p1 == leftest_and_lowest) return 1;//leftmost is the min means p1 < p2
        if(p2 == leftest_and_lowest) return 0;//leftmost is the min means p1 > p2
        
        double arccos_angle1 = (CP_Vector2D(MAX_INT, leftest_and_lowest.m_y) - leftest_and_lowest) * (p1 - leftest_and_lowest) / (CP_Vector2D(MAX_INT, leftest_and_lowest.m_y) - leftest_and_lowest).mf_getLength() / (p1 - leftest_and_lowest).mf_getLength();
        double arccos_angle2 = (CP_Vector2D(MAX_INT, leftest_and_lowest.m_y) - leftest_and_lowest) * (p2 - leftest_and_lowest) / (CP_Vector2D(MAX_INT, leftest_and_lowest.m_y) - leftest_and_lowest).mf_getLength() / (p2 - leftest_and_lowest).mf_getLength();
        //cos  decreasing
        return arccos_angle1 > arccos_angle2 ;
    }
    //compareEdgeByAngle, increasing 
    int compareEdgeByAngle(ExtremeEdge &p1, ExtremeEdge &p2)
    {
        double arccos_angle1 = (p1.end.m_x - p1.start.m_x) / (p1.end - p1.start).mf_getLength();
        double arccos_angle2 = (p2.end.m_x - p2.start.m_x) / (p2.end - p2.start).mf_getLength();
        if(p1.end.m_y - p1.start.m_y >= 0 && p2.end.m_y - p2.start.m_y >= 0 )//[0,pi], monotone decreasing
            return arccos_angle1 >= arccos_angle2;
        if(p1.end.m_y - p1.start.m_y < 0 && p2.end.m_y - p2.start.m_y < 0)//[pi, 2*pi] monotone increasing
            return arccos_angle1 < arccos_angle2;
        if((p1.end.m_y - p1.start.m_y) * (p2.end.m_y - p2.start.m_y) < 0)//one [0,pi],the other[pi,2pi]
            return p1.end.m_y - p1.start.m_y > 0 ;//if true, theta(p1) < theta(p2) ,otherwise theta(p1) > theta(p2)
        //remove the warning message.
        return 1;
    }

    void for_debug()
    {
        {
            leftest_and_lowest = CP_Vector2D(0, 0);
            CP_Vector2D p1 = CP_Vector2D(1, 0);
            CP_Vector2D p2 = CP_Vector2D(0, 1);
            CP_Vector2D p3 = CP_Vector2D(1, 1);
            int t = compareByAngle(CP_Vector2D(1,0), CP_Vector2D(0,1));
            assert(t);//true,means cos(p1) > cos(p2), theta(p1) < thita(p2);

            assert(to_left(p2, leftest_and_lowest, p1));
            assert(in_triangle(leftest_and_lowest, CP_Vector2D(-1, 0), p2, p1));
        }

        {
            leftest_and_lowest = CP_Vector2D(568, 160);
            CP_Vector2D p1 = CP_Vector2D(258, 407);
            CP_Vector2D p2 = CP_Vector2D(294, 275);
            
            int t = compareByAngle(p1, p2);
            assert(t == 1);//theta p1 < p2, cos(p1) > cos(p2)
          
        }
    }
    //O(n^4)
    vector<CP_Vector2D> cal_extreme_points()
    {
        vector<CP_Vector2D> extreme_points;
        vector<MyPoint> my_points;
        unsigned int size = points.size();
        //inital all points to extreme
        for (unsigned int i = 0; i < size; i++)
            my_points.push_back(MyPoint(points[i], true));

        for (unsigned int p = 0; p < size; p++)
        {
            for (unsigned int q = p+1 ; q < size; q++)
            {
                for (unsigned int r = q+1; r < size; r++)
                {
                    for (unsigned int k = 0; k < size; k++)
                    {
                        if (k == p || k == q || k == r || !my_points[k].isExtreme) 
                            continue;
                        if (in_triangle(my_points[k].point, my_points[p].point, my_points[q].point, my_points[r].point))
                            my_points[k].isExtreme = false;
                    }
                }
            }
        }

        for (unsigned int i = 0; i < size; i++)
        {
            if(my_points[i].isExtreme)
                extreme_points.push_back(my_points[i].point);
        }
        //sort the extrme_points
        sort(extreme_points.begin(), extreme_points.end(), compare);
        leftest_and_lowest = extreme_points[0];
        //sort the extrme_points by left_and_lowest polar
        sort(extreme_points.begin()+1, extreme_points.end(), compareByAngle);
        return extreme_points;
    }

    //O(n^3)
    vector<CP_Vector2D> cal_extreme_edges()
    {
        vector<ExtremeEdge> edgs;
        
        unsigned int size = points.size();

        for (unsigned int p = 0; p < size; p++)
        {
            for (unsigned int q = p+1 ; q < size; q++)
            {
                bool all_left = true, all_right = true;
                for (unsigned int k = 0; k < size; k++)
                {
                    if (k == p || k == q || p == q)
                        continue;
                    if(to_left(points[k], points[p], points[q]))
                        all_left = false;
                    else
                        all_right = false;
                }
                //let the direction be all the same
                if (all_left)
                    edgs.push_back(ExtremeEdge(points[p], points[q]));
                if (all_right)
                    edgs.push_back(ExtremeEdge(points[q], points[p]));
            }
        }
        
        vector<CP_Vector2D> convexPoints;

        ////O(n^2), should better sort the edges by angle.
        //convexPoints.push_back(edgs[0].end);
        //while(convexPoints.size() != edgs.size())
        //{
        //    for (unsigned int i = 0; i < edgs.size(); i++)
        //    {
        //        if (convexPoints[convexPoints.size() - 1] == edgs[i].start)
        //            convexPoints.push_back(edgs[i].end);
        //    }
        //}

        //O(nlogn)
        sort(edgs.begin(), edgs.end(), compareEdgeByAngle);
        for (unsigned int i=0; i < edgs.size(); i++)
        {
            convexPoints.push_back(edgs[i].start);
        }

        return convexPoints;
    }

    //O(n^2),in fact O(n*h),h = convex hull.size, output sensitive
    vector<CP_Vector2D> giftwraping()
    {
        int ltl = 0;//find the lowest-and-leftmost point
        for (unsigned int i = 1; i < points.size(); i++)
        {
            if(points[i].m_y < points[ltl].m_y || 
               (points[i].m_y == points[ltl].m_y && points[i].m_x < points[ltl].m_x))
               ltl = i;
        }
        vector<CP_Vector2D> convexPoints;
        convexPoints.push_back(points[ltl]);
        for(int p = ltl; ;)
        {
            int q = -1;
            for(unsigned int k = 0; k < points.size(); k++)
            {
                if( (k != p) &&  (q < 0 || !to_left(points[k], points[p], points[q])))
                    q = k;//update q if k lies in right of pq
            }
            if(ltl == q)
                break;// has find out the convex hull circle
            //find q
            convexPoints.push_back(points[q]);
            p = q;
        }
         return convexPoints;
    }

    //O(nlgn)
    vector<CP_Vector2D> graham_scan()
    {
        int ltl = 0;//find the lowest-and-leftmost point
        for (unsigned int i = 1; i < points.size(); i++)
        {
            if(points[i].m_y < points[ltl].m_y || 
                (points[i].m_y == points[ltl].m_y && points[i].m_x < points[ltl].m_x))
                ltl = i;
        }

        leftest_and_lowest = points[ltl];
        //compareByAngle should use the ltl point
        sort(points.begin(), points.end(), compareByAngle);

        vector<CP_Vector2D> convexPoints;
        stack<CP_Vector2D> S, T;
        S.push(points[0]);S.push(points[1]);
        for (int i = points.size()-1; i >= 2;i--)
        {
            T.push(points[i]);
        }
        CP_Vector2D s0, s1;
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
        unsigned int size = S.size();
        for(int i = 0; i < size; i++)
        {
            convexPoints.push_back(S.top());
            S.pop();//size is change, cannot use i < S.size()
        }
        return convexPoints;
    }
   

    // COpenGLDemoView 构造/析构

    COpenGLDemoView::COpenGLDemoView()
    {
        // TODO: 在此处添加构造代码

    }

    COpenGLDemoView::~COpenGLDemoView()
    {
    }

    BOOL COpenGLDemoView::PreCreateWindow(CREATESTRUCT& cs)
    {
        // TODO: 在此处通过修改
        //  CREATESTRUCT cs 来修改窗口类或样式

        return CView::PreCreateWindow(cs);
    }
 
    void drawResult(vector<CP_Vector2D> &result)
    {
        if (result.size() == 0)
            return;
        glBegin(GL_LINE_LOOP);
        for (unsigned int i=0; i < result.size(); i++)
        {
            glColor3d(1, 0, 0);
            CP_Vector2D p = result[i];
            glVertex3d(p.m_x, p.m_y, 0);
        }
        glEnd();
    }

    void drawPoints()
    {
        glClearColor(1,1,1,1);
        glClear(GL_COLOR_BUFFER_BIT);
        glPointSize(3);
        glBegin(GL_POINTS);
        for (unsigned int i=0; i < points.size(); i++)
        {
            glColor3d(1, 0, 0);
            CP_Vector2D p = points.at(i);
            glVertex3d(p.m_x, p.m_y, 0);
        }
        glEnd();
    }

    void COpenGLDemoView::OnDraw(CDC* pDC)
    {
        COpenGLDemoDoc* pDoc = GetDocument();
        ASSERT_VALID(pDoc);
        if (!pDoc)
            return;
        wglMakeCurrent(pDC->m_hDC, m_hRC);

        glClearColor(0,0,0,0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        drawPoints();
        /*
        if (method == ExtreamPoint_Method)
            drawResult(convexHullResult);
        if (method == ExtremeEdge_Method)
            drawResult(convexHullResult);
        */
        drawResult(convexHullResult);

        SwapBuffers(pDC->m_hDC);
        wglMakeCurrent(NULL, NULL);

    }


    // COpenGLDemoView 打印


    void COpenGLDemoView::OnFilePrintPreview()
    {
#ifndef SHARED_HANDLERS
        AFXPrintPreview(this);
#endif
    }


    BOOL COpenGLDemoView::OnPreparePrinting(CPrintInfo* pInfo)
    {
        // 默认准备
        return DoPreparePrinting(pInfo);
    }

    void COpenGLDemoView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
    {
        // TODO: 添加额外的打印前进行的初始化过程
    }

    void COpenGLDemoView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
    {
        // TODO: 添加打印后进行的清理过程
    }

    void COpenGLDemoView::OnRButtonUp(UINT nflags, CPoint point)
    {
        CView::OnRButtonUp(nflags, point);
    }
    /*
    void COpenGLDemoView::OnContextMenu(CWnd*   pWnd  , CPoint point)
    {
    #ifndef SHARED_HANDLERS
    theApp.GetContextMenuManager()->ShowPopupMenu(IDR_POPUP_EDIT, point.x, point.y, this, TRUE);
    #endif
    }
    */

    // COpenGLDemoView 诊断

#ifdef _DEBUG
    void COpenGLDemoView::AssertValid() const
    {
        CView::AssertValid();
    }

    void COpenGLDemoView::Dump(CDumpContext& dc) const
    {
        CView::Dump(dc);
    }

    COpenGLDemoDoc* COpenGLDemoView::GetDocument() const // 非调试版本是内联的
    {
        ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(COpenGLDemoDoc)));
        return (COpenGLDemoDoc*)m_pDocument;
    }
#endif //_DEBUG


    // COpenGLDemoView 消息处理程序


    int COpenGLDemoView::OnCreate(LPCREATESTRUCT lpCreateStruct)
    {
        if (CView::OnCreate(lpCreateStruct) == -1)
            return -1;

        // The PIXELFORMATDESCRIPTOR structure describes
        //		the pixel format of a drawing surface.
        PIXELFORMATDESCRIPTOR pfd =
        { 
            sizeof(PIXELFORMATDESCRIPTOR),  //  size of this pfd 
            1,                     			// version number 
            PFD_DRAW_TO_WINDOW |   	// support window 
            PFD_SUPPORT_OPENGL |   	// support OpenGL 
            PFD_DOUBLEBUFFER,	// double buffered
            PFD_TYPE_RGBA,
            24,                    	// 24-bit color depth 
            0, 0, 0, 0, 0, 0,      // color bits ignored 
            0,                     	// no alpha buffer 
            0,                     	// shift bit ignored 
            0,                     	// no accumulation buffer 
            0, 0, 0, 0,            	// accum bits ignored 
            32,                    	// 32-bit z-buffer (depth)
            0,                     	// no stencil buffer 
            0,                     	// no auxiliary buffer 
            PFD_MAIN_PLANE,        // main layer 
            0,                     	// reserved 
            0, 0, 0                	// layer masks ignored 
        }; 
        CClientDC dc(this);
        // Get the best available match of pixel format for the device context
        // In other words, if this computer doesn't support features that I
        // asked for, try to get the next best thing.  i.e. 16-bit color mode
        // instead of 24-bit color mode.
        int pixelFormat = ChoosePixelFormat(dc.m_hDC, &pfd);

        // Set the pixel format to the best pixel format I can get (see above)
        // and if that operation fails, bring up a message box that tells the
        // user the error.
        if (!SetPixelFormat(dc.m_hDC, pixelFormat, &pfd))
        {
            MessageBox(_T("Error: Unable to Set Pixel Format in CGLTemplate1View::OnCreate( ) "),
                _T("Application Error"), MB_ICONERROR);
        }
        // Creates an OpenGL rendering context so that OpenGL knows how to draw
        // to this view. You can't use OpenGL in MFC without using the handle
        // that this function returns
        m_hRC = wglCreateContext(dc.m_hDC);

        for_debug();
        return 0;
    }


    void COpenGLDemoView::OnDestroy()
    {
        CView::OnDestroy();
        // Set : a specified OpenGL rendering context ==> NULL
        // Set : current rendering context ==> NULL
        wglMakeCurrent(NULL, NULL);

        // Delete the handle to an OpenGL rendering context 
        wglDeleteContext(m_hRC);
        m_hRC=NULL;

    }


    void COpenGLDemoView::OnSize(UINT nType, int cx, int cy)
    {

        CView::OnSize(nType, cx, cy);
        CClientDC dc(this);
        wglMakeCurrent(dc.m_hDC, m_hRC);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

        gluOrtho2D(0, cx, 0, cy);

        glMatrixMode(GL_MODELVIEW);
        glViewport(0, 0, cx, cy);
        wglMakeCurrent(NULL, NULL);

    }


    BOOL COpenGLDemoView::OnEraseBkgnd(CDC* pDC)
    {

        //return CView::OnEraseBkgnd(pDC);
        return TRUE;
    }

    void COpenGLDemoView::OnMouseMove(UINT nFlags, CPoint point)
    {
        /*
        if(!firstClick)
        return;

        CRect rect;
        GetClientRect(&rect);
        int x = point.x;
        int y = rect.Height() - point.y;

        hoverPoint = CP_Vector2D(x, y);
        */
        CView::OnMouseMove(nFlags, point);
    }

    void COpenGLDemoView::OnLButtonUp(UINT nFlags, CPoint point)
    {

        {
            CRect rect;
            GetClientRect(&rect);
            int x = point.x;
            int y = rect.Height() - point.y;

            CP_Vector2D p = CP_Vector2D(x, y);
            points.push_back(p);
        }
        Invalidate(TRUE);

        CView::OnLButtonUp(nFlags, point);
    }

    void COpenGLDemoView::OnLButtonDblClk(UINT nFlags, CPoint point)
    {
        /*
        MessageBox(L"清空控制点，重画");
        ctrlPoints.clear();
        isReady = false;
        Invalidate(FALSE);
        */
       
        CView::OnLButtonDblClk(nFlags,point);
    }

    void COpenGLDemoView::OnNsqr4()
    {
        convexHullResult.clear();
        convexHullResult = cal_extreme_points();
        Invalidate(FALSE);
    }


    void COpenGLDemoView::OnNsqr3()
    {
        convexHullResult.clear();
        convexHullResult = cal_extreme_edges();
        Invalidate(FALSE);
    }


    void COpenGLDemoView::OnGiftwrapping()
    {
        convexHullResult.clear();
        convexHullResult = giftwraping();
        Invalidate(FALSE);
    }


    void COpenGLDemoView::OnGrahamscan()
    {
        convexHullResult.clear();
        convexHullResult = graham_scan();
        Invalidate(FALSE);
    }
