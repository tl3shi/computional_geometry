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
#include <math.h>

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
    END_MESSAGE_MAP()

    const int MAX_INT = 65536;
    class MyPoint
    {
    public:
        CP_Vector2D point;
        bool isExtreme;
        MyPoint(CP_Vector2D &p):point(p),isExtreme(false)
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
    CP_Vector2D leftest_and_lowest;

    int compareByAngle(CP_Vector2D &p1, CP_Vector2D &p2)
    {
        double angle1 = (CP_Vector2D(MAX_INT, leftest_and_lowest.m_y) - leftest_and_lowest) * (p1 - leftest_and_lowest) / (CP_Vector2D(MAX_INT, leftest_and_lowest.m_y) - leftest_and_lowest).mf_getLength() / (p1 - leftest_and_lowest).mf_getLength();
        double angle2 = (CP_Vector2D(MAX_INT, leftest_and_lowest.m_y) - leftest_and_lowest) * (p2 - leftest_and_lowest) / (CP_Vector2D(MAX_INT, leftest_and_lowest.m_y) - leftest_and_lowest).mf_getLength() / (p2 - leftest_and_lowest).mf_getLength();
        
        return angle1 > angle2 ;
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
    vector<CP_Vector2D> cal_extreme_points(vector<MyPoint> &points)
    {
        vector<CP_Vector2D> extreme_points;

        unsigned int size = points.size();
        //inital all points to extreme
        for (unsigned int i = 0; i < size; i++)
            points[i].isExtreme = true;

        for (unsigned int p = 0; p < size; p++)
        {
            for (unsigned int q = p+1 ; q < size; q++)
            {
                for (unsigned int r = q+1; r < size; r++)
                {
                    for (unsigned int k = 0; k < size; k++)
                    {
                        if (k == p || k == q || k == r || !points[k].isExtreme) 
                            continue;
                        if (in_triangle(points[k].point, points[p].point, points[q].point, points[r].point))
                            points[k].isExtreme = false;
                    }
                }
            }
        }

        for (unsigned int i = 0; i < size; i++)
        {
            if(points[i].isExtreme)
                extreme_points.push_back(points[i].point);
        }
        //sort the extrme_points
        sort(extreme_points.begin(), extreme_points.end(), compare);
        leftest_and_lowest = extreme_points[0];
        //sort the extrme_points by left_and_lowest polar
        sort(extreme_points.begin()+1, extreme_points.end(), compareByAngle);
        return extreme_points;
    }


    vector<MyPoint> points;

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
    vector<CP_Vector2D> convexHullResult;
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
            CP_Vector2D p = points.at(i).point;
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
            points.push_back(MyPoint(p));
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
        convexHullResult = cal_extreme_points(points);
        Invalidate(FALSE);
    }
