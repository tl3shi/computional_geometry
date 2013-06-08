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
#include "voronoi_diagram.h"
#include <gl/glut.h>
#include <iostream>
#include <cstdlib>
#include <fstream>


bool debug = true;
using namespace std;

vector<Point> points;
VoronoiDiagram * result = NULL;

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
    ON_COMMAND(ID_32771, &COpenGLDemoView::OnDevideConquer)
END_MESSAGE_MAP()

void initLights();
void drawString(int x, int y, char* str);
void drawString(int x, int y, CString cstr);

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
void writetofile()
{
    ofstream outf("f:/pointsout.txt");
    streambuf *default_buf=cout.rdbuf();   
    cout.rdbuf( outf.rdbuf() );   
    for (unsigned int i = 0; i < points.size(); i++)
    {
        Point p = points.at(i);
        cout << p.x() << "," << p.y() << endl;
    }
    cout.rdbuf(default_buf);
}
vector<Point> readfromfile(char* filename)
{
    vector<Point> ps;
    ifstream f(filename);
    while (!f.eof())
    {
        int x = -1, y = -1;
        char comma;
        f >> x >> comma >> y;
        if(x == -1 && y == -1) continue;
        ps.push_back(Point(x, y));
    }
    f.close();
    return ps;
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
        Point p = points.at(i);
        glVertex3d(p.x(), p.y(), 0);
    }
    glEnd();

    CString str;
    for (unsigned int i=0; i < points.size(); i++)
    {
        glColor3d(1, 0, 0);
        Point p = points.at(i);
        str.Format(L"%d", i);
        drawString(p.x(), p.y()+5, str);
    }
}

void drawString(int x, int y, char* str)
{
    glColor3d(0.0, 0.0, 0.0);
    int n = strlen(str);  
    //设置要在屏幕上显示字符的起始位置 
    glRasterPos2i(x,y);  
    //逐个显示字符串中的每个字符  
    for (int i = 0; i < n; i++)  
        glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *(str+i)); 
}
void drawString(int x, int y, CString cstr)
{
    char * str = (LPSTR)(LPCTSTR)cstr;
    glColor3d(0.0, 0.0, 0.0);
    int n = strlen(str);  
    //设置要在屏幕上显示字符的起始位置 
    glRasterPos2i(x,y);  
    //逐个显示字符串中的每个字符  
    for (int i = 0; i < n; i++)  
        glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *(str+i)); 
}

const int INFINITE_LENGTH  = 2000;

//if the half edge's top infomation is all right, may be the following function is OK
void drawResult()
{
    Point infinitePoint = Point(DBL_MAX, DBL_MAX);
    glLineWidth(1);
    glBegin(GL_LINES);
    glColor3d(1, 0, 0);
    if(result->sites.size() == 1)//one site, the whole face
        return;
    if(result->sites.size() == 2)//two sites,draw a line
    {
        Halfedge * e  = result->halfedges.at(0);
        Point * p = e->midPoint();
        Point oneP, anotherP ;
        oneP = (*p) + (-*(e->direction())) * INFINITE_LENGTH;
        anotherP = (*p) - (-*(e->direction())) * INFINITE_LENGTH;
        glVertex3d(oneP.x(), oneP.y(), 0);
        glVertex3d(anotherP.x(), anotherP.y(), 0);
        glEnd();
        return;
    }


    ofstream outf("f:/pointsout.debug.txt");
    streambuf *default_buf=cout.rdbuf();   
    cout.rdbuf( outf.rdbuf() );   
     
    for (unsigned int i = 0; i < result->halfedges.size(); i++)
    {
        Halfedge * edge = result->halfedges.at(i);
        
        if (edge->nextEdge() != NULL)// edge point to regular point
        {
            Point p = edge->oriVertex()->p;
            if(p == infinitePoint)//edge from infinite
            {
                Point np = edge->twinEdge()->oriVertex()->p;
                glVertex3d(np.x(), np.y(), 0);

                p = np + (-*(edge->direction())) * INFINITE_LENGTH;
                glVertex3d(p.x(), p.y(), 0);

                cout<< i<<":"<< p.x() << "," << p.y() <<"," << np.x() << ","<< np.y() << endl;
                continue;
            }
            //regular edge
            glVertex3d(p.x(), p.y(), 0);
        
            
            Point np = edge->twinEdge()->oriVertex()->p;
            glVertex3d(np.x(), np.y(), 0);
            cout<< i<<":"<< p.x() << "," << p.y() <<"," << np.x() << ","<< np.y() << endl;
        }
        else//ray, point to infinite
        {
            Point p = edge->oriVertex()->p;
            if (p == infinitePoint)
            {
                Point np = edge->twinEdge()->oriVertex()->p;
                if (np == infinitePoint)// a line from infinite to infinite 
                {
                    p = Point(0,0) - (*(edge->direction())) * INFINITE_LENGTH;
                    glVertex3d(p.x(), p.y(), 0);
                  

                    p = Point(0,0) + (*(edge->direction())) * INFINITE_LENGTH;
                    glVertex3d(p.x(), p.y(), 0);
                    cout<< i<<":"<<  p.x() << "," << p.y() <<"," << np.x() << "," << np.y() << endl;
                    continue;
                }
            }
            glVertex3d(p.x(), p.y(), 0);
            Point np = p + (*(edge->direction())) * INFINITE_LENGTH;
            glVertex3d(np.x(), np.y(), 0);

            cout<< i<<":"<<  p.x() << "," << p.y() <<"," << np.x() << "," << np.y() << endl;
        }
    }
    glEnd();
    cout.rdbuf(default_buf);   

    //draw string
    return;
    CString str;
    for (unsigned int i = 0; i < result->halfedges.size(); i++)
    {
        Halfedge * edge = result->halfedges.at(i);
        Point* p = edge->midPoint();
        str.Format(L"%d", i);
        if(edge->twinEdge()->hasDraw)
            drawString(p->x(), p->y()-5, str);
        else
            drawString(p->x(), p->y(), str);
        edge->hasDraw = true;
    }
}

void drawResultFUCK()
{
    Point infinitePoint = Point(DBL_MAX, DBL_MAX);
    glLineWidth(1);
    glBegin(GL_LINES);
    glColor3d(1, 0, 0);
    if(result->sites.size() == 1)//one site, the whole face
        return;
    if(result->sites.size() == 2)//two sites,draw a line
    {
        Halfedge * e  = result->halfedges.at(0);
        Point * p = e->midPoint();
        Point oneP, anotherP ;
        oneP = (*p) + (-*(e->direction())) * INFINITE_LENGTH;
        anotherP = (*p) - (-*(e->direction())) * INFINITE_LENGTH;
        glVertex3d(oneP.x(), oneP.y(), 0);
        glVertex3d(anotherP.x(), anotherP.y(), 0);
        glEnd();
        return;
    }


    ofstream outf("f:/pointsout.debug.txt");
    streambuf *default_buf=cout.rdbuf();   
    cout.rdbuf( outf.rdbuf() );   

    for (unsigned int i = 0; i < result->halfedges.size(); i++)
    {
        Halfedge * edge = result->halfedges.at(i);

        if (edge->nextEdge() != NULL)// edge point to regular point
        {
            Point p = edge->oriVertex()->p;
            if(p == infinitePoint)//edge from infinite
            {
                Point np = edge->twinEdge()->oriVertex()->p;
                glVertex3d(np.x(), np.y(), 0);
                
                //fuck begin
                if(edge->endVertex != NULL && edge->endVertex->p != infinitePoint)
                {
                    p = edge->endVertex->p;
                    glVertex3d(p.x(), p.y(), 0);
                    continue;
                }
                //fuck end

                p = np + (-*(edge->direction())) * INFINITE_LENGTH;
                glVertex3d(p.x(), p.y(), 0);

                cout<< i<<":"<< p.x() << "," << p.y() <<"," << np.x() << ","<< np.y() << endl;
                continue;
            }
            //regular edge
            glVertex3d(p.x(), p.y(), 0);
            Point np = edge->twinEdge()->oriVertex()->p;
            glVertex3d(np.x(), np.y(), 0);
            cout<< i<<":"<< p.x() << "," << p.y() <<"," << np.x() << ","<< np.y() << endl;
        }
        else//ray, point to infinite
        {
            Point p = edge->oriVertex()->p;
            if (p == infinitePoint)
            {
                Point np = edge->twinEdge()->oriVertex()->p;
                if (np == infinitePoint)// a line from infinite to infinite 
                {
                    p = Point(0,0) - (*(edge->direction())) * INFINITE_LENGTH;
                    glVertex3d(p.x(), p.y(), 0);
                    p = Point(0,0) + (*(edge->direction())) * INFINITE_LENGTH;
                    glVertex3d(p.x(), p.y(), 0);
                    cout<< i<<":"<<  p.x() << "," << p.y() <<"," << np.x() << "," << np.y() << endl;
                    continue;
                }
            }
            glVertex3d(p.x(), p.y(), 0);

            //fuck begin
            if(edge->endVertex != NULL && edge->endVertex->p != infinitePoint)
            {
                p = edge->endVertex->p;
                glVertex3d(p.x(), p.y(), 0);
                continue;
            }
            //fuck end

            Point np = p + (*(edge->direction())) * INFINITE_LENGTH;
            glVertex3d(np.x(), np.y(), 0);

            cout<< i<<":"<<  p.x() << "," << p.y() <<"," << np.x() << "," << np.y() << endl;
        }
    }
    glEnd();
    cout.rdbuf(default_buf);   
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
    drawString(0, 0 , ("tanglei-begin"));
    if (result != NULL)
    {
        drawResult();
    }
    drawString(100, 100 , ("tanglei-end"));

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
    if(debug)
    {
        points = readfromfile("f:/points-5.txt");
        OnDevideConquer();
    }
	//MessageBox(L"鼠标左键选择控制点位置\r\n右键生成画bezier曲线\r\n双击左键清空控制点");
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

	    points.push_back(Point(x, y));
    }
	Invalidate(TRUE);
	
	CView::OnLButtonUp(nFlags, point);
}

void COpenGLDemoView::OnLButtonDblClk(UINT nFlags, CPoint point)
{
	
	MessageBox(L"清空控制点，重画");
	points.clear();
    result->halfedges.clear();
	Invalidate(FALSE);
	
	CView::OnLButtonDblClk(nFlags,point);
}



void COpenGLDemoView::OnDevideConquer()
{
    /*vector<Point> testPoints;
    testPoints.push_back(Point(200,200)+Point(0,0));
    testPoints.push_back(Point(200,200)+Point(40,0));
    testPoints.push_back(Point(200,200)+Point(10,50));
    points = testPoints;*/
    if(points.size() > 0 )
        writetofile();
  
    VoronoiDiagram * vd = new VoronoiDiagram();
    result = vd->DevideConquerConstruction(points);
    Invalidate(TRUE);
}
