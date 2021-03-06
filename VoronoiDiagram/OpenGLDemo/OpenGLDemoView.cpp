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
#include  <io.h>
#include <set>
set<Point> pointset;

Point infinitePoint = Point(DBL_MAX, DBL_MAX);

bool debug = false;
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
    ON_COMMAND(ID_32772, &COpenGLDemoView::OnIncrementalConstruction)
    ON_COMMAND(ID_FILE_SAVE_AS, &COpenGLDemoView::OnFileSaveAs)
    ON_COMMAND(ID_FILE_OPEN, &COpenGLDemoView::OnFileOpen)
    ON_COMMAND(ID_FILE_SAVE, &COpenGLDemoView::OnFileSave)
END_MESSAGE_MAP()

void initLights();
void drawString(int x, int y, char* str);
void drawString(int x, int y, CString cstr);
void processPoints();

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
void writetofile(string filename="pointsout.txt")
{
    //ofstream outf("pointsout.txt");
    ofstream outf(filename);
    streambuf *default_buf=cout.rdbuf();   
    cout.rdbuf( outf.rdbuf() );   
    for (unsigned int i = 0; i < points.size(); i++)
    {
        Point p = points.at(i);
        cout << p.x() << "," << p.y() << ";" << i <<endl;
    }
    cout.rdbuf(default_buf);
}
vector<Point> readfromfile(string filename)
{
    vector<Point> ps;
    ifstream f(filename);
    if(!f.good())
        return ps;
    int idindex=-1;
    while (!f.eof())
    {
        int x = -1, y = -1;
        char comma;
        f >> x >> comma >> y >> comma >> idindex;
        if(x == -1 && y == -1)
        {
            cout << "Make sure that the input file format is OK." << endl;
            break;
        }
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

    char * str = new char[3];
    for (unsigned int i=0; i < points.size(); i++)
    {
        glColor3d(1, 0, 0);
        Point p = points.at(i);
        sprintf(str, "%d", i);
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
    int n = cstr.GetLength();//strlen(cstr);  // the force converter may ensercurity. todo
    //设置要在屏幕上显示字符的起始位置 
    glRasterPos2i(x,y);  
    //逐个显示字符串中的每个字符  
    for (int i = 0; i < n; i++)  
        glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *(str+i)); 
}

//if the half edge's top infomation is all right, may be the following function is OK
void drawResult()
{
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


    ofstream outf("pointsout.debug.txt");
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
                if(np == infinitePoint)
                {
                    np = (-*(edge->direction())) * INFINITE_LENGTH;
                }
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
                //edge from infinite to a fix point
                glVertex3d(np.x(), np.y(), 0);
                p = np - (*(edge->direction())) * INFINITE_LENGTH;
                glVertex3d(p.x(), p.y(), 0);
                continue;
            }
            glVertex3d(p.x(), p.y(), 0);
            Point np = p + (*(edge->direction())) * INFINITE_LENGTH;
            glVertex3d(np.x(), np.y(), 0);

            cout<< i<<":"<<  p.x() << "," << p.y() <<"," << np.x() << "," << np.y() << endl;
        }
    }
    glEnd();
    

    //draw string
    if(!debug)
        return;
    char* str = new char[3];
    for (unsigned int i = 0; i < result->halfedges.size(); i++)
    {
        Halfedge * edge = result->halfedges.at(i);
        Point* p = edge->midPoint();
        sprintf(str, "%d", i);
        cout << i << p->toString()<< str << endl;
        if(edge->twinEdge()->hasDraw)
        {
            //cout << i << p->toString()<< str << endl;
            drawString(p->x(), p->y()-10, str);
            //do not draw....
        }else
        {
            drawString(p->x(), p->y(), str);
        }
        edge->hasDraw = true;
    }
    cout.rdbuf(default_buf);  
}

void drawResultForIncrementalConstruction()
{
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


    ofstream outf("pointsout.debug.txt");
    streambuf *default_buf=cout.rdbuf();   
    cout.rdbuf( outf.rdbuf() );   

    for (unsigned int i = 0; i < result->halfedges.size(); i++)
    {
        Halfedge * edge = result->halfedges.at(i);

        if ( edge->nextEdge()->oriVertex()->x() != DBL_MAX && edge->oriVertex()->x() != DBL_MAX )// edge point to regular point
        {
            Point p = edge->oriVertex()->p;
            if(p == infinitePoint)//edge from infinite
            {
                Point np = edge->twinEdge()->oriVertex()->p;
                if(np == infinitePoint)
                {
                    np = (-*(edge->direction())) * INFINITE_LENGTH;
                }
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
                //edge from infinite to a fix point
                glVertex3d(np.x(), np.y(), 0);
                p = np - (*(edge->direction())) * INFINITE_LENGTH;
                glVertex3d(p.x(), p.y(), 0);
                continue;
            }
            glVertex3d(p.x(), p.y(), 0);
            Point np = p + (*(edge->direction())) * INFINITE_LENGTH;
            glVertex3d(np.x(), np.y(), 0);

            cout<< i<<":"<<  p.x() << "," << p.y() <<"," << np.x() << "," << np.y() << endl;
        }
    }
    glEnd();


    //draw string
    if(!debug)
        return;
    char* str = new char[3];
    for (unsigned int i = 0; i < result->halfedges.size(); i++)
    {
        Halfedge * edge = result->halfedges.at(i);
        Point* p = edge->midPoint();
        sprintf(str, "%d", i);
        cout << i << p->toString()<< str << endl;
        if(edge->twinEdge()->hasDraw)
        {
            //cout << i << p->toString()<< str << endl;
            drawString(p->x(), p->y()-10, str);
            //do not draw....
        }else
        {
            drawString(p->x(), p->y(), str);
        }
        edge->hasDraw = true;
    }
    cout.rdbuf(default_buf);  
}

void drawResultFUCK()
{
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
                
                ////fuck begin
                //if(edge->endVertex != NULL && edge->endVertex->p != infinitePoint)
                //{
                //    p = edge->endVertex->p;
                //    glVertex3d(p.x(), p.y(), 0);
                //    continue;
                //}
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

            ////fuck begin
            //if(edge->endVertex != NULL && edge->endVertex->p != infinitePoint)
            //{
            //    p = edge->endVertex->p;
            //    glVertex3d(p.x(), p.y(), 0);
            //    continue;
            //}
            ////fuck end

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
    
    processPoints();

    drawPoints();
    //drawString(0, 0 , ("tanglei-begin"));
    if (result != NULL)
    {
        if ( result->method == 0 ) {
            drawResultForIncrementalConstruction();
        } else if ( result->method == 1 ) {
            drawResult();
        }
    }
    //drawString(100, 100 , ("tanglei-end"));

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
        points = readfromfile("pointsout.txt");
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

	    //points.push_back(Point(x, y));
        pointset.insert(Point(x, y));
    }
	Invalidate(TRUE);
	
	CView::OnLButtonUp(nFlags, point);
}

void COpenGLDemoView::OnLButtonDblClk(UINT nFlags, CPoint point)
{
	
	MessageBox(L"清空控制点，重画");
	points.clear();
    pointset.clear();
    result->halfedges.clear();

	Invalidate(FALSE);
	
	CView::OnLButtonDblClk(nFlags,point);
}
const double POINT_TOLERANCE = 0.1; 
bool pointComp(const Point &p0, const Point &p1)
{
    if (abs(p0.x()-p1.x()) < POINT_TOLERANCE
        && abs(p1.y()-p0.y()) < POINT_TOLERANCE)
        return true;
    return false;
}

void processPoints()
{
    if(points.size() == 0 && pointset.size() == 0) return;
    if(pointset.size() == 0)//may read from file...for test..
    {
        for (unsigned int i = 0; i < points.size(); i++)
        {
            pointset.insert(points[i]);
        }
    }
    points.resize(pointset.size());
    std::copy(pointset.begin(), pointset.end(), points.begin());
}
vector<Point *> getPoints(vector<Point> pts)
{
    vector<Point *> re;
    for (unsigned int i =0; i < pts.size(); i++ )
    {
        re.push_back(new Point(pts[i].x(),pts[i].y()));
    }
    return re;
}

void COpenGLDemoView::OnDevideConquer()
{
    /*vector<Point> testPoints;
    testPoints.push_back(Point(200,200)+Point(0,0));
    testPoints.push_back(Point(200,200)+Point(40,0));
    testPoints.push_back(Point(200,200)+Point(10,50));
    points = testPoints;*/

    //preprocess the points data, unique.

    processPoints();


    if(debug && points.size() > 0 )
        writetofile();

    VoronoiDiagram * vd = new VoronoiDiagram();
    result = vd->DevideConquerConstruction(points);
    //vector<Point *> pointerPoints = getPoints(points);
    //result = vd->IncrementalConstruction(pointerPoints);
    result->method = 1;
    Invalidate(TRUE);

    if(debug && points.size() > 0 )
        writetofile();

    //ofstream outf("halfedges.address.txt");
    //streambuf *default_buf=cout.rdbuf();   
    //cout.rdbuf( outf.rdbuf() );   
    //for (unsigned int i = 0; i < points.size(); i++)
    //{
    //    Point p = points.at(i);
    //    cout << p.x() << "," << p.y() << endl;
    //}

    //for (unsigned int i = 0; i < result->halfedges.size(); i++)
    //{
    //    Halfedge *e = result->halfedges[i];
    //    if(e->prevEdge() != NULL)
    //    {
    //        //the address should the same
    //        cout << e->oriVertex()  << ",";
    //        cout << e->prevEdge()->oriVertex()  << endl;
    //    }
    //}
    //cout.rdbuf(default_buf);

    


    ///if(true)// if not throw exception above, should write the sorted points 
    //    return;
    //if(points.size() > 0 )
    //    writetofile();

}


void COpenGLDemoView::OnIncrementalConstruction()
{

    processPoints();


    if(debug && points.size() > 0 )
        writetofile();

    VoronoiDiagram * vd = new VoronoiDiagram();
    //result = vd->DevideConquerConstruction(points);
    vector<Point *> pointerPoints = getPoints(points);
    try
    {
        result = vd->IncrementalConstruction(pointerPoints);
    }
    catch (CException* e)
    {

    }
    
    
    result->method = 0;
    Invalidate(TRUE);

    if(debug && points.size() > 0 )
        writetofile();


}

TCHAR* CString2TCHAR(CString &str) 
{ 
    int iLen = str.GetLength(); 
    TCHAR* szRs = new TCHAR[iLen]; 
    lstrcpy(szRs, str.GetBuffer(iLen)); 
    str.ReleaseBuffer(); 
    return szRs; 
} 

char* CString2char(CString &str) 
{ 
    int len = str.GetLength(); 
    char* chRtn = (char*)malloc((len*2+1)*sizeof(char));//CString的长度中汉字算一个长度 
    memset(chRtn, 0, 2*len+1); 
    USES_CONVERSION; 
    strcpy((LPSTR)chRtn,OLE2A(str.LockBuffer())); 
    return chRtn; 
} 



void COpenGLDemoView::OnFileSaveAs()
{
    CString strExt = L".txt";                                // extention 
    CString strFilePath;
    CString strFilter;

    char * filename ;
    strFilter.Format(L"Text Files (*.txt)|*.txt|All Files (*.*)|*.*||");
    CFileDialog dlg(FALSE, NULL, L"", NULL, strFilter);
    if (!(dlg.DoModal() == IDOK))
        return;
    else
    {
        strFilePath = dlg.GetPathName();
        filename = CString2char(strFilePath);
        if (strFilePath.Find(strExt) == -1)//查找扩展名，如果没有输入则自动加
        {
            strFilePath += strExt;
        }
        if (access(filename, 0) == 0 )//should include  <io.h>
        {
            CString strQuery;
            strQuery.Format(L"%s has exists, REPLACE it?", strFilePath);
            if (IDNO == ::MessageBox(m_hWnd, strQuery, L"File replace QA", MB_ICONQUESTION | MB_YESNO) )
            {
                return;
            }
        }
    }
    //const char * filename ="testfile.txt";
  
    writetofile(filename);
    CString tip;
    tip.Format(L"%s Save OK.", strFilePath);
    MessageBox(tip);
}


void COpenGLDemoView::OnFileOpen()
{
    CString strExt = _T(".txt");                                // 扩展名
    CString strFilePath;
    CString strFilter;
    char*  filename;
    strFilter.Format(L"Text Files (*txt)|*txt|All Files (*.*)|*.*||");
    CFileDialog dlg(true, NULL, L"", NULL, strFilter);
    if (!(dlg.DoModal() == IDOK))
        return;
    else
    {
        strFilePath = dlg.GetPathName();
        filename = CString2char(strFilePath);
        if (strFilePath.Find(strExt) == -1)//查找扩展名，如果没有输入则自动加
        {
            strFilePath += strExt;
        }
        if (access(filename, 0) == 0 )
        {
        }else
        {
            CString strQuery;
            strQuery.Format(L"%s does Not exist", strFilePath);
            MessageBox(strQuery);
        }
    }
   
    //const char * filename ="testfile.txt";
    points = readfromfile(filename);
    pointset.clear();

    Invalidate(TRUE);
}


void COpenGLDemoView::OnFileSave()
{
    COpenGLDemoView::OnFileSaveAs();
}
