// CP_PointVector.cpp ʵ����CP_Point2D��CP_Point3D��CP_Vector2D��CP_Vector3D
#include "stdafx.h"
#include "CP_PointVector.h"

#include <math.h>

// ////////////////////////////////////////////////////////////////////////////
// ʵ����CP_Point2D��ʼ
CP_Point2D::CP_Point2D(double newx, double newy):m_x(newx), m_y(newy)
{
} // ��CP_Point2D���캯������
// ʵ����CP_Point2D����
// ////////////////////////////////////////////////////////////////////////////

// ////////////////////////////////////////////////////////////////////////////
// ʵ����CP_Point3D��ʼ
CP_Point3D::CP_Point3D(double newx, double newy, double newz):m_x(newx), m_y(newy), m_z(newz)
{
} // ��CP_Point3D���캯������
// ʵ����CP_Point3D����
// ////////////////////////////////////////////////////////////////////////////

// ////////////////////////////////////////////////////////////////////////////
// ʵ����CP_Point3D��ʼ
CP_Vector2D::CP_Vector2D (double newx, double newy):m_x(newx), m_y(newy)
{
} // ��CP_Vector2D���캯������

CP_Vector2D& CP_Vector2D::operator += (const CP_Vector2D& v)
{ 
    m_x += v.m_x;
    m_y += v.m_y; 
    return *this; 
} //��Ա����operator +=����

CP_Vector2D& CP_Vector2D::operator -= (const CP_Vector2D& v)
{ 
    m_x -= v.m_x;
    m_y -= v.m_y; 
    return *this; 
} //��Ա����operator -=����

CP_Vector2D& CP_Vector2D::operator *= (double num)
{
    m_x *= num;
    m_y *= num; 
    return *this; 
} //��Ա����operator *=����

CP_Vector2D& CP_Vector2D::operator /= (double num)
{
    m_x /= num;  // ע������û�д������Ϊ0������
    m_y /= num; 
    return *this;
} //��Ա����operator /=����

double CP_Vector2D::operator ^(const CP_Vector2D& v)
{
    return( m_x* v.m_y- m_y* v.m_x );
} //��Ա����operator ^����

CP_Vector2D CP_Vector2D::operator - () const
{ 
    return CP_Vector2D (-m_x, -m_y); 
} //��Ա����operator -����

double CP_Vector2D::mf_getLength( )  const                             
{ 
    return sqrt(m_x*m_x + m_y*m_y); 
} //��Ա����mf_getLength����

CP_Vector2D CP_Vector2D::mf_getPerpendicularVector() const
{
    return CP_Vector2D(-m_y, m_x);
} //��Ա����mf_getPerpendicularVector����

void CP_Vector2D::mf_normalize( )
{
    double a = mf_getLength( );
    (*this) /= a; // ע��: ����û�д�������Ϊ0�����
} //��Ա����mf_normalize����

void CP_Vector2D::mf_setValue(double newx, double newy)
{
    m_x=newx;
    m_y=newy;
} //��Ա����mf_setValue����

// ʵ����CP_Vector2D����
// ////////////////////////////////////////////////////////////////////////////


// ////////////////////////////////////////////////////////////////////////////
// ʵ����CP_Vector3D��ʼ
CP_Vector3D::CP_Vector3D (double newx, double newy, double newz):m_x(newx), m_y(newy), m_z(newz)
{
} // ��CP_Vector3D���캯������

CP_Vector3D& CP_Vector3D::operator += (const CP_Vector3D& v)
{ 
    m_x += v.m_x;
    m_y += v.m_y;
    m_z += v.m_z;  
    return *this; 
} //��Ա����operator +=����

CP_Vector3D& CP_Vector3D::operator -= (const CP_Vector3D& v)
{
    m_x -= v.m_x;
    m_y -= v.m_y;
    m_z -= v.m_z; 
    return *this; 
} //��Ա����operator -=����

CP_Vector3D& CP_Vector3D::operator *= (double num)
{ 
    m_x *= num;
    m_y *= num;
    m_z *= num; 
    return *this; 
} //��Ա����operator *=����

CP_Vector3D& CP_Vector3D::operator /= (double num)
{
    num = 1.0/num;
    m_x *= num;
    m_y *= num;
    m_z *= num;
    return *this;
} //��Ա����operator /=����

CP_Vector3D& CP_Vector3D::operator ^= (const CP_Vector3D& v)
{ 
    double a =   m_y * v.m_z - m_z * v.m_y;
    double b = - m_x * v.m_z + m_z * v.m_x;
    double c =   m_x * v.m_y - m_y * v.m_x;

    m_x = a;
    m_y = b;
    m_z = c;
    return *this;
} //��Ա����operator ^=����

CP_Vector3D CP_Vector3D::operator - ( ) const
{ 
    return CP_Vector3D (-m_x, -m_y, -m_z); 
} //��Ա����operator -����

double CP_Vector3D::mf_getLength( )  const                             
{ 
    return sqrt(m_x*m_x + m_y*m_y + m_z*m_z); 
} //��Ա����mf_getLength����

CP_Vector3D CP_Vector3D::mf_getPerpendicularVector( ) const
{
    CP_Vector3D vecReturn;
    if( fabs(m_y)<fabs(m_z))
    {
        vecReturn.m_x=m_z;
        vecReturn.m_y=0.0;
        vecReturn.m_z=-m_x;
    }
    else
    {
        vecReturn.m_x=-m_y;
        vecReturn.m_y=m_x;
        vecReturn.m_z=0.0;
    }
    return vecReturn;
} //��Ա����mf_getPerpendicularVector����

void CP_Vector3D::mf_normalize( )
{
    double a = mf_getLength( );
    (*this) /= a; // ע��: ����û�д������Ϊ0�����
} //��Ա����mf_normalize����

void CP_Vector3D::mf_setValue(double newx, double newy, double newz)
{
    m_x=newx;
    m_y=newy;
    m_z=newz;
} //��Ա����mf_setValue����

// ʵ����CP_Vector3D����
// ////////////////////////////////////////////////////////////////////////////

CP_Point3D operator + (const CP_Point3D& pt, const CP_Vector3D& v)
{
    return CP_Point3D (pt.m_x + v.m_x, pt.m_y + v.m_y, pt.m_z + v.m_z); 
} //����operator +����

CP_Point3D operator - (const CP_Point3D& pt, const CP_Vector3D& v)
{
    return CP_Point3D (pt.m_x - v.m_x, pt.m_y - v.m_y, pt.m_z - v.m_z); 
} //����operator -����

CP_Vector3D operator - (const CP_Point3D& p, const CP_Point3D& q)
{
    return CP_Vector3D (p.m_x - q.m_x, p.m_y - q.m_y, p.m_z - q.m_z); 
} //����operator -����

CP_Point2D operator + (const CP_Point2D& p, const CP_Vector2D& v)
{
    return CP_Point2D (p.m_x + v.m_x, p.m_y + v.m_y); 
} //����operator +����

CP_Point2D operator - (const CP_Point2D& pt, const CP_Vector2D& v)
{
    return CP_Point2D (pt.m_x - v.m_x, pt.m_y - v.m_y); 
} //����operator -����

//��д��������ȣ�������ȼ����
bool operator == (const CP_Point2D& u, const CP_Point2D& v)
{
    if (u.m_x == v.m_x && u.m_y == v.m_y) 
        return true;
    return false;
}

CP_Vector2D operator - (const CP_Point2D& p, const CP_Point2D& q)
{
    return CP_Vector2D (p.m_x - q.m_x, p.m_y - q.m_y); 
} //����operator -����


CP_Vector2D operator + (const CP_Vector2D& u, const CP_Vector2D& v)
{
    return CP_Vector2D (u.m_x + v.m_x, u.m_y + v.m_y); 
} //����operator +����

CP_Vector2D operator - (const CP_Vector2D& u, const CP_Vector2D& v)
{
    return CP_Vector2D (u.m_x - v.m_x, u.m_y - v.m_y);
} //����operator -����

// ���
double  operator * (const CP_Vector2D& u, const CP_Vector2D& v)
{
    return u.m_x * v.m_x + u.m_y * v.m_y;
} //����operator *����

double operator ^ (const CP_Vector2D& u, const CP_Vector2D& v)
{
    return(u.m_x* v.m_y - u.m_y* v.m_x );
} //����operator ^����


CP_Vector2D operator * (const CP_Vector2D& v, double num)
{
    return CP_Vector2D (v.m_x * num, v.m_y * num);
} //����operator *����

CP_Vector2D operator / (const CP_Vector2D& v, double num)
{
    return CP_Vector2D (v.m_x / num, v.m_y / num); // ע��: ����û�д������Ϊ0�����
} //����operator /����

CP_Vector3D operator + (const CP_Vector3D& u, const CP_Vector3D& v)
{
   return CP_Vector3D(u.m_x + v.m_x, u.m_y + v.m_y, u.m_z + v.m_z);
} //����operator +����

CP_Vector3D operator - (const CP_Vector3D& u, const CP_Vector3D& v)
{
    return CP_Vector3D (u.m_x - v.m_x, u.m_y - v.m_y, u.m_z - v.m_z);
} //����operator -����

// ���
double operator * (const CP_Vector3D& u, const CP_Vector3D& v)
{
    return (u.m_x * v.m_x+u.m_y * v.m_y+ u.m_z * v.m_z);
} //����operator *����

// ���
CP_Vector3D operator ^ (const CP_Vector3D& u, const CP_Vector3D& v)
{
    return CP_Vector3D(u.m_y * v.m_z - u.m_z*v.m_y, 
                       -u.m_x*v.m_z+u.m_z*v.m_x,
                       u.m_x*v.m_y-u.m_y*v.m_x
                      );
} //����operator ^����

CP_Vector3D operator * (const CP_Vector3D& v, double num)
{
    return CP_Vector3D (v.m_x * num, v.m_y * num, v.m_z * num);
} //����operator *����

CP_Vector3D operator / (const CP_Vector3D& v, double num)
{
    num = 1.0/num; // ע��: ����û�д������Ϊ0�����
    return CP_Vector3D (v.m_x * num, v.m_y * num, v.m_z * num);
} //����operator /����


