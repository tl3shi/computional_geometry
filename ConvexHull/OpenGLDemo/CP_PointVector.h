// CP_PointVector.h: ������CP_Point2D��CP_Point3D��CP_Vector2D��CP_Vector3D
#ifndef CP_POINTVECTOR_H
#define CP_POINTVECTOR_H

#define PI2         6.28318530717958647692
#define PI          3.14159265358979323846
#define PI_2        1.57079632679489661923

class CP_Point2D
{
public:
    double    m_x, m_y;
public:
    //���캯��
    CP_Point2D (double newx=0.0, double newy=0.0);
};

class CP_Point3D
{
public:
    double    m_x, m_y, m_z;
public:
    //���캯��
    CP_Point3D (double newx=0.0, double newy=0.0, double newz=0.0);
};

class CP_Vector2D
{
public:
    double    m_x, m_y;

public:
    CP_Vector2D (double newx=0.0, double newy=0.0);

    // ��ֵ����
    CP_Vector2D& operator += (const CP_Vector2D& v);
    CP_Vector2D& operator -= (const CP_Vector2D& v);
    CP_Vector2D& operator *= (double num);
    CP_Vector2D& operator /= (double num);
    double operator ^(const CP_Vector2D& v);

    //��Ŀ��
    CP_Vector2D operator - ( ) const;

    double mf_getLength( ) const; // ȡ����
    CP_Vector2D mf_getPerpendicularVector( ) const; //�õ�һ����ֱ������

    void mf_normalize( ); // ��λ��
    void mf_setValue(double newx=0.0, double newy=0.0);
};

class CP_Vector3D
{
public:
    double    m_x, m_y, m_z;

public:
    CP_Vector3D (double newx=0.0, double newy=0.0, double newz=0.0);

    //��ֵ����
    CP_Vector3D& operator += (const CP_Vector3D& v);
    CP_Vector3D& operator -= (const CP_Vector3D& v);
    CP_Vector3D& operator *= (double num);
    CP_Vector3D& operator /= (double num);
    CP_Vector3D& operator ^= (const CP_Vector3D& v);

    //��Ŀ��
    CP_Vector3D operator - () const;

    double mf_getLength ( ) const; // ȡ����
    CP_Vector3D mf_getPerpendicularVector( ) const; //�õ�һ����ֱ������

    void mf_normalize( ); // ��λ��
    void mf_setValue(double newx=0.0, double newy=0.0,double newz=0.0);
};

extern CP_Point2D operator + (const CP_Point2D& pt, const CP_Vector2D& v);
extern bool operator == (const CP_Point2D& u, const CP_Point2D& v);
extern CP_Point2D operator - (const CP_Point2D& pt, const CP_Vector2D& v);
extern CP_Vector2D operator - (const CP_Point2D& p, const CP_Point2D& q);

extern CP_Point3D operator + (const CP_Point3D& pt, const CP_Vector3D& v);
extern CP_Point3D operator - (const CP_Point3D& pt, const CP_Vector3D& v);
extern CP_Vector3D operator - (const CP_Point3D& p, const CP_Point3D& q);

extern CP_Vector2D operator + (const CP_Vector2D& u, const CP_Vector2D& v); 
extern CP_Vector2D operator - (const CP_Vector2D& u, const CP_Vector2D& v); 
extern double  operator * (const CP_Vector2D& u, const CP_Vector2D& v); // ���
extern CP_Vector2D operator * (const CP_Vector2D& v, double num);
extern CP_Vector2D operator / (const CP_Vector2D& v, double num); 

extern CP_Vector3D operator + (const CP_Vector3D& u, const CP_Vector3D& v);
extern CP_Vector3D operator - (const CP_Vector3D& u, const CP_Vector3D& v);
extern double operator * (const CP_Vector3D& u, const CP_Vector3D& v); // ���
extern CP_Vector3D operator ^ (const CP_Vector3D& u, const CP_Vector3D& v); // ���
extern CP_Vector3D operator * (const CP_Vector3D& v, double num);
extern CP_Vector3D operator / (const CP_Vector3D& v, double num);

#endif