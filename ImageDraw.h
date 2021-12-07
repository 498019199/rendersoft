//
// Created by zhangbei on 2021/10/24.
//
#ifndef TINYRENDER_IMAGEDRAW_H
#define TINYRENDER_IMAGEDRAW_H
#include "TGA_Image.h"
#include "Math/Math.h"
#include <cmath>
struct IShader
{
    virtual float4 vertex(int iface, int nthvert) = 0;
    virtual bool fragment(float3 bar, TGAColor &color) = 0;
};

float4x4 viewport(int x, int y, int w, int h)
{
    float4x4 Viewport =  float4x4::Identity();
    Viewport(0,3) = x+w/2.f;
    Viewport(1,3) = y+h/2.f;
    Viewport(2,3) = 255.f/2.f;
    Viewport(0,0) = w/2.f;
    Viewport(1,1) = h/2.f;
    Viewport(2,2) = 255.f/2.f;
    return Viewport;
}

float4x4 projection(float coeff) {
    float4x4 Projection =  float4x4::Identity();
    Projection(3,2) = coeff;
    return Projection;
}

void Image_DDADrawLine(int x0,int y0, int x1, int y1, TGA_Image& image, TGAColor color)
{
    float slope = (y1 - y0) / static_cast<float>(x1 - x0);
    float y = static_cast<float>(y0);
    for (int i = x0; i < x1; ++i)
    {
        y += slope;
        image.set(i, std::lround(y), color);
    }
}

// https://www.cs.montana.edu/courses/spring2009/425/dslectures/Bresenham.pdf
void BresehamDrawLine(int x0,int y0, int x1, int y1, TGA_Image& image, TGAColor color)
{
    if (x0 == x1 && y0 == y1)
    {
        return;
    }

    bool Step = false;
    int dx = std::abs(x1 - x0);
    int dy = std::abs(y1 - y0);
    if (dx < dy)
    {
        Step = true;
        std::swap(x0, y0);
        std::swap(x1, y1);
        std::swap(dx, dy);
    }
    if (x0>x1)
    {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }

    int incrE = 2* dy;
    int incrNE = 2 * (dy - dx);
    int d = incrE - dx;

    int tx = (x1 - x0) > 0 ? 1 : -1;
    int ty = (y1 - y0) > 0 ? 1 : -1;
    for(int x = x0,y = y0; x != x1; x+=tx)
    {
        if (d < 0)
        {
            d += incrE;
        }
        else
        {
            y+=ty;
            d += incrNE;
        }

        if (Step)
        {
            image.set(y, x, color);
        }
        else
        {
            image.set(x, y, color);
        }
    }
}

void BresehamDrawLine(int2 t0, int2 t1, TGA_Image &image, TGAColor color)
{
    BresehamDrawLine(t0.x(), t0.y(), t1.x(), t1.y(), image, color);
}

//http://www.sunshine2k.de/coding/java/TriangleRasterization/TriangleRasterization.html
void BarycenticDrawTriangle(int2 t0, int2 t1, int2 t2, TGA_Image& image, TGAColor color)
{

    /* 包围盒 */
    int maxX = std::max(t0.x(), std::max(t1.x(), t2.x()));
    int minX = std::min(t0.x(), std::min(t1.x(), t2.x()));
    int maxY = std::max(t0.y(), std::max(t1.y(), t2.y()));
    int minY = std::min(t0.y(), std::min(t1.y(), t2.y()));
//http://www.sunshine2k.de/coding/java/PointInTriangle/PointInTriangle.html
// https://gdbooks.gitbooks.io/3dcollisions/content/Chapter4/point_in_triangle.html
    // 重心计算
    int2 vs1(t1.x() - t0.x(), t1.y() - t0.y());
    int2 vs2(t2.x() - t0.x(), t2.y() - t0.y());
    for (int x = minX; x <= maxX; x++)
    {
        for (int y = minY; y <= maxY; y++)
        {
            int2 q = int2(x - t0.x(), y - t0.y());
            float dev = (float)MathLib::Cross(vs1, vs2);
            float s = (float)MathLib::Cross(q, vs2) / dev;
            float t = (float)MathLib::Cross(vs1, q) / dev;
            if ( (s >= 0) && (t >= 0) && (s + t <= 1))
            {
                /* inside triangle */
                image.set(x, y, color); // attention, due to int casts t0.y+i != A.y
            }
        }
    }
}

void BarycenticDrawTriangle2(float3 Pts[3], float *zbuffer, TGA_Image& image, TGAColor color, int width)
{
    /* 包围盒 */
    float maxX = std::max(Pts[0].x(), std::max(Pts[1].x(), Pts[2].x()));
    float minX = std::min(Pts[0].x(), std::min(Pts[1].x(), Pts[2].x()));
    float maxY = std::max(Pts[0].y(), std::max(Pts[1].y(), Pts[2].y()));
    float minY = std::min(Pts[0].y(), std::min(Pts[1].y(), Pts[2].y()));

    // 重心计算
    for (float x = minX; x <= maxX; x++)
    {
        for (float y = minY; y <= maxY; y++)
        {
            float3 P(x, y, 0.f);
            float3 BaryCentric = MathLib::ComputeBaryCentric2D(P, Pts[0], Pts[1], Pts[2]);
            if (BaryCentric.x()<0 || BaryCentric.y()<0 || BaryCentric.z()<0)
            {
                continue;
            }

            for (int i=0; i<3; i++)
            {
                P.z() += Pts[i][2] * BaryCentric[i];
            }
            int zIndex = int(P.x()+P.y()*width);
            if (zbuffer[zIndex]<P.z())
            {
                zbuffer[zIndex] = P.z();
                image.set(P.x(), P.y(), color);
            }
        }
    }
}

// https://github.com/ssloy/tinyrenderer/wiki/Lesson-1:-Bresenham%E2%80%99s-Line-Drawing-Algorithm
void TETS_DrawLine(int x0,int y0, int x1, int y1, TGA_Image& image, TGAColor color)
{
    bool steep = false;
    if (std::abs(x0-x1)<std::abs(y0-y1))
    {
        std::swap(x0, y0);
        std::swap(x1, y1);
        steep = true;
    }
    if (x0>x1)
    {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }
    int dx = x1-x0;
    int dy = y1-y0;
    int derror2 = std::abs(dy)*2;
    int error2 = 0;
    int y = y0;
    for (int x=x0; x<=x1; x++)
    {
        if (steep)
        {
            image.set(y, x, color);
        }
        else
        {
            image.set(x, y, color);
        }

        error2 += derror2;
        if (error2 > dx)
        {
            y += (y1>y0?1:-1);
            error2 -= dx*2;
        }
    }
}

float3 world2screen(float3 v, float width, float height)
{
    return float3(int((v.x()+1.)*width/2.+.5), int((v.y()+1.)*height/2.+.5), v.z());
}

void  TETS_DrawTriangle1(int2 t0, int2 t1, int2 t2, TGA_Image &image, TGAColor color)
{
    BresehamDrawLine(t0, t1, image, color);
    BresehamDrawLine(t1, t2, image, color);
    BresehamDrawLine(t2, t0, image, color);
}

void  TETS_DrawTriangle2(int2 t0, int2 t1, int2 t2, TGA_Image &image, TGAColor color)
{
    // sort the vertices, t0, t1, t2 lower−to−upper (bubblesort yay!)
    if (t0.y() > t1.y()) std::swap(t0, t1);
    if (t0.y() > t2.y()) std::swap(t0, t2);
    if (t1.y() > t2.y()) std::swap(t1, t2);

    BresehamDrawLine(t0, t1, image, color);
    BresehamDrawLine(t1, t2, image, color);
    BresehamDrawLine(t2, t0, image, color);
}

// 扫面线法，把不规则三角形切割为上直角三角，下直角三角。一行一行填充像素点
void  TETS_DrawTriangle3(float2 t0, float2 t1, float2 t2, TGA_Image &image, TGAColor color)
{
    // 按照y的值排序
    if (t0.y() > t1.y()) std::swap(t0, t1);
    if (t0.y() > t2.y()) std::swap(t0, t2);
    if (t1.y() > t2.y()) std::swap(t1, t2);

    // 上三角形，按y行x方向以次填充
    float total_height = t2.y() - t0.y();
    for (float y = t0.y(); y <= t1.y(); y++)
    {
        float segment_height = t1.y() - t0.y() + 1;
        float invslope1 = (float)(y - t0.y()) / total_height;
        float invslope2 = (float)(y - t0.y()) / segment_height; // be careful with divisions by zero
        float2 A = t0 + (t2 - t0)*invslope1;
        float2 B = t0 + (t1 - t0)*invslope2;
        if (A.x() > B.x()) std::swap(A, B);
        for (float j = A.x(); j <= B.x(); j++)
        {
            image.set(j, y, color); // attention, due to int casts t0.y+i != A.y
        }
    }

    // 下三角形
    for (float y = t1.y(); y <= t2.y(); y++)
    {
        int segment_height = t2.y() - t1.y() + 1;
        float invslope1 = (float)(y - t0.y()) / total_height;
        float invslope2 = (float)(y - t1.y()) / segment_height; // be careful with divisions by zero
        int2 A = t0 + (t2 - t0)*invslope1;
        int2 B = t1 + (t2 - t1)*invslope2;
        if (A.x() > B.x()) std::swap(A, B);
        for (int j = A.x(); j <= B.x(); j++)
        {
            image.set(j, y, color); // attention, due to int casts t0.y+i != A.y
        }
    }
}

void TETS_DrawTriangle4(float2 t0, float2 t1, float2 t2, TGA_Image& image, TGAColor color)
{
    if (t0.y() == t1.y() && t0.y() == t2.y()) return; // i dont care about degenerate triangles
    if (t0.y() > t1.y()) std::swap(t0, t1);
    if (t0.y() > t2.y()) std::swap(t0, t2);
    if (t1.y() > t2.y()) std::swap(t1, t2);

    int total_height = t2.y() - t0.y();
    for (float i = 0; i < total_height; i++)
    {
        bool second_half = i > t1.y() - t0.y() || t1.y() == t0.y();
        float segment_height = second_half ? t2.y() - t1.y() : t1.y() - t0.y();
        float alpha = (float)i / total_height;
        float beta = (float)(i - (second_half ? t1.y() - t0.y() : 0)) / segment_height; // be careful: with above conditions no division by zero here
        float2 A = t0 + (t2 - t0) * alpha;
        float2 B = second_half ? t1 + (t2 - t1) * beta : t0 + (t1 - t0) * beta;
        if (A.x() > B.x()) std::swap(A, B);
        for (float j = A.x(); j <= B.x(); j++)
        {
            image.set(j, t0.y() + i, color); // attention, due to int casts t0.y()+i != A.y()
        }
    }
}

float3 barycentric(int2* pts, float2 P)
{
    float3 u = MathLib::Cross(float3(float(pts[2][0] - pts[0][0]), float(pts[1][0] - pts[0][0]), float(pts[0][0] - P[0])),
                              float3(float(pts[2][1] - pts[0][1]), float(pts[1][1] - pts[0][1]), float(pts[0][1] - P[1])));
    /* `pts` and `P` has integer value as coordinates
       so `abs(u[2])` < 1 means `u[2]` is 0, that means
       triangle is degenerate, in this case return something with negative coordinates */
    if (std::abs(u[2]) < 1) return
                float3(-1, 1, 1);
    return float3(1.f - (u.x() + u.y()) / u.z(), u.y() / u.z(), u.x() / u.z());
}

void  TETS_DrawTriangle5(int2* pts, TGA_Image& image, TGAColor color)
{
    int2 bboxmin(image.get_width() - 1, image.get_height() - 1);
    int2 bboxmax(0, 0);
    int2 clamp(image.get_width() - 1, image.get_height() - 1);
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            bboxmin[j] = std::max<float>(0, std::min<float>(bboxmin[j], pts[i][j]));
            bboxmax[j] = std::min<float>(clamp[j], std::max<float>(bboxmax[j], pts[i][j]));
        }
    }

    int2 P;
    for (P.x() = bboxmin.x(); P.x() <= bboxmax.x(); P.x()++)
    {
        for (P.y() = bboxmin.y(); P.y() <= bboxmax.y(); P.y()++)
        {
            float3 bc_screen = barycentric(pts, P);
            if (bc_screen.x() < 0 || bc_screen.y() < 0 || bc_screen.z() < 0)
                continue;
            image.set(P.x(), P.y(), color);
        }
    }
}

void  TETS_DrawTriangle6(int2 p1, int2 p2, int2 p3, TGA_Image& image, TGAColor color)
{
    int2 pts[] = {p1, p2, p3};
    TETS_DrawTriangle5(pts, image, color);
}

float3 barycentric(float3 A, float3 B, float3 C, float3 P)
{
    float3 s[2];
    for (int i=2; i--; )
    {
        s[i][0] = C[i]-A[i];
        s[i][1] = B[i]-A[i];
        s[i][2] = A[i]-P[i];
    }
    float3 u = MathLib::Cross(s[0], s[1]);
    if (std::abs(u[2])>1e-2) // dont forget that u[2] is integer. If it is zero then triangle ABC is degenerate
        return float3(1.f-(u.x()+u.y())/u.z(), u.y()/u.z(), u.x()/u.z());
    return float3(-1,1,1); // in this case generate negative coordinates, it will be thrown away by the rasterizator
}

void TETS_DrawTriangle7(float3 *pts, float *zbuffer, TGA_Image &image, TGAColor color, int width)
{
    float2 bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max());
    float2 bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    float2 clamp(image.get_width()-1, image.get_height()-1);
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<2; j++)
        {
            bboxmin[j] = std::max(0.f,      std::min(bboxmin[j], pts[i][j]));
            bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
        }
    }
    float3 P;
    for (P.x()=bboxmin.x(); P.x()<=bboxmax.x(); P.x()++)
    {
        for (P.y()=bboxmin.y(); P.y()<=bboxmax.y(); P.y()++)
        {
            float3 bc_screen  = barycentric(pts[0], pts[1], pts[2], P);
            if (bc_screen.x()<0 || bc_screen.y()<0 || bc_screen.z()<0) continue;
            P.z() = 0;
            for (int i=0; i<3; i++) P.z() += pts[i][2]*bc_screen[i];
            int _1 = int(P.x()+P.y()*width);
            if (zbuffer[int(P.x()+P.y()*width)]<P.z())
            {
                zbuffer[int(P.x()+P.y()*width)] = P.z();
                image.set(P.x(), P.y(), color);
            }
        }
    }
}

float2 proj(float4 p)
{
    return float2(p.x(), p.y());
}

// uAB+wAC=AP
float3 barycentric(float2 A, float2 B, float2 C, float2 P) {
    float3 s[2];
    for (int i=2; i--; ) {
        s[i][0] = C[i]-A[i];//AC
        s[i][1] = B[i]-A[i];//AB
        s[i][2] = A[i]-P[i];//PA
    }

    //[u,v,1]和[AB,AC,PA]对应的x和y向量都垂直，所以叉乘
    float3 u = MathLib::Cross(s[0], s[1]);
    if (std::abs(u[2])>1e-2) // dont forget that u[2] is integer. If it is zero then triangle ABC is degenerate
        return float3(1.f-(u.x()+u.y())/u.z(), u.y()/u.z(), u.x()/u.z());
    return float3(-1,1,1); // in this case generate negative coordinates, it will be thrown away by the rasterizator
}

void TETS_DrawTriangle8(float4 *pts, IShader &shader, TGA_Image &image, TGA_Image &zbuffer)
{
    float2 bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max());
    float2 bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            bboxmin[j] = std::min(bboxmin[j], pts[i][j]/pts[i][3]);
            bboxmax[j] = std::max(bboxmax[j], pts[i][j]/pts[i][3]);
        }
    }

    int2 P;
    TGAColor color;
    for (P.x()=bboxmin.x(); P.x()<=bboxmax.x(); P.x()++)
    {
        for (P.y()=bboxmin.y(); P.y()<=bboxmax.y(); P.y()++)
        {
            float2 A = proj(pts[0]/pts[0][3]);
            float2 B =proj(pts[1]/pts[1][3]);
            float2 C =proj(pts[2]/pts[2][3]);
            float2 P1(P.x(), P.y());
            float3 c = barycentric(A, B, C, P1);
            float z = pts[0][2]*c.x() + pts[1][2]*c.y() + pts[2][2]*c.z();
            float w = pts[0][3]*c.x() + pts[1][3]*c.y() + pts[2][3]*c.z();
            int frag_depth = std::max(0, std::min(255, int(z/w+.5)));
            if (c.x()<0 || c.y()<0 || c.z()<0 || zbuffer.get(P.x(), P.y())[0]>frag_depth) continue;
            bool discard = shader.fragment(c, color);
            if (!discard)
            {
                zbuffer.set(P.x(), P.y(), TGAColor(frag_depth));
                image.set(P.x(), P.y(), color);
            }
        }
    }
}
#endif //TINYRENDER_IMAGEDRAW_H
