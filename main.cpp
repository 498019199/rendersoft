#include <filesystem>
#include "softrender/ImageDraw.h."
#include "softrender/Render/Model.h"

const int width  = 800;
const int height = 800;

// 光照方向
float3  LightDir(1,1,1);
// 相机
float3  v3Eye(0,-1,3);
float3  v3Center(0,0,0);
float3  v3Up(0,1,0);
Model* g_Model = nullptr;
float4x4 g_Viewport = float4x4::Identity();
float4x4 g_Projection = float4x4::Identity();
float4x4 g_ModelView = float4x4::Identity();

struct GouraudShader : public IShader
{
    // written by vertex shader, read by fragment shader
    float3 varying_intensity;

    virtual float4 vertex(int iface, int nthvert)
    {
        float4 gl_Vertex = g_Model->vert2(iface, nthvert); // read the vertex from .obj file
        auto mvp = g_Viewport * g_Projection * g_ModelView;     // transform it to screen coordinates
        gl_Vertex = MathLib::MatrixMulVector(mvp, gl_Vertex);
        varying_intensity[nthvert] = std::max(0.f, MathLib::Dot(g_Model->normal(iface, nthvert), LightDir)); // get diffuse lighting intensity
        return gl_Vertex;
    }

    virtual bool fragment(float3 bar, TGAColor &color)
    {
        // no, we do not discard this pixel
        float intensity = MathLib::Dot(varying_intensity, bar);   // interpolate intensity for the current pixel
        color = TGAColor(255, 255, 255) * intensity; // well duh
        return false;                              // no, we do not discard this pixel
    }
};

struct Shader :public IShader
{
    // written by vertex shader, read by fragment shader
    float3 varying_intensity;
    // same as above
    std::array<float2, 3> varying_uv;

    virtual float4 vertex(int iface, int nthvert)
    {
        varying_uv[nthvert] = g_Model->uv(iface, nthvert);
        varying_intensity[nthvert] = std::max(0.f,
                                              MathLib::Dot(g_Model->normal(iface, nthvert),LightDir)); // get diffuse lighting intensity
        float4 gl_Vertex = g_Model->vert2(iface, nthvert); // read the vertex from .obj file
        auto mvp = g_Viewport * g_Projection * g_ModelView;
        return MathLib::MatrixMulVector(mvp, gl_Vertex); // transform it to screen coordinates
    }

    virtual bool fragment(float3 bar, TGAColor &color)
    {
        float intensity = MathLib::Dot(varying_intensity, bar);   // interpolate intensity for the current pixel
        float2 uv = float2(MathLib::Dot(float3(varying_uv[0].x(), varying_uv[1].x(), varying_uv[2].x()), bar),
                           MathLib::Dot(float3(varying_uv[0].y(), varying_uv[1].y(), varying_uv[2].y()), bar));                 // interpolate uv for the current pixel
        color = g_Model->diffuse(uv)*intensity;      // well duh
        return false;                              // no, we do not discard this pixel
    }
};
int main(int argc, char** argv)
{
    auto strFilePath = std::filesystem::current_path();
    auto NewFilePath = strFilePath.parent_path() += "//Res//obj//african_head.obj";
    g_Model = new Model(NewFilePath.string().c_str());

    TGA_Image image(width, height, TGA_Image::RGB);
/*
//      https://github.com/ssloy/tinyrenderer/wiki/Lesson-1:-Bresenham%E2%80%99s-Line-Drawing-Algorithm
//    for (int i = 0; i < model->nfaces(); ++i)
//    {
//        int3 indexs = model->face(i);
//        for (int j = 0; j < 3; ++j)
//        {
//            float3 v0 = model->vert(indexs[j]);
//            float3 v1 = model->vert(indexs[(j+1)%3]);
//            int x0 = (v0.x()+1.)*width/2.;
//            int y0 = (v0.y()+1.)*height/2.;
//            int x1 = (v1.x()+1.)*width/2.;
//            int y1 = (v1.y()+1.)*height/2.;
//            BresehamDrawLine(x0, y0, x1, y1, image, white);
//        }
//    }

//https://github.com/ssloy/tinyrenderer/wiki/Lesson-2:-Triangle-rasterization-and-back-face-culling
//    int2 t0[3] = {int2(10, 70),   int2(50, 160),  int2(70, 80)};
//    int2 t1[3] = {int2(180, 50),  int2(150, 1),   int2(70, 180)};
//    int2 t2[3] = {int2(180, 150), int2(120, 160), int2(130, 180)};
//    BarycenticDrawTriangle(t0[0], t0[1], t0[2], image, red);
//    BarycenticDrawTriangle(t1[0], t1[1], t1[2], image, white);
//    BarycenticDrawTriangle(t2[0], t2[1], t2[2], image, green);
//    for (int i=0; i<model->nfaces(); i++)
//    {
//        int3 face = model->face(i);
//        int2 screen_coords[3];
//        for (int j=0; j<3; j++)
//        {
//            float3 world_coords = model->vert(face[j]);
//            //屏幕坐标    (-1,-1)映射为(0,0)  （1,1）映射为(width,height)
//            screen_coords[j] = int2((world_coords.x()+1.)*width/2., (world_coords.y()+1.)*height/2.);
//        }
//        BarycenticDrawTriangle(screen_coords[0], screen_coords[1], screen_coords[2], image, TGAColor(rand()%255, rand()%255, rand()%255, 255));
//    }

// error
//https://github.com/ssloy/tinyrenderer/wiki/Lesson-3:-Hidden-faces-removal-(z-buffer)
//    float3 light(0, 0, -1);
//    float *zbuffer = new float[width*height];
//    float SetValue = -std::numeric_limits<float>::max();
//    memset(zbuffer, SetValue, sizeof(float) * width*height);
//    for (int i=0; i<model->nfaces(); i++)
//    {
//        int3 face = model->face(i);
//        float3 world_coords[3];
//        float3 screen_coords[3];
//        for (int j=0; j<3; j++)
//        {
//            world_coords[j] = model->vert(face[j]);
//            screen_coords[j] = world2screen(model->vert(face[j]), width, height);
//        }
//
//        // 用世界坐标计算法向量
//        float3 norm = MathLib::Cross(world_coords[2] - world_coords[0], world_coords[1] - world_coords[0]);
//        norm = MathLib::Normalize(norm);
//        // 光照强度=法向量*光照方向   即法向量和光照方向重合时，亮度最高
//        float intensity = MathLib::Dot(norm, light);
//        if (intensity > 0)
//        {
//            auto color = TGAColor(intensity*255,intensity*255,intensity*255,255);
//            // 有点问题
//            BarycenticDrawTriangle2(screen_coords, zbuffer, image, color, width);
//            //TETS_DrawTriangle7(screen_coords, zbuffer, image, color, width);
//        }
//    }
*/
    g_ModelView = MathLib::LookAtRH(v3Center, v3Eye, v3Up);
    g_Viewport = viewport(width/8, height/8, width*3/4, height*3/4);
    g_Projection = projection(-1.f/ MathLib::Length(v3Eye-v3Center));
    LightDir = MathLib::Normalize(LightDir);
    TGA_Image zbuffer(width, height, TGA_Image::GRAYSCALE);

    Shader shader;
    for (int i=0; i<g_Model->nfaces(); i++)
    {
        float4 screen_coords[3];
        for (int j=0; j<3; j++)
        {
            screen_coords[j] = shader.vertex(i, j);
        }
        TETS_DrawTriangle8(screen_coords, shader, image, zbuffer);
    }

    //auto ProjectionMat = MathLib::PerspectiveFovLH();
    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    strFilePath += "//output.tga";
    image.write_tga_file(strFilePath.string().c_str());
    return 0;
}

