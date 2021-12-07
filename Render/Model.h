//
// Created by zhangbei on 2021/10/24.
//

#ifndef TINYRENDER_MODEL_H
#define TINYRENDER_MODEL_H

#include <vector>
#include "../Math/Vector.h"
#include "../Math/MathDefine.h"
#include "../TGA_Image.h"
#include "../Math/Math.h"

class Model
{
public:
	Model(const char* FileName);

	~Model();

	int nverts();
	int nfaces();

    float2 uv(int iface, int nvert);
    float3 normal(float2 uvf);
    float3 normal(int iface, int nthvert);

	float3 vert(int i);
    float3 vert(int iface, int nthvert);
    float4 vert2(int iface, int nthvert);

    std::vector<int> face(int idx);
    TGAColor diffuse(float2 uv);
    float specular(float2 uvf);
private:
    void load_texture(std::string filename, const char *suffix, TGA_Image &img);
private:
    std::vector<float3> vets;
    std::vector<std::vector<int3> > faces; // attention, this Vec3i means vertex/uv/normal
    std::vector<float3> norms;
    std::vector<float2> uvs;

    TGA_Image diffusemap_;
    TGA_Image normalmap_;
    TGA_Image specularmap_;
};


#endif //TINYRENDER_MODEL_H
