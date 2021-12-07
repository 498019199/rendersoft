//
// Created by zhangbei on 2021/10/24.
//

#include "Model.h"
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>

Model::Model(const char *FileName)
		:vets(),faces()
{
	std::ifstream in;
	in.open (FileName, std::ifstream::in);
	if (in.fail())
	{
		return;
	}

	std::string line;
	while (!in.eof())
	{
		std::getline(in, line);
		std::istringstream iss(line.c_str());
		char trash;
		if (!line.compare(0, 2, "v "))
		{
			iss >> trash;
			float3 v;
			for (int i=0;i<3;i++) iss >> v[i];
			vets.push_back(v);
		}
        else if (!line.compare(0, 3, "vn "))
        {
            iss >> trash >> trash;
            float3 n;
            for (int i = 0; i < 3; i++) iss >> n[i];
            norms.push_back(n);
        }
         else if (!line.compare(0, 3, "vt "))
         {
            iss >> trash >> trash;
            float2 uv;
            for (int i=0;i<2;i++) iss >> uv[i];
            uvs.push_back(uv);
        }
        else if (!line.compare(0, 2, "f "))
		{
            std::vector<int3> f;
			int index = 0;
			int3 tmp;
			int32_t itrash, idx;
			iss >> trash;
            while (iss >> tmp[0] >> trash >> tmp[1] >> trash >> tmp[2])
			{
                for (int i=0; i<3; i++) tmp[i]--; // in wavefront obj all indices start at 1, not zero
                f.push_back(tmp);
			}
			faces.push_back(f);
		}
	}

	std::cerr << "# v# " << vets.size() << " f# "  << faces.size() << std::endl;
	load_texture(FileName, "_diffuse.tga", diffusemap_);
	load_texture(FileName, "_nm.tga",      normalmap_);
	load_texture(FileName, "_spec.tga",    specularmap_);
}

void Model::load_texture(std::string filename, const char *suffix, TGA_Image &img)
{
    std::string texfile(filename);
    size_t dot = texfile.find_last_of(".");
    if (dot!=std::string::npos)
    {
        texfile = texfile.substr(0,dot) + std::string(suffix);
        std::cerr << "texture file " << texfile << " loading " << (img.read_tga_file(texfile.c_str()) ? "ok" : "failed") << std::endl;
        img.flip_vertically();
    }
}

int Model::nverts()
{
    return static_cast<int>(vets.size());
}

int Model::nfaces()
{
    return static_cast<int>(faces.size());
}

float3 Model::normal(float2 uvf)
{
    int2 uv(uvf[0]*normalmap_.get_width(), uvf[1]*normalmap_.get_height());
    TGAColor c = normalmap_.get(uv[0], uv[1]);
    float3 res;
    for (int i=0; i<3; i++)
        res[2-i] = (float)c[i]/255.f*2.f - 1.f;
    return res;
}

float3 Model::normal(int iface, int nthvert)
{
    int idx = faces[iface][nthvert][2];
    return MathLib::Normalize(norms[idx]);
}

float3 Model::vert(int i)
{
    return vets[i];
}

float3 Model::vert(int iface, int nthvert)
{
    return vets[faces[iface][nthvert][0]];
}

std::vector<int> Model::face(int idx)
{
    std::vector<int> face;
    for (int i=0; i<(int)faces[idx].size(); i++)
        face.push_back(faces[idx][i][0]);
    return face;
}

TGAColor Model::diffuse(float2 uvf)
{
    int2 uv(uvf[0]*diffusemap_.get_width(), uvf[1]*diffusemap_.get_height());
    return diffusemap_.get(uv.x(), uv.y());
}

float Model::specular(float2 uvf)
{
    int2 uv(uvf[0]*specularmap_.get_width(), uvf[1]*specularmap_.get_height());
    return specularmap_.get(uv[0], uv[1])[0]/1.f;
}

float2 Model::uv(int iface, int nthvert)
{
    return uvs[faces[iface][nthvert][1]];
}

float4 Model::vert2(int iface, int nthvert)
{
    float3 tmp = vert(iface, nthvert);
    return float4(tmp.x(), tmp.y(), tmp.z(), 1);
}
