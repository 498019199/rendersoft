//
// Created by zhangbei on 2021/10/25.
//

#ifndef TINYRENDER_MATHDEFINE_H
#define TINYRENDER_MATHDEFINE_H
#pragma once
#include <cstdint>

template<typename T, int SIZE>
class Vector_T;
typedef Vector_T<int32_t, 1> int1;
typedef Vector_T<int32_t, 2> int2;
typedef Vector_T<int32_t, 3> int3;
typedef Vector_T<int32_t, 4> int4;
typedef Vector_T<uint32_t, 1> uint1;
typedef Vector_T<uint32_t, 2> uint2;
typedef Vector_T<uint32_t, 3> uint3;
typedef Vector_T<uint32_t, 4> uint4;
typedef Vector_T<float, 1> float1;
typedef Vector_T<float, 2> float2;
typedef Vector_T<float, 3> float3;
typedef Vector_T<float, 4> float4;

template <typename T>
class Matrix4_T;
typedef Matrix4_T<float> float4x4;

template <typename T>
class Quaternion_T;
typedef Quaternion_T<float> Quaternion;

template<typename T>
class Color_T;
typedef Color_T<float> Color;
#endif //TINYRENDER_MATHDEFINE_H
