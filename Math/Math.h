//
// Created by admin on 2021/10/25.
//

#ifndef UNTITLED_MATH_H
#define UNTITLED_MATH_H
#include "MathDefine.h"
#include <limits>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <iterator>

namespace MathLib
{
// tool math
///////////////////////////////////////////////////////////////////////////////
    const float PI = 3.141592f;
    const float PI2 = 6.283185f;
    const float PIdiv2 = 1.570796f;
    const float SMALL_NUMBER		(1.e-8f);
    float const DEG90 = 1.570796f;			// 90 度
    float const DEG270 = -1.570796f;		// 270 度
    float const DEG45 = 0.7853981f;			// 45 度
    float const DEG5 = 0.0872664f;			// 5 度
    float const DEG10 = 0.1745329f;			// 10 度
    float const DEG20 = 0.3490658f;			// 20 度
    float const DEG30 = 0.5235987f;			// 30 度
    float const DEG60 = 1.047197f;			// 60 度
    float const DEG120 = 2.094395f;			// 120 度
    float const DEG40 = 0.6981317f;			// 40 度
    float const DEG80 = 1.396263f;			// 80 度
    float const DEG140 = 2.443460f;			// 140 度
    float const DEG160 = 2.792526f;			// 160 度
    float const DEG2RAD = 0.01745329f;			// 角度化弧度因数
    float const RAD2DEG = 57.29577f;			// 弧度化角度因数

    // 角度化弧度
	template <typename T>
	inline T deg2rad(const T& x) noexcept
	{
		return static_cast<T>(x * DEG2RAD);
	}
    // 弧度化角度
	template <typename T>
	inline T rad2deg(const T& x) noexcept
	{
		return static_cast<T>(x * RAD2DEG);
	}


	float Sin(float angle) noexcept;
	float Cos(float angle) noexcept;
	float Tan(float angle) noexcept;
	float Asin(float angle) noexcept;
	float Acos(float angle) noexcept;
	float Atan2(float x, float y) noexcept;
	void SinCos(float angle, float& fs, float& fc) noexcept;
	float Log(float x) noexcept;
	float Log10(float x) noexcept;

    // 平方
	template <typename T>
		inline T Sqrt(T x) noexcept
	{
		return static_cast<T>(std::sqrt(x));
	}
	template<typename T> inline 
		T Sqr(const T& data) noexcept
	{
		return data * data;
	}
    // 立方
	template<typename T> inline 
		T Cube(const T& data) noexcept
	{
		return Sqr(data) * data;
	}

    // 求绝对值
	template <typename T>
	inline T Abs(const T& lhs) noexcept
	{
		return (lhs > T(0) ? lhs : -lhs);
	}

    // 相等
	template <typename T>
	inline bool Equal(const T& lhs, const T& rhs) noexcept
	{
		return lhs == rhs;
	}
	template<>
	inline bool Equal<double>(const double& lhs, const double& rhs) noexcept
	{
		return (Abs<double>(lhs - rhs)
			<= std::numeric_limits<double>::epsilon());
	}
	template<>
	inline bool Equal<float>(const float& lhs, const float& rhs) noexcept
	{
		return (Abs<double>(lhs - rhs)
			<= std::numeric_limits<float>::epsilon());
	}

	//  平方根倒数速算法
	float RecipSqrt(float number) noexcept;

    // 限制 val 在 low 和 high 之间
	template <typename T>
	inline const T& Clamp(const T& val, const T& low, const T& high) noexcept;

    // 线性插值
    template <typename T>
    T Lerp(const T& lhs, const T& rhs, float s) noexcept;

    // 取出最大
    template <typename T>
    T Maximize(const T& lhs, const T& rhs) noexcept;

    // 取出最小
    template <typename T>
    T Minimize(const T& lhs, const T& rhs) noexcept;
// Color///////////////////////////////////////////////////////////////////////////////
	template <typename T>
	Color_T<T> Modulate(const Color_T<T>& lhs, const Color_T<T>& rhs) noexcept;

//向量//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 叉积
	template<typename T>
	T Cross(const Vector_T<T, 2> & lhs, const Vector_T<T, 2> & rhs) noexcept;
	template<typename T>
	Vector_T<T, 3> Cross(const Vector_T<T, 3> & lhs, const Vector_T<T, 3> & rhs) noexcept;
	template<typename T>
	Vector_T<T, 4> Cross(const Vector_T<T, 4> & lhs, const Vector_T<T, 4> & rhs) noexcept;

	// 点积
	template<typename T>
	typename T::value_type Dot(const T & lhs, const T & rhs) noexcept;

	// 长度的平方
	template<typename T>
	typename T::value_type LengthSq(const T & rhs) noexcept;

	// 求模，向量的长
	template<typename T>
	typename T::value_type Length(const T & rhs) noexcept;

	// 向量标准化
	template<typename T>
	T Normalize(const T & rhs) noexcept;

	// 获取重心坐标
    float3 ComputeBaryCentric2D(const float3& Point, const float3& A, const float3& B, const float3& C);
//// 矩阵//////////////////////////////////////////////////////////////////////////////////////////////////////////
	template<typename T>
	Matrix4_T<T> Mul(const Matrix4_T<T>& rhs, const Matrix4_T<T>& lhs) noexcept;

	// 矩阵转置
	template<typename T>
	Matrix4_T<T> Transpose(const Matrix4_T<T>& mat) noexcept;

	// 缩放
	template<typename T>
	Matrix4_T<T> MatrixScale(const T x, const T y, const T z) noexcept;
	template<typename T>
	Matrix4_T<T> MatrixScale(const Vector_T<T, 3>& n) noexcept;

	// 旋转x轴
	template<typename T>
	Matrix4_T<T> MatrixRotateX(const T & x) noexcept;

	// 旋转y轴
	template<typename T>
	Matrix4_T<T> MatrixRotateY(const T & y) noexcept;

	// 旋转z轴
	template<typename T>
	Matrix4_T<T> MatrixRotateZ(const T & z) noexcept;

	// 旋转
	template<typename T>
	Matrix4_T<T> MatrixRotate(const Vector_T<T, 3>& n, const T rotate) noexcept;

	//行列式值
	template<typename T>
	T Determinant(const Matrix4_T<T>& mat) noexcept;

	// 矩阵的逆
	template<typename T>
	Matrix4_T<T> Inverse(const Matrix4_T<T>& mat) noexcept;

	// 矩阵平移
	template<typename T>
	Matrix4_T<T> MatrixMove(const T x, const T y, const T z) noexcept;
	template<typename T>
	Matrix4_T<T> MatrixMove(const Vector_T<T, 3>& n) noexcept;

	// 获取观察矩阵
    template<typename T>
    Matrix4_T<T> LookAtRH(const Vector_T<T, 3>& vPostion, const Vector_T<T, 3>& vLookAt, const Vector_T<T, 3>& vUp) noexcept;
    // 正交投影
    template<typename T>
    Matrix4_T<T> OrthoOffCenterRH(T left,T right, T bottom, T top, T near, T far) noexcept;

	// 获取透视投影矩阵
    template<typename T>
    Matrix4_T<T> PerspectiveOffCenterRH(T left,T right, T bottom, T top, T near, T far) noexcept;
	template<typename T>
	Matrix4_T<T>  PerspectiveFovLH(T fFov, T fNearPlane, T fFarPlane, T fAspect) noexcept;

	// 行向量乘矩阵
	template<typename T>
	Vector_T<T, 4> MatrixMulVector(const Vector_T<T, 4>& vec, const Matrix4_T<T>& mat) noexcept;
    // 列向量乘矩阵
    template<typename T>
    Vector_T<T, 4> MatrixMulVector(const Matrix4_T<T>& mat, const Vector_T<T, 4>& vec) noexcept;

	template <typename T>
	Matrix4_T<T> ToMatrix(const Quaternion_T<T>& quat) noexcept;
//// 四元数//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 四元数乘法
	template <typename T>
	Quaternion_T<T> Mul(Quaternion_T<T> const& lhs, Quaternion_T<T> const& rhs) noexcept;

    // 四元数共轭
	template <typename T>
	Quaternion_T<T> Conjugate(Quaternion_T<T> const& rhs) noexcept;

	//四元数取逆
	template <typename T>
	Quaternion_T<T> Inverse(const Quaternion_T<T>& rhs) noexcept;

	// 矩阵转四元数
	template <typename T>
	Quaternion_T<T> ToQuaternion(const Matrix4_T<T>& mat) noexcept;

	// 向量转四元数
	template <typename T>
	Quaternion_T<T>  ToQuaternion(const Vector_T<T, 3>& rhs);

    // 四元数的对数
	template <typename T>
	Quaternion_T<T>  Exp(const Quaternion_T<T>& rhs);
	template <typename T>
	Quaternion_T<T>  In(const Quaternion_T<T>& quat);

	// 四元数旋转
	template <typename T>
	Quaternion_T<T> QuaternionRotate(const Vector_T<T, 3>& n, const T rotate);

    // 球面线性插值
	template <typename T>
	Quaternion_T<T>  Slerp(const Quaternion_T<T>& p1, const Quaternion_T<T>& p2, T ft);
	template <typename T>
	Quaternion Squad(const Quaternion_T<T>& q1, const Quaternion_T<T>& a,
		const Quaternion_T<T>& b, const Quaternion_T<T>& c, T ft);

	// 四元数旋转
	template <typename T>
	void ToYawPitchRoll(T& yaw, T& pitch, T& roll, Quaternion_T<T> const& quat) noexcept;
}
#include "Vector.h"
#include "Matrix.h"
#include "Quaternation.h"
#include "Color.h"
#endif //UNTITLED_MATH_H
