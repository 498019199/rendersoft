#ifndef _STX_MATH_COLOR_
#define _STX_MATH_COLOR_
#pragma once
#include "MathDefine.h"
#include "Vector.h"

///////////////////////////////////////////////////////////////////////////////
template <typename T>
class Color_T final
{
public:
	enum { elem_num = 4 };

	typedef T value_type;

	typedef typename Vector_T<T, elem_num>::pointer pointer;
	typedef typename Vector_T<T, elem_num>::const_pointer const_pointer;

	typedef typename Vector_T<T, elem_num>::reference reference;
	typedef typename Vector_T<T, elem_num>::const_reference const_reference;

	typedef typename Vector_T<T, elem_num>::iterator iterator;
	typedef typename Vector_T<T, elem_num>::const_iterator const_iterator;

public:
	constexpr Color_T() noexcept
	{
	}
	explicit constexpr Color_T(T const * rhs) noexcept
		: col(rhs)
	{
	}
	Color_T(Color_T const & rhs) noexcept;
	Color_T(Color_T&& rhs) noexcept;
	constexpr Color_T(T r, T g, T b, T a) noexcept
		: col(r, g, b, a)
	{
	}
	explicit Color_T(uint32_t dw) noexcept;

	iterator begin() noexcept
	{
		return col.begin();
	}
	constexpr const_iterator begin() const noexcept
	{
		return col.begin();
	}
	iterator end() noexcept
	{
		return col.end();
	}
	constexpr const_iterator end() const noexcept
	{
		return col.end();
	}
	reference operator[](size_t index) noexcept
	{
		return col[index];
	}
	constexpr const_reference operator[](size_t index) const noexcept
	{
		return col[index];
	}

	reference r() noexcept
	{
		return col[0];
	}
	constexpr const_reference r() const noexcept
	{
		return col[0];
	}
	reference g() noexcept
	{
		return col[1];
	}
	constexpr const_reference g() const noexcept
	{
		return col[1];
	}
	reference b() noexcept
	{
		return col[2];
	}
	constexpr const_reference b() const noexcept
	{
		return col[2];
	}
	reference a() noexcept
	{
		return col[3];
	}
	constexpr const_reference a() const noexcept
	{
		return col[3];
	}

	void RGBA(uint8_t& R, uint8_t& G, uint8_t& B, uint8_t& A) const noexcept;
	uint32_t ARGB() const noexcept;
	uint32_t ABGR() const noexcept;
	static const Color_T Zero() noexcept;


	Color_T& operator+=(Color_T<T> const & rhs) noexcept;
	Color_T& operator-=(Color_T<T> const & rhs) noexcept;
	Color_T& operator*=(T rhs) noexcept;
	Color_T& operator*=(Color_T<T> const & rhs) noexcept;
	Color_T& operator/=(T rhs) noexcept;

	Color_T& operator=(Color_T const & rhs) noexcept;
	Color_T& operator=(Color_T&& rhs) noexcept;

	Color_T const operator+() const noexcept;
	Color_T const operator-() const noexcept;

	bool operator==(Color_T<T> const & rhs) const noexcept;

private:
	Vector_T<T, elem_num> col;
};

template<typename U>
inline Color_T<U> operator*(const Color_T<U>& lhs, const U rhs) noexcept
{
	return Color_T<U>(lhs).operator*=(rhs);
}
template<typename U>
inline Color_T<U> operator/(const Color_T<U>& lhs, const U rhs) noexcept
{
	return Color_T<U>(lhs).operator/=(rhs);
}

template<typename U>
inline Color_T<U> operator+(const Color_T<U>& lhs, const Color_T<U>& rhs) noexcept
{
	return Color_T<U>(lhs).operator+=(rhs);
}
template<typename U>
inline Color_T<U> operator-(const Color_T<U>& lhs, const Color_T<U>& rhs) noexcept
{
	return Color_T<U>(lhs).operator-=(rhs);
}
template<typename U>
inline Color_T<U> operator*(const Color_T<U>& lhs, const Color_T<U>& rhs) noexcept
{
	return Color_T<U>(lhs).operator*=(rhs);
}

template<typename U>
inline Color_T<U> operator==(const Color_T<U>& lhs, const Color_T<U>& rhs) noexcept
{
	return lhs.operator==(rhs);
}
#endif _STX_MATH_COLOR_