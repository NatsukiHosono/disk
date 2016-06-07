#pragma once

namespace math{
	const double pi = 4.0 * atan(1.0);
	const double NaN = +0.0 / 0.0;
	template <typename type> inline type plus(const type a){
		return (a > 0) ? a : 0;
	}
	template <typename type> inline type square(const type a){
		return a * a;
	}
	template <typename type> inline type cube(const type a){
		return a * a * a;
	}
	template <typename type> inline type abs(const type a){
		return (a > 0) ? a : -a;
	}
	template <typename type> inline type min(const type a, const type b){
		return (a < b) ? a : b;
	}
	template <typename type> inline type max(const type a, const type b){
		return (a > b) ? a : b;
	}
	template <typename type> inline type powdim(const type a, const int Ndim){
		return pow(a, Ndim);
	}
	template <typename type> inline type rootdim(const type a, const int Ndim){
		return pow(a, 1.0 / (double)(Ndim));
	}
	template <typename type> inline type sign(const type a){
		return (a < 0) ? -1.0 : 1.0;
	}
}

