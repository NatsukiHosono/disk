#pragma once
#include "vector.h"

template <typename type> struct mat33{//3x3 matrix
	type xx, xy, xz, yx, yy, yz, zx, zy, zz;
	mat33(void) : xx(0.0), xy(0.0), xz(0.0), yx(0.0), yy(0.0), yz(0.0), zx(0.0), zy(0.0), zz(0.0){
	}
	mat33(type _a) : xx(_a), xy(_a), xz(_a), yx(_a), yy(_a), yz(_a), zx(_a), zy(_a), zz(_a){
	}
	~mat33(){
	}
	mat33& operator = (const type a){
		xx = yy = zz = xy = xz = yx = yz = zx = zy = a;
		return *this;
	}
	mat33 operator + () const{
		return *this;
	}
	void PrintData(void){
		printf("%e %e %e\n", xx, xy, xz);
		printf("%e %e %e\n", yx, yy, yz);
		printf("%e %e %e\n", zx, zy, zz);
	}
	//operator -
	mat33 operator - () const{
		mat33<double> r;
		r.xx = xx;
		r.xy = xy;
		r.xz = xz;
		r.yx = yx;
		r.yy = yy;
		r.yz = yz;
		r.zx = zx;
		r.zy = zy;
		r.zz = zz;
		return r;
	}
	//operator +=
	mat33 operator += (const mat33 a){
		xx += a.xx;
		xy += a.xy;
		xz += a.xz;
		yx += a.yx;
		yy += a.yy;
		yz += a.yz;
		zx += a.zx;
		zy += a.zy;
		zz += a.zz;
		return *this;
	}
	//operator +=
	mat33 operator -= (const mat33 a){
		xx -= a.xx;
		xy -= a.xy;
		xz -= a.xz;
		yx -= a.yx;
		yy -= a.yy;
		yz -= a.yz;
		zx -= a.zx;
		zy -= a.zy;
		zz -= a.zz;
		return *this;
	}
	//operator *=
	mat33 operator *= (const type a){
		xx *= a;
		xy *= a;
		xz *= a;
		yx *= a;
		yy *= a;
		yz *= a;
		zx *= a;
		zy *= a;
		zz *= a;
		return *this;
	}
	//operator /=
	mat33 operator /= (const type a){
		xx /= a;
		xy /= a;
		xz /= a;
		yx /= a;
		yy /= a;
		yz /= a;
		zx /= a;
		zy /= a;
		zz /= a;
		return *this;
	}
	//operator +
	inline mat33 operator + (const mat33& a) const{
		mat33 r;
		r.xx = xx + a.xx;
		r.xy = xy + a.xy;
		r.xz = xz + a.xz;
		r.yx = yx + a.yx;
		r.yy = yy + a.yy;
		r.yz = yz + a.yz;
		r.zx = zx + a.zx;
		r.zy = zy + a.zy;
		r.zz = zz + a.zz;
		return r;
	}
	//operator -
	inline mat33 operator - (const mat33& a) const{
		mat33 r;
		r.xx = xx - a.xx;
		r.xy = xy - a.xy;
		r.xz = xz - a.xz;
		r.yx = yx - a.yx;
		r.yy = yy - a.yy;
		r.yz = yz - a.yz;
		r.zx = zx - a.zx;
		r.zy = zy - a.zy;
		r.zz = zz - a.zz;
		return r;
	}
	//operator * scalar
	inline mat33 operator * (const type& a) const{
		mat33 r;
		r.xx = a * xx;
		r.xy = a * xy;
		r.xz = a * xz;
		r.yx = a * yx;
		r.yy = a * yy;
		r.yz = a * yz;
		r.zx = a * zx;
		r.zy = a * zy;
		r.zz = a * zz;
		return r;
	}
	//operator * vec3
	inline vec3<type> operator * (const vec3<type>& a) const{
		vec3<type> r;
		r.x = xx * a.x + xy * a.y + xz * a.z;
		r.y = yx * a.x + yy * a.y + yz * a.z;
		r.z = zx * a.x + zy * a.y + zz * a.z;
		return r;
	}
	//operator * mat33
	inline mat33<type> operator * (const mat33<type>& a) const{
		mat33<type> r;
		r.xx = xx * a.xx + xy * a.yx + xz * a.zx;
		r.xy = xx * a.xy + xy * a.yy + xz * a.zy;
		r.xz = xx * a.xz + xy * a.yz + xz * a.zz;

		r.yx = yx * a.xx + yy * a.yx + yz * a.zx;
		r.yy = yx * a.xy + yy * a.yy + yz * a.zy;
		r.yz = yx * a.xz + yy * a.yz + yz * a.zz;

		r.zx = zx * a.xx + zy * a.yx + zz * a.zx;
		r.zy = zx * a.xy + zy * a.yy + zz * a.zy;
		r.zz = zx * a.xz + zy * a.yz + zz * a.zz;
		return r;
	}
	//operator / scalar
	inline mat33 operator / (const type& a) const{
		mat33 r;
		r.xx = xx / a;
		r.xy = xy / a;
		r.xz = xz / a;
		r.yx = yx / a;
		r.yy = yy / a;
		r.yz = yz / a;
		r.zx = zx / a;
		r.zy = zy / a;
		r.zz = zz / a;
		return r;
	}
	inline type det() const{
		return xx * yy * zz + yx * zy * xz + zx * xy * yz - xx * zy * yz - zx * yy * xz - yx * xy * zz;
	}
	inline mat33 inv(const int Ndim = 3) const{
		mat33 r;
		r.xx = (yy * zz - yz * zy);
		r.xy = (xz * zy - xy * zz);
		r.xz = (xy * yz - xz * yy);
		r.yx = (yz * zx - yx * zz);
		r.yy = (xx * zz - xz * zx);
		r.yz = (xz * yx - xx * yz);
		r.zx = (yx * zy - yy * zx);
		r.zy = (xy * zx - xx * zy);
		r.zz = (xx * yy - xy * yx);
		return r / ((*this).det() + 1.0e-16);
	}
	inline type trace() const{
		return xx + yy + zz;
	}
	inline friend mat33 operator * (const type& a, const mat33<type>& b){
		return b * a;
	}
	static mat33 I(){
		mat33<type> r;
		r.xx = r.yy = r.zz = 1.0;
		return r;
	}
};

template <typename type> inline mat33<type> transpose(const mat33<type>& a){
	mat33<type> r;
	r.xx = a.xx;
	r.xy = a.yx;
	r.xz = a.zx;
	r.yx = a.xy;
	r.yy = a.yy;
	r.yz = a.zy;
	r.zx = a.xz;
	r.zy = a.yz;
	r.zz = a.zz;
	return r;
}

template <typename type> inline mat33<type> tensor(const vec3<type> &a, const vec3<type> &b){
	mat33<type> r;
	r.xx = a.x * b.x;
	r.xy = a.x * b.y;
	r.xz = a.x * b.z;
	r.yx = a.y * b.x;
	r.yy = a.y * b.y;
	r.yz = a.y * b.z;
	r.zx = a.z * b.x;
	r.zy = a.z * b.y;
	r.zz = a.z * b.z;
	return r;
}

