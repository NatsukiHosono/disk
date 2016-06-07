#pragma once
//vector.h
//Last modefied 14/12/11

template <typename type> struct vec3{
	type x , y , z;
	//
	vec3(void) : x(0.0), y(0.0), z(0.0){//
	}
	vec3(type _x) : x(_x), y(0.0), z(0.0){//
	}
	vec3(type _x , type _y) : x(_x), y(_y), z(0.0){
	}
	vec3(type _x , type _y , type _z) : x(_x), y(_y), z(_z){
	}
	//operator =
	vec3& operator=(const vec3& a){
		this->x = a.x;
		this->y = a.y;
		this->z = a.z;
		return *this;
	}
	vec3& operator=(const type a){
		this->x = a;
		this->y = a;
		this->z = a;
		return *this;
	}
	//operator +vec3
	vec3 operator+() const{
		return *this;
	}
	//operator -vec3
	vec3 operator-() const {
		return vec3(-x, -y, -z);
	}
	//operator +=
	vec3 operator+=(const vec3& a){
		this->x += a.x;
		this->y += a.y;
		this->z += a.z;
		return *this;
	}
	//operator -=
	vec3 operator-=(const vec3& a){
		this->x -= a.x;
		this->y -= a.y;
		this->z -= a.z;
		return *this;
	}
	//operator *=
	vec3 operator*=(const type& a){
		this->x *= a;
		this->y *= a;
		this->z *= a;
		return *this;
	}
	//operator /=
	vec3 operator/=(const type& a){
		this->x /= a;
		this->y /= a;
		this->z /= a;
		return *this;
	}
	//operator +
	inline vec3 operator+(const vec3& a)const{
		return vec3(this->x + a.x , this->y + a.y , this->z + a.z);
	}
	//operator -
	inline vec3 operator-(const vec3& a)const{
		return vec3(this->x - a.x , this->y - a.y , this->z - a.z);
	}
	//operator *
	inline vec3 operator*(const type& a)const{
		return vec3(this->x * a , this->y * a , this->z * a);
	}
	//operator /
	inline vec3 operator/(const type& a)const{
		return vec3(this->x / a , this->y / a , this->z / a);
	}
	//
	~vec3(void){
	}
};

template <typename type> inline vec3<type> operator*(type a , const vec3<type>& y){
	return vec3<type> (a * y.x , a * y.y , a * y.z);
}

template <typename type> inline type abs(const vec3<type>& a){
	return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

template <typename type> inline type abs2(const vec3<type>& a){
	return (a.x * a.x + a.y * a.y + a.z * a.z);
}

template <typename type> inline type inner(const vec3<type> &a , const vec3<type> &b){
	return (a.x * b.x + a.y * b.y + a.z * b.z);
}

template <typename type> inline vec3<type> outer(const vec3<type> &a , const vec3<type> &b){
	return vec3<type>(a.y * b.z - a.z * b.y , a.z * b.x - a.x * b.z , a.x * b.y - a.y * b.x);
}
