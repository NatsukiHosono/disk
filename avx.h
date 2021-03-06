#pragma once
#ifdef __AVX__
//#warning AVX is available

typedef double   v4df __attribute__((vector_size(32)));
typedef float    v8sf __attribute__((vector_size(32)));
typedef long int v4di __attribute__((vector_size(32)));

//typedef double v4df[4] alignas(32);

/*
using v4df = double   __attribute__((vector_size(32)));
using v8sf = float    __attribute__((vector_size(32)));
using v4di = long int __attribute__((vector_size(32)));
*/

union AVX{
	//4 = 256 bit / sizeof(double) bit;
	v4df v;
	double x[4];
	//default constructer
	AVX(){
	}
	//destructer
	~AVX(){
	}
	//copy constructer
	AVX(const double& a){
		v = (v4df){a, a, a, a};
	}
	AVX(const v4df& a){
		v = a;
	}
	//
	double& operator [](const short int i){
		return x[i];
	}
	double operator [](const short int i) const{
		return x[i];
	}
	AVX operator = (const double& a){
		v = (v4df){a, a, a, a};
		return *this;
	}
	AVX operator = (const v4df& a){
		v = a;
		return *this;
	}
	//mul
	AVX operator * (const v4df a) const{
		return a * v;
	}
	AVX operator * (const AVX a) const{
		return a * v;
	}
	AVX operator * (const double a) const{
		return (v4df){a, a, a, a} * v;
	}
	friend AVX operator * (const double a, const AVX b){
		return (v4df){a, a, a, a} * b.v;
	}
	//div
	AVX operator / (const v4df a) const{
		return v / a;
	}
	AVX operator / (const AVX a) const{
		return v / a.v;
	}
	friend AVX operator / (const v4df a, const AVX b){
		return a / b.v;
	}

	//plus
	AVX operator + (const v4df a) const{
		return a + v;
	}
	AVX operator + (const AVX a) const{
		return a.v + v;
	}
	AVX operator + (const double a) const{
		return v + (v4df){a, a, a, a};
	}
	friend AVX operator + (const double a, const AVX b){
		return b.v + (v4df){a, a, a, a};
	}
	//minus
	AVX operator - (const v4df a) const{
		return a - v;
	}
	AVX operator - (const AVX a) const{
		return a.v - v;
	}
	AVX operator - (const double a) const{
		return v - (v4df){a, a, a, a};
	}
	friend AVX operator - (const double a, const AVX b){
		return (v4df){a, a, a, a} - b.v;
	}

	AVX operator += (const v4df a){
		return v += a;
	}
	AVX operator += (const AVX a){
		return v += a.v;
	}
	AVX invsqrt(){
		/*
		v4df y0 = __builtin_ia32_cvtps2pd256(__builtin_ia32_rsqrtps(__builtin_ia32_cvtpd2ps256(x2)));
		v4df y1 = (v4df){-0.5, -0.5, -0.5, -0.5} * (x2 * y0 * y0 - (v4df){3.0, 3.0, 3.0, 3.0}) * y0;
		return y1;
		*/
		#if 1
			return (v4df){1.0, 1.0, 1.0, 1.0} / __builtin_ia32_sqrtpd256(v);
		#else
			const v4df y0 = __builtin_ia32_cvtps2pd256(__builtin_ia32_rsqrtps(__builtin_ia32_cvtpd2ps256(v)));
			const v4df c0 = (v4df){ 1.0    ,  1.0    ,  1.0    ,  1.0    };
			const v4df c1 = (v4df){ 0.5    ,  0.5    ,  0.5    ,  0.5    };
			const v4df c2 = (v4df){ 3./  8.,  3./  8.,  3./  8.,  3./  8.};
			const v4df c3 = (v4df){ 5./ 16.,  5./ 16.,  5./ 16.,  5./ 16.};
			const v4df c4 = (v4df){35./128., 35./128., 35./128., 35./128.};

			const v4df h  = c0 - x2 * y0 * y0;
			const v4df y1 = y0 * (c0 + h * (c1 + h * (c2 + h * (c3 + h * (c4)))));
			//const v4df y1 = y0 * (c0 + h * (c1 + h * (c2 + h * (c3))));
			return y1;
		#endif
	}
	AVX sqrt(){
		return __builtin_ia32_sqrtpd256(v);
	}
	AVX max(const double& a) const{
		return __builtin_ia32_maxpd256(v, (v4df){a, a, a, a});
	}
	AVX pow(const unsigned int N) const{
		v4df ret = (v4df){1.0, 1.0, 1.0, 1.0};
		for(int i = 0 ; i < N ; ++ i){
			ret *= v;
		}
		return ret;
	}
};

inline AVX operator * (const v4df& a, const AVX& y){
	return a * y.v;
}
#endif
