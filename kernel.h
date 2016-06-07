#pragma once

template <typename type> inline type pow4(const type x){
	return x * x * x * x;
}

template <typename type> inline type pow5(const type x){
	return x * x * x * x * x;
}

template <typename type> inline type pow6(const type x){
	return x * x * x * x * x * x;
}

template <typename type> inline type pow7(const type x){
	return x * x * x * x * x * x * x;
}

template <typename type> inline type pow8(const type x){
	const type x2 = x * x;
	const type x4 = x2 * x2;
	return x4 * x4;
}

/*
template <typename type> struct kernel_t{
	type width, Cnorm;//kernel width and normalized constant
	const int Ndim;
	kernel_t(const int _Ndim) : Ndim(_Ndim){
		if(Ndim == 1){
			width = 2.;
			Cnorm = 8.0 / 3.0;
		}else if(Ndim == 2){
			width = 3.;
			Cnorm = 80.0 / (7.0 * math::pi);
		}else{
			width = 2.;
			Cnorm = 16.0 / math::pi;
		}
	}
	//W
	type operator() (const vec3<type> dr, const type h) const{
		const type H = width * h;
		const type s = abs(dr) / H;
		type r_value;
		r_value = math::plus(math::cube(1.0 - s)) - 4.0 * math::plus(math::cube(0.5 - s));
		r_value *= Cnorm / math::powdim(H, Ndim);
		return r_value;
	}
};
*/

/*
template <typename type> struct kernel_t{
	type width, Cnorm;//kernel width and normalized constant
	const int Ndim;
	kernel_t(int _Ndim) : Ndim(_Ndim){
		if(Ndim == 1){
			width = 2.;
			Cnorm = 3125.0 / 768.0;
		}else if(Ndim == 2){
			width = 2.;
			Cnorm = 46875.0 / (2398.0 * math::pi);
		}else{
			width = 2.;
			Cnorm = 15625.0 / (512.0 * math::pi);
		}
	}
	//W
	type W(vec3<type> dr, type h) const{
		const type H = width * h;
		const type s = abs(dr) / H;
		type r_value;
		r_value = pow4(math::plus(1.0 - s)) - 5.0 * pow4(math::plus(0.6 - s)) + 10.0 * pow4(math::plus(0.2 - s));
		r_value *= Cnorm / math::powdim(H, Ndim);
		return r_value;
	}
};
*/

/*
template <typename type> struct kernel_t{
	type width, Cnorm;//kernel width and normalized constant
	const int Ndim;
	kernel_t(const int _Ndim) : Ndim(_Ndim){
		if(Ndim == 1){
			width = 2.0;
			Cnorm = 6.075;
		}else if(Ndim == 2){
			width = 3.0;
			Cnorm = 15309.0 / (478.0 * math::pi);
		}else{
			width = 2.;
			Cnorm = 2187.0 / (40.0 * math::pi);
		}
	}
	//W
	type operator ()(const vec3<type> dr, const type h) const{
		const type H = width * h;
		const type s = abs(dr) / H;
		type r_value;
		r_value = pow5(math::plus(1.0 - s)) - 6.0 * pow5(math::plus(2.0 / 3.0 - s)) + 15.0 * pow5(math::plus(1.0 / 3.0 - s));
		r_value *= Cnorm / math::powdim(H, Ndim);
		return r_value;
	}
};
*/

/*
//Wendland C2
template <typename type> struct kernel_t{
	type width, Cnorm;//kernel width and normalized constant
	const int Ndim;
	kernel_t(int _Ndim) : Ndim(_Ndim){
		if(Ndim == 1){
			width = 1.620185;
			Cnorm = 1.25;
		}else if(Ndim == 2){
			width = 2.5;
			Cnorm = 7.0 / math::pi;
		}else{
			width = sqrt(15./4.);
			Cnorm = 21.0 / 2.0 / math::pi;
		}
	}
	//W
	type W(vec3<type> dr, type h) const{
		const type H = width * h;
		const type s = abs(dr) / H;
		type r_value;
		if(Ndim == 1){
			r_value = math::cube(math::plus(1.0 - s)) * (1.0 + 3.0 * s);
		}else{
			r_value = pow4(math::plus(1.0 - s)) * (1.0 + 4.0 * s);
		}
		r_value *= Cnorm / math::powdim(H, Ndim);
		return r_value;
	}
};
*/


//Wendland C6
template <typename type> struct kernel_t{
	type width, Cnorm;//kernel width and normalized constant
	const int Ndim;
	kernel_t(const int _Ndim) : Ndim(_Ndim){
		if(Ndim == 1){
			width = 2.5;
			Cnorm = 55./32.;
		}else if(Ndim == 2){
			width = 3.5;//
			Cnorm = 78. / (7. * math::pi);
		}else{
			width = 5.0;
			Cnorm = 1365. / (64. * math::pi);
		}
	}
	//W
	type operator() (const vec3<type> dr, const type h) const{
		const type H = width * h;
		const type s = abs(dr) / H;
		type r_value;
		if(Ndim == 1){
			r_value = pow7(math::plus(1. - s)) * (1.0 + 7.0 * s + 19.0 * s * s + 21.0 * s * s * s);
		}else{
			r_value = pow8(math::plus(1. - s)) * (1.0 + s * (8.0 + s * (25.0 + 32.0 * s)));
		}
		r_value *= Cnorm / math::powdim(H, Ndim);
		return r_value;
	}
	
	#ifdef __AVX__
	AVX operator() (const AVX dr, const AVX h) const{
		const AVX H = (v4df){width, width, width, width} * h;
		const AVX s = dr / H;
		AVX r_value;
		if(Ndim == 1){
			r_value = ((1. - s).max(0.0)).pow(7) * (1.0 + 7.0 * s + 19.0 * s * s + 21.0 * s * s * s);
		}else{
			r_value = ((1. - s).max(0.0)).pow(8) * (1.0 + s * (8.0 + s * (25.0 + 32.0 * s)));
		}
		return r_value * v4df{Cnorm, Cnorm, Cnorm, Cnorm} / H.pow(Ndim);
	}
	#endif
	
	vec3<type> gradW(const vec3<type> dr, const type h) const{
		const type H = width * h;
		const type s = abs(dr) / H;
		type r_value;
		if(Ndim == 1){
			r_value = pow6(math::plus(1.0 - s)) * (math::plus(1.0 - s) * (7.0 + 38.0 * s + 63.0 * s * s) - 7.0 * (1.0 + 7.0 * s + 19.0 * s * s + 21.0 * s * s * s));
		}else{
			r_value = pow7(math::plus(1.0 - s)) * (math::plus(1.0 - s) * (8.0 + 50.0 * s + 96.0 * s * s) - 8.0 * (1.0 + 8.0 * s + 25.0 * s * s + 32.0 * s * s * s));
		}
		r_value *= Cnorm / math::powdim(H, Ndim);
		return dr * r_value / (abs(dr) * H  + 1.0e-6 * h);
	}
};

