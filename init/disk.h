#pragma once

namespace PARAM{
	const double T_end = 2.0 * math::pi * 10;
	const double OUTPUT_TIME_INTERVAL = T_end / 100.0;
	const short int Ndim = 2;
}

template <typename Ptcl> void SetupIC(vector<Ptcl>& ptcl, boundary_t<double>& boundary){
	const int N = 64;
	const double dx = 1./(double)(N);
	const double dens = 1.0;
	const double pres = 1.0e-6;
	const double r_out = 2.0;
	const double r_in  = 0.5;
	std::size_t cnt = 0;
	for(double x = -r_out ; x <= r_out ; x += dx){
		for(double y = -r_out ; y <= r_out ; y += dx){
			const double r = sqrt(x * x + y * y);
			if(r < r_in || r > r_out) continue;
			++ cnt;
		}
	}
	for(double x = -r_out ; x <= r_out ; x += dx){
		for(double y = -r_out ; y <= r_out ; y += dx){
			const double r = sqrt(x * x + y * y);
			if(r < r_in || r > r_out) continue;
			const vec3<double> pos = vec3<double>(x, y, 0.0);

			const double eng = pres / (1.4 - 1.0) / dens;
			const double mass = math::pi * (r_out * r_out - r_in * r_in) * dens /(double)(cnt);
			const double theta = atan2(y, x);
			//const vec3<double> vel = vec3<double>(-sin(theta), cos(theta), 0.0) * sqrt(r) / sqrt(r * r + 0.25 * 0.25);
			const vec3<double> vel = vec3<double>(-sin(theta), cos(theta), 0.0) / sqrt(r);
			Ptcl ith;
			ith.createParticle(pos, vel, mass, dens, eng);
			ith.setEoS(&Diatomic);
			ptcl.push_back(ith);
		}
	}
	return;
}

template <typename real> void ieSPH<real>::setSource(){
	
	const double r3 = math::cube(abs(r));
	const double eps3 = (r3 < 0.25 * 0.25 * 0.25) ? 0.25 * 0.25 * 0.25 : 0;
	ext_a = - 1.0 / (r3 + eps3) * r;
	/*
	const double r2 = abs2(r);
	ext_a = - 1.0 / (r2 + 0.25 * 0.25) * r / sqrt(r2);
	*/

	a += ext_a;
	teng_dot += inner(v, ext_a);
}

