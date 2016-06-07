#pragma once

namespace PARAM{
	const double T_end = 2.0 * math::pi * 10.0;
	const double OUTPUT_TIME_INTERVAL = 2.0 * math::pi / 10.0;
	const short int Ndim = 2;
}

template <typename Ptcl> void SetupIC(vector<Ptcl>& ptcl, boundary_t<double>& boundary){
	const int N = 128;
	const double dens = 1.0;
	const double pres = 1.0e-6;
	const double r_out = 2.0;
	const double r_in  = 0.5;
	std::size_t cnt = 0;

	const double dx = r_out / (double)(N);
	const int Np = 6;
	for(int i = floor(r_in / dx) ; i <= N ; ++ i){
		for(int j = 0 ; j < i * Np ; ++ j){
			++ cnt;
		}
	}
	
	for(int i = floor(r_in / dx) ; i <= N ; ++ i){
		const double r_ = i * dx;
		const double angle = (double)(rand()) / (double)(RAND_MAX) * 0.0;
		//
		for(int j = 0 ; j < i * Np ; ++ j){
			const double phi_ = j * 2.0 * math::pi / (double)(i * Np);
			const vec3<double> pos = vec3<double>(r_ * cos(phi_ + angle), r_ * sin(phi_ + angle), 0.0);
			const double eng = pres / (1.4 - 1.0) / dens;
			const double mass = math::pi * (r_out * r_out - r_in * r_in) * dens /(double)(cnt);
			const double theta = atan2(pos.y, pos.x);
			const vec3<double> vel = vec3<double>(-sin(theta), cos(theta), 0.0) / sqrt(abs(pos));
			Ptcl ith;
			ith.createParticle(pos, vel, mass, dens, eng);
			ith.setEoS(&Diatomic);
			ith.tag = i;
			ptcl.push_back(ith);
		}
	}
}

template <typename real> void ieSPH<real>::setSource(){
	const double r3 = math::cube(abs(r));
	if(r3 < 0.253 * 0.253 * 0.253) exit(0);
	const double eps3 = (r3 < 0.25 * 0.25 * 0.25) ? 0.25 * 0.25 * 0.25 : 0;
	ext_a = - 1.0 / (r3 + eps3) * r;
	
	ext_a -= 1.0 * inner(r, v) / abs(r) / dt * r / abs(r);

	a += ext_a;
}

