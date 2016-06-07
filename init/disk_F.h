#pragma once

namespace PARAM{
	const double T_end = 2.0 * math::pi * 10;
	const double OUTPUT_TIME_INTERVAL = T_end / 100.0;
	const short int Ndim = 2;
}

template <typename Ptcl> void SetupIC(vector<Ptcl>& ptcl, boundary_t<double>& boundary){
	const int N = 65536/2;
	const double dens = 1.0;
	const double pres = 1.0e-6;
	const double r_out = 2.0;
	const double r_in  = 0.5;

	double angle = 0.0;
	for(int i = 0 ; i < N ; ++ i){
		const double m = (double)(i + 0.5) / (double)(N);
		const double r = sqrt(3.75 * m + 0.25);
		//const double angle = (4.0 * math::pi / (3.0 + sqrt(5.0))) * i;
		//const double angle = sqrt(math::pi * (sqrt(5.0) - 1.0)) * i;
		//
		const vec3<double> pos = r * vec3<double>(cos(angle), sin(angle), 0.0);
		const double eng = pres / (1.4 - 1.0) / dens;
		const double mass = math::pi * (r_out * r_out - r_in * r_in) * dens /(double)(N);
		const vec3<double> vel = vec3<double>(-sin(angle), cos(angle), 0.0) / sqrt(r);
		Ptcl ith;
		ith.createParticle(pos, vel, mass, dens, eng);
		ith.setEoS(&Diatomic);
		ptcl.push_back(ith);
		
		const double m_next = (double)(i + 1.5) / (double)(N);
		const double r_next = sqrt(3.75 * m_next + 0.25);
		angle += sqrt(2.0 * math::pi * (1.0 - r / r_next));
		
	}
	return;
}

template <typename real> void ieSPH<real>::setSource(){
	const double r3 = math::cube(abs(r));
	const double eps3 = (r3 < 0.25 * 0.25 * 0.25) ? 0.25 * 0.25 * 0.25 : 0;
	ext_a = - 1.0 / (r3 + eps3) * r;
	a += ext_a;
	teng_dot += inner(v, ext_a);
}

