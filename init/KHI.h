#pragma once

namespace PARAM{
	const double T_end = 0.71 * 3.0;
	const short int Ndim = 2;
	const double OUTPUT_TIME_INTERVAL = T_end / 100.0;
}

template <typename Ptcl> void SetupIC(vector<Ptcl>& ptcl, boundary_t<double>& boundary){
	boundary.setXRange(0.0, 1.0);
	boundary.setYRange(0.0, 1.0);
	
	const int N = 512;
	const double dx = 1./(double)(N);
	const double pres = 2.5;
	const double hcr = 1.4;
	for(double x = 0 ; x < 1.0 ; x += dx){
		for(double y = 0 ; y < 1.0 ; y += dx){
			const vec3<double> pos = vec3<double>(x, y, 0.0);
			double dens;
			vec3<double> vel;
			if(0.25 < y && y < 0.75){
				dens = 1.0;
				vel.x = -0.5;
			}else{
				dens = 2.0;
				vel.x = +0.5;
			}
			double sigma2 = 0.05 * 0.05 / 2.0;
			vel.y = 0.1 * sin(2.0 * 2.0 * math::pi * x) * exp(- (y - 0.25) * (y - 0.25) / (2.0 * sigma2)) + 0.01 * sin(6.0 * 2.0 * math::pi * x) * exp(- (y - 0.75) * (y - 0.75) / (2.0 * sigma2));
			const double mass = dens / (double)(N * N);
			const double eng  = pres / (hcr - 1.0) / dens;

			Ptcl ith;
			ith.createParticle(pos, vel, mass, dens, eng);
			ith.setEoS(&Diatomic);
			ptcl.push_back(ith);
		}
	}
	return;
}

template <typename real> void ieSPH<real>::setSource(){

}

