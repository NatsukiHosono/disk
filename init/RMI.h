#pragma once

namespace PARAM{
	const double T_end = 1.0;
	const short int Ndim = 2;
	const double OUTPUT_TIME_INTERVAL = T_end / 50.0;
}

template <typename Ptcl> void SetupIC(vector<Ptcl>& ptcl, boundary_t<double>& boundary){
	boundary.setXRange(0.0, 0.5);
	boundary.setYRange(0.0, 1.0);
	
	const int N = 512;
	const double dx = 1./(double)(N);
	const double hcr = 1.4;
	for(double x = 0 ; x < 0.5 ; x += dx){
		for(double y = 0 ; y < 1.0 ; y += dx){
			const vec3<double> pos = vec3<double>(x, y, 0.0);
			double dens;
			vec3<double> vel;
			double eng;
			if(0.25 < y && y < 0.35){
				const double pres = 1.0;
				dens = 4.0;
				eng  = pres / (hcr - 1.0) / dens;
			}else if(y < 0.5 + 0.05 * cos(2.0 * math::pi * x / 0.5)){
				const double pres = 0.1;
				dens = 1.0;
				eng  = pres / (hcr - 1.0) / dens;
			}else{
				const double pres = 0.1;
				dens = 2.0;
				eng  = pres / (hcr - 1.0) / dens;
			}
			const double mass = dens / (double)(N * N);
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

