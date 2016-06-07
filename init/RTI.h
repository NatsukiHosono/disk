#pragma once

namespace PARAM{
	const double T_end = 5.0;
	const short int Ndim = 2;
	const double OUTPUT_TIME_INTERVAL = T_end / 250.0;
}

template <typename Ptcl> void SetupIC(vector<Ptcl>& ptcl, boundary_t<double>& boundary){
	boundary.setXRange(0.0, 0.25);
	//boundary.setYRange(0.0, 1.0);
	
	const int N = 1024;
	const double dx = 1./(double)(N);
	const double pres = 10.0/7.0;
	const double hcr = 1.4;
	for(double x = 0 ; x < 0.25 ; x += dx){
		for(double y = 0 ; y < 1.0 ; y += dx){
			const vec3<double> pos = vec3<double>(x, y, 0.0);
			double dens;
			vec3<double> vel;
			if(y < 0.5){
				//dens = 1.0 * pow(1 + (hcr - 1.0) / hcr * 1.0 * 0.5 * (y - 0.5) / pres, 1.0 / (hcr - 1.0));
				dens = 1.0;
			}else{
				//dens = 2.0 * pow(1 + (hcr - 1.0) / hcr * 2.0 * 0.5 * (y - 0.5) / pres, 1.0 / (hcr - 1.0));
				dens = 2.0;
			}
			if(0.3 < y && y < 0.7) vel.y = 0.1 * (1.0 + cos(8.0 * math::pi * (x + 0.25))) * (1.0 + cos(5.0 * math::pi * (y - 0.5))) * 0.25;
			const double mass = dens / (double)(N * N);
			const double eng  = (pres - 0.5 * dens * (y - 0.5)) / (hcr - 1.0) / dens;
			Ptcl ith;
			if(y < 0.1 || y > 0.9) ith.type = FREEZE;
			ith.createParticle(pos, vel, mass, dens, eng);
			ith.setEoS(&Diatomic);
			ptcl.push_back(ith);
		}
	}
	return;
}

template <typename real> void ieSPH<real>::setSource(){
	if(type != HYDRO) return;
	ext_a.y = - 0.5;
	a += ext_a;
}

