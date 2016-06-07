#pragma once

namespace PARAM{
	const double T_end = 0.71 * 3.0;
	const short int Ndim = 2;
	const double OUTPUT_TIME_INTERVAL = T_end / 50.0;
}

template <typename Ptcl> void SetupIC(vector<Ptcl>& ptcl, boundary_t<double>& boundary){
	boundary.setXRange(0.0, 1.0);
	boundary.setYRange(0.0, 1.0);
	
	const int N = 512;
	const double dx = 1./(double)(N);
	const double dens_low = 1.0;
	const double dens_hi  = 2.0;
	const double pres = 2.5;
	const double delta = 0.025;
	const double hcr = 5.0 / 3.0;
	for(double x = 0 ; x < 1.0 ; x += dx){
		for(double y = 0 ; y < 1.0 ; y += dx){
			const vec3<double> pos = vec3<double>(x, y, 0.0);
			double dens;
			vec3<double> vel;
			if(y < 0.25){
				dens  = dens_hi - 0.5 * (dens_hi - dens_low) * exp((y - 0.25) / delta);
				vel.x = -0.5 + 0.5 * exp((y - 0.25) / delta);
			}else if(y < 0.5){
				dens  = dens_low + 0.5 * (dens_hi - dens_low) * exp(-(y - 0.25) / delta);
				vel.x = 0.5 - 0.5 * exp(-(y - 0.25) / delta);
			}else if(y < 0.75){
				dens  = dens_low + 0.5 * (dens_hi - dens_low) * exp((y - 0.75) / delta);
				vel.x = 0.5 - 0.5 * exp((y - 0.75) / delta);
			}else{
				dens  = dens_hi - 0.5 * (dens_hi - dens_low) * exp(-(y - 0.75) / delta);
				vel.x = -0.5 + 0.5 * exp(-(y - 0.75) / delta);
			}
			vel.y = 0.01 * sin(4.0 * math::pi * x);
			const double mass = dens / (double)(N * N);
			const double eng  = pres / (hcr - 1.0) / dens;


			Ptcl ith;
			ith.createParticle(pos, vel, mass, dens, eng);
			ith.setEoS(&Monoatomic);
			ptcl.push_back(ith);
		}
	}
	return;
}

template <typename real> void ieSPH<real>::setSource(){

}

