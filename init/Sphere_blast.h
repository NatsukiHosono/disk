#pragma once

namespace PARAM{
	const double T_end = 1.5;
	const short int Ndim = 2;
	const double OUTPUT_TIME_INTERVAL = T_end / 50.0;
}

template <typename Ptcl> void SetupIC(vector<Ptcl>& ptcl, boundary_t<double>& boundary){
	const double x_max = 1.0;
	const double y_max = 1.5;
	boundary.setXRange(0, x_max);
	boundary.setYRange(0, y_max);
	
	const int N = 512;
	const double dx = 1./(double)(N);
	const double hcr = 5./3.;
	const double dens = 1.0;
	for(double x = 0 ; x < x_max ; x += dx){
		for(double y = 0 ; y < y_max ; y += dx){
			double pres;
			vec3<double> vel;
			const vec3<double> pos = vec3<double>(x, y, 0.0);
			const double R = sqrt((x - 0.5 * x_max) * (x - 0.5 * x_max) + (y - 0.5 * y_max) * (y - 0.5 * y_max));
			if(R < 0.1){
				pres = 10.0;
			}else{
				pres = 0.1;
			}
			const double eng = pres / ((hcr - 1.0) * dens);
			const double mass = 1.0 / (N * N);

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

