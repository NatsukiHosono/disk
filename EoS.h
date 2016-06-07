#pragma once

namespace EoS{
	//Enumerate
	enum Type{
		Monoatomic,
		Diatomic,
	};
	//////////////////
	//abstract class
	//////////////////
	template <typename type> class EoS_t{
		public:
		EoS_t(){
			return;
		}
		~EoS_t(){
			return;
		}
		virtual type Pressure  (const type dens, const type eng) const = 0;
		virtual type SoundSpeed(const type dens, const type eng) const = 0;
		virtual type HeatCapacityRatio() const = 0;
	};
	//////////////////
	//EoSs
	//////////////////
	template <typename type> class IdealGas : public EoS_t<type>{
		const type hcr;//heat capacity ratio;
		public:
		IdealGas(const type _hcr) : hcr(_hcr){
		}
		inline type Pressure(const type dens, const type eng) const{
			//return math::max((hcr - 1.0) * dens * eng, 0.0);
			return (hcr - 1.0) * dens * eng;
		}
		inline type SoundSpeed(const type dens, const type eng) const{
			//return sqrt(math::max(hcr * (hcr - 1.0) * eng, 0.0));
			//return sqrt(math::abs(hcr * (hcr - 1.0) * eng));
			return sqrt(hcr * (hcr - 1.0) * eng);
		}
		inline type HeatCapacityRatio() const{
			return hcr;
		}
	};
}

static EoS::IdealGas<double> Diatomic(1.4);
static EoS::IdealGas<double> Monoatomic(5./3.);

