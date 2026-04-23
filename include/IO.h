#pragma once
#include <istream>
#include <ostream>
#include <cmath>      // log10, pow, floor
#include <sstream>    // std::stringstream
#include <iomanip>    // std::setw, std::setprecision, std::setfill
#include <string>     
#include <cstdlib>
#include "Event.h"



//   N
//   pdg px py pz E m x y z t Thydro cvx cvy cvz ipx ipy ipz wt
namespace IO {
  bool readEventText(std::istream& in, Event& ev);
  void writeEventText(std::ostream& out, const Event& ev);

  Event makeToyEvent(int event_id = 0);
}


// Fortran style float format // this is copied from Lido by Weiyao Ke
class ff {
public:
    ff(double x): value(x) {}
    const double value;
	friend std::ostream & operator<< (std::ostream & stream, const ff & x) {
		// So that the log does not scream
		if (x.value == 0.) {
		    stream << "0.000000E+00";
		    return stream;
		}
		int exponent = floor(log10(std::abs(x.value)));
		double base = x.value / pow(10, exponent);
		// Transform here
		base /= 10;
		exponent += 1;
		if (base >= 0){
			std::stringstream buff;
			buff << std::setw(8) << std::setprecision(6);
			buff << std::fixed << base;
			stream << std::setw(8) << buff.str();
		}
		else{
			std::stringstream buff;
			buff << std::setw(8) << std::setprecision(6);
			buff << std::fixed << base;
			//Checking if we have a leading minus sign
			std::string newbase = "-" + buff.str().substr(2, buff.str().size()-1);
			stream << std::setw(8) << newbase;
		}

		if (exponent >= 0) stream << "E+" << std::setw(2) << std::setfill('0') << exponent;
		else stream << "E-" << std::setw(2) << std::setfill('0') << std::abs(exponent);
		return stream;
	}
};
