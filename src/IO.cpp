#include "IO.h"
#include <limits>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <iostream>

bool IO::readEventText(std::istream& in, Event& ev)
{
  ev.particles.clear();

  std::string line;

  // ---skip OSCAR header ---
  std::getline(in, line); // OSC1997A
  //std::cout << "[debug] line1 = " << line << std::endl; //debug
  std::getline(in, line); // final_id_p_x
  std::getline(in, line); // system line

  std::getline(in, line);
  std::istringstream iss(line);

  int dummy = 0;
  int N = 0;
  iss >> dummy >> N;  

  if (!iss || N <= 0) return false;

  ev.event_id += 1;
  ev.particles.reserve(N);

  for (int i = 0; i < N; ++i) {
    Particle p;
    int idx;
    double wt_dummy;//for UrQMD
    
    if (!(in >> idx >> p.pdg
              >> p.px >> p.py >> p.pz >> p.E >> p.m
              >> p.x  >> p.y  >> p.z  >> p.t
              >> p.Thydro
              >> p.cvx >> p.cvy >> p.cvz
              >> p.ipx >> p.ipy >> p.ipz
	      >> p.iE  //for UrQMD
              >> p.wt
	      >> wt_dummy //for UrQMD
	      )) {
//	    std::cout << "idx= " << idx << "  " << "p.wt= " << p.wt << std::endl; 
//	    std::cout << "[debug] failed when reading particle i = " << i << std::endl;
//  std::cout << "[debug] in.fail() = " << in.fail()
//            << ", in.eof() = " << in.eof()
//            << ", in.bad() = " << in.bad() << std::endl;
      return false;
    }

    ev.particles.push_back(p);
  }

  return true;
}



void IO::writeEventText(std::ostream& out, const Event& ev)
{
  const int Nparticles = static_cast<int>(ev.particles.size());

  // -------- OSCAR header --------
  out << "OSC1997A\n";
  out << "final_id_p_x\n";
  out << "    lbt  1.0alpha   208    82   208    82   aacm  0.1380E+04        1\n";
  out << "        1  " << std::setw(10) << Nparticles
      << "    0.001    0.001    1    1       1\n";

  // -------- particle list --------
  int i = 0;
  for (const auto& p : ev.particles) {

    out << std::setw(10) << std::setfill(' ') << i << "  "
        << std::setw(10) << std::setfill(' ') << p.pdg << "  "

        << ff(p.px) << "  "
        << ff(p.py) << "  "
        << ff(p.pz) << "  "
        << ff(p.E ) << "  "
        << ff(p.m ) << "  "

        << ff(p.x) << "  "
        << ff(p.y) << "  "
        << ff(p.z) << "  "
        << ff(p.t) << "  "

        << ff(p.Thydro) << "  "

        << ff(p.cvx) << "  "
        << ff(p.cvy) << "  "
        << ff(p.cvz) << "  "

        << ff(p.ipx) << "  "
        << ff(p.ipy) << "  "
        << ff(p.ipz) << "  "

	<< ff(p.iE) << "  "  // for UrQMD

        << ff(p.wt)  
	<< "  "
	<< ff(0.0)  //for UrQMD
	<< "\n";

    ++i;
  }
}

 

/////////////////////////
/*
bool IO::readEventText(std::istream& in, Event& ev) {
  ev.particles.clear();

  int N = 0;
  if (!(in >> N)) return false;
  if (N < 0) return false;

  ev.event_id += 1;
  ev.particles.reserve(static_cast<size_t>(N));

  for (int i = 0; i < N; ++i) {
    Particle p;
    if (!(in>> i  >>p.pdg
            >> p.px >> p.py >> p.pz >> p.E >> p.m
            >> p.x  >> p.y  >> p.z  >> p.t
            >> p.Thydro
            >> p.cvx >> p.cvy >> p.cvz
	    >> p.ipx >> p.ipy >> p.ipz
//	    >> p.ix >>p.iy >> p.iz
            >> p.wt)) {
      return false;
    }
    ev.particles.push_back(p);
  }
  return true;
}*/



/*
void IO::writeEventText(std::ostream& out, const Event& ev) {
  out << ev.particles.size() << "\n";
  for (const auto& p : ev.particles) {
    out
      << p.pdg << " "
      << p.px << " " << p.py << " " << p.pz << " " << p.E << " " << p.m << " "
      << p.x  << " " << p.y  << " " << p.z  << " " << p.t << " "
      << p.Thydro << " "
      << p.cvx << " " << p.cvy << " " << p.cvz << " "
      << p.ipx << " " << p.ipy << " " << p.ipz << " "
//    << p.ix << " " << p.iy << " " << p.iz << " "
      << p.wt << "\n";
  }
}*/

Event IO::makeToyEvent(int event_id) {
  Event ev;
  ev.event_id = event_id;

  // One charm quark + one pion as spectator
  Particle c;
  c.pdg = 4;
  c.px = 5.0; c.py = 1.0; c.pz = 2.0;
  c.m  = 1.8;
  c.E  = std::sqrt(c.px*c.px + c.py*c.py + c.pz*c.pz + c.m*c.m);
  c.Thydro = 0.16;
  c.cvx = 0.1; c.cvy = 0.0; c.cvz = 0.0;
  c.wt = 1.0;

  Particle pi;
  pi.pdg = 211;
  pi.px = 0.2; pi.py = 0.1; pi.pz = 0.0;
  pi.m = 0.139;
  pi.E = std::sqrt(pi.px*pi.px + pi.py*pi.py + pi.pz*pi.pz + pi.m*pi.m);
  pi.Thydro = 0.16;
  pi.wt = 1.0;

  ev.particles.push_back(c);
  ev.particles.push_back(pi);
  return ev;
}

