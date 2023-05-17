#ifndef litchi_paramStruct
#define litchi_paramStruct

#include <cmath>
#include <string>

//TODO doxygen
///Helper struct containing all necessary minkmap generation parameters
struct paramStruct{
        uint rankA{0}, rankB{0}, curvIndex{0}, numt{1}, Nside{0}, smooth{0}, NsideOut{0};
        double mint{0}, maxt{1}, maskThresh{0.9}, smoothRad{0};
        bool linThresh{true}, forceOutname{false}, sequence{false};
        std::string function{"trace"}, maskname{""};
};


#endif
