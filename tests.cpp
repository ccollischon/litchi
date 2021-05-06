#include <cassert>
#include <vector>
#include <iostream>

#include "healpix_cxx/healpix_map.h"
#include "healpix_cxx/healpix_map_fitsio.h"

#include "tensor2D.hpp"
#include "tensorOperations.hpp"
#include "tensorFamily.hpp"
#include "minkTensorIntegrand.hpp"
#include "litchi_pulp.hpp"
#include "litchi_peel.hpp"
#include "litchi_eat.hpp"

using namespace std;
int main ()
{
    
    const double pi = 3.14159265358979;
    
    //Test midpoint, geometry helpers
    pointing A(pi/2,0);
    pointing B(pi/4,0);
    pointing C(pi/2,pi/4);
    pointing D(pi/2,pi/2);
    
    assert( abs(midpoint(A,D).theta-pi/2)<1e-13 && abs(midpoint(A,D).phi-pi/4)<1e-13 && "midpoint broken!");
    assert( abs(midpoint(C,C).theta-C.theta)<1e-13 && abs(midpoint(C,C).phi-C.phi)<1e-13 && "midpoint broken when entering same vector twice!");
    assert( abs(midpoint(A,B).theta*(1/pi)-0.375)<1e-13 && abs(midpoint(A,B).phi)<1e-13 && "midpoint broken!");
    
    
    flipPointing(A);
    assert( abs(A.theta-pi/2) < 1e-12 && abs(A.phi-pi) < 1e-12 && "flipPointing broken!" );
    pointing n1 = getN_cartesian(B,A);
    minkTensorIntegrand asdf(0,2,0,A,B);
    
    assert( abs(sphereArea(A, B, C) - sphereArea(A,C,B)) < 1e-13 && abs(sphereArea(A, B, C) - sphereArea(B,C,A)) <1e-13 && abs(sphereArea(A,B,C)-sphereArea(B,A,C))<1e-13 && "sphereArea not commutative!");
    assert( abs(sphereArea(A,D,pointing(0,0))/pi - 0.5) < 1e-13 && "sphereArea broken!" );
    vec3 axis(1,0,0);
    vec3 input(0,1,0);
    vec3 testvec = rotateAroundAxis(axis, input, pi/2);
    assert((abs(testvec.x)<1e-14 && abs(testvec.y)<1e-14 && abs(testvec.z-1)<1e-14) && "rotateAroundAxis() broken!");
    
    
    Healpix_Map<double> map = read_Healpix_map_from_fits<double>("../litchi/COM_CMB_IQU-smica_2048_R3.00_hm1.fits", 1, 2);
    //Test ispolar
    minkmapSphere testmink(map);
    normalHealpixInterface interface(testmink);
    fix_arr<int, 4> neighbors; 
    fix_arr<double, 4> weight; 
    pointing north(0,0);
    map.get_interpol(north,neighbors,weight);
    for(uint i=0;i<neighbors.size();i++)
    {
        assert(interface.ispolar(neighbors[i])==1 && "ispolar RING north broken!");
    }
    pointing south(pi,0);
    map.get_interpol(south,neighbors,weight);
    for(uint i=0;i<neighbors.size();i++)
    {
        assert(interface.ispolar(neighbors[i])==2 && "ispolar RING south broken!");
    }
    pointing random(4.76,1.345);
    random.normalize();
    map.get_interpol(random,neighbors,weight);
    for(uint i=0;i<neighbors.size();i++)
    {
        assert(!interface.ispolar(neighbors[i]) && "ispolar RING broken (false-positive)!");
    }
    
    
    map.get_interpol(north,neighbors,weight);
    for(uint i=0;i<neighbors.size();i++)
    {
        assert(interface.ispolar(neighbors[i])==1 && "ispolar NEST north broken!");
    }
    map.get_interpol(south,neighbors,weight);
    for(uint i=0;i<neighbors.size();i++)
    {
        assert(interface.ispolar(neighbors[i])==2 && "ispolar NEST south broken!");
    }
    random.normalize();
    map.get_interpol(random,neighbors,weight);
    for(uint i=0;i<neighbors.size();i++)
    {
        assert(!interface.ispolar(neighbors[i]) && "ispolar NEST broken (false-positive)!");
    }
    
    //Test more pointing-related functions
    pointing r(0.5, 0);
    pointing n(0.5, pi);
    
    pointing interpolation  = interpPointing(r, 1.5, n, 0.5, 1);
    assert(abs(interpolation.theta)<1e-12 && "interpPointing() broken!");
    interpolation  = interpPointing(n, 0.5, r, 1.5, 1);
    assert(abs(interpolation.theta)<1e-12 && "interpPointing() broken!");
    
    //Test tensor functionality
    pointing nr(0.45*pi, 0.5*pi);
    minkTensorIntegrand mytensor(0,2,0,r,n);
    minkTensorIntegrand mytensor2(0,2,0,nr,nr);
    minkTensorIntegrand mytensor3(1,1,0,r,nr);
    tensor2D testtensor(mytensor);
    assert( abs(mytensor.accessElement({1,1})-(pi*pi)) < 1e-13 && "minkTensorIntegrand::accessElement() broken!" );
    assert( abs(mytensor2.accessElement({0,0})-(0.45*pi*0.45*pi)) < 1e-13 && "minkTensorIntegrand::accessElement() broken!" );
    assert( abs(mytensor3.accessElement({1,0})-(0.5* (0.5*0.5*pi + 0*0.45*pi))) < 1e-13 && "minkTensorIntegrand::accessElement() broken!" );
    assert( abs(testtensor.accessElement({0,1}) - (0.5*pi)) < 1e-13 && "tensor2D::accessElement() or assign() broken!");
    
    testtensor = mytensor2;
    tensor2D zahlenfriedhof (mytensor+mytensor2);
    assert( abs(zahlenfriedhof.accessElement({1,1})-mytensor.accessElement({1,1})-mytensor2.accessElement({1,1}))<1e-13 && "minkTensorSum broken!" );
    assert( abs(zahlenfriedhof.accessElement({0,0})-mytensor.accessElement({0,0})-mytensor2.accessElement({0,0}))<1e-13 && "minkTensorSum broken!" );
    zahlenfriedhof = mytensor + 3*mytensor2;
    assert( abs(zahlenfriedhof.accessElement({1,1})-mytensor.accessElement({1,1})-3*mytensor2.accessElement({1,1}))<1e-13 && "minkTensorTimes or minkTensorSum broken!" );
    
    
    int drei = 3;
    minkTensorIntegrand mytensor4(3,4,0,r,n);
    minkTensorIntegrand mytensor5(3,4,0,nr,nr);
    auto newtensor = drei*(3*(mytensor4 + mytensor5));
    assert( abs( newtensor.accessElement({1,1,0,1,1,0,0}) -(9*(mytensor4.accessElement({1,1,0,1,1,0,0}) + mytensor5.accessElement({1,1,0,1,1,0,0}))) ) < 1e-13 && "'auto newtensor = drei*(3*(mytensor4 + mytensor5));' does not work! " );
    
    
    tensor2D testTensor(2,1,0);
    testTensor.writeElement({0,0,1},5);
    assert(abs(testTensor.accessElement({0,0,1})-5)<1e-14 && "testTensor::writeElement broken!");
    
    zahlenfriedhof = mytensor + mytensor2;
    cout << "trace additivity check deactivated" << '\n';
    //assert( abs(trace(zahlenfriedhof) - (trace(mytensor) +trace(mytensor2))) < 1e-13 && "trace not additive!" );
    assert( abs(trace(zahlenfriedhof) - (zahlenfriedhof.accessElement({1,1})*pow(sin(zahlenfriedhof.r.theta),2) + zahlenfriedhof.accessElement({0,0}))) < 1e-13 && "trace not working! (check normalization?)" );
    
    
    //Test actual map, trace of W^(0,2)_1 should be equal to boundary length
    Healpix_Map<double> degradedMap(128, map.Scheme(), SET_NSIDE);
    degradedMap.Import_degrade(map);
    map = degradedMap;
    
    std::vector<double> thresholds = makeIntervals_lin(0, 1e-5, 3);
    std::vector<minkmapSphere> mapsScalar;
    std::vector<minkmapSphere> mapsTensor;
    for(double thresh : thresholds)
    {
        mapsScalar.push_back(minkmapSphere(map, 0, 0, 1, thresh));
        mapsTensor.push_back(minkmapSphere(map, 0, 2, 1, thresh));
    }
    int pixnum = 68; //1933
    minkmapStack sumOfMapsS(mapsScalar);
    minkmapStack sumOfMapsT(mapsTensor);
    
    auto pix1S = mapsScalar.at(0).at(pixnum);
    auto pix1T = mapsTensor.at(0).at(pixnum);
    assert( abs(pix1S-trace(pix1T))/max(double(pix1S), trace(pix1T) )<1e-6 && "Trace of (0,2,1) and boundary different at thresh 0" );
    
    auto pix2S = mapsScalar.at(1).at(pixnum);
    auto pix2T = mapsTensor.at(1).at(pixnum);
    assert( abs(pix2S-trace(pix2T))/max(double(pix2S), trace(pix2T) )<2e-5 && "Trace of (0,2,1) and boundary different at thresh > 0" );
    
    auto pix3S = sumOfMapsS.at(pixnum);
    auto pix3T = sumOfMapsT.at(pixnum);
    assert( abs(pix3S-trace(pix3T))/max(double(pix3S), trace(pix3T) )<2e-5 && "Trace of (0,2,1) and boundary different in sum of maps" );
    
    
    const normalHealpixInterface interfaceS(sumOfMapsS);
    const normalHealpixInterface interfaceT(sumOfMapsT);
    auto pix4S = interfaceS.at(pixnum);
    auto pix4T = interfaceT.at(pixnum);
    auto pix4Tfail = sumOfMapsT.at(67);
    
    cout << pix4S << " " << trace(pix4T) << " " << abs(pix4S-trace(pix4T))/max(double(pix4S), trace(pix4T) ) << endl;
    assert( abs(pix4S-trace(pix4T))/max(double(pix4S), trace(pix4T) )<1e-3 && "Trace of (0,2,1) and boundary different in Interface" );
    
    return 0;
}
