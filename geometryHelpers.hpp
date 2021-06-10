#ifndef geometryhelpers
#define geometryhelpers


#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <cmath>

#include "healpix_cxx/pointing.h"
#include "healpix_cxx/rotmatrix.h"


vec3 rotateAroundAxis(vec3 axis, const vec3& input, double angle)
{
    axis.Normalize();
    rotmatrix theRotation;
    theRotation.Make_Axis_Rotation_Transform(axis, angle);
    return theRotation.Transform(input);
}


double arclength(pointing A, pointing B) //Distance between A and B on Sphere
{
    if( (std::abs(A.theta-B.theta)<1e-14) && (std::abs(A.phi-B.phi)<1e-14) )  return 0.;
    const double pi = 3.14159265359;
    double sigma = acos( sin(pi/2-A.theta)*sin(pi/2-B.theta) + cos(pi/2-A.theta)*cos(pi/2-B.theta)*cos(A.phi-B.phi));
    return sigma;
}

pointing midpoint(pointing A, pointing B) //rotate starting from A by half arclength
{
    double sigma = arclength(A, B);
    if(sigma < 1e-12) return A;
    
    vec3 rotAxis = crossprod(A.to_vec3(),B.to_vec3());
    vec3 outVec = rotateAroundAxis(rotAxis, A, sigma/2);
    return pointing(outVec);
}

pointing interpPointing(pointing A, double valA, pointing B, double valB, double thresh) // Do the Mantz et al 2008 interpolation between 2 lattice points
{
    const double pi = 3.14159265359;
    double sigma = acos( sin(pi/2-A.theta)*sin(pi/2-B.theta) + cos(pi/2-A.theta)*cos(pi/2-B.theta)*cos(A.phi-B.phi));//Angle between A and B
    double quotient = (thresh - valA)/(valB - valA); //Mantz et al 2008 style interpolation
    double angle = quotient*sigma;
    
    vec3 rotAxis = crossprod(A.to_vec3(),B.to_vec3());
    vec3 Avec = A.to_vec3();
    vec3 rotatedVec = rotateAroundAxis(rotAxis, Avec, angle);
    pointing retpoint(rotatedVec);
    retpoint.normalize();
    return retpoint;
}



void normalizeVectorOnSphere(pointing& input, double theta)
{
    double lengthsquared = pow(input.theta,2) + pow(input.phi*sin(theta),2); //EDIT
    //double lengthsquared = pow(input.theta,2) + pow(input.phi,2);
    input.theta /= sqrt(lengthsquared);
    input.phi /= sqrt(lengthsquared);
}

//Calculate vecor normal to contour given by two points
pointing getN_cartesian(pointing A, pointing B)//A, B should be given such that AxB points away from body
{
    const double pi = 3.14159265359;
    A.normalize();
    B.normalize();
    if(std::abs(A.theta)<1e-15 || (std::abs(A.theta-pi))<1e-15)
    {
        std::cerr << "Error: Trying to get normal vector very close to pole: (theta,phi) = " << A << " , this is undefined" << std::endl;
        throw std::invalid_argument( "getN_cartesian: n close to pole" );
    }
    double delTheta = A.theta-B.theta; //Direction along curve in cartesian approximation
    double delPhi  = A.phi-B.phi;
    
    double n1 = delPhi; //elements of normalized vector normal to (deltheta,delPhi)
    double n2 = -delTheta/(sin(A.theta)*sin(A.theta));
    
    pointing n(n1,n2);
    normalizeVectorOnSphere(n,A.theta);
    return n;
}


double sphereArea(pointing A, pointing B, pointing C) //Area of spherical triangle = sum of angles minus pi
{
    const double pi = 3.14159265359;
    vec3 AcrossB = crossprod(A.to_vec3(),B.to_vec3());
    vec3 AcrossC = crossprod(A.to_vec3(),C.to_vec3());
    AcrossB.Normalize();
    AcrossC.Normalize();
    double angA = acos(dotprod(AcrossB,AcrossC));
    
    AcrossB.Flip();
    vec3 BcrossA = AcrossB;
    vec3 BcrossC = crossprod(B.to_vec3(),C.to_vec3());
    BcrossC.Normalize();
    double angB = acos(dotprod(BcrossA,BcrossC));
    
    AcrossC.Flip();
    vec3 CcrossA = AcrossC;
    BcrossC.Flip();
    vec3 CcrossB = BcrossC;
    double angC = acos(dotprod(CcrossA,CcrossB));
    double area = angA+angB+angC - pi;
    return area;
}

void flipPointing(pointing& input) //make pointing point in opposite direction
{
    const double pi = 3.14159265359;
    input.theta = pi-input.theta;
    input.phi = pi+input.phi;
    input.normalize();
}

pointing parallelTransport(pointing start, pointing stop, pointing initialVector) //Transport along geodesic
{
    start.normalize();
    stop.normalize();
    
    //Schild's ladder
    pointing initShort(initialVector.theta*0.001, initialVector.phi*0.001);
    pointing fromStartAlongInit(start.theta + initShort.theta, start.phi + initShort.phi); //move small distance along init
    pointing midpointToAim = midpoint(fromStartAlongInit, stop); //find midpoint between pointing above and endpoint of transport
    
    //Now find point at twice distance from start to midpoint along geodesic, the geodesic through stop along transported vector passes through here
    double dist = 2*arclength(start,midpointToAim);
    auto rotaxis = crossprod(start.to_vec3(),midpointToAim.to_vec3());
    pointing fromStopAlongFinal(rotateAroundAxis(rotaxis, start.to_vec3(), dist));
    
    pointing finalVector(fromStopAlongFinal.theta-stop.theta, fromStopAlongFinal.phi-stop.phi); //Take difference to return to tengent space
    normalizeVectorOnSphere(finalVector,stop.theta);
    
    
    return finalVector;
}



#endif
