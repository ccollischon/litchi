#ifndef geometryhelpers
#define geometryhelpers


#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <cmath>

#include "healpix_cxx/pointing.h"
#include "healpix_cxx/rotmatrix.h"

/** \file geometryHelpers.hpp
 * \brief Helpful functions for geometry on the sphere
 */

/** 
 * Rotates input around axis by angle
 * \param axis Rotation axis. Gets normalized before rotation
 * \param input Vector to be rotated
 * \param angle Rotation angle. Since axis ist normalized, input will be rotated by exactly this angle
 * \return The rotated vector
 */
vec3 rotateAroundAxis(vec3 axis, const vec3& input, double angle)
{
    axis.Normalize();
    rotmatrix theRotation;
    theRotation.Make_Axis_Rotation_Transform(axis, angle);
    return theRotation.Transform(input);
}

/**
 * Angular distance between A and B
 */
double arclength(pointing A, pointing B) //Distance between A and B on Sphere
{
    if( (std::abs(A.theta-B.theta)<1e-14) && (std::abs(A.phi-B.phi)<1e-14) )  return 0.;
    const double pi = 3.14159265359;
    double sigma = acos( sin(pi/2-A.theta)*sin(pi/2-B.theta) + cos(pi/2-A.theta)*cos(pi/2-B.theta)*cos(A.phi-B.phi));
    return sigma;
}

/**
 * Midpoint between A and B. Only calculated if distance larger than 1e-12
 */
pointing midpoint(pointing A, pointing B) //rotate starting from A by half arclength
{
    double sigma = arclength(A, B);
    if(sigma < 1e-12) return A;
    
    vec3 rotAxis = crossprod(A.to_vec3(),B.to_vec3());
    vec3 outVec = rotateAroundAxis(rotAxis, A, sigma/2);
    return pointing(outVec);
}

/**
 * Mantz et al 2008 interpolation between 2 lattice points. If the field value at one point is close to the threshold the interpolated point will be close to it and vice versa
 * \param A Position on sphere of one point
 * \param valA Field value at A
 * \param B Position on sphere of other point
 * \param valB Field value at B
 * \param thresh Threshold applied to the field for which the interpolated position is wanted
 * \return A point bewteen A and B whose exact position depends on the field values at those points
 */
pointing interpPointing(pointing A, double valA, pointing B, double valB, double thresh) // Do the Mantz et al 2008 interpolation between 2 lattice points
{
    const double pi = 3.1415926535897962643;
    double sigma = arclength(A,B);//Angle between A and B
    double quotient = (thresh - valA)/(valB - valA); //Mantz et al 2008 style interpolation
    //double fun = sin(quotient*(pi/2))*sin(quotient*(pi/2)); //nonlinearity via sin^2(x), scaled to fit into [0,1]
    //double fun = asin(sqrt(quotient))/(pi/2);
    //double fun = 0.5;
    double angle = quotient*sigma;
    
    vec3 rotAxis = crossprod(A.to_vec3(),B.to_vec3());
    vec3 Avec = A.to_vec3();
    vec3 rotatedVec = rotateAroundAxis(rotAxis, Avec, angle);
    pointing retpoint(rotatedVec);
    retpoint.normalize();
    if(retpoint.phi-2*pi > 0) {
        std::cout << retpoint.phi << '\n';
        retpoint.phi -= 2*pi;
    }
    return retpoint;
}


/**
 * Normalizes a vector on the sphere in tangent space at colatitude theta
 * \param input Tangent space vector to be normalized (the two coordinate components of the pointing class are used as vector components)
 * \param theta Colatitude of the vector's position
 */
void normalizeVectorOnSphere(pointing& input, double theta)
{
    double lengthsquared = pow(input.theta,2) + pow(input.phi*sin(theta),2); //EDIT
    //double lengthsquared = pow(input.theta,2) + pow(input.phi,2);
    input.theta /= sqrt(lengthsquared);
    input.phi /= sqrt(lengthsquared);
}

/**
 * Calculate vector normal to contour given by two close points. Vector connecting the two positions is calculated in cartesian approximation, so their distance should not be too large
 * \param A One position
 * \param B Other position
 * \return Tangent space vector at position A normal to line from A to B
 */
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

/**
 * Calculate vector normal to contour given by two close points. Vector connecting the two positions is calculated by rotating the first point
 * \param A One position
 * \param B Other position
 * \return Tangent space vector at position A normal to line from A to B
 */
pointing getN_rotation(pointing A, pointing B)//A, B should be given such that AxB points away from body
{
    if(A.phi<0 || A.theta>3.14159) A.normalize(); //A with negative phi would cause problems when constructing n
    
    vec3 AcrossB = crossprod(A.to_vec3(),B.to_vec3());
    vec3 rotAxis = crossprod(A.to_vec3(), AcrossB); //rotating A around this axis moves it along the desired normal vector
    
    pointing Arot(rotateAroundAxis(rotAxis,A.to_vec3(),0.001));
    if(Arot.phi-A.phi > 5){
        Arot.phi = Arot.phi - 2*3.14159;
    }
    pointing n(Arot.theta-A.theta, Arot.phi-A.phi);
    
    normalizeVectorOnSphere(n,A.theta);
    return n;
}

/**
 * Area on sphere in triangle defined by the three given corners
 */
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

/**
 * Make pointing point in opposite direction on sphere
 */
void flipPointing(pointing& input)
{
    const double pi = 3.14159265359;
    input.theta = pi-input.theta;
    input.phi = pi+input.phi;
    input.normalize();
}

/**
 * Parallel transport tangent vector on Sphere along geodesic from start to stop. Uses Schild's ladder procedure.
 * \param start Starting position
 * \param stop Final position
 * \param initialVector Tangent space vector to be transported (the two coordinate components of the pointing class are used as vector components)
 * \return Transported vector at position stop
 */
pointing parallelTransport(pointing start, pointing stop, pointing initialVector) //Transport along geodesic
{
    start.normalize();
    stop.normalize();/*
    //std::cout << depth << "\n" ;
    for(int i=depth; (arclength(start,stop)/sin(start.theta)>1e-2) && i<10; ++i) //If two points are further apart use more steps
    {
        pointing newstart = midpoint(start,stop);
        initialVector = parallelTransport(start,newstart,initialVector,i+1);
        start = newstart;
    }*/
    
    //Schild's ladder
    pointing initShort(initialVector.theta*0.001, initialVector.phi*0.001);
    pointing fromStartAlongInit(start.theta + initShort.theta, start.phi + initShort.phi); //move small distance along init
    pointing midpointToAim = midpoint(fromStartAlongInit, stop); //find midpoint between pointing above and endpoint of transport
    
    //Now find point at twice distance from start to midpoint along geodesic, the geodesic through stop along transported vector passes through here
    double dist = 2*arclength(start,midpointToAim);
    auto rotaxis = crossprod(start.to_vec3(),midpointToAim.to_vec3());
    pointing fromStopAlongFinal(rotateAroundAxis(rotaxis, start.to_vec3(), dist));
    
    pointing finalVector(fromStopAlongFinal.theta-stop.theta, fromStopAlongFinal.phi-stop.phi); //Take difference to return to tangent space
    normalizeVectorOnSphere(finalVector,stop.theta);
    
    
    return finalVector;
}



#endif
