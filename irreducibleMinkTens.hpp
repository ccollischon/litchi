/*
 * This file is part of litchi, a lightweight C++ library
 * for Minkowski analysis
 * 
 * Copyright (C) 2021-2024 Caroline Collischon <caroline.collischon@fau.de>
 * 
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */


#ifndef irreducible
#define irreducible

#include "healpix_cxx/pointing.h"

#include <complex>
#include <cmath>

struct irreducibleMinkTens {
    const int l{};//, m{};
    pointing r{0,0};
    pointing n{1,0};
    
    irreducibleMinkTens(int lNew, pointing rNew, pointing nNew) : l{lNew}, r{std::move(rNew)}, n{std::move(nNew)}
    {
    }
    
    //~ std::complex<double> accessElement() const // Y*lm(n)
    //~ {
        //~ double thetapart = std::sph_legendre(l_,m_,n_.theta);
        //~ return thetapart*std::exp(std::complex<double>(0.0, -1.0) * (double)m_ * n_.phi);
    //~ }
    
    std::complex<double> accessElement() const // exp(i l phi(n))
    {
        double phi = giveAngle(n,r);
        return std::exp(std::complex<double>(0.0, +1.0) * (double)l * phi);
    }
    
    
};



#endif
