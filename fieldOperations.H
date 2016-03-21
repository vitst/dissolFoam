/* 
 * File:   fieldOperations.H
 * Author: vstar
 *
 * Created on October 16, 2015, 1:05 PM
 */

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    Auxiliary class for the dissolFoam solver.

Description
    All the operations with fields i.e. rescaling, boundary conditions.
  
Contributors/Copyright:
    Vitaliy Starchenko (2015)

 
\*---------------------------------------------------------------------------*/


#ifndef FIELDOPERATIONS_H
#define	FIELDOPERATIONS_H


#include "fvCFD.H"
  
class fieldOperations{

private:
  const argList args_;
  const label inletID_;
  
  // aux. parameter in order to calculate inlet flow rate at t0 only once
  bool inletFlowRateT0calculated;
  // aux. parameter where we store inlet flow rate ar t0
  scalar inletFlowRateT0;
  
  // functions
  scalar getInletFlowRateT0(const surfaceScalarField&);
  
public:
  
  // Constractors
  fieldOperations(const argList&, const label);
  
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // functions
  scalar getInletFlowRate(const surfaceScalarField&, bool);

  vectorField getWallPointMotion(const fvMesh&, const volScalarField&, const scalar, const label);
  
  scalar getInletAreaT0();
};


#endif	/* FIELDOPERATIONS_H */

