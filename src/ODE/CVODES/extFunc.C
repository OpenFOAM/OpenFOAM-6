/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "extFunc.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extFunc::extFunc(label n, scalarField& y, const ODESystem& ode)
:
    FuncEval(),
    n_(n),
    y_(y),
    odes_(ode)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::extFunc::eval
(
    doublereal t, 
    doublereal* y,
    doublereal* ydot, 
    doublereal* p
)
{
    scalarField y0(n_);
    scalarField dy(n_);
    forAll(y0, i)
    {
        y0[i] = y[i];
    }
    
    odes_.derivatives(t, y0, dy);
    
    forAll(dy, i)
    {
        ydot[i] = dy[i];
    }
}


void Foam::extFunc::getState
(
    doublereal* y
)
{
    forAll(y_, i)
    {        
        y[i] = y_[i];
    }
}

// ************************************************************************* //
