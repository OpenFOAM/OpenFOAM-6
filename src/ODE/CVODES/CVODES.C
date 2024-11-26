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

#include "CVODES.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(CVODES, 0);

    addToRunTimeSelectionTable(ODESolver, CVODES, dictionary);

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CVODES::CVODES(const ODESystem& ode, const dictionary& dict)
:
    ODESolver(ode, dict),
    integ_(Cantera::newIntegrator("CVODE")),
    y_(n_)
{
    // use backward differencing, with a full Jacobian computed
    // numerically, and use a Newton linear iterator
    integ_->setMethod(Cantera::BDF_Method);
    integ_->setProblemType(Cantera::DENSE + Cantera::NOJAC);
    integ_->setIterator(Cantera::Newton_Iter);
    
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
bool Foam::CVODES::resize()
{
    if (ODESolver::resize())
    {
        resizeField(y_);
        init_ = false;

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::CVODES::solve
(
    const scalar xStart,
    const scalar xEnd,
    scalarField& y,
    scalar& dxEst
) const
{
    y_ = y;
    
    extFunc ext_(n_, y_, odes_);
    time_ = xStart;

    if (!init_)
    {        
        doublereal CvodeAbsTol[n_];
        for(label i=0; i<n_; i++)
        {
            CvodeAbsTol[i] = absTol_[i];
        }

        integ_->setTolerances(relTol_[0], n_, CvodeAbsTol);
        integ_->setSensitivityTolerances(senRelTol_, senAbsTol_);
        integ_->setMaxStepSize(xEnd - xStart);
        integ_->setMaxErrTestFails(maxErrTestFails_); 
        integ_->initialize(time_, ext_);
        integrator_init_ = true;
        init_ = true;
    }
    else if (!integrator_init_)
    {
        integ_->reinitialize(time_, ext_);
        integrator_init_ = true;
    }
    
    integ_->integrate(xEnd);
    time_ = xEnd;
    forAll (y_, i)
    {
        y[i] = integ_->solution()[i];
    }
    
    integrator_init_ = false;    
    dxEst = 1; //make it big
}


// ************************************************************************* //
