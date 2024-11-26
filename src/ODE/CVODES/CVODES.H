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

Class
    Foam::CVODES

Description


SourceFiles
    CVODES.C

\*---------------------------------------------------------------------------*/

#ifndef CVODES_H
#define CVODES_H


#include "extFunc.H"
#include "ODESolver.H"
#include "cantera/numerics/CVodesIntegrator.h"
#include "cantera/numerics/Integrator.h"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class CVODES Declaration
\*---------------------------------------------------------------------------*/

class CVODES
:
    public ODESolver
{
    // Private data

        // pointer to access cantera Integrators
        std::unique_ptr<Cantera::Integrator> integ_;
        // current solution
        mutable scalarField y_;        
        // Sensitivity Tolerances
        scalar senRelTol_ = 1.0e-4;
        scalar senAbsTol_ = 1.0e-6;
        int maxErrTestFails_ = 0;
        
        // current time
        mutable scalar time_ = 0.0;
                
        // if already initialized, true
        mutable bool init_ = false;
        
        // if already reinitialized, true
        mutable bool integrator_init_ = false;

public:

    //- Runtime type information
    TypeName("CVODES");


    // Constructors

        //- Construct from ODESystem
        CVODES(const ODESystem& ode, const dictionary& dict);


    //- Destructor
    virtual ~CVODES()
    {}


    // Member Functions    
        
        //- Resize the ODE solver
        virtual bool resize();

        //- Solve the ODE system from xStart to xEnd, update the state
        //  and return an estimate for the next step in dxTry
        virtual void solve
        (
            const scalar xStart,
            const scalar xEnd,
            scalarField& y,
            scalar& dxEst
        ) const;        
        
        

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
