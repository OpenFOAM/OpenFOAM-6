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
    Foam::extFunc

Description


SourceFiles
    extFunc.C

\*---------------------------------------------------------------------------*/

#ifndef extFunc_H
#define extFunc_H

#include "cantera/numerics/FuncEval.h"
#include "ODESystem.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class CVODES Declaration
\*---------------------------------------------------------------------------*/

class extFunc
:
    public Cantera::FuncEval
{
    // Private data

        label n_;
        scalarField& y_;        
        const ODESystem& odes_;
public:

    // Constructors

        //- Construct from ODESystem
        extFunc(label n, scalarField& y, const ODESystem& ode);


    //- Destructor
    virtual ~extFunc()
    {}


    // Member Functions    
        virtual size_t neq()
        {
            return n_;
        }    
        virtual void getState(double* y);
        virtual void eval
        (
            doublereal t, 
            doublereal* y,
            doublereal* ydot, 
            doublereal* p
        );    

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
