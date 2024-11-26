/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::reactingMixture

Description
    Foam::reactingMixture

SourceFiles
    reactingMixture.C

\*---------------------------------------------------------------------------*/

#ifndef reactingMixture_H
#define reactingMixture_H

#include "speciesTable.H"
#include "chemistryReader.H"
#include "multiComponentMixture.H"
#include "fvc.H"
#include "cantera/transport.h"
#include "cantera/transport/MultiTransport.h"
#include "cantera/transport/MixTransport.h"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class reactingMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class reactingMixture
:
    public speciesTable,
    public autoPtr<chemistryReader<ThermoType>>,
    public multiComponentMixture<ThermoType>,
    public PtrList<Reaction<ThermoType>>
{
    // Private member data
    
        //- cantera object
        mutable Cantera::IdealGasMix gas_;
        Cantera::Transport* tran_;
        Cantera::MixTransport* tranMix_;
        Cantera::MultiTransport* tranMulti_;        

        //- number of species
        int nsp_;
        
        //- if mixtureAveraged transport model used
        bool mixtureAveraged_ = false;

        //- if multiComponent transport model used        
        bool multiComponent_ = false;
        
        //- if mixtureAveraged transport model used
        bool constantLewis_ = false;

        //- if multiComponent transport model used        
        bool unityLewis_ = false;        
        
        //- Table of species composition
        speciesCompositionTable speciesComposition_;
        


    // Private Member Functions

        //- Disallow default bitwise copy construct
        reactingMixture(const reactingMixture&);

        //- Disallow default bitwise assignment
        void operator=(const reactingMixture&);


public:

    //- The type of thermo package this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary, mesh and phase name
        reactingMixture(const dictionary&, const fvMesh&, const word&);


    //- Destructor
    virtual ~reactingMixture()
    {}


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "reactingMixture<" + ThermoType::typeName() + '>';
        }

        //- Read dictionary
        void read(const dictionary&);

        label size() const
        {
            return PtrList<Reaction<ThermoType>>::size();
        }

        Reaction<ThermoType>& operator[](const label i)
        {
            return PtrList<Reaction<ThermoType>>::operator[](i);
        }

        const Reaction<ThermoType>& operator[](const label i) const
        {
            return PtrList<Reaction<ThermoType>>::operator[](i);
        }

        //- Table of species composition
        const speciesCompositionTable& specieComposition() const
        {
            return speciesComposition_;
        }
        
        //- species mass diffusivity
        inline PtrList<volScalarField> D
        (
            const volScalarField& p, 
            const volScalarField& T, 
            const PtrList<volScalarField>& Y,
            const volScalarField& rho        
        ) const;
        
        inline PtrList<volScalarField> Dmass
        (
            const volScalarField& p, 
            const volScalarField& T, 
            const PtrList<volScalarField>& Y
        ) const;        
        
        //- species thermal diffusivity
        inline PtrList<volScalarField> DT
        (
            const volScalarField& p, 
            const volScalarField& T, 
            const PtrList<volScalarField>& Y
        ) const;
        
        inline PtrList<volVectorField> Dm
        (
            const volScalarField& p, 
            const volScalarField& T, 
            const PtrList<volScalarField>& Y,
            const volScalarField& rho    
        ) const;
        
/*        //- mixture viscosity
        inline scalar mixVisc
        (
            const scalar p,
            const scalar T,
            const scalarField& Y
        ) const;

        //- mixture thermal conductivity
        inline scalar mixThermalCond
        (
            const scalar p,
            const scalar T,
            const scalarField& Y
        ) const;
*/

        //- binary diffusivity coefficients
        inline PtrList<volScalarField> bD
        (
            const volScalarField& p, 
            const volScalarField& T, 
            const PtrList<volScalarField>& Y
        ) const;    

        //- species fluxes
        inline PtrList<volVectorField> speFlux
        (
            const volScalarField& p,
            const volScalarField& T,
            const PtrList<volScalarField>& c
        ) const;            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "reactingMixtureI.H"

#ifdef NoRepository
    #include "reactingMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
