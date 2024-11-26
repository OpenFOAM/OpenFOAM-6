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
    Foam::multiComponentMixture

Description
    Foam::multiComponentMixture

SourceFiles
    multiComponentMixture.C

\*---------------------------------------------------------------------------*/

#ifndef multiComponentMixture_H
#define multiComponentMixture_H

#include "basicSpecieMixture.H"
#include "HashPtrTable.H"
#include "cantera/IdealGasMix.h"
#include "cantera/transport.h"
#include "cantera/transport/MultiTransport.h"
#include "cantera/transport/MixTransport.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class multiComponentMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class multiComponentMixture
:
    public basicSpecieMixture
{
    // Private data
    
        //- cantera object
        Cantera::IdealGasMix gas_;
        Cantera::Transport* tran_;
        Cantera::MixTransport* tranMix_;
        Cantera::MultiTransport* tranMulti_;
        
        //- Species data
        PtrList<ThermoType> speciesData_;

        //- Temporary storage for the cell/face mixture thermo data
        mutable ThermoType mixture_;

        //- Temporary storage for the volume weighted
        //  cell/face mixture thermo data
        mutable ThermoType mixtureVol_;


    // Private Member Functions

        //- Construct the species data from the given dictionary and return the
        //  data for the first specie to initialise the mixture thermo data
        const ThermoType& constructSpeciesData(const dictionary& thermoDict);

        //- Correct the mass fractions to sum to 1
        void correctMassFractions();

        //- Construct as copy (not implemented)
        multiComponentMixture(const multiComponentMixture<ThermoType>&);      
    
protected:
        //- if mixtureAveraged transport model used
        bool MmixtureAveraged_ = false;

        //- if multiComponent transport model used        
        bool MmultiComponent_ = false;

        //- if mixtureAveraged transport model used
        bool MconstantLewis_ = false;

        //- if multiComponent transport model used        
        bool MunityLewis_ = false;
        
        //- Lewis number
        scalarField LewisNumber_;
public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary, specie names, thermo database,
        //  mesh and phase name
        multiComponentMixture
        (
            const dictionary&,
            const wordList& specieNames,
            const HashPtrTable<ThermoType>& thermoData,
            const fvMesh&,
            const word&,
            Cantera::IdealGasMix&
        );

        //- Construct from dictionary, mesh and phase name
        multiComponentMixture(const dictionary&, const fvMesh&, const word&);


    //- Destructor
    virtual ~multiComponentMixture()
    {}


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "multiComponentMixture<" + ThermoType::typeName() + '>';
        }

        const ThermoType& cellMixture(const label celli) const;

        const ThermoType& patchFaceMixture
        (
            const label patchi,
            const label facei
        ) const;

        const ThermoType& cellVolMixture
        (
            const scalar p,
            const scalar T,
            const label celli
        ) const;

        const ThermoType& patchFaceVolMixture
        (
            const scalar p,
            const scalar T,
            const label patchi,
            const label facei
        ) const;
        
        scalar muCellMixture
        (
            const scalar p,
            const scalar T,
            const label celli
        );

        scalar alphahCellMixture
        (
            const scalar p,
            const scalar T,
            const label celli
        );

        scalar muPatchFaceMixture
        (
            const scalar p,
            const scalar T,
            const label patchi,
            const label facei            
        );

        scalar alphahPatchFaceMixture
        (
            const scalar p,
            const scalar T,
            const label patchi,
            const label facei
        );        

        //- Return the raw specie thermodynamic data
        const PtrList<ThermoType>& speciesData() const
        {
            return speciesData_;
        }
        
        Cantera::IdealGasMix canteraGas()
        {
            return gas_;
        }

        //- Read dictionary
        void read(const dictionary&);

        //- Return thermo based on index
        inline const ThermoType& getLocalThermo(const label speciei) const
        {
            return speciesData_[speciei];
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "multiComponentMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
