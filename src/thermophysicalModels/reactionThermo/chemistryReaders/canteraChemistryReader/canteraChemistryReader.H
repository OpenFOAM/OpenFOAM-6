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
    Foam::canteraChemistryReader

Description
    Cantera chemistry reader for Cantera format

SourceFiles
    canteraChemistryReader.C

\*---------------------------------------------------------------------------*/

#ifndef canteraChemistryReader_H
#define canteraChemistryReader_H

#include "chemistryReader.H"
#include "fileName.H"
#include "typeInfo.H"
#include "Switch.H"
#include "HashPtrTable.H"
#include "ReactionList.H"
#include "DynamicList.H"
#include "labelList.H"
#include "speciesTable.H"
#include "atomicWeights.H"
#include "reactionTypes.H"
#include "cantera/IdealGasMix.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class canteraChemistryReader Declaration
\*---------------------------------------------------------------------------*/

class canteraChemistryReader
:
    public chemistryReader<gasHThermoPhysics>
{
        string ctiFileString;
        
            string ctiFileID;
        
          Cantera::IdealGasMix gas_;
        
        //- List of elements
        DynamicList<word> elementNames_;

        //- Element indices
        HashTable<label> elementIndices_;

        //- Table of species
        speciesTable& speciesTable_;

        //- Table of species composition
        speciesCompositionTable speciesComposition_;     

        //- Table of the thermodynamic data given in the foamChemistry file
        HashPtrTable<gasHThermoPhysics> speciesThermo_;

        //- List of the reactions
        ReactionList<gasHThermoPhysics> reactions_;
    
    // Private Member Functions
    
        //- Set the species list
        speciesTable& setSpecies
        (
            Cantera::IdealGasMix gas, 
            speciesTable& species
        );
        
        // Set species thermo data
        void setSpeciesThermo
        (
            Cantera::IdealGasMix gas, 
            const dictionary& transportDict
        );
        
        // Set species thermo data
        void setSpeciesComposition(Cantera::IdealGasMix gas);

        // Set reaction data
        void setReactions(Cantera::IdealGasMix gas);
    
        //initialize cantera
        void canteraRead(Cantera::IdealGasMix gas);
        
        enum reactionType
        {
            irreversible,
            reversible,
            unknownReactionType
        };
        
        enum reactionRateType
        {
            Arrhenius,
            thirdBodyArrhenius,
            unimolecularFallOff,
            chemicallyActivatedBimolecular,
            LandauTeller,
            Janev,
            powerSeries,
            unknownReactionRateType
        };
        
        template<class ReactionRateType>
        void addReactionType
        (
            const reactionType rType,
            DynamicList<gasHReaction::specieCoeffs>& lhs,
            DynamicList<gasHReaction::specieCoeffs>& rhs,
            const ReactionRateType& rr
        );

        template<template<class, class> class PressureDependencyType>
        void addPressureReactionType
        (
            const reactionType rType,
            DynamicList<gasHReaction::specieCoeffs>& lhs,
            DynamicList<gasHReaction::specieCoeffs>& rhs,
            const scalarList& efficiencies,
            const scalarList& k0Coeffs,
            const scalarList& kInfCoeffs,
            const scalarList& fallOffCoeffs,
            int fofType
        );        
        
public:

    //- Runtime type information
    TypeName("canteraChemistryReader");

    // Constructors

    canteraChemistryReader
    (
        const dictionary& thermoDict,
        speciesTable& species
    );

    //- Destructor
    virtual ~canteraChemistryReader()
    {}
    
    // member functions
            Cantera::IdealGasMix& canteraGas()
            {
                return gas_;
            }
    
        //- List of elements from cantera
        const wordList& elementNames() const
        {
            return elementNames_;
        }

        //- Element indices cantera
        const HashTable<label>& elementIndices() const
        {
            return elementIndices_;
        }
        
        // Return the species table from cantera
        const speciesTable& species() const
        {
            return speciesTable_;
        } 
        
        //- Table of species composition from cantera
        const speciesCompositionTable& specieComposition() const
        {
            return speciesComposition_;
        }

        //- Table of the thermodynamic data from cantera
        const HashPtrTable<gasHThermoPhysics>& speciesThermo() const
        {
            return speciesThermo_;
        }

        //- List of the reactions from cantera
        const ReactionList<gasHThermoPhysics>& reactions() const
        {
            return reactions_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


#endif

// ************************************************************************* //
