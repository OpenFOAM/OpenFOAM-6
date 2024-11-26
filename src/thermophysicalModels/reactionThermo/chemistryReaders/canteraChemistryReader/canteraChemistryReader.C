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

\*---------------------------------------------------------------------------*/

#include "canteraChemistryReader.H"
#include "IFstream.H"
#include "atomicWeights.H"
#include "ReactionProxy.H"
#include "IrreversibleReaction.H"
#include "ReversibleReaction.H"
#include "NonEquilibriumReversibleReaction.H"
#include "ArrheniusReactionRate.H"
#include "thirdBodyArrheniusReactionRate.H"
#include "FallOffReactionRate.H"
#include "ChemicallyActivatedReactionRate.H"
#include "LindemannFallOffFunction.H"
#include "TroeFallOffFunction.H"
#include "SRIFallOffFunction.H"
#include "LandauTellerReactionRate.H"
#include "JanevReactionRate.H"
#include "powerSeriesReactionRate.H"
#include "addToRunTimeSelectionTable.H"
#include "cantera/IdealGasMix.h" 
#include "cantera/kinetics/GasKinetics.h"
#include "cantera/kinetics/reaction_defs.h"

namespace Foam
{
    addChemistryReaderType(canteraChemistryReader, gasHThermoPhysics);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Obtain species table from cantera
Foam::speciesTable& Foam::canteraChemistryReader::setSpecies
(
    Cantera::IdealGasMix gas,
    speciesTable& species   
)
{
    label nsp = gas.nSpecies();
    wordList s(nsp);    
    for(label i=0; i<nsp; i++)
    {
        s[i] = gas.speciesName(i);
    }
    species.transfer(s);
//    Info<< species << endl;
    return species;    
}

// Obtain species thermo data from cantera
void Foam::canteraChemistryReader::setSpeciesThermo
(
    Cantera::IdealGasMix gas,
    const dictionary& transportDict    
)
{
    word currentSpecieName;
    scalar molecularWeight;
    int type;
    doublereal c[15];
    doublereal minTemp, maxTemp, refPressure;
    scalar currentLowT, currentHighT, currentCommonT;
    gasHThermoPhysics::coeffArray highCpCoeffs(scalarList(7));
    gasHThermoPhysics::coeffArray lowCpCoeffs(scalarList(7));
    label nsp = gas.nSpecies();
    
    // get a reference to the species thermo property manager
    Cantera::MultiSpeciesThermo& sp = gas.speciesThermo();
    
    for(label i=0; i<nsp; i++)
    {
        currentSpecieName = gas.speciesName(i);
        molecularWeight = gas.molecularWeight(i);
        
        // get the NASA coefficients in array c[1-14];
        // mid point temperature is c[0];        
        sp.reportParams(i, type, c, minTemp, maxTemp, refPressure);
        currentLowT = minTemp;
        currentHighT = maxTemp;
        currentCommonT = c[0];
        for (label j=0; j<7; j++)
        {
            highCpCoeffs[j] = c[j+1];
            lowCpCoeffs[j] = c[j+8];            
        }
        
        /* Note that the transportDict is provided in the same way as CHEMKIN reader,
           but the data is not used in calculating transport data. The viscosity, 
           diffusion coefficients et al. are computed in the cantera transport class */
        speciesThermo_.insert
        (
            currentSpecieName,
            new gasHThermoPhysics
            (
                janafThermo<perfectGas<specie>>
                (
                    specie
                    (
                        currentSpecieName,
                        1.0,
                        molecularWeight
                    ),
                    currentLowT,
                    currentHighT,
                    currentCommonT,
                    highCpCoeffs,
                    lowCpCoeffs,
                    true
                ),
                transportDict.subDict(currentSpecieName)
            )
        );
    }
//    Info << speciesThermo_ << endl;
    
}

// Obtain species composition from cantera
void Foam::canteraChemistryReader::setSpeciesComposition
(
    Cantera::IdealGasMix gas  
)
{
    label ne = gas.nElements();
    wordList e(ne);    
    
    DynamicList<word> elementNames_;
    HashTable<label> elementIndices_;
    label currentElementIndex(0);
    
    for(label i=0; i<ne; i++)
    {
        e[i] = gas.elementName(i);
    }
    
    forAll(e, ei)
    {
        elementIndices_.insert(e[ei], currentElementIndex++);
        elementNames_.append(e[ei]);
    }

    // Loop through all species to retrieve
    // the species composition from cantera
    forAll(speciesTable_, si)
    {
        label elementArraySize(0);
        forAll(e, ei)
        {
            if (gas.nAtoms(si, ei) > 0)
            {
                elementArraySize += 1;
            }

        }

        List<specieElement> currentComposition(elementArraySize);
        label elementNotZero(0);
        for(label eni=0; eni<ne; eni++)
        {            
            if (gas.nAtoms(si, eni) > 0)
                {                    
                    currentComposition[elementNotZero].name() = 
                    e[eni];
                    currentComposition[elementNotZero].nAtoms() =
                    gas.nAtoms(si, eni);
                    elementNotZero += 1;
                }
        }

        // Add current specie composition to the hash table
        speciesCompositionTable::iterator specieCompositionIter
        (
            speciesComposition_.find(speciesTable_[si])
        );

        if (specieCompositionIter != speciesComposition_.end())
        {
            speciesComposition_.erase(specieCompositionIter);
        }

        speciesComposition_.insert(speciesTable_[si], currentComposition);        
    }    
//    Info<< speciesComposition_ << endl;    
}

// obtain reactions from cantera
void Foam::canteraChemistryReader::setReactions
(
    Cantera::IdealGasMix gas  
)
{
    DynamicList<gasHReaction::specieCoeffs> lhs;
    DynamicList<gasHReaction::specieCoeffs> rhs;
    gasHReaction::specieCoeffs currentSpecieCoeff;
    
    scalarList efficiencies_;        
    label nsp = gas.nSpecies();
    efficiencies_.setSize(nsp);
    efficiencies_ = 1.0;
        
    reactionType rType = unknownReactionType;

    scalarList ArrheniusCoeffs(3);
    scalarList k0Coeffs(3);
    scalarList kInfCoeffs(3); 
    scalarList TroeCoeffs(4);
    scalarList SRICoeffs(5);
    double Troeparams[4];
    double SRIparams[5];

    Cantera::ElementaryReaction R1;
    Cantera::ThreeBodyReaction R2;
    Cantera::FalloffReaction R3;
    Cantera::Troe TroeR;
    Cantera::SRI SRIR;
    int falloffType;     
    
    // in cantera, the non-equilibrium reaction (REV / / in chemkin files)
    // is counted and considered as a new reaction.
    label nr = gas.nReactions();
    for (label i=0; i<nr; i++)
    {
        if (gas.isReversible(i))
    {
        rType = reversible;
    }
    else
    {
        rType = irreversible;
    }
        
    // R is a shared_prt;
    auto R = gas.reaction(i);
    switch (R->reaction_type)
    {
        case Cantera::ELEMENTARY_RXN:        

            R1 = dynamic_cast<Cantera::ElementaryReaction&>(*R);
            
        ArrheniusCoeffs[0] = R1.rate.preExponentialFactor(); // kmol, s, m
        ArrheniusCoeffs[1] = R1.rate.temperatureExponent(); 
        ArrheniusCoeffs[2] = R1.rate.activationEnergy_R(); // K 
            
            for (const auto& sp : R1.reactants) 
        {
            label isp = gas.speciesIndex(sp.first);
            currentSpecieCoeff.index = isp;
            currentSpecieCoeff.stoichCoeff = sp.second;
            currentSpecieCoeff.exponent = sp.second;
            lhs.append(currentSpecieCoeff);
        }
            
        for (const auto& sp : R1.products) 
        {
        label isp = gas.speciesIndex(sp.first);
        currentSpecieCoeff.index = isp;
        currentSpecieCoeff.stoichCoeff = sp.second;
        currentSpecieCoeff.exponent = sp.second;
        rhs.append(currentSpecieCoeff);
        }
        
            addReactionType
        (
                rType,
                lhs, rhs,
                ArrheniusReactionRate
                (
                    ArrheniusCoeffs[0],
                    ArrheniusCoeffs[1],
                    ArrheniusCoeffs[2]
                )
            );
        break;
            
        case Cantera::THREE_BODY_RXN:        
        
        R2 = dynamic_cast<Cantera::ThreeBodyReaction&>(*R);
            
        ArrheniusCoeffs[0] = R2.rate.preExponentialFactor(); // kmol, s, m
            ArrheniusCoeffs[1] = R2.rate.temperatureExponent(); 
            ArrheniusCoeffs[2] = R2.rate.activationEnergy_R(); // K 

            for (const auto& sp : R2.third_body.efficiencies) 
            {
                label isp = gas.speciesIndex(sp.first);
                efficiencies_[isp] = sp.second;
            }
            
            for (const auto& sp : R2.reactants) 
            { 
        label isp = gas.speciesIndex(sp.first);
        currentSpecieCoeff.index = isp;
        currentSpecieCoeff.stoichCoeff = sp.second;
        currentSpecieCoeff.exponent = sp.second;
        lhs.append(currentSpecieCoeff);
        }
            
        for (const auto& sp : R2.products) 
        {
        label isp = gas.speciesIndex(sp.first);
        currentSpecieCoeff.index = isp;
        currentSpecieCoeff.stoichCoeff = sp.second;
        currentSpecieCoeff.exponent = sp.second;
        rhs.append(currentSpecieCoeff);
        }
        
            addReactionType
            (
                rType,
                lhs, rhs,
                thirdBodyArrheniusReactionRate
                (
                    ArrheniusCoeffs[0],
                    ArrheniusCoeffs[1],
                    ArrheniusCoeffs[2],
                    thirdBodyEfficiencies(speciesTable_, efficiencies_)
                )
            );
            break;            
            
        case Cantera::FALLOFF_RXN: 
            case Cantera::CHEMACT_RXN:        
            R3 = dynamic_cast<Cantera::FalloffReaction&>(*R);
            
        k0Coeffs[0] = R3.low_rate.preExponentialFactor(); // kmol, s, m
        k0Coeffs[1] = R3.low_rate.temperatureExponent(); 
        k0Coeffs[2] = R3.low_rate.activationEnergy_R(); // K 
            
            kInfCoeffs[0] = R3.high_rate.preExponentialFactor(); // kmol, s, m
        kInfCoeffs[1] = R3.high_rate.temperatureExponent(); 
        kInfCoeffs[2] = R3.high_rate.activationEnergy_R(); // K         

        for (const auto& sp : R3.third_body.efficiencies) 
        {
        label isp = gas.speciesIndex(sp.first);
        efficiencies_[isp] = sp.second;
        }
            
        for (const auto& sp : R3.reactants) 
        {
        label isp = gas.speciesIndex(sp.first);
        currentSpecieCoeff.index = isp;
        currentSpecieCoeff.stoichCoeff = sp.second;
        currentSpecieCoeff.exponent = sp.second;
        lhs.append(currentSpecieCoeff);
        }
            
        for (const auto& sp : R3.products) 
        {
        label isp = gas.speciesIndex(sp.first);
        currentSpecieCoeff.index = isp;
        currentSpecieCoeff.stoichCoeff = sp.second;
        currentSpecieCoeff.exponent = sp.second;
        rhs.append(currentSpecieCoeff);
        }            
        
            // simple fall-off, troe fall-off and SRI fall-off in cantera.            
        // auto falloffPtr = R3.falloff;
        falloffType = R3.falloff->getType();
            
        switch(falloffType)
        {    
                case Cantera::SIMPLE_FALLOFF:
                if (R->reaction_type == 4 ) //FALLOFF_RXN, hard-coded here 
                {
                    addPressureReactionType<FallOffReactionRate>
                    (
                rType,
            lhs,
            rhs,
            efficiencies_,
            k0Coeffs,
            kInfCoeffs,
            k0Coeffs, // dummy values, not used in function addPressurereactionType
            Cantera::SIMPLE_FALLOFF
                    );
                }
                else
                {
                    addPressureReactionType<ChemicallyActivatedReactionRate>
                    (
                rType,
            lhs,
            rhs,
            efficiencies_,
            k0Coeffs,
            kInfCoeffs,
            k0Coeffs, // dummy values, not used in function addPressurereactionType
            Cantera::SIMPLE_FALLOFF
                    );                    
                }
        break;
                
                case Cantera::TROE_FALLOFF:
                TroeR = dynamic_cast<Cantera::Troe&>(*R3.falloff);
        TroeR.getParameters(Troeparams);
        TroeCoeffs[0] = Troeparams[0];
        TroeCoeffs[1] = Troeparams[1];
        TroeCoeffs[2] = Troeparams[2];
        TroeCoeffs[3] = Troeparams[3];
                if (R->reaction_type == 4)   //FALLOFF_RXN 
                {
                    addPressureReactionType<FallOffReactionRate>
                    (
                rType,
            lhs,
            rhs,
            efficiencies_,
            k0Coeffs,
            kInfCoeffs,
            TroeCoeffs,
            Cantera::TROE_FALLOFF
                    );
                }
                else 
                {
                    addPressureReactionType<ChemicallyActivatedReactionRate>
                    (
                rType,
            lhs,
            rhs,
            efficiencies_,
            k0Coeffs,
            kInfCoeffs,
            TroeCoeffs,
            Cantera::TROE_FALLOFF
                    );
                }
                break;

                case Cantera::SRI_FALLOFF:
                SRIR = dynamic_cast<Cantera::SRI&>(*R3.falloff);
                SRIR.getParameters(SRIparams);
        SRICoeffs[0] = SRIparams[0];
        SRICoeffs[1] = SRIparams[1];
        SRICoeffs[2] = SRIparams[2];
        SRICoeffs[3] = SRIparams[3];
        SRICoeffs[4] = SRIparams[4];
                if (R->reaction_type == 4)   // FALLOFF_RXN 
                {                
                    addPressureReactionType<FallOffReactionRate>
                    (
                 rType,
                lhs,
                rhs,
                efficiencies_,
                k0Coeffs,
                kInfCoeffs,
                SRICoeffs,
                Cantera::SRI_FALLOFF
                    );
                }
                else
                {
                    addPressureReactionType<ChemicallyActivatedReactionRate>
                    (
                 rType,
                lhs,
                rhs,
                efficiencies_,
                k0Coeffs,
                kInfCoeffs,
                SRICoeffs,
                Cantera::SRI_FALLOFF
                    );                 
                }

                break;                    
                
                default:
        FatalErrorInFunction
        << "Fall off reaction type not supported: Reaction: "
        << i << exit(FatalError);                                                                    
        }
            
        break;
            default:
        FatalErrorInFunction
        << "reaction type not supported: Reaction "
        << i << exit(FatalError); 
    }
        
        // reinitialize parameters
    rType = unknownReactionType;
    lhs.clear();
    rhs.clear();
    efficiencies_ = 1.0;
    }
//    Info << reactions_ << endl;         
}

template<class ReactionRateType>
void Foam::canteraChemistryReader::addReactionType
(
    const reactionType rType,
    DynamicList<gasHReaction::specieCoeffs>& lhs,
    DynamicList<gasHReaction::specieCoeffs>& rhs,
    const ReactionRateType& rr
)
{
    switch (rType)
    {
        case irreversible:
        {
            reactions_.append
            (
                new IrreversibleReaction
                <Reaction, gasHThermoPhysics, ReactionRateType>
                (
                    ReactionProxy<gasHThermoPhysics>
                    (
                        speciesTable_,
                        lhs.shrink(),
                        rhs.shrink(),
                        speciesThermo_
                    ),
                    rr
                )
            );
        }
        break;

        case reversible:
        {
            reactions_.append
            (
                new ReversibleReaction
                <Reaction, gasHThermoPhysics, ReactionRateType>
                (
                    ReactionProxy<gasHThermoPhysics>
                    (
                        speciesTable_,
                        lhs.shrink(),
                        rhs.shrink(),
                        speciesThermo_
                    ),
                    rr
                )
            );
        }
        break;

        default:
            FatalErrorInFunction
                << "Unknown reaction type " << rType
                << exit(FatalError);
         
    }
}

template<template<class, class> class PressureDependencyType>
void Foam::canteraChemistryReader::addPressureReactionType
(
    const reactionType rType,
    DynamicList<gasHReaction::specieCoeffs>& lhs,
    DynamicList<gasHReaction::specieCoeffs>& rhs,
    const scalarList& efficiencies,
    const scalarList& k0Coeffs,
    const scalarList& kInfCoeffs,
    const scalarList& fallOffCoeffs,
    int fofType
)
{
    scalarList TroeCoeffs(4);
    scalarList SRICoeffs(5);
    switch (fofType)
    {
        case Cantera::SIMPLE_FALLOFF:
        {
            addReactionType
            (
                rType,
                lhs, rhs,
                PressureDependencyType
                    <ArrheniusReactionRate, LindemannFallOffFunction>
                (
                    ArrheniusReactionRate
                    (
                        k0Coeffs[0],
                        k0Coeffs[1],
                        k0Coeffs[2]
                    ),
                    ArrheniusReactionRate
                    (
                        kInfCoeffs[0],
                        kInfCoeffs[1],
                        kInfCoeffs[2]
                    ),
                    LindemannFallOffFunction(),
                    thirdBodyEfficiencies(speciesTable_, efficiencies)
                )
            );
            break;
        }
        case Cantera::TROE_FALLOFF:
        {
            TroeCoeffs = fallOffCoeffs;

            addReactionType
            (
                rType,
                lhs, rhs,
                PressureDependencyType
                    <ArrheniusReactionRate, TroeFallOffFunction>
                (
                    ArrheniusReactionRate
                    (
                        k0Coeffs[0],
                        k0Coeffs[1],
                        k0Coeffs[2]
                    ),
                    ArrheniusReactionRate
                    (
                        kInfCoeffs[0],
                        kInfCoeffs[1],
                        kInfCoeffs[2]
                    ),
                    TroeFallOffFunction
                    (
                        TroeCoeffs[0],
                        TroeCoeffs[1],
                        TroeCoeffs[2],
                        TroeCoeffs[3]
                    ),
                    thirdBodyEfficiencies(speciesTable_, efficiencies)
                )
            );
            break;
        }
        case Cantera::SRI_FALLOFF:
        {
            SRICoeffs = fallOffCoeffs;

            addReactionType
            (
                rType,
                lhs, rhs,
                PressureDependencyType
                    <ArrheniusReactionRate, SRIFallOffFunction>
                (
                    ArrheniusReactionRate
                    (
                        k0Coeffs[0],
                        k0Coeffs[1],
                        k0Coeffs[2]
                    ),
                    ArrheniusReactionRate
                    (
                        kInfCoeffs[0],
                        kInfCoeffs[1],
                        kInfCoeffs[2]
                    ),
                    SRIFallOffFunction
                    (
                        SRICoeffs[0],
                        SRICoeffs[1],
                        SRICoeffs[2],
                        SRICoeffs[3],
                        SRICoeffs[4]
                    ),
                    thirdBodyEfficiencies(speciesTable_, efficiencies)
                )
            );
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Fall-off function type "
                << "not implemented"
                << exit(FatalError);
        }
    }
}

void Foam::canteraChemistryReader::canteraRead
(
    Cantera::IdealGasMix gas
)
{
    int nsp = gas.nSpecies();
    Info<< "Reading cantera files, number of species: " << nsp << endl;
}

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::canteraChemistryReader::canteraChemistryReader
(
    const dictionary& thermoDict,
    speciesTable& species
) 
: 
    ctiFileString(thermoDict.lookup("canteraChemistryFile")),
    ctiFileID(thermoDict.lookup("canteraFileID")),
    gas_(ctiFileString, ctiFileID),
    speciesTable_(setSpecies(gas_,species)),
    reactions_(speciesTable_, speciesThermo_)       
{    
    canteraRead(gas_);
    fileName transportFile
    (
        fileName(thermoDict.lookup("canteraTransportFile")).expand()
    );
    dictionary transportDict;
    transportDict.read(IFstream(transportFile)());
    setSpeciesThermo(gas_, transportDict);
    setSpeciesComposition(gas_);
    setReactions(gas_);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
