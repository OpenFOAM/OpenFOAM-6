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

#include "reactingMixture.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::reactingMixture<ThermoType>::reactingMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName 
)
:
    speciesTable(),
    autoPtr<chemistryReader<ThermoType>>
    (
        chemistryReader<ThermoType>::New(thermoDict, *this)
    ),
    multiComponentMixture<ThermoType>
    (
        thermoDict,
        *this,
        autoPtr<chemistryReader<ThermoType>>::operator()().speciesThermo(),
        mesh,
        phaseName,
        autoPtr<chemistryReader<ThermoType>>::operator()().canteraGas()
    ),
    PtrList<Reaction<ThermoType>>
    (
        autoPtr<chemistryReader<ThermoType>>::operator()().reactions()
    ),
    gas_(autoPtr<chemistryReader<ThermoType>>::operator()().canteraGas()),    
    speciesComposition_
    (
        autoPtr<chemistryReader<ThermoType>>::operator()().specieComposition()
    )
    
{
    autoPtr<chemistryReader<ThermoType>>::clear();
    nsp_ = speciesComposition_.size();
    
/*    if (thermoDict.found("transportModel"))
    {
        const word transportModel(thermoDict.lookup("transportModel"));
        int log = 0;    
        if (transportModel == "mixtureAveraged") 
        {
            mixtureAveraged_ = true;
            tran_ = Cantera::newTransportMgr("Mix", &gas_, log=0);
        }
        if (transportModel == "multiComponent") 
        {
            multiComponent_ = true;
            tran_ = Cantera::newTransportMgr("Multi", &gas_, log=0);
        }
        
        tranMix_ = dynamic_cast<Cantera::MixTransport*>(tran_);
        tranMulti_ = dynamic_cast<Cantera::MultiTransport*>(tran_);
    }    
*/    
       mixtureAveraged_ = this->MmixtureAveraged_;
    multiComponent_ = this->MmultiComponent_;
    constantLewis_ = this->MconstantLewis_;
        unityLewis_ = this->MunityLewis_;
    if (mixtureAveraged_ || multiComponent_)
    {
        int log = 0;
        if (mixtureAveraged_) 
        {
            tran_ = Cantera::newTransportMgr("Mix", &gas_, log=0);
        }
        if (multiComponent_) 
        {
            tran_ = Cantera::newTransportMgr("Multi", &gas_, log=0);
        }    
        if (constantLewis_) 
        {
            tran_ = Cantera::newTransportMgr("Mix", &gas_, log=0);
        }        
        tranMix_ = dynamic_cast<Cantera::MixTransport*>(tran_);
        tranMulti_ = dynamic_cast<Cantera::MultiTransport*>(tran_);
    }
    Info << "Dezhi in reactingMixture: " << gas_.nSpecies() << endl; 
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::reactingMixture<ThermoType>::read(const dictionary& thermoDict)
{}


// ************************************************************************* //
