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

#include "canteraChemistryModel.H"
#include "reactingMixture.H"
//#include "canteraMixture.H"
#include "UniformField.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::canteraChemistryModel<ReactionThermo, ThermoType>::canteraChemistryModel
(
    ReactionThermo& thermo
)
:
    BasicChemistryModel<ReactionThermo>(thermo),
    ODESystem(),
    reactions_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>(this->thermo())
    ),    
    gas_
    (
        dynamic_cast<reactingMixture<ThermoType>&>
            (this->thermo()).canteraGas()
    ),

    Y_(this->thermo().composition().Y()),
    specieThermo_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>
            (this->thermo()).speciesData()
    ),
    
    nSpecie_(Y_.size()),
    nReaction_(gas_.nReactions()),
    Treact_
    (
        BasicChemistryModel<ReactionThermo>::template lookupOrDefault<scalar>
        (
            "Treact",
            0
        )
    ),
    RR_(nSpecie_),
    c_(nSpecie_),
    dcdt_(nSpecie_)
{
    // Info << gas_.nSpecies() << endl;
    // Create the fields for the chemistry sources
    forAll(RR_, fieldi)
    {
        RR_.set
        (
            fieldi,
            new volScalarField::Internal
            (
                IOobject
                (
                    "RR." + Y_[fieldi].name(),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                thermo.p().mesh(),
                dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0)
            )
        );
    }

    Info<< "StandardChemistryModel: Number of species = " << nSpecie_
        << " and reactions = " << nReaction_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::canteraChemistryModel<ReactionThermo, ThermoType>::
~canteraChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
void Foam::canteraChemistryModel<ReactionThermo, ThermoType>::omega
(
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalarField& dcdt
) const
{
    dcdt = Zero;
    
    doublereal X[nSpecie_];
    doublereal wdot[nSpecie_];
    
    for(label i=0; i<nSpecie_; i++)
    {
        X[i] = c[i];
    }    

    gas_.setState_TPX(T, p, X);
    gas_.getNetProductionRates(wdot);
    
    for(label i=0; i<nSpecie_; i++)
    {
        dcdt[i] = wdot[i];
    }        
}


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::canteraChemistryModel<ReactionThermo, ThermoType>::omegaI
(
    const label index,
    const scalarField& c,
    const scalar T,
    const scalar p
) const
{
    doublereal X[nSpecie_];
    doublereal rdot[nSpecie_];
    
    for(label i=0; i<nSpecie_; i++)
    {
        X[i] = c[i];
    }    

    gas_.setState_TPX(T, p, X);
    gas_.getNetRatesOfProgress(rdot);
    
    scalar r = rdot[index];
    
    return(r);
}

// This is a dummy function for EulerImplicit.C
// CanterChemistryModel should not be used with the EulerImplicit solver.
template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::canteraChemistryModel<ReactionThermo, ThermoType>::omegaI
(
    const label index,
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    return(0.0);
}

template<class ReactionThermo, class ThermoType>
void Foam::canteraChemistryModel<ReactionThermo, ThermoType>::derivatives
(
    const scalar time,
    const scalarField& c,
    scalarField& dcdt
) const
{
    const scalar T = c[nSpecie_];
    const scalar p = c[nSpecie_ + 1];

    forAll(c_, i)
    {
        c_[i] = max(c[i], 0);
    }

    omega(c_, T, p, dcdt);

    // Constant pressure
    // dT/dt = ...
    scalar rho = 0;
    scalar cSum = 0;
    for (label i=0; i < nSpecie_; i++)
    {
        const scalar W = specieThermo_[i].W();
        cSum += c_[i];
        rho += W*c_[i];
    }
    scalar cp = 0;
    for (label i=0; i<nSpecie_; i++)
    {
        cp += c_[i]*specieThermo_[i].cp(p, T);
    }
    cp /= rho;

    scalar dT = 0;
    for (label i=0; i < nSpecie_; i++)
    {
        const scalar hi = specieThermo_[i].ha(p, T);
        dT += hi*dcdt[i];
    }
    dT /= rho*cp;

    dcdt[nSpecie_] = -dT;

    // dp/dt = ...
    dcdt[nSpecie_ + 1] = 0;
}


//Cantera is not able to return analytical Jacobian at this time.
// Jacobian is evaluated by numerical perturbation here.
//if analytical Jacobian is needed, please use StandardChemistryModel.
template<class ReactionThermo, class ThermoType>
void Foam::canteraChemistryModel<ReactionThermo, ThermoType>::jacobian
(
    const scalar t,
    const scalarField& c,
    scalarField& dcdt,
    scalarSquareMatrix& J
) const
{    
    scalarField cc(nSpecie_+2, 0.0);
    scalarField dcdtP(nSpecie_+2, 0.0);
    const scalar T = c[nSpecie_];
    const scalar p = c[nSpecie_ + 1];
    
    J = Zero;
    dcdt = Zero;
    
    forAll (cc, i)
    {
    cc[i] = c[i];
        cc[i] = max(cc[i], 0);        
    }   
    derivatives(t, cc, dcdt);  
    
    for (label i=0; i<nSpecie_+1; i++)
    {
    scalar c_save = cc[i];
    scalar dc = 1.0e-15 + c_save*1.0e-9;
    cc[i] = c_save + dc;
    dc = cc[i] - c_save;
        
    derivatives(t, cc, dcdtP);
        
        // species terms, including species on temperature        
        for (label j=0; j<nSpecie_; j++)
        {
        J(j, i) = (dcdtP[j] - dcdt[j])/dc;
        }   

    cc[i] = c_save;
    }
    
    scalarField hi(nSpecie_);
    scalarField cpi(nSpecie_);
    scalarField cpiP(nSpecie_);
    for (label i = 0; i < nSpecie_; i++)
    {
        hi[i] = specieThermo_[i].ha(p, T);
        cpi[i] = specieThermo_[i].cp(p, T);
    }
    scalar cpMean = 0;
    for (label i=0; i<nSpecie_; i++)
    {
        cpMean += c_[i]*cpi[i]; // J/(m3.K)
    }
    scalar dTdt = 0.0;
    for (label i=0; i<nSpecie_; i++)
    {
        dTdt += hi[i]*dcdt[i]; // J/(m3.s)
    }
    dTdt /= -cpMean; // K/s 

    // Temperature derivative on species 
    for (label i = 0; i < nSpecie_; i++)
    {
        J(nSpecie_, i) = 0;
        for (label j = 0; j < nSpecie_; j++)
        {
            J(nSpecie_, i) += hi[j]*J(j, i);
        }
        J(nSpecie_, i) += cpi[i]*dTdt; // J/(mol.s)
        J(nSpecie_, i) /= -cpMean;    // K/s / (mol/m3)
    }

    // purturbation to calculate dcpdTMean

    scalar T_save = T;
    scalar dT = 1.0e-15 + T_save*1.0e-9;
    T_save = T_save + dT;
    dT = T_save - T;
    scalar dcpdTMean = 0;
    for (label i=0; i<nSpecie_; i++)
    {        
    cpiP[i] = specieThermo_[i].cp(p, T_save);
    dcpdTMean += c_[i]*(cpiP[i]-cpi[i])/dT; 
    }

    // ddT of dTdt
    J(nSpecie_, nSpecie_) = 0;
    for (label i = 0; i < nSpecie_; i++)
    {
        J(nSpecie_, nSpecie_) += cpi[i]*dcdt[i] + hi[i]*J(i, nSpecie_);
    }
    J(nSpecie_, nSpecie_) += dTdt*dcpdTMean;
    J(nSpecie_, nSpecie_) /= -cpMean;
// This last term is commented out by Dezhi
// based on his derivation of the temperature source term (MacArt et al. JCP NGA paper)
//    J(nSpecie_, nSpecie_) += dTdt/T;
    
}


template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::canteraChemistryModel<ReactionThermo, ThermoType>::tc() const
{
    tmp<volScalarField> ttc
    (
        new volScalarField
        (
            IOobject
            (
                "tc",
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimTime, small),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );

    scalarField& tc = ttc.ref();

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();
    
    doublereal fwdRate[nReaction_];
    doublereal X[nSpecie_];

    if (this->chemistry_)
    {
        forAll(rho, celli)
        {
            const scalar rhoi = rho[celli];
            const scalar Ti = T[celli];
            const scalar pi = p[celli];

            scalar cSum = 0;

            for (label i=0; i<nSpecie_; i++)
            {
                c_[i] = rhoi*Y_[i][celli]/specieThermo_[i].W();
                cSum += c_[i];
            }
            for(label i=0; i<nSpecie_; i++)
            {
                X[i] = c_[i];
            }            
            gas_.setState_TPX(Ti, pi, X);
            gas_.getFwdRatesOfProgress(fwdRate);


            for (label i=0; i<nReaction_; i++)
            {
                auto R = gas_.reaction(i);
                for (const auto& sp : R->products)
                {
                    tc[celli] += sp.second*fwdRate[i];
                }
                    
            }

            tc[celli] = nReaction_*cSum/tc[celli];
        }
    }

    ttc.ref().correctBoundaryConditions();

    return ttc;
}


template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::canteraChemistryModel<ReactionThermo, ThermoType>::Qdot() const
{
    tmp<volScalarField> tQdot
    (
        new volScalarField
        (
            IOobject
            (
                "Qdot",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0)
        )
    );

    if (this->chemistry_)
    {
        scalarField& Qdot = tQdot.ref();

        forAll(Y_, i)
        {
            forAll(Qdot, celli)
            {
                const scalar hi = specieThermo_[i].Hc();
                Qdot[celli] -= hi*RR_[i][celli];
            }
        }
    }

    return tQdot;
}


template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::canteraChemistryModel<ReactionThermo, ThermoType>::calculateRR
(
    const label ri,
    const label si
) const
{
    tmp<volScalarField::Internal> tRR
    (
        new volScalarField::Internal
        (
            IOobject
            (
                "RR",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0)
        )
    );

    volScalarField::Internal& RR = tRR.ref();

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();
    
    doublereal netRate[nReaction_];
    doublereal X[nSpecie_];
    
    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        const scalar Ti = T[celli];
        const scalar pi = p[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            const scalar Yi = Y_[i][celli];
            c_[i] = rhoi*Yi/specieThermo_[i].W();
        }

        for(label i=0; i<nSpecie_; i++)
        {
            X[i] = c_[i];
        }
            
        gas_.setState_TPX(Ti, pi, X);
        gas_.getNetRatesOfProgress(netRate);

        auto R = gas_.reaction(ri);
        for (const auto& sp : R->reactants)
        {
            if (si == static_cast<int>(gas_.speciesIndex(sp.first)))
            {
                RR[celli] -= sp.second*netRate[ri];
            }
            
        }            
        for (const auto& sp : R->products)
        {
            if (si == static_cast<int>(gas_.speciesIndex(sp.first)))
            {
                RR[celli] += sp.second*netRate[ri];
            }
        }                    
        
        RR[celli] *= specieThermo_[si].W();
    }

    return tRR;
}


template<class ReactionThermo, class ThermoType>
void Foam::canteraChemistryModel<ReactionThermo, ThermoType>::calculate()
{
    if (!this->chemistry_)
    {
        return;
    }

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        const scalar Ti = T[celli];
        const scalar pi = p[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            const scalar Yi = Y_[i][celli];
            c_[i] = rhoi*Yi/specieThermo_[i].W();
        }

        omega(c_, Ti, pi, dcdt_);

        for (label i=0; i<nSpecie_; i++)
        {
            RR_[i][celli] = dcdt_[i]*specieThermo_[i].W();
        }
    }
}


template<class ReactionThermo, class ThermoType>
template<class DeltaTType>
Foam::scalar Foam::canteraChemistryModel<ReactionThermo, ThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    BasicChemistryModel<ReactionThermo>::correct();

    scalar deltaTMin = great;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    scalarField c0(nSpecie_);

    forAll(rho, celli)
    {
        scalar Ti = T[celli];

        if (Ti > Treact_)
        {
            const scalar rhoi = rho[celli];
            scalar pi = p[celli];

            for (label i=0; i<nSpecie_; i++)
            {
                c_[i] = rhoi*Y_[i][celli]/specieThermo_[i].W();
                c0[i] = c_[i];
            }

            // Initialise time progress
            scalar timeLeft = deltaT[celli];

            // Calculate the chemical source terms
            while (timeLeft > small)
            {
                scalar dt = timeLeft;
                this->solve(c_, Ti, pi, dt, this->deltaTChem_[celli]);
                timeLeft -= dt;
            }

            deltaTMin = min(this->deltaTChem_[celli], deltaTMin);

            this->deltaTChem_[celli] =
                min(this->deltaTChem_[celli], this->deltaTChemMax_);

            for (label i=0; i<nSpecie_; i++)
            {
                RR_[i][celli] =
                    (c_[i] - c0[i])*specieThermo_[i].W()/deltaT[celli];
                
//                Info << RR_[i][celli] << " In Solve" << endl;
            }
        }
        else
        {
            for (label i=0; i<nSpecie_; i++)
            {
                RR_[i][celli] = 0;
            }
        }
    }

    return deltaTMin;
}


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::canteraChemistryModel<ReactionThermo, ThermoType>::solve
(
    const scalar deltaT
)
{
    // Don't allow the time-step to change more than a factor of 2
    return min
    (
        this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT)),
        2*deltaT
    );
}


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::canteraChemistryModel<ReactionThermo, ThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


// ************************************************************************* //
