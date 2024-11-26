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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
// mixture diffusivity (m2/s)  12.180
template<class ThermoType>
inline Foam::PtrList<Foam::volScalarField> Foam::reactingMixture<ThermoType>::D
(
    const volScalarField& p, 
    const volScalarField& T, 
    const PtrList<Foam::volScalarField>& Y,
        const volScalarField& rho
) const
{        
    PtrList<volScalarField> DNew(nsp_);
    
    forAll(Y, i)
    {
        DNew.set
        (
            i,        
            new volScalarField
            (
                IOobject
                (
                    "D",
                    p.mesh().time().timeName(),
                    p.mesh()
                ),
                p.mesh(),
                dimensionedScalar
                (
                    "D", 
                    dimensionSet(0,2,-1,0,0,0,0), 
                    Zero
                )
            )        
        );                
    }
    
    Cantera::vector_fp YNew(nsp_, 0.0);
    Cantera::vector_fp mixDiffs(nsp_, 0.0);    

    // for internal field
    forAll(p.mesh().cells(), celli)
    {
        scalar ptmp = p[celli];
        scalar Ttmp = T[celli];
        scalar rhotmp = rho[celli];
    
           forAll (Y, isp)
           {
            YNew[isp] = Y[isp][celli];
           }
           
        gas_.setState_TPY(Ttmp, ptmp, YNew.data());
          if (constantLewis_ || unityLewis_)
          {
            forAll (Y, isp)
            {
                if ((this->LewisNumber_[isp]) < 1.0e-20)
                {
                    FatalErrorInFunction
                    << "LewisNumber of species " << isp
                    << " is zero, check your input"<< exit(FatalError);
                    }
                mixDiffs[isp] = tranMix_->thermalConductivity()/gas_.cp_mass();
                mixDiffs[isp] = mixDiffs[isp]/(rhotmp*(this->LewisNumber_[isp]));
            }
        }
        else
        {
           tran_->getMixDiffCoeffs(mixDiffs.data());
        }
            
                       
        forAll (Y, isp)
        {
                volScalarField& tD = DNew[isp];
                tD[celli] = mixDiffs[isp];
        }          
    }    
        
    // for boundary fields    
    forAll(p.boundaryField(), patchi)
    {
        forAll(p.boundaryField()[patchi], facei)
        {
            scalar ptmp = p.boundaryField()[patchi][facei];
            scalar Ttmp = T.boundaryField()[patchi][facei];
                        scalar rhotemp = rho.boundaryField()[patchi][facei];
            forAll (Y, isp)
            {
                YNew[isp] = Y[isp].boundaryField()[patchi][facei];
            }
            
            gas_.setState_TPY(Ttmp, ptmp, YNew.data());
            
            if (constantLewis_ || unityLewis_)
            {
                forAll (Y, isp)
                {
                    if ((this->LewisNumber_[isp]) < 1.0e-20)
                    {
                        FatalErrorInFunction
                        << "LewisNumber of species " << isp
                        << " is zero, check your input"<< exit(FatalError);
                    }
                mixDiffs[isp] = tranMix_->thermalConductivity()/gas_.cp_mass();
                mixDiffs[isp] = mixDiffs[isp]/((this->LewisNumber_[isp])*rhotemp);
                        }           
            }
            else
            {
                tran_->getMixDiffCoeffs(mixDiffs.data());
            }
            
            forAll (Y, isp)
            {
                volScalarField& tD = DNew[isp];
        volScalarField::Boundary& Bf = tD.boundaryFieldRef();
        Bf[patchi][facei] = mixDiffs[isp];
            }            
        }
        
    }    
                
    return DNew;        
}


// Eq. 12.178 in reacting flow 
template<class ThermoType>
inline Foam::PtrList<Foam::volScalarField> Foam::reactingMixture<ThermoType>::Dmass
(
    const volScalarField& p, 
    const volScalarField& T, 
    const PtrList<Foam::volScalarField>& Y
) const
{        
    PtrList<volScalarField> DNew(nsp_);
    
    forAll(Y, i)
    {
        DNew.set
        (
            i,        
            new volScalarField
            (
                IOobject
                (
                    "Dmass",
                    p.mesh().time().timeName(),
                    p.mesh()
                ),
                p.mesh(),
                dimensionedScalar
                (
                    "Dmass", 
                    dimensionSet(0,2,-1,0,0,0,0), 
                    Zero
                )
            )        
        );                
    }
    
    Cantera::vector_fp YNew(nsp_, 0.0);
    Cantera::vector_fp mixDiffs(nsp_, 0.0);    

    // for internal field
    forAll(p.mesh().cells(), celli)
    {
        scalar ptmp = p[celli];
        scalar Ttmp = T[celli];
    
        forAll (Y, isp)
        {
            YNew[isp] = Y[isp][celli];
        }
// not working, when Y = 1, put D = 0, need more revision    
        gas_.setState_TPY(Ttmp, ptmp, YNew.data());
        tran_->getMixDiffCoeffsMass(mixDiffs.data());
                   
        forAll (Y, isp)
        {
            volScalarField& tD = DNew[isp];
            tD[celli] = mixDiffs[isp];
        }          
    }    
        
    // for boundary fields    
    forAll(p.boundaryField(), patchi)
    {
        forAll(p.boundaryField()[patchi], facei)
        {
            scalar ptmp = p.boundaryField()[patchi][facei];
            scalar Ttmp = T.boundaryField()[patchi][facei];

            forAll (Y, isp)
            {
                YNew[isp] = Y[isp].boundaryField()[patchi][facei];
            }
            
            gas_.setState_TPY(Ttmp, ptmp, YNew.data());
            
            tran_->getMixDiffCoeffsMass(mixDiffs.data());
            
            forAll (Y, isp)
            {
                volScalarField& tD = DNew[isp];
                volScalarField::Boundary& Bf = tD.boundaryFieldRef();
                Bf[patchi][facei] = mixDiffs[isp];
            }            
        }
        
    }    
                
    return DNew;        
}

// soret effect coefficient
template<class ThermoType>
inline Foam::PtrList<Foam::volScalarField> Foam::reactingMixture<ThermoType>::DT
(
    const volScalarField& p, 
    const volScalarField& T, 
    const PtrList<Foam::volScalarField>& Y
) const
{        
    PtrList<volScalarField> DTNew(nsp_);
    
    forAll(Y, i)
    {
        DTNew.set
        (
            i,        
            new volScalarField
            (
                IOobject
                (
                    "DT",
                    p.mesh().time().timeName(),
                    p.mesh()
                ),
                p.mesh(),
                dimensionedScalar
                (
                    "DT", 
                    dimensionSet(1,-1,-1,0,0,0,0), 
                    Zero
                )
            )        
        );                
    }
    
    Cantera::vector_fp YNew(nsp_, 0.0);
    Cantera::vector_fp thermDiff(nsp_, 0.0);    

    // for internal field
    forAll(p.mesh().cells(), celli)
    {
        scalar ptmp = p[celli];
        scalar Ttmp = T[celli];
    
        forAll (Y, isp)
        {
            YNew[isp] = Y[isp][celli];
        }
    
        gas_.setState_TPY(Ttmp, ptmp, YNew.data());
        
        if (mixtureAveraged_)
        {
            tranMix_->getThermalDiffCoeffs(thermDiff.data());
        }
    
        if (multiComponent_)
        {
            tranMulti_->getThermalDiffCoeffs(thermDiff.data());
        }

                   
        forAll (Y, isp)
        {
            volScalarField& tDT = DTNew[isp];
            tDT[celli] = thermDiff[isp];
        }          
    }    
        
    // for boundary fields    
    forAll(p.boundaryField(), patchi)
    {
        forAll(p.boundaryField()[patchi], facei)
        {
            scalar ptmp = p.boundaryField()[patchi][facei];
            scalar Ttmp = T.boundaryField()[patchi][facei];

            forAll (Y, isp)
            {
                YNew[isp] = Y[isp].boundaryField()[patchi][facei];
            }
            
            gas_.setState_TPY(Ttmp, ptmp, YNew.data());
            
            if (mixtureAveraged_)
            {
                tranMix_->getThermalDiffCoeffs(thermDiff.data());
            }
    
            if (multiComponent_)
            {
                tranMulti_->getThermalDiffCoeffs(thermDiff.data());
            }
            
            forAll (Y, isp)
            {
                volScalarField& tDT = DTNew[isp];
                volScalarField::Boundary& Bf = tDT.boundaryFieldRef();
                Bf[patchi][facei] = thermDiff[isp];
            }            
        }
        
    }    
                
    return DTNew;        
}

// binary species diffusion coefficients  m2/s
// bD[nsp_*j + i] = binaryD[i][j]
template<class ThermoType>
inline Foam::PtrList<Foam::volScalarField> Foam::reactingMixture<ThermoType>::bD
(
    const volScalarField& p, 
    const volScalarField& T, 
    const PtrList<Foam::volScalarField>& Y
) const
{        
    PtrList<volScalarField> bDNew(nsp_*nsp_);
    
    for(label i=0; i<nsp_*nsp_; i++)
    {
        bDNew.set
        (
            i,        
            new volScalarField
            (
                IOobject
                (
                    "bD",
                    p.mesh().time().timeName(),
                    p.mesh()
                ),
                p.mesh(),
                dimensionedScalar
                (
                    "bD", 
                    dimensionSet(0,2,-1,0,0,0,0), 
                    Zero
                )
            )        
        );                
    }
    
    doublereal YNew[nsp_];
    doublereal Bdiff[nsp_*nsp_];    

    // for internal field
    forAll(p.mesh().cells(), celli)
    {
        scalar ptmp = p[celli];
        scalar Ttmp = T[celli];
    
        forAll (Y, isp)
        {
            YNew[isp] = Y[isp][celli];
        }
        
        gas_.setState_TPY(Ttmp, ptmp, YNew);    
        tran_->getMultiDiffCoeffs(nsp_, Bdiff);

/*      for (label i=0; i<nsp_; i++)
        {
            for (label j=0; j<nsp_; j++)
            {
                binaryD[i][j] = Bdiff(nsp_*j + i);
            }
        }
*/        
        for(label isp=0; isp<nsp_*nsp_; isp++)
        {
            volScalarField& tbD = bDNew[isp];
            tbD[celli] = Bdiff[isp];
        }          
    }    
        
    // for boundary fields    
    forAll(p.boundaryField(), patchi)
    {
        forAll(p.boundaryField()[patchi], facei)
        {
            scalar ptmp = p.boundaryField()[patchi][facei];
            scalar Ttmp = T.boundaryField()[patchi][facei];

            forAll (Y, isp)
            {
                YNew[isp] = Y[isp].boundaryField()[patchi][facei];
            }
            
            gas_.setState_TPY(Ttmp, ptmp, YNew);            
            tran_->getBinaryDiffCoeffs(nsp_, Bdiff);
            
            for(label isp=0; isp<nsp_*nsp_; isp++)
            {
                volScalarField& tbD = bDNew[isp];
                volScalarField::Boundary& Bf = tbD.boundaryFieldRef();
                Bf[patchi][facei] =  Bdiff[isp];
            }             
        
        }
        
    }    
                
    return bDNew;        
}

// mixture diffusivity D_km^'
// see Eq. 12.169 and Eq. 12.180 in the book 
// "chemically reacting flow".
template<class ThermoType>
inline Foam::PtrList<Foam::volVectorField> 
       Foam::reactingMixture<ThermoType>::Dm
(
    const volScalarField& p, 
    const volScalarField& T, 
    const PtrList<Foam::volScalarField>& Y,
    const volScalarField& rho    
) const
{            
    
    // VNew: diffusion velocity multiply mass fraction     
    PtrList<volVectorField> VNew(nsp_);
    PtrList<volVectorField> DmNew(nsp_);    
    volVectorField gradT = fvc::grad(T).ref();
    PtrList<volScalarField> X(nsp_);    
    PtrList<volVectorField> gradX(nsp_);
    PtrList<volVectorField> gradY(nsp_);
    
    volScalarField Wm(0.0*Y[0]);  // mixture molecular weight
    volScalarField moleT(0.0*Y[0]);
    
    forAll (Y, isp)
    {
        moleT += Y[isp]/this->W(isp);
    }
    Wm = 1.0/moleT;    
    
    forAll(Y, i)
    {
        VNew.set
        (
            i,        
            new volVectorField
            (
                IOobject
                (
                    "V",
                    p.mesh().time().timeName(),
                    p.mesh()
                ),
                p.mesh(),
                dimensionedVector
                (
                    "V", 
                    dimensionSet(0,1,-1,0,0,0,0), 
                    Zero
                )
            )        
        );

        DmNew.set
        (
            i,        
            new volVectorField
            (
                IOobject
                (
                    "Dm",
                    p.mesh().time().timeName(),
                    p.mesh()
                ),
                p.mesh(),
                dimensionedVector
                (
                    "Dm", 
                    dimensionSet(0,1,-1,0,0,0,0), 
                    Zero
                )
            )        
        );
        
        X.set
        (
            i,
            Y[i]*Wm/this->W(i)
        );
    }
    
    forAll(Y, i)
    {
        gradX.set
        (
            i,
            fvc::grad(X[i]).ref()
        );
        
        gradY.set
        (
            i,
            fvc::grad(Y[i]).ref()
        );        
    }    
            
    volVectorField sumV
    (
        IOobject
        (
            "sumV",
            p.mesh().time().timeName(),
            p.mesh()
        ),
        p.mesh(),
        dimensionedVector
        (
            "sumV", 
            dimensionSet(0,1,-1,0,0,0,0), 
            Zero
        )
    );     
    
    
    if (mixtureAveraged_ || constantLewis_ || unityLewis_)
    {
        forAll(Y, i) 
        {
            VNew[i] = -gradX[i]*D(p, T, Y, rho)[i]*this->W(i)/Wm;
        }                
    }
    
    if (multiComponent_)
    {
        forAll(Y, i)
        {
            sumV *= scalar(0); //put sumV to zero
            forAll(Y, j)
            {
                sumV += this->W(j)*bD(p, T, Y)[nsp_*j + i]*gradX[j]; 
            }
            VNew[i] = sumV*this->W(i)/(Wm*Wm);            
        }
    }
    
    // add thermal diffusion
    // note in mixture averaged, DT = 0    
    forAll(Y, i)
    {
        VNew[i] -= DT(p, T, Y)[i]*gradT/(rho*T);
    }
    
    // add velocity correction, needed especially for mixture averaged 
    volVectorField sumCV(0.0*VNew[0]);
    forAll(Y, i)
    {
        sumCV += VNew[i]; // correct it to make it zero
    }
    forAll(Y, i)
    {
        VNew[i] -= sumCV*Y[i];
    }
    
/*    dimensionedScalar smallNum
    (
        "smallNum",
        gradY[0].dimension(),
        1.0e-20
    ); */
    forAll(Y, i)
    {
        DmNew[i] = VNew[i];
    }
    return DmNew;        
}
// species fluxes in a specific dimension//kg m-2 s-1
template<class ThermoType>
inline Foam::PtrList<Foam::volVectorField> Foam::reactingMixture<ThermoType>::speFlux
(
    const volScalarField& p, 
    const volScalarField& T, 
    const PtrList<Foam::volScalarField>& c  // molar fraction, not concentration
) const
{        
        volVectorField gradT = fvc::grad(T).ref();    
    PtrList<volVectorField> gradC(nsp_);        
    PtrList<volVectorField> speFluxNew(nsp_);
    
    forAll(c, i)
    {
        speFluxNew.set
        (
        i,        
        new volVectorField
        (
            IOobject
            (
                "speFlux",
                p.mesh().time().timeName(),
                p.mesh()
            ),
            p.mesh(),
            dimensionedVector
            (
                "speFlux", 
                dimensionSet(1,-2,-1,0,0,0,0), 
                Zero
            )
        )
        
    );

        gradC.set
        (
            i,
            fvc::grad(c[i]).ref()
        );                
    }
    
    Cantera::vector_fp grad_T(3, 0.0);
    Cantera::Array2D grad_X(nsp_, 3, 0.0);    
    Cantera::Array2D fluxes(nsp_, 3, 0.0);
    Cantera::vector_fp X(nsp_, 0.0);
    

    // for internal field
    forAll(p.mesh().cells(), celli)
    {
        scalar ptmp = p[celli];
        scalar Ttmp = T[celli];
        
        grad_T[0] = gradT[celli].x();
        grad_T[1] = gradT[celli].y();
        grad_T[2] = gradT[celli].z();
    
        forAll (c, isp)
        {
            grad_X(isp, 0) = gradC[isp][celli].x();
            grad_X(isp, 1) = gradC[isp][celli].y();            
            grad_X(isp, 2) = gradC[isp][celli].z();
            X[isp] = c[isp][celli];
        }
    
        gas_.setState_TPX(Ttmp, ptmp, X.data());
    
        if (mixtureAveraged_)
        {
            tranMix_->getSpeciesFluxes
            (
                3, grad_T.data(), nsp_, grad_X.ptrColumn(0), 
                nsp_, fluxes.ptrColumn(0)
            );
        }
    
        if (multiComponent_)
        {
            tranMulti_->getSpeciesFluxes
            (
                3, grad_T.data(), nsp_, grad_X.ptrColumn(0), 
                nsp_, fluxes.ptrColumn(0)
            );
        }
                   
        forAll (c, isp)
        {
            volVectorField& tSF = speFluxNew[isp];
            tSF[celli].x() = fluxes(isp, 0);
            tSF[celli].y() = fluxes(isp, 1);
            tSF[celli].z() = fluxes(isp, 2);
        }          
    }    
    
    
    // for boundary fields    
    forAll(p.boundaryField(), patchi)
    {
        forAll(p.boundaryField()[patchi], facei)
        {
            scalar ptmp = p.boundaryField()[patchi][facei];
            scalar Ttmp = T.boundaryField()[patchi][facei];
            grad_T[0] = gradT.boundaryField()[patchi][facei].x();
            grad_T[1] = gradT.boundaryField()[patchi][facei].y();
            grad_T[2] = gradT.boundaryField()[patchi][facei].z();

            forAll (c, isp)
            {
                grad_X(isp, 0) = gradC[isp].boundaryField()[patchi][facei].x();
                grad_X(isp, 1) = gradC[isp].boundaryField()[patchi][facei].y();    
                grad_X(isp, 2) = gradC[isp].boundaryField()[patchi][facei].z();
                X[isp] = c[isp].boundaryField()[patchi][facei];
            }
            
            gas_.setState_TPX(Ttmp, ptmp, X.data());
            
            if (mixtureAveraged_)
            {
                tranMix_->getSpeciesFluxes
                (
                    3, grad_T.data(), nsp_, grad_X.ptrColumn(0), 
                    nsp_, fluxes.ptrColumn(0)
                );
            }
            
            if (multiComponent_)
            {
                tranMulti_->getSpeciesFluxes
                (
                    3, grad_T.data(), nsp_, grad_X.ptrColumn(0), 
                    nsp_, fluxes.ptrColumn(0)
                );
            }
            forAll (c, isp)
            {
                volVectorField& tSF = speFluxNew[isp];
                volVectorField::Boundary& Bf = tSF.boundaryFieldRef();
                vector flux(fluxes(isp,0), fluxes(isp,1),fluxes(isp,2));
                Bf[patchi][facei] = flux;
            }            
        }
        
    }    
                
    return speFluxNew;    
}


// ************************************************************************* //
