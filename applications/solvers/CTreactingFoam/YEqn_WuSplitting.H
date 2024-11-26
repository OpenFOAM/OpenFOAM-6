tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);


PtrList<volVectorField> Dm =
dynamic_cast<const reactingMixture<gasHThermoPhysics>&>
(composition).Dm(p, T, Y, rho);

PtrList<volScalarField> hsi(Y.size());
PtrList<surfaceScalarField> J(Y.size());
PtrList<volScalarField> D(Y.size());
PtrList<volScalarField> tempY(Y.size()), YSave(Y.size());
forAll(Y, i)
{
    hsi.set
    (
        i,        
        new volScalarField
        (
            IOobject
            (
                "hsi",
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar
            (
                "hsi", 
                dimEnergy/dimMass, 
                Zero
            )
        )
    );
    
    tempY.set
    (
        i,        
        new volScalarField
        (
            IOobject
            (
                "tempY",
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            scalar(0.0)
        )
    );
    
    YSave.set(i, Y[i]);
    
    D.set
    (
        i, 
        mag(Dm[i])/(dimensionedScalar("SMALL", dimensionSet(0,-1,0,0,0), Foam::SMALL)
        +mag(fvc::grad(Y[i]).ref()))        
    );

    J.set
    (
        i,         
        linearInterpolate(Dm[i]*rho) & mesh.Sf()
    );    
}

forAll (Y, i)
{
    volScalarField& tHsi = hsi[i];    
    forAll(tHsi, celli)
    {
        tHsi[celli] = composition.Hs(i, p[celli], T[celli]);
    }
    volScalarField::Boundary& Bf = tHsi.boundaryFieldRef();    
    forAll(Bf, patchi)
    {
        forAll(Bf[patchi], facei)    
        {
            Bf[patchi][facei] = 
            composition.Hs
            (
                i,
                p.boundaryField()[patchi][facei],
                T.boundaryField()[patchi][facei]
            );
        }
    }
}


{
    tempY = Y;
    volScalarField Yt(0.0*Y[0]);
    
    // do simple balanced splitting
    // reference: Computer Physics Communications 243 (2019) 81–96
    // evaluate cn term, which will be constant during the splitting integration 

    // -----------------------------step 1, integrate R-cn by dt------------------------------------
    // cn is included in the chemistry model
    // but we do not necessarily need cn to be integrated in the chemistry model 
    // because we can use explicit for Transport (Cn) and then S, it will be no splitting error   
    scalar dtSave = runTime.deltaT().value();
    reaction->correct();
    Qdot = reaction->Qdot();
     
    // --------Substep 1a, integrate the transport term explicitly by dt--------
    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvcDiv(phi, Yi)
              - fvc::laplacian(rho*D[i], Yi)
              == fvOptions(rho, Yi)
            );

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;

        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
    //YSave[inertIndex] = Y[inertIndex];


    // update Rho oldTime field, this could be done with a smarter way: write 
    //    fvm::ddt(rho, Yi)-fvc::ddt(rho, Yi) == reaction->R(Yi)
    // for the transport Equation, here we still use the "update oldTime field" way
    volScalarField& tRho = rho.oldTime(); 
    volScalarField::Boundary& Bf = tRho.boundaryFieldRef();    
    forAll (tRho, celli)
    {
        tRho[celli] = rho[celli];
    }
 
    forAll(Bf, patchi)
    {
        forAll(Bf[patchi], facei)   
        {
            Bf[patchi][facei] = rho.boundaryField()[patchi][facei];
        }
    } 
    
    //update Y.oldTime(), which will be used in fvm::ddt
    forAll (Y, i)
    {
        volScalarField& tYi = Y[i].oldTime();
        
        forAll(tYi, celli)
        {
            tYi[celli] = Y[i][celli];
        }
        volScalarField::Boundary& Bf = tYi.boundaryFieldRef(); 
        forAll(Bf, patchi)
        {
            forAll(Bf[patchi], facei)   
            {
                Bf[patchi][facei] = Y[i].boundaryField()[patchi][facei];
            }
        }
    } 
    // --------Substep 1b, integrate the chemical source term by dt--------    
    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];
            fvScalarMatrix YiEqn
            (
               fvm::ddt(rho, Yi)
             ==
               reaction->R(Yi)
            );
            YiEqn.relax();

            YiEqn.solve(mesh.solver("Yi"));

            Yi.max(0.0);
            Yt += Yi;           
        }
    }
    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);     
    
    // ----------------------- step 2, integrate T+cn by dt/2 ----------------------------
    //update Y.oldTime(), which will be used in fvm::ddt
    forAll (Y, i)
    {
        volScalarField& tYi = Y[i].oldTime();
        
        forAll(tYi, celli)
        {
            tYi[celli] = Y[i][celli];
        }
        volScalarField::Boundary& Bf = tYi.boundaryFieldRef(); 
        forAll(Bf, patchi)
        {
            forAll(Bf[patchi], facei)   
            {
                Bf[patchi][facei] = Y[i].boundaryField()[patchi][facei];
            }
        }
    }     
    Yt = 0.0;
    runTime.setDeltaT(dtSave/2.0);
    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
                + mvConvection->fvmDiv(phi, Yi)
                - fvm::laplacian(rho*D[i], Yi)
             ==
                mvConvection->fvcDiv(phi, YSave[i])
                - fvc::laplacian(rho*D[i], YSave[i])
              + fvOptions(rho, Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
        }
    }
    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);     

    runTime.setDeltaT(dtSave);   
    // update oldY here for CoDAC
    // oldY = tempY;
}


