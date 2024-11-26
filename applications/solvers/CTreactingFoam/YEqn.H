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


// Added by Taaresh - Mix Avg Diffusivity (m2/s)

PtrList<volScalarField> D =
dynamic_cast<const reactingMixture<gasHThermoPhysics>&>
(composition).D(p, T, Y, rho);

PtrList<volScalarField> hsi(Y.size());
PtrList<surfaceScalarField> J(Y.size());
PtrList<volScalarField> tempY(Y.size());
PtrList<volVectorField> YV(Y.size());


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

    // product of mass fraction and diffusion vel of each specie

    YV.set
    (
    i,
    new volVectorField
    (
        IOobject
        (
        "YV",
        mesh.time().timeName(),
        mesh
        ),
        mesh,
        dimensionedVector("YV", dimVelocity, vector(0.0, 0.0, 0.0))
    )
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

// add turbulence effect on D, assuming constant turbulence Lewis number
// turbulence thermal diffusivity alphat
volScalarField alphat = turbulence->alphat();
//volScalarField nut = turbulence->nut();
//Info << rho*nut/(dimensionedScalar("SMALL", dimensionSet(1,-1,-1,0,0), Foam::SMALL)+alphat) << endl;
const scalar Let(thermo.lookupOrDefault<scalar>("Let", 1.0));
volScalarField Dt = alphat/(rho*Let);
forAll(Y, i)
{
    D[i] = D[i] + Dt;
}

{

    reaction->correct();
    Qdot = reaction->Qdot();
    volScalarField Yt(0.0*Y[0]);

    Vc *= 0.0;

// Compute Correction Velocity

    forAll(Y, i)
    {
        Vc += D[i]*fvc::grad(Y[i]) + (D[i]*Y[i]/Wmix)*fvc::grad(Wmix);
    }

    phiVc = linearInterpolate(rho*Vc) & mesh.Sf();


    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
         == fvm::laplacian(rho*D[i], Yi) + fvc::laplacian(rho*D[i]*Yi/Wmix, Wmix)
        - fvm::div(phiVc, Yi, "div(phiVc,Yi_h)") + 
                reaction->R(Yi)
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

    volVectorField YV_balance(0.0*YV[0]);

    forAll(YV, i)
    {
    // Calculate diff velocity for each species (including inertSpecies) for use in EEqn.H

    YV[i] = (-D[i]*fvc::grad(Y[i])) - (D[i]*Y[i]/Wmix)*fvc::grad(Wmix);
    
        if (i != inertIndex && composition.active(i))
        {
        YV_balance += YV[i]; 
        }
    }
    
    YV[inertIndex] = -YV_balance; // since sum(Y[i]*V[i]) over all species = 0 for global mass conservation

    forAll(Y, i)
    {

    // Species Diffusion Flux term used in EEqn.H

        J.set
        (
            i,         
            linearInterpolate(YV[i]*rho) & mesh.Sf() 
        );    
    }

}


D_O2 = D[O2_index];
D_N2 = D[N2_index];
D_O = D[O_index];

Info << "Max O2 diff = " << max(D[O2_index]) << endl; 

