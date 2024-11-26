#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;
autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));

psiReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

// initilize omega, diff, omegaT and diffT for use in CEMA
//PtrList<volScalarField>& omega = composition.omega();
//PtrList<volScalarField>& diff = composition.diff();
//volScalarField& omegaT = composition.omegaT();
//volScalarField& diffT = composition.diffT();

// initilize oldY for CoDAC, it will be large numbers to aviod similarity initially
//PtrList<volScalarField>& oldY = composition.oldY(); 
//forAll(Y, i)
//{
//    oldY[i] = 1000.0;
//}
const word inertSpecie(thermo.lookup("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species() << exit(FatalIOError);
}

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volVectorField Vc
(
    IOobject
    (
        "Vc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("Vc", dimVelocity, vector(0,0,0))
);

volScalarField& p = thermo.p();

#include "compressibleCreatePhi.H"

pressureControl pressureControl(p, rho, pimple.dict(), false);

mesh.setFluxRequired(p.name());

Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating reaction model\n" << endl;
autoPtr<CombustionModel<psiReactionThermo>> reaction
(
    CombustionModel<psiReactionThermo>::New(thermo, turbulence())
);

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

volScalarField Qdot
(
    IOobject
    (
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Qdot", dimEnergy/dimVolume/dimTime, 0.0)
);

// add mixture fraction field
volScalarField Zmix 
(
    IOobject
    (
        "Zmix",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    scalar(0.0)
);

// add progress variable field
volScalarField Prog 
(
    IOobject
    (
        "Prog",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    scalar(0.0)
);
#include "createMRF.H"
#include "createFvOptions.H"

label O_index, O2_index, N2_index;

forAll(Y,i)
{
    if (Y[i].member() == "O")
    {
        O_index = i;
        Info << "O_index = " << O_index << endl;
        Info << Y[O_index].member() << endl;
    }
    else if (Y[i].member() == "O2")
    {
        O2_index = i;
        Info << "O2_index = " << O2_index << endl;
        Info << Y[O2_index].member() << endl;
    }
    else if (Y[i].member() == "N2")
    {
        N2_index = i;
        Info << "N2 index = " << N2_index <<  endl;
    }

}

volScalarField D_O2
{
    IOobject
    (
        "D_O2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("D_O2", dimArea/dimTime, 0.0)
};

volScalarField D_N2
{
    IOobject
    (
        "D_N2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("D_N2", dimArea/dimTime, 0.0)
};


volScalarField D_O
{
    IOobject
    (
        "D_O",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("D_O", dimArea/dimTime, 0.0)
};



volScalarField gradY_O2
{
    IOobject
    (
        "gradY_O2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("gradY_O2", dimless/dimLength, 0.0)
};

OFstream maxDataFile
(
        "max_data_myrf.dat"
);

volScalarField Wmix
{
    IOobject
    (
        "Wmix",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.W()
};


surfaceScalarField phiVc
(
    IOobject
    (
        "phiVc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearInterpolate(rho*Vc) & mesh.Sf()
);




