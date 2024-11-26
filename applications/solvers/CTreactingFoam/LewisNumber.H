/*OFstream LewisNumber
(
    "LewisNumber.dat"
);

LewisNumber
    << "species" << "      Lewis Number" << endl;

forAll(mesh.cells(), celli)
{
    if ( abs( T[celli]-max(T).value() ) < 1.0e-20 )
    {
        forAll (Y, k)
        {
            scalar Diff =  dynamic_cast<const reactingMixture<gasHThermoPhysics>&>
                              (composition).D(p, T, Y)[k][celli];
            if (Diff > 1.0e-20) 
            {
                scalar Lewis = thermo.alpha()[celli]/Diff;
                LewisNumber
                << k << tab << tab << Lewis << endl;    
            }

        }
    }
}*/
