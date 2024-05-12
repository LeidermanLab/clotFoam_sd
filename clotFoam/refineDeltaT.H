maxCoRefine = maxCo / 2.;

Info << "\nTime = " << runTime.timeName() 
<< "\n!!!!! max(Theta_T) = "<< max(Theta_T).value()
<< "\n!!! Refinement of deltaT required !!! (# of refinements = " << refinements << ")"
<< "\nCurrent maxCo = " << maxCo << ", current deltaT = "<< runTime.deltaTValue()
<< endl;

if ( refinements > maxRefinements )
{
    
    Info<< "\n!!! Maximum # of deltaT Refinements Reached !!!" 
    << "\n!!!!! max(Theta_T) = "<< max(Theta_T).value()
    << ", # of refinements = " << refinements - 1 << endl;
    FatalError<< nl << exit(FatalError);
}

//restore the previous time state
runTime.TimeState::operator=(tSCurrent);
    
// Refine deltaT by decreasing maxCo and update time
maxCo = max(minCo, maxCo * .9); // was 0.5 * prevMaxCo
#include "CourantNo.H"
#include "setDeltaT.H"
runTime++;

Info << "New time = " << runTime.time().value() 
<< ", New maxCo = " << maxCo << ", new deltaT = "<< runTime.deltaTValue()
<< endl;

// Restore U, p, and platelets to old time
U = U.oldTime();
p = p.oldTime();
forAll(Plt.field,k)
{
    Plt.field[k] = Plt.field[k].oldTime();
}

