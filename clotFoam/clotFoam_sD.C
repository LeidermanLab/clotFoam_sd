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
    for more deta_Uils.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    clotFoam_sD

Description
    This solver simulates blood clotting in any type of domain.
    The solver can be broken down into 3 major compents:
    1) Fluid:  Transient solver for incompressible, laminar flow of Newtonian 
       fluid with an additional Darcy term

        du/dt = - grad(p') - div[u*grad(u) - nu*grad(u)] - nu*alpha(theta_U_B)*u,
        div(u) = 0,

    2) Platelet Aggregation: hindered transport of 7 platelet species with 
       activation by ADP, thrombin (e2), and shear rate
       
         dPmu/dt = - div[W(theta_U_T)*(u*Pmu - Dp*grad(Pmu))] 
                   + Rmu(Pmu,Pma,Pbvu,Pseu,Psea),
         dPma/dt = - div[W(theta_U_T)*(u*Pma - Dp*grad(Pma))] 
                   + Rma(Pmu,Pma,Pbvu,Pbva,Pseu,Psea),
       dPbvu/dt  = Rbvu(Pmu,Pma,Pbvu,Pbva,Pseu,Psea),
       dPbva/dt  = Rbva(Pmu,Pma,Pbvu,Pbva,Pseu,Psea),
        dPseu/dt = Rseu(Pmu,Pma,Pbvu,Pbva,Pseu,Psea),
        dPsea/dt = Rsea(Pmu,Pma,Pbvu,Pbva,Pseu,Psea),
       d[ADP]/dt = - div[u*[ADP] - Dp*grad([ADP])] + sigma_release(Pbvu,Pbva)

    3) Coagulation: Thrombin (E2) generation via enzymatic reactions
    
       S1 + E0 <=> C0 -> E0 + E1,
       S1 + P1 <=> S1b
       E1 + P1 <=> E1b,
       S2 + P2 <=> S2b,
       E2 + P2 <=> E2b,
       S2b + E1b <=> C1 -> E1b + E2b
       S1b + E2b <=> C2 -> E1b + E2b

    Notes: 
        - The reaction zone must be defined as a patch called "injuryWalls".
          This can be done in blockMesh or with the topoSet tool.
        - The reactive boundary conditions for the fluidPhase species are
          specified in the $FOAM_CASE/0 directory for that species. 

Author: David Montgomery 
        PhD Candidate at Colorado School of Mines 2023
        with help from Dr. Federico Municchi and Dr. Karin Leiderman
\*---------------------------------------------------------------------------*/

// Classes from OpenFOAM
#include "fvCFD.H"
#include "pisoControl.H"
#include "mathematicalConstants.H"

// Classes/structures for managing various species
#include "plateletConstants.H"
#include "chemConstants.H"
#include "Species_baseClass.H"
#include "Species_platelet.H"
#include "Species_seBound.H"
#include "Species_fluidPhase.H"
#include "Species_pltBound.H"

// RK4 Solver
#include "odeSolver.H"

// Rate of ADP release bell function R(tau)
double R_ADP(const double& tau)
{
    return std::exp(-1.*std::pow( tau - 3.0, 2.0)) 
        / std::sqrt(constant::mathematical::pi);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    // Create the constants and fields for simulation
    #include "createConstants.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // Set up the sigma release source term for ADP
    #include "initSigmaReleaseADP.H"

    // Set necessary pointers for each species
    #include "setSpeciesPointers.H"
    
    // Calculate initial Theta_T, Theta_U, Theta_A
    Plt.updateFractions();

    // Define parameters for controling CFL and refining when needed
    scalar maxCoControlDict = 1.;
    scalar maxCoRefine = 1.;
    scalar minCo = 0.20; //was .25
    label refinements = 0;
    label maxRefinements = 20;
    scalar& threshold = pltConst.thresholdTheta;

    //--- Start time loop
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {   
        if (runTime.write())
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;
            
            /*
            // If deltaT has been refined, slowly increase it 
            // Note: this can lead to too many refinements)
            if (maxCoRefine < maxCoControlDict)
            {
                maxCoRefine = max(minCo, min(maxCoRefine + 0.1, maxCoControlDict));
            }
            */
        }

        // Variable time step control variables and adjustments
        const bool adjustTimeStep =
            runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

        maxCoControlDict = 
            runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);
        
        scalar maxCo = min(maxCoControlDict, maxCoRefine);

        scalar maxDeltaT =
            runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);
        #include "CourantNo.H"
        #include "setDeltaT.H"

        // Make a backup copy of curren TimeState (in case refinement is needed)
        TimeState tSCurrent(runTime);
        
        // Update the time using deltaT from CourantNo.H
        runTime++;

        // Solve the Navier-Stokes-Brinkman Equations
        #include "solveFluids.H"

        // Calculate shearRate (used in shear-dependent fxns for Plt reactions)
        shearRate = Foam::sqrt(2.0) * mag(symm( fvc::grad(U) )) ;

        // Transport the platelets dp/dt = - div(W*J)
        #include "plateletTransport.H" 

        // Check if Theta_T > 1, if so refine deltaT and try again
        if (max(Theta_T).value() >= threshold && maxCo > minCo)
        {
            refinements++;
            threshold = min(threshold + 0.005, pltConst.thresholdThetaMax);
            #include "refineDeltaT.H"
            #include "solveFluids.H"
            shearRate = Foam::sqrt(2.0) * mag(symm( fvc::grad(U) )) ;
            #include "plateletTransport.H" 
        }
        
        // Solve the reaction equations
        h_rxn = runTime.deltaT()/M_rxn; // update the reaction time-step size
        
        #include "plateletReactions.H"
        Plt.updateFractions();
        
        // Check if Theta_T > 1, if so refine h_rxn and try again
        if (max(Theta_T).value() >= threshold)
        {
            Info << "\n!!! Platelet Reaction Refinement !!!" << nl 
            << "Time = " << runTime.time().value()  
            << " New m_rxn = " << M_rxn*2 
            << ", max(Theta_T) = "<< max(Theta_T).value() << endl;
            
            h_rxn = runTime.deltaT()/M_rxn/2.; 

            // Restore to value after transport, but before reactions
            forAll(Plt.field,k)
            {
                Plt.field[k] = Plt.fieldOldTime[k]; 
            } 
            #include "plateletReactions.H"
        }

        if (coagReactionsOn)
        {
            #include "fluidPhaseChemTransport.H"
            #include "chemReactions.H"
        }

        // Transport ADP and update sigma_release
        #include "ADP.H" 

        runTime.write();

        if (runTime.write())
        {
            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << "\n max(Theta_T) = "<< max(Theta_T).value()
                << ",  max(shearRate) = "<< max(shearRate).value() <<" 1/s"<< nl << endl;
        }

        // Check if solution is diverging
        #include "isSolutionDiverging.H"
      
    } // End of time loop

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
