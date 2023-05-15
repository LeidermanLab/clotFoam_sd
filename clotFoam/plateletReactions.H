/**************************************************************************
    React the Platelets M times
***************************************************************************/
{    
    //Info<< "Platelet Reaction Solver: M_react = " << M_rxn << 
    //        ", dt_react = " << h_rxn.value() << endl;
    
    // Solve for substance eta_U
    eta_U = max(0*Theta_U, min(1.+0*Theta_U,Theta_U));
    eta_U.correctBoundaryConditions(); 
    solve(
        fvm::Sp(1.0,eta_U)
        - fvc::Sp(1.0,eta_U)
        -fvm::laplacian(pltConst.l_D_eta*pltConst.l_D_eta/4.0, eta_U)
    ); 
    eta_U.correctBoundaryConditions(); 
    eta_U = max(0*eta_A, min(1.+0*eta_U,eta_U));

    // Solve for substance eta_A
    eta_A = max(0*Theta_A, min(1.+0*Theta_A,Theta_A));
    eta_A.correctBoundaryConditions(); 
    solve(
        fvm::Sp(1.0,eta_A)
        - fvc::Sp(1.0,eta_A)
        -fvm::laplacian(pltConst.l_D_eta*pltConst.l_D_eta/4.0, eta_A)
    ); 
    eta_A.correctBoundaryConditions(); 
    eta_A = max(0*eta_A, min(1.+0*eta_A,eta_A));

    // Solve for platelet reactions M_rxn times
    for (int m = 0; m < M_rxn; m++ )
    {      
        //- RK4 Solver
        odeSolver(h_rxn, ptrToPlt, num_pltObjs);
    }

    forAll(Plt.field,i)
    {
        // ignore small oscilations due to numerical methods
        Plt.field[i] = max(0.*Plt.field[i] ,Plt.field[i]); 
    }
}

// Update the BC's for the mobile Platelet Species
for(int i=0; i < pltConst.num_plts_mobile; i++)
{
    Plt.field[i].correctBoundaryConditions();
}
