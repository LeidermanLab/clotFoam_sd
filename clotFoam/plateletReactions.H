/**************************************************************************
    React the Platelets M times
***************************************************************************/
{    
    //Info<< "Platelet Reaction Solver: M_react = " << M_rxn << 
    //        ", dt_react = " << h_rxn.value() << endl;
    
    // Solve for substance eta_B
    eta_B = max(0*Theta_B, min(1.+0*Theta_B,Theta_B));
    eta_B.correctBoundaryConditions(); 
    solve(
        fvm::Sp(1.0,eta_B)
        - fvc::Sp(1.0,eta_B)
        -fvm::laplacian(pltConst.l_D_eta*pltConst.l_D_eta/4.0, eta_B)
    ); 
    eta_B.correctBoundaryConditions(); 
    eta_B = max(0*eta_F, min(1.+0*eta_B,eta_B));

    // Solve for substance eta_F
    eta_F = max(0*Theta_F, min(1.+0*Theta_F,Theta_F));
    eta_F.correctBoundaryConditions(); 
    solve(
        fvm::Sp(1.0,eta_F)
        - fvc::Sp(1.0,eta_F)
        -fvm::laplacian(pltConst.l_D_eta*pltConst.l_D_eta/4.0, eta_F)
    ); 
    eta_F.correctBoundaryConditions(); 
    eta_F = max(0*eta_F, min(1.+0*eta_F,eta_F));

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
