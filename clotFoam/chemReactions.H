/**************************************************************************
    React the Chemicals M times
***************************************************************************/

// Solve for chemical reactions M_rxn times
for (int m = 0; m < M_rxn; m++ )
{
    //- RK4 Solver
    odeSolver(h_rxn, ptrToChems, num_chemObjs);

    // Ignore oscillations in numerical solutions for all chemicals
    for(int i = 0; i < num_chemObjs; i++)
    {
        forAll(ptrToChems[i]->field, k)
        {
            ptrToChems[i]->field[k] = max
                                        (
                                            0.*ptrToChems[i]->field[k], 
                                            ptrToChems[i]->field[k]
                                        ); // set chems > 0
        }
    }
}

// Update the BC's for the fluidPhase Species
forAll(fluidPhase.field, i)
{
    fluidPhase.field[i].correctBoundaryConditions();
}




