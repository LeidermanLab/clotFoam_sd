/**************************************************************************
    Fluid Fields
***************************************************************************/
// Pressure Field
Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Velocity Field
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

// Shear rate field 
Info<< "Reading field shearRate\n" << endl;
volScalarField shearRate
(
    IOobject
    (
        "shearRate",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    Foam::sqrt(2.0) * mag(symm( fvc::grad(U) ))
);

#include "createPhi.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

/**************************************************************************
    Fields for Platelet Reactions
***************************************************************************/
// Virtual substance secreted from all Bound Platelets
Info<< "Reading field eta_B\n" << endl;
volScalarField eta_B
(
    IOobject
    (
        "eta_B",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Virtual substance secreted from eligible fibrinogen bound Platelets
Info<< "Reading field eta_F\n" << endl;
volScalarField eta_F
(
    IOobject
    (
        "eta_F",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    0.0*eta_B
);

// Hadh Field (defines adhesion region)
Info<< "Reading field Hadh\n" << endl;
volScalarField Hadh
(
    IOobject
    (
        "Hadh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

/**************************************************************************
    Platelet Fraction Fields
***************************************************************************/
// Total Platelet Fraction
Info<< "Reading field Theta_T\n" << endl;
volScalarField Theta_T
(
    IOobject
    (
        "Theta_T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Bound Platelet Fraction
Info<< "Reading field Theta_B\n" << endl;
volScalarField Theta_B
(
    IOobject
    (
        "Theta_B",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    0.0*Theta_T
);

// Fibrinogen Bound Platelet Fraction
Info<< "Reading field Theta_F\n" << endl;
volScalarField Theta_F
(
    IOobject
    (
        "Theta_F",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    0.0*Theta_T
);



/**************************************************************************
    Fields for ADP and ADP source term
***************************************************************************/
// ADP Field
Info<< "Reading field ADP\n" << endl;
volScalarField ADP
(
    IOobject
    (
        "ADP",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimMoles/dimVolume, 0.),
    eta_B.boundaryField().types()
);

// sigma_release Field
Info<< "Reading field sigma_release\n" << endl;
volScalarField sigma_release
(
    IOobject
    (
        "sigma_release",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimMoles/dimVolume/dimTime, 0.),
    eta_B.boundaryField().types()
);



/**************************************************************************
    Platelet Object 
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
1) Create derived Species object Plt
2) Set the ptrList sizes
3) Read in IC/BC for the following fields:
        Plt.field[0] -> Mobile Unactivated
        Plt.field[1] -> Mobile Activated
        Plt.field[2] -> vWF Bound Unactivated
        Plt.field[2] -> vWF Bound Unactivated
        Plt.field[3] -> SEBound Activated
4) Initialize the PtrList/fields needed for RK solver
***************************************************************************/

// Object for for platelets
Species_platelet Plt
(
    pltConst.num_plts,
    pltConst,
    Hadh, 
    eta_B,
    eta_F,  
    ADP, 
    Theta_T, 
    Theta_B,
    Theta_F, 
    shearRate
);

// Set ALL ptrList Sizes
Plt.setPtrListSize();

// Mobile Unactivated Platelet Field
Plt.field.set
(
    0,
    new volScalarField
    (   
        IOobject
        (   
            "plt_Pmu",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,                                      
            IOobject::AUTO_WRITE
        ),  
        mesh
    )   
);

// Mobile Activated Platelet Field
Plt.field.set
(    
    1,
    new volScalarField
    (
        IOobject
        (
            "plt_Pma",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);

// vWF Bound Unactivated Platelet Field
Info<< "Reading field Pbvu\n" << endl;
Plt.field.set
(
    2,
    new volScalarField
    (
        IOobject
        (
            "plt_Pbvu",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);

// vWF Bound Activated Platelet Field
Info<< "Reading field Pbva\n" << endl;
Plt.field.set
(
    3,
    new volScalarField
    (
        IOobject
        (
            "plt_Pbva",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);

// Fibrinogen Bound Activated Platelet Field
Info<< "Reading field Pbfa\n" << endl;
Plt.field.set
(
    4,
    new volScalarField
    (
        IOobject
        (
            "plt_Pbfa",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);

// Subendothelial-Bound Unactivated Platelet Field
Info<< "Reading field Pseu\n" << endl;
Plt.field.set
(
    5, 
    new volScalarField
    (
        IOobject
        (
            "plt_Pseu",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);

// Subendothelial-Bound Anactivated Platelet Field
Info<< "Reading field Psea\n" << endl;
Plt.field.set
(
    6, 
    new volScalarField
    (
        IOobject
        (
            "plt_Psea",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);

// Set all necessary fields for RK method = 0*Plt.field
Plt.initializeRKFields();

/**************************************************************************
    SE-Bound Object
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
1) Create derived Species object seBound
2) Set the ptrList sizes
3) Read in IC/BC for the following fields:
        seBound.field[0] -> E0
        seBound.field[1] -> C0
4) Initialize the PtrList/fields needed for RK solver
5) Set initial conditions for seBound species
***************************************************************************/
Info<< "Creating species object for seBound Chemicals\n" << endl;

// Object for for seBound chemicals
Species_seBound seBound
(
    num_seBound,
    injuryPatchID,
    chemConst
);

// Set ALL ptrList Sizes
seBound.setPtrListSize();

// E0 (TF:VIIa)
seBound.field.set
(
    0,
    new volScalarField
    (   
        IOobject
        (   
            "seBound_E0",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,                                      
            IOobject::AUTO_WRITE
        ),  
        mesh
    )   
);

// C0 
seBound.field.set
(    
    1,
    new volScalarField
    (
        IOobject
        (
            "seBound_C0",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);

// Set all necessary fields for RK method = 0 with proper dimensions
seBound.initializeRKFields();

// Set the initial concentration for the seBound species
seBound.setInitialConcentration();

/**************************************************************************
    Fluid-Phase Object
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
1) Create derived Species object fluidPhase
2) Set the ptrList sizes
3) Read in IC/BC for the following fields:
        fluidPhase.field[0] -> S1
        fluidPhase.field[1] -> S2
        fluidPhase.field[2] -> E1
        fluidPhase.field[3] -> E2
4) Initialize the PtrList/fields needed for RK solver
***************************************************************************/

Info<< "Creating species object for fluidPhase Chemicals\n" << endl;
// Object for for fluidPhase chemicals
Species_fluidPhase fluidPhase
(
    num_fluidPhase,
    chemConst
);

// Set ALL ptrList Sizes
fluidPhase.setPtrListSize();

// S1 (FV and FX)
fluidPhase.field.set
(
    0,
    new volScalarField
    (   
        IOobject
        (   
            "fluidPhase_S1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,                                      
            IOobject::AUTO_WRITE
        ),  
        mesh
    )   
);

// S2 (Prothrombin)
fluidPhase.field.set
(
    1,
    new volScalarField
    (   
        IOobject
        (   
            "fluidPhase_S2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,                                      
            IOobject::AUTO_WRITE
        ),  
        mesh
    )   
);

// E1 (FXa)
fluidPhase.field.set
(    
    2,
    new volScalarField
    (
        IOobject
        (
            "fluidPhase_E1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);

// E2 (thrombin)
fluidPhase.field.set
(    
    3,
    new volScalarField
    (
        IOobject
        (
            "fluidPhase_E2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);

// Set all necessary fields for RK method = 0 with proper dimensions
fluidPhase.initializeRKFields();

/**************************************************************************
    Platelet-Bound Object
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
1) Create derived Species object fluidPhase
2) Set the ptrList sizes
3) Read in IC/BC for the following fields:
        pltBound.field[0] -> S1b
        pltBound.field[1] -> S2b
        pltBound.field[2] -> E1b
        pltBound.field[3] -> E2b
        pltBound.field[4] -> C1
        pltBound.field[5] -> C2
4) Initialize the PtrList/fields needed for RK solver
***************************************************************************/

Info<< "Creating species object for pltBound Chemicals\n" << endl;
Species_pltBound pltBound
(
    num_pltBound,
    chemConst
);

// Set ALL ptrList Sizes
pltBound.setPtrListSize();

// S1b (FV)
pltBound.field.set
(
    0,
    new volScalarField
    (   
        IOobject
        (   
            "pltBound_S1b",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,                                      
            IOobject::AUTO_WRITE
        ),  
        mesh
    )   
);

// S2b (Prothrombin)
pltBound.field.set
(
    1,
    new volScalarField
    (   
        IOobject
        (   
            "pltBound_S2b",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,                                      
            IOobject::AUTO_WRITE
        ),  
        mesh
    )   
);

// E1b (Prothrombinase)
pltBound.field.set
(    
    2,
    new volScalarField
    (
        IOobject
        (
            "pltBound_E1b",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);

// E2b (Bound Thrombin)
pltBound.field.set
(    
    3,
    new volScalarField
    (
        IOobject
        (
            "pltBound_E2b",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);

// C1 (S2b:E1b)
pltBound.field.set
(    
    4,
    new volScalarField
    (
        IOobject
        (
            "pltBound_C1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);

// C2 (S1b:E2b)
pltBound.field.set
(    
    5,
    new volScalarField
    (
        IOobject
        (
            "pltBound_C2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);

// Set all necessary fields for RK method = 0 with proper dimensions
pltBound.initializeRKFields();