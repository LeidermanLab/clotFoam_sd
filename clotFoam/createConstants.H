Info<< "\nReading inputParameters\n" << endl;

IOdictionary inputParameters
(
    IOobject
    (
        "inputParameters",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Patch identifier of the injury boundary
label injuryPatchID = mesh.boundaryMesh().findPatchID("injuryWalls");

/**************************************************************************
    Fluid Scalars and Parameters
***************************************************************************/

// Dynamic Viscosity
dimensionedScalar nu
(
    "nu",
    dimViscosity,
    inputParameters.lookup("nu")
);

// Carman-Kozeny Consant
dimensionedScalar C_CK
(
    "C_CK",
    dimless/dimArea,
    inputParameters.lookup("C_CK")
);

/**************************************************************************
    Constants for Platelets
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 See plateletConstants.H. The constants for platelets are created in a 
 structure so they can be easily passed to various classes.  This makes 
 the code easier to modify for more complex platelet/coagulation models.
 ***************************************************************************/

// Create an instance of plateletConstants for passing to various classes
Info<< "\nCreating an instance of plateletConstants\n" << endl;
plateletConstants pltConst(inputParameters);

/**************************************************************************
    Constants for ADP
***************************************************************************/
Info<< "\nReading constants for ADP and sigma_release\n" << endl;
// ADP Diffusion Constant
dimensionedScalar D_ADP
(
    "D_ADP",
    dimViscosity,
    inputParameters.lookup("D_ADP")
);

// Total ADP Released (\hat A)
dimensionedScalar ADP_rel
(
    "ADP_rel",
    dimMoles,
    inputParameters.lookup("ADP_rel")
);

// How often the sigma_release function is calculated
scalar sigma_dt
(
    readScalar(inputParameters.lookup("sigma_dt"))
);

// The length of time sigma_release is calculated (0 -> sigma_Tf)
scalar sigma_Tf
(
    readScalar(inputParameters.lookup("sigma_Tf"))
);

/**************************************************************************
    Params for Chemical Reactions (not Kinetic Rates)
***************************************************************************/
// Avogadro's Constant
dimensionedScalar N_Avo
(
    "N_Avo", 
    dimless/dimMoles,
    inputParameters.lookup("N_Avo")
);

// Number of sebound species
label num_seBound
(
    readLabel(inputParameters.lookup("num_seBound"))
);

// Number of sebound species
label num_fluidPhase
(
    readLabel(inputParameters.lookup("num_fluidPhase"))
);

// Number of platelet-bound chem species
label num_pltBound
(
    readLabel(inputParameters.lookup("num_pltBound"))
);

// Number of time steps for reactions, so that h_rxn = dt / M_rxn
label M_rxn
(   
    readLabel(inputParameters.lookup("M_rxn"))
);

// Time step for reaction solver
dimensionedScalar h_rxn //time step for reactions
(   "h_rxn", 
    dimTime,
    1.0
);

// Diffusion Constant for all chemicals
dimensionedScalar Dc
(
    "Dc",
    dimViscosity,
    inputParameters.lookup("Dc")
);

// Number of binding sites on platelet surface for specific chems
scalar N_1 = readScalar(inputParameters.lookup("N_1"));
scalar N_2 = readScalar(inputParameters.lookup("N_2"));

/**************************************************************************
    Kinetic Rates for Chemical Reactions 
***************************************************************************/
// Kinetic Rates for C0
dimensionedScalar kC0_on
(
    "kC0_on",
    dimVolume/dimMoles/dimTime,
    inputParameters.lookup("kC0_on")
); 
dimensionedScalar kC0_off
(
    "kC0_off",
    dimless/dimTime,
    inputParameters.lookup("kC0_off")
); 

dimensionedScalar kC0_cat
(
    "kC0_cat",
    dimless/dimTime,
    inputParameters.lookup("kC0_cat")
); 

// Kinetic Rates for C1
dimensionedScalar kC1_on
(
    "kC1_on",
    dimVolume/dimMoles/dimTime,
    inputParameters.lookup("kC1_on")
); 
dimensionedScalar kC1_off
(
    "kC1_off",
    dimless/dimTime,
    inputParameters.lookup("kC1_off")
); 

dimensionedScalar kC1_cat
(
    "kC1_cat",
    dimless/dimTime,
    inputParameters.lookup("kC1_cat")
); 

// Kinetic Rates for C2
dimensionedScalar kC2_on
(
    "kC2_on",
    dimVolume/dimMoles/dimTime,
    inputParameters.lookup("kC2_on")
); 
dimensionedScalar kC2_off
(
    "kC2_off",
    dimless/dimTime,
    inputParameters.lookup("kC2_off")
); 

dimensionedScalar kC2_cat
(
    "kC2_cat",
    dimless/dimTime,
    inputParameters.lookup("kC2_cat")
); 

// Kinetic Rates for S1
dimensionedScalar kS1_on
(
    "kS1_on",
    dimVolume/dimMoles/dimTime,
    inputParameters.lookup("kS1_on")
); 
dimensionedScalar kS1_off
(
    "kS1_off",
    dimless/dimTime,
    inputParameters.lookup("kS1_off")
); 

// Kinetic Rates for S2
dimensionedScalar kS2_on
(
    "kS2_on",
    dimVolume/dimMoles/dimTime,
    inputParameters.lookup("kS2_on")
); 
dimensionedScalar kS2_off
(
    "kS2_off",
    dimless/dimTime,
    inputParameters.lookup("kS2_off")
); 

// Kinetic Rates for E1
dimensionedScalar kE1_on
(
    "kE1_on",
    dimVolume/dimMoles/dimTime,
    inputParameters.lookup("kE1_on")
); 
dimensionedScalar kE1_off
(
    "kE1_off",
    dimless/dimTime,
    inputParameters.lookup("kE1_off")
); 

// Kinetic Rates for E2
dimensionedScalar kE2_on
(
    "kE2_on",
    dimVolume/dimMoles/dimTime,
    inputParameters.lookup("kE2_on")
); 
dimensionedScalar kE2_off
(
    "kE2_off",
    dimless/dimTime,
    inputParameters.lookup("kE2_off")
); 