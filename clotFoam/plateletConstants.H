struct plateletConstants
{
    // Reference to the inputParameters dict
    IOdictionary& inputParameters;

    // Total number of platelets and number of mobile platelets
    label num_plts
    {
        readLabel(inputParameters.lookup("num_plts"))
    };
    label num_plts_mobile
    {
        readLabel(inputParameters.lookup("num_plts_mobile"))
    };

    // Platelet diameter length
    dimensionedScalar Pdiam
    {
        "Pdiam",
        dimLength,
        inputParameters.lookup("Pdiam")
    };

    // Maximum Packing Limit of Platelets
    dimensionedScalar Pmax
    {
        "Pmax",
        dimless/dimVolume,
        inputParameters.lookup("Pmax")
    };

    // Platelet Diffusion Constant
    dimensionedScalar Dp
    {
        "Dp",
        dimViscosity,
        inputParameters.lookup("Dp")
    };

    // Platelet adhesion binding on-rate via collagen
    dimensionedScalar kAdhCol
    {
        "kAdhCol",
        dimVolume/dimTime,
        inputParameters.lookup("kAdhCol")
    };

    // Platelet activation rate by collagen
    dimensionedScalar kActCol
    {
        "kActCol",
        dimless/dimTime,
        inputParameters.lookupOrDefault("kActCol",kAdhCol.value()*Pmax.value())
    };

    // Platelet adhesion binding on-rate via vWF
    dimensionedScalar kAdhVwfOn
    {
        "kAdhVwfOn",
        dimVolume/dimTime,
        inputParameters.lookup("kAdhVwfOn")
    };

    // Platelet adhesion binding off-rate via vWF
    dimensionedScalar kAdhVwfOff
    {
        "kAdhVwfOff",
        dimless/dimTime,
        inputParameters.lookup("kAdhVwfOff")
    };

    // Platelet cohesion binding on-rate via vWF
    dimensionedScalar kCohVwfOnPmax
    {
        "kCohVwfOnPmax",
        dimless/dimTime,
        inputParameters.lookup("kCohVwfOnPmax")
    };

    // Platelet cohesion binding off-rate via vWF
    dimensionedScalar kCohVwfOff
    {
        "kCohVwfOff",
        dimless/dimTime,
        inputParameters.lookup("kCohVwfOff")
    };

    // Platelet cohesion binding on-rate via fibrinogen
    dimensionedScalar kCohFbgOnPmax
    {
        "kCohFbgOnPmax",
        dimless/dimTime,
        inputParameters.lookup("kCohFbgOnPmax")
    };

    // Platelet activation via vWF/shear
    dimensionedScalar kActVwf
    {
        "kActVwf",
        dimless/dimTime,
        inputParameters.lookup("kActVwf")
    };

    // Platelet activation rate by ADP k_{act}^{ADP}
    dimensionedScalar kActADP
    {
        "kActADP",
        dimless/dimTime,
        inputParameters.lookup("kActADP")
    };

    // Critical concentration of ADP [ADP]^\ast for platelet activation
    dimensionedScalar ADP_ast
    {
        "ADP_ast",
        dimMoles/dimVol,
        inputParameters.lookup("ADP_ast")
    };

    // Degree of the Hill function for ADP
    dimensionedScalar ADP_deg
    {
        "ADP_deg",
        dimless,
        inputParameters.lookupOrDefault("ADP_deg",1.0)
    };

    // Platelet activation rate by thrombin k_{act}^{E2}
    dimensionedScalar kActE2
    {
        "kActE2",
        dimless/dimTime,
        inputParameters.lookup("kActE2")
    };

    // Critcal concentration of thrombin for platelet activation
    dimensionedScalar E2_ast
    {
        "E2_ast",
        dimMoles/dimVolume,
        inputParameters.lookup("E2_ast")
    };

    // Degree of the Hill function for thrombin
    dimensionedScalar E2_deg
    {
        "E2_deg",
        dimless,
        inputParameters.lookupOrDefault("E2_deg",1.0)
    };

    // Critical concentration of virtual substance eta, eta_ast
    scalar eta_ast
    {
        readScalar(inputParameters.lookup("eta_ast"))
    };

    // Threshold for virtual substance eta, eta_t
    scalar eta_t
    {
        readScalar(inputParameters.lookup("eta_t"))
    };

    // Lenght of diffusion for virtual substance eta, l_D_eta
    dimensionedScalar l_D_eta
    {
        "l_D_eta",
        dimLength,
        inputParameters.lookup("l_D_eta")
    };

    // Normalization constant for binding affinity function g(eta), g0
    scalar g0
    {
        ( Foam::pow(eta_ast,3.0) 
                + Foam::pow(1.0 - eta_t,3.0) ) 
               / Foam::pow(1.0-eta_t, 3.0)
    };

    // Critical shear rate for sAdhOn (platelet adhesion on-rate)
    dimensionedScalar sAdhOn_C
    {
        "sAdhOn_C",
        dimless/dimTime,
        inputParameters.lookup("sAdhOn_C")
    };

    // Steepness parameter for sAdhOn function    
    scalar sigmaAdhOn
    {
        readScalar(inputParameters.lookup("sigmaAdhOn"))
    };

    // Critical shear rate for sAdhOff (platelet adhesion off-rate)
    dimensionedScalar sAdhOff_C
    {
        "sAdhOff_C",
        dimless/dimTime,
        inputParameters.lookup("sAdhOff_C")
    };

    // Shape parameters for sAdhOff function a1
    scalar sAdhOff_a1
    {
        readScalar(inputParameters.lookup("sAdhOff_a1"))
    };

    // Shape parameters for sAdhOff function a2
    scalar sAdhOff_a2
    {
        readScalar(inputParameters.lookup("sAdhOff_a2"))
    };

    // Shape parameters for sAdhOff function a3
    scalar sAdhOff_a3
    {
        readScalar(inputParameters.lookup("sAdhOff_a3"))
    };

    // Critical shear rate for sCohOn (platelet cohesion on-rate)
    dimensionedScalar sCohOn_C
    {
        "sCohOn_C",
        dimless/dimTime,
        inputParameters.lookup("sCohOn_C")
    };

    // Steepness parameter for sCohOn function    
    scalar sigmaCohOn
    {
        readScalar(inputParameters.lookup("sigmaCohOn"))
    };

    // Critical shear rate for sCohOff (platelet adhesion on-rate)
    dimensionedScalar sCohOff_C
    {
        "sCohOff_C",
        dimless/dimTime,
        inputParameters.lookup("sCohOff_C")
    };

    // Shape parameters for sCohOff function a1
    scalar sCohOff_a1
    {
        readScalar(inputParameters.lookup("sCohOff_a1"))
    };

    // Shape parameters for sCohOff function a2
    scalar sCohOff_a2
    {
        readScalar(inputParameters.lookup("sCohOff_a2"))
    };

    // Shape parameters for sCohOff function a3
    scalar sCohOff_a3
    {
        readScalar(inputParameters.lookup("sCohOff_a3"))
    };


    // Critical shear rate for sAct (platelet activation by shear)
    dimensionedScalar sAct_C
    {
        "sAct_C",
        dimless/dimTime,
        inputParameters.lookup("sAct_C")
    };

    
    // Steepness parameter for sAct function    
    scalar sigmaAct
    {
        readScalar(inputParameters.lookup("sigmaAct"))
    };

    /***********************************************************************
        Constructor
    ************************************************************************/
    plateletConstants(IOdictionary& inputParameters_)
    :
    inputParameters(inputParameters_)
    {}
};
