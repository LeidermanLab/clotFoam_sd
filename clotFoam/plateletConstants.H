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

    // Threshold value for when to refine deltaT. Refine if Theta_T > thresholdTheta
    scalar thresholdTheta = inputParameters.lookupOrDefault("thresholdTheta",1.0025);

    // Max threshold for ThetaT where an error is thrown
    scalar thresholdThetaMax = inputParameters.lookupOrDefault("thresholdThetaMax",1.015);
   
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

    // Value of low shear rate = 300
    dimensionedScalar lowShearRate
    {
        "lowShearRate",
        dimless/dimTime,
        inputParameters.lookup("lowShearRate")
    };

    // Value of low shear rate = 1500
    dimensionedScalar highShearRate
    {
        "highShearRate",
        dimless/dimTime,
        inputParameters.lookup("highShearRate")
    };

    // Value of max shear rate = 2500 (default)
    dimensionedScalar maxShearRate
    {
        "maxShearRate",
        dimless/dimTime,
        inputParameters.lookupOrDefault("maxShearRate",2500.0)
    };

    // Unit shear rate (1 [1/s]), used to non-dimensionalize shearRate
    dimensionedScalar one
    {
        "one",
        dimless/dimTime,
        1.0
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
        inputParameters.lookup("kActCol")
    };

    // Platelet adhesion binding on-rate via vWF at low shear
    dimensionedScalar kAdhVwfOn300
    {
        "kAdhVwfOn300",
        dimVolume/dimTime,
        inputParameters.lookup("kAdhVwfOn300")
    };

    // Platelet adhesion binding on-rate via vWF at high shear
    dimensionedScalar kAdhVwfOn1500
    {
        "kAdhVwfOn1500",
        dimVolume/dimTime,
        inputParameters.lookup("kAdhVwfOn1500")
    };

    // Platelet adhesion binding off-rate via vWF at low shear
    dimensionedScalar kAdhVwfOff300
    {
        "kAdhVwfOff300",
        dimless/dimTime,
        inputParameters.lookup("kAdhVwfOff300")
    };

    // Platelet adhesion binding off-rate via vWF at high shear
    dimensionedScalar kAdhVwfOff1500
    {
        "kAdhVwfOff1500",
        dimless/dimTime,
        inputParameters.lookup("kAdhVwfOff1500")
    };

    // Platelet cohesion binding on-rate via vWF at low shear
    dimensionedScalar kCohVwfOnPmax300
    {
        "kCohVwfOnPmax300",
        dimless/dimTime,
        inputParameters.lookup("kCohVwfOnPmax300")
    };

    // Platelet cohesion binding on-rate via vWF at high shear
    dimensionedScalar kCohVwfOnPmax1500
    {
        "kCohVwfOnPmax1500",
        dimless/dimTime,
        inputParameters.lookup("kCohVwfOnPmax1500")
    };

    // Platelet cohesion binding off-rate via vWF at low shear
    dimensionedScalar kCohVwfOff300
    {
        "kCohVwfOff300",
        dimless/dimTime,
        inputParameters.lookup("kCohVwfOff300")
    };

    // Platelet cohesion binding off-rate via vWF at high shear
    dimensionedScalar kCohVwfOff1500
    {
        "kCohVwfOff1500",
        dimless/dimTime,
        inputParameters.lookup("kCohVwfOff1500")
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
    
    // Right shifting of the Hill function for ADP 
    dimensionedScalar ADP_shift
    {
        "ADP_shift",
        dimMoles/dimVolume,
        inputParameters.lookupOrDefault("ADP_shift",0.0)
    };

    // Platelet activation rate by thrombin k_{act}^{E2}
    dimensionedScalar kActE2
    {
        "kActE2",
        dimless/dimTime,
        inputParameters.lookup("kActE2")
    };

    // Critical concentration of thrombin for platelet activation
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

    // Right shifting of the Hill function for E2 
    dimensionedScalar E2_shift
    {
        "E2_shift",
        dimMoles/dimVolume,
        inputParameters.lookupOrDefault("E2_shift",0.0)
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

    // Length of diffusion for virtual substance eta, l_D_eta
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
