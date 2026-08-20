class Species_platelet : public Species 
{    
    private:
        const plateletConstants& c;
        const volScalarField& Hadh;
        const volScalarField& eta_U;
        const volScalarField& eta_A;
        const volScalarField& cohOnOff;
        const volScalarField& ADP;
        volScalarField& Theta_T;
        volScalarField& Theta_B;
        volScalarField& Theta_U;
        volScalarField& Theta_A;
        volScalarField& shearRate;
        Species *fluidPhase {nullptr};

    public:
        //- Constructor
        Species_platelet
        (
            const int numSpecies_,
            const plateletConstants& c_,
            const volScalarField& Hadh_,
            const volScalarField& eta_U_,
            const volScalarField& eta_A_,
            const volScalarField& cohOnOff_,
            const volScalarField& ADP_,
            volScalarField& Theta_T_,
            volScalarField& Theta_B_,
            volScalarField& Theta_U_,
            volScalarField& Theta_A_,
            volScalarField& shearRate_
            
        )
        :
            Species(numSpecies_),
            c(c_),
            Hadh(Hadh_), 
            eta_U(eta_U_), 
            eta_A(eta_A_), 
            cohOnOff(cohOnOff_), 
            ADP(ADP_), 
            Theta_T(Theta_T_), 
            Theta_B(Theta_B_), 
            Theta_U(Theta_U_),
            Theta_A(Theta_A_),
            shearRate(shearRate_),
            fluidPhase(nullptr)
        {}
        
        //- Member Functions
        void setPtrs(Species& fluidPhase_)
        {
            fluidPhase = &fluidPhase_;
        }
        
        // update the platelet fractions
        void updateFractions()
        {
            // Platelets that are bound and unactivated
            // Theta_U = (Pbvu + Pseu) / Pmax
            Theta_U = (field[2] + field[5]) / c.Pmax;

            // Platelets that are bound and activated
            // Theta_A = (Pbva + Pbfa + Psea) / Pmax
            Theta_A = (field[3] + field[4] + field[6]) / c.Pmax;

            // All bound platelets (used in Brinkman term)
            // Theta_B = Theta_U + Theta_A
            Theta_B = Theta_U + Theta_A;

            // All platelets (used in hindered transport function)
            // Theta_T = Theta_B + mobile species
            Theta_T = Theta_B + (field[0] + field[1]) / c.Pmax;
        }

        // Ignore small oscilations due to numerical methods
        void limitPlatelets()
        {
            forAll(field,i)
            {
                field[i] = max
                            (
                                0.*field[i],
                                min
                                (
                                    field[i], 
                                    0.*field[i] +  c.Pmax
                                )
                            ); 
            }
        }

        // Activation hill function for any concentration C, given rate kC_Pla & C_ast
        volScalarField A
        (
            const dimensionedScalar& kC_Pla, // Kinitec Rate Constant
            const dimensionedScalar& C_ast, // Critical Value
            const dimensionedScalar& d, // degree of Hill function
            const dimensionedScalar& s, // parameter for shifting right
            const volScalarField& C // Concentration
        )
        {
            return kC_Pla * max(0.*pow(C,d), pow(C-s,d) ) / (pow(C_ast,d) + max(0.*pow(C,d),pow(C-s,d)) );
        }

        // Platelet binding affinity function g(eta)  
        volScalarField g
        (
            const volScalarField& eta
        )
        {
            return  c.g0 * max(
                            0*eta,
                            min
                            (
                                1. + 0*eta,
                                Foam::pow(eta - c.eta_t,3.0) /
                                ( Foam::pow(c.eta_ast,3.0) 
                                  + Foam::pow(eta - c.eta_t,3.0) )
                            )
                            );
        }

        // Shear-dependent linear function for vWF off rates
        volScalarField linearScaling
        (
            const dimensionedScalar& K_L,
            const dimensionedScalar& K_H
        )
        {
            const volScalarField& S = shearRate;
            const dimensionedScalar& S_L = c.lowShearRate; 
            const dimensionedScalar& S_H = c.highShearRate; 
            const dimensionedScalar& S_M = c.maxShearRate;
            // 1st order Lagrange interpolating polynomial to 2k then exponential but flatish to 6000 1/s (Edited by Jake Grdadolnik 9/16/25)
            return 
	    min
            (
               K_L * (S-S_H)/(S_L-S_H) + K_H * (S-S_L)/(S_H-S_L),
               (K_H/K_H.value() * (0*S+S_L)/(S_L))*Foam::exp( Foam::log(K_L.value() * (2000.0-S_H.value())/(S_L.value()-S_H.value()) + K_H.value() * (2000.0-S_L.value())/(S_H.value()-S_L.value())) + Foam::pow(mag(S/1.0-2000.0)/4000.0,4.6)*Foam::log(1.02))
	    );
	}

        // Thresholded shear-dependent linear function for vWF on rates
        volScalarField minLinearScaling
        (
            const dimensionedScalar& K_L,
            const dimensionedScalar& K_H
        )
        {
            const volScalarField& S = shearRate;
            const dimensionedScalar& S_L = c.lowShearRate; 
            const dimensionedScalar& S_H = c.highShearRate; 
            const dimensionedScalar& S_M = c.maxShearRate;

            // 1st order Lagrange interpolating polynomial then hyperbolic, then flat (Edited by Jake Grdadolnik 1/29/25)
            return 
            max
            (
             Foam::tanh(0.0005*(S/1.0-S_H.value()))*(K_L * (.8*S_M-S_H+0*S)/(S_L-S_H) + K_H * (.8*S_M-S_L+0*S)/(S_H-S_L)),
	     min
             (
                K_L * (S-S_H)/(S_L-S_H) + K_H * (S-S_L)/(S_H-S_L),
                K_L * (.8*S_M-S_H+0*S)/(S_L-S_H) + K_H * (.8*S_M-S_L+0*S)/(S_H-S_L)
             )
	    ) ;
        }

        // Shear-dependent activation on-rate function S_{act}
        volScalarField sAct()
        {
            return 0.5 * ( 1.  
                + tanh(c.sigmaAct * (shearRate - c.sAct_C) / c.sAct_C));
        }

        // Function that edits the k's for RK4 method.  This assumes that 
        // all k's are defined as 
        //      k[i] = R_i(input) 
        // where R represents the reactions in the ADR equation 
        //      d/dt(C_i) = - div(u C_i - D grad(C_i)) + R_i(C1,C2,...)
        void updateKs()
        {
            /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
            Key for variable references (note mobile platetlets come first)
            Plt[0] -> Pmu,   fluidPhase[0] -> S1,   
            Plt[1] -> Pma,   fluidPhase[1] -> S2,   
            Plt[2] -> Pbvu,  fluidPhase[2] -> E1,   
            Plt[3] -> Pbva,  fluidPhase[3] -> E2 (thrombin)
            Plt[4] -> Pbfa,
            Plt[5] -> Pseu,
            Plt[6] -> Psea,

            input refers to the function input from RK4, 
            e.g. input[0] = Pmu_0 + h/2 * k1

            Output is k[i] = R_i(input)
            \* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

            // Quick references to various fields
            const volScalarField& Pmu = input[0];
            const volScalarField& Pma = input[1];
            const volScalarField& Pbvu = input[2];
            const volScalarField& Pbva = input[3];
            const volScalarField& Pbfa = input[4];
            const volScalarField& Pseu = input[5];
            const volScalarField& Psea = input[6];
            const volScalarField& E2 = fluidPhase->field[3];

            // Temporary volScalarFields
            const volScalarField K_AdhOnVwf = minLinearScaling(c.kAdhVwfOn300,c.kAdhVwfOn1500);
            const volScalarField K_AdhOffVwf = linearScaling(c.kAdhVwfOff300,c.kAdhVwfOff1500);
            const volScalarField K_CohOnVwfPmax = minLinearScaling(c.kCohVwfOnPmax300,c.kCohVwfOnPmax1500);
            const volScalarField K_CohOff = linearScaling(c.kCohVwfOff300,c.kCohVwfOff1500);
            const volScalarField K_ActVwf = c.kActVwf * sAct();
            const volScalarField A_ADP = A(c.kActADP,c.ADP_ast,c.ADP_deg,c.ADP_shift,ADP);
            const volScalarField A_E2 = A(c.kActE2,c.E2_ast,c.E2_deg,c.E2_shift,E2); 
            volScalarField gEtaU = g(eta_U); 
            volScalarField gEtaA = g(eta_A); 
            const volScalarField Hadh_x_seBoundPlts
            (
                Hadh * 
                    max
                    (
                        0*Pseu, 
                        min
                        (
                            c.Pmax - Pseu - Psea - Pbvu - Pbva - Pbfa,
                            c.Pmax + 0*Pseu
                        )
                    )
            );
            
            // Pmu: Mobile Unactivated 
            k[0] = - K_AdhOnVwf * Hadh_x_seBoundPlts * Pmu 
                    + K_AdhOffVwf * Pseu
                    - K_CohOnVwfPmax * (gEtaU + gEtaA)* Pmu
                    + K_CohOff * Pbvu
                    - (A_ADP + A_E2)*Pmu;
            
            // Pma: Mobile Activated 
            k[1] = - K_AdhOnVwf * Hadh_x_seBoundPlts * Pma  
                   - K_CohOnVwfPmax * (gEtaU + gEtaA) * Pma
                   - K_CohOnVwfPmax * gEtaA * Pma
                   + K_CohOff * Pbva
                   + (A_ADP + A_E2)*Pmu;
            
            // Pbvu: vWF Bound Unctivated 
            k[2] = - K_AdhOnVwf * Hadh_x_seBoundPlts * Pbvu
                   + K_CohOnVwfPmax * (gEtaU + gEtaA)* Pmu
                   - K_CohOff * Pbvu
                   - (A_ADP + A_E2) * Pbvu
                   - K_ActVwf * Pbvu;
            
            // Pbva: vWF Bound Activated 
            k[3] = - c.kAdhCol * Hadh_x_seBoundPlts * Pbva
                   + K_CohOnVwfPmax * (gEtaU + gEtaA) * Pma
                   - K_CohOff * Pbva
                   - K_CohOnVwfPmax * gEtaA * Pbva
                   + (A_ADP + A_E2) * Pbvu
                   + K_ActVwf * Pbvu;
            
            // Pbfa: Fibrinogen Bound Activated
            k[4] = - c.kAdhCol * Hadh_x_seBoundPlts * Pbfa
                   + K_CohOnVwfPmax * gEtaA * Pma 
                   + K_CohOnVwfPmax * gEtaA * Pbva;
                   
            // Pseu: SE-Bound Unactivated
            k[5] =  K_AdhOnVwf * Hadh_x_seBoundPlts * (Pmu + Pbvu) 
                    - K_AdhOffVwf * Pseu
                    - (A_ADP + A_E2) * Pseu
                    - K_ActVwf * Pseu
                    - c.kActCol * Pseu;

            // Psea: SE-Bound Activated
            k[6] = K_AdhOnVwf * Hadh_x_seBoundPlts * Pma  
                   + c.kAdhCol * Hadh_x_seBoundPlts * (Pbva + Pbfa)
                   + (A_ADP + A_E2) * Pseu
                   + K_ActVwf * Pseu
                   + c.kActCol * Pseu;

        }
};
