class Species_platelet : public Species 
{    
    private:
        const plateletConstants& c;
        const volScalarField& Hadh;
        const volScalarField& eta_U;
        const volScalarField& eta_A;
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
            // Theta_T = Theta_B + mobile speices
            Theta_T = Theta_B + (field[0] + field[1]) / c.Pmax;
        }

        // Activation hill function for any concentration C, given rate kC_Pla & C_ast
        volScalarField A
        (
            const dimensionedScalar& kC_Pla, // Kinitec Rate Constant
            const dimensionedScalar& C_ast, // Critical Value
            const volScalarField& C // Concentration
        )
        {
            return kC_Pla * max(0.*C, C) / (C_ast + max(0.*C,C) );
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

        // Shear-dependent adhesion on-rate function S_{adh+}
        volScalarField sAdhOn()
        {
            // For testing, we make the sAdhOn function 1 everywhere
            return 0.*Theta_B + 1.;
        }

        volScalarField sAdhOff()
        {
            // For testing, we make the sAdhOff function 1 everywhere
            return 0.*Theta_B  + 1.;
        }

        // Shear-dependent cohesion on-rate function S_{coh+}
        volScalarField sCohOn()
        {
            // For testing, we make the sCohOn function 1 everywhere
            return 0.*Theta_B  + 1.;
        }

        // Shear-dependent cohesion off-rate function S_{coh-}
        volScalarField sCohOff()
        {
            // For testing, we make the sCohOff function 1 everywhere
            return 0.*Theta_B  + 1.;
        }

        /*
        volScalarField sAdhOn()
        {
            return 0.5 * ( 1. 
                   + tanh( c.sigmaAdhOn * (shearRate - c.sAdhOn_C) / c.sAdhOn_C ));
        }
        

        // Shear-dependent cohesion on-rate function S_{coh+}
        volScalarField sCohOn()
        {
            return 0.5 * ( 1. 
                   + tanh(c.sigmaCohOn * (shearRate - c.sCohOn_C) / c.sCohOn_C));
        }

        
        
        // Shear-dependent adhesion off-rate function S_{adh-}
        volScalarField sAdhOff()
        {
            return ( exp(c.sAdhOff_a1 * shearRate / c.sAdhOff_C)
                   + exp(c.sAdhOff_a2 * shearRate / c.sAdhOff_C)
                   + exp(c.sAdhOff_a3 * shearRate / c.sAdhOff_C) ) / 3.;
        }
        

        // Shear-dependent cohesion off-rate function S_{coh-}
        volScalarField sCohOff()
        {
            return ( exp(c.sCohOff_a1 * shearRate / c.sCohOff_C)
                   + exp(c.sCohOff_a2 * shearRate / c.sCohOff_C)
                   + exp(c.sCohOff_a3 * shearRate / c.sCohOff_C) ) / 3.;
        }
        */

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
            const volScalarField S_AdhOn = sAdhOn();
            const volScalarField S_AdhOff = sAdhOff();
            const volScalarField S_CohOn = sCohOn();
            const volScalarField S_CohOff = sCohOff();
            const volScalarField S_Act = sAct();
            const volScalarField A_ADP = A(c.kActADP,c.ADP_ast,ADP);
            const volScalarField A_E2 = A(c.kActE2,c.E2_ast,E2); 
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
            k[0] = - c.kAdhVwfOn * S_AdhOn * Hadh_x_seBoundPlts * Pmu 
                    + c.kAdhVwfOff * S_AdhOff * Pseu
                    - c.kCohVwfOnPmax * S_CohOn * (gEtaU + gEtaA)* Pmu
                    + c.kCohVwfOff * S_CohOff * Pbvu
                    - (A_ADP + A_E2)*Pmu;
            
            // Pma: Mobile Activated 
            k[1] = - c.kAdhVwfOn * S_AdhOn * Hadh_x_seBoundPlts * Pma  
                   - c.kCohVwfOnPmax * S_CohOn * (gEtaU + gEtaA) * Pma
                   - c.kCohFbgOnPmax * S_CohOn * gEtaA * Pma
                   + c.kCohVwfOff * S_CohOff * Pbva
                   + (A_ADP + A_E2)*Pmu;
            
            // Pbvu: vWF Bound Unctivated 
            k[2] = - c.kAdhCol * Hadh_x_seBoundPlts * Pbvu
                   + c.kCohVwfOnPmax * S_CohOn * (gEtaU + gEtaA)* Pmu
                   - c.kCohVwfOff * S_CohOff * Pbvu
                   - (A_ADP + A_E2) * Pbvu
                   - c.kActVwf * S_Act * Pbvu;
            
            // Pbva: vWF Bound Activated 
            k[3] = - c.kAdhCol * Hadh_x_seBoundPlts * Pbva
                   + c.kCohVwfOnPmax * S_CohOn * (gEtaU + gEtaA) * Pma
                   - c.kCohVwfOff * S_CohOff * Pbva
                   - c.kCohFbgOnPmax * S_CohOn * gEtaA * Pbva
                   + (A_ADP + A_E2) * Pbvu
                   + c.kActVwf * S_Act * Pbvu;
            
            // Pbfa: Fibrinogen Bound Activated
            k[4] = - c.kAdhCol * Hadh_x_seBoundPlts * Pbfa
                   + c.kCohFbgOnPmax * S_CohOn * gEtaA * (Pma + Pbva);
                   
            // Pseu: SE-Bound Unactivated
            k[5] =  c.kAdhVwfOn * S_AdhOn * Hadh_x_seBoundPlts * Pmu 
                    - c.kAdhVwfOff * S_AdhOff * Pseu
                    - (A_ADP + A_E2) * Pseu
                    - c.kActVwf * S_Act * Pseu
                    - c.kActCol * Pseu;

            // Psea: SE-Bound Activated
            k[6] = c.kAdhVwfOn * S_AdhOn * Hadh_x_seBoundPlts * Pma  
                   + c.kAdhCol * Hadh_x_seBoundPlts * (Pbvu + Pbva + Pbfa)
                   + (A_ADP + A_E2) * Pseu
                   + c.kActVwf * S_Act * Pseu
                   + c.kActCol * Pseu;

        }
};
