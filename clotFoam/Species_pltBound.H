class Species_pltBound : public Species
{
    private:
        const dimensionedScalar& N_Avo;
        const dimensionedScalar& kC1_on;
        const dimensionedScalar& kC1_off; 
        const dimensionedScalar& kC1_cat;
        const dimensionedScalar& kC2_on;
        const dimensionedScalar& kC2_off; 
        const dimensionedScalar& kC2_cat;
        const dimensionedScalar& kE1_on;
        const dimensionedScalar& kE1_off; 
        const dimensionedScalar& kE2_on;
        const dimensionedScalar& kE2_off; 
        const dimensionedScalar& kS1_on;
        const dimensionedScalar& kS1_off;
        const dimensionedScalar& kS2_on;
        const dimensionedScalar& kS2_off; 
        const scalar& N_1; 
        const scalar& N_2; 
        Species *Plt {nullptr};
        Species *fluidPhase {nullptr};
        

    public: 
        //- Constructor
        Species_pltBound
        (
            const int& numSpecies_,
            const dimensionedScalar& N_Avo_,
            const dimensionedScalar& kC1_on_, 
            const dimensionedScalar& kC1_off_, 
            const dimensionedScalar& kC1_cat_,
            const dimensionedScalar& kC2_on_, 
            const dimensionedScalar& kC2_off_, 
            const dimensionedScalar& kC2_cat_,
            const dimensionedScalar& kE1_on_, 
            const dimensionedScalar& kE1_off_,
            const dimensionedScalar& kE2_on_, 
            const dimensionedScalar& kE2_off_,
            const dimensionedScalar& kS1_on_, 
            const dimensionedScalar& kS1_off_,
            const dimensionedScalar& kS2_on_, 
            const dimensionedScalar& kS2_off_,
            const scalar& N_1_, 
            const scalar& N_2_
        ) 
        : 
            Species(numSpecies_), 
            N_Avo(N_Avo_),
            kC1_on(kC1_on_), kC1_off(kC1_off_), kC1_cat(kC1_cat_),
            kC2_on(kC2_on_), kC2_off(kC2_off_), kC2_cat(kC2_cat_),
            kE1_on(kE1_on_), kE1_off(kE1_off_),
            kE2_on(kE2_on_), kE2_off(kE2_off_),
            kS1_on(kS1_on_), kS1_off(kS1_off_),
            kS2_on(kS2_on_), kS2_off(kS2_off_),
            N_1(N_1_), N_2(N_2_)
        {}

        void setPtrs(Species& Plt_, Species& fluidPhase_)
        {
            Plt = &Plt_;
            fluidPhase = &fluidPhase_;
        }

        // Function that edits the k's for RK4 method.
        void updateKs()
        {
            /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
            This assumes that all k's are defined as 
                k[i] = R_i(input) 
            where R represents the reactions in the ADR equation 
                d/dt(C_i) = - div(u C_i - D grad(C_i)) + R_i(C1,C2,...)
            
            Key for variable references:
            pltBound[0] = S1b,     fluidPhase[0] -> S1,     Plt[0] -> Pmu,
            pltBound[1] = S2b,     fluidPhase[1] -> S2,     Plt[1] -> Pma,
            pltBound[2] = E1b,     fluidPhase[2] -> E1,     Plt[2] -> Pba,
            pltBound[3] = E2b,     fluidPhase[3] -> E2,     Plt[3] -> Pbse.
            pltBound[4] = C1,
            pltBound[5] = C2,
            
            input refers to the function input from RK4, 
            e.g. input[0] = S2_b + h/2 * k1

            Output is k[i] = R_i(input)
            \* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

            // pltBound Reference Vars
            const volScalarField& S1b = input[0];
            const volScalarField& S2b = input[1];
            const volScalarField& E1b = input[2];
            const volScalarField& E2b = input[3];
            const volScalarField& C1 = input[4];
            const volScalarField& C2 = input[5];

             // fluidPhase Reference Vars
            const volScalarField& S1 = fluidPhase->input[0];
            const volScalarField& S2 = fluidPhase->input[1];
            const volScalarField& E1 = fluidPhase->input[2];
            const volScalarField& E2 = fluidPhase->input[3];

            // Platelet Reference Vars
            const volScalarField& Pba = Plt->input[2];
            const volScalarField& Pbse = Plt->input[3];

            // S1b
            k[0] = kS1_on * (N_1 * (Pba + Pbse) / N_Avo
                   - (S1b + E1b + C1 + C2)) * S1
                   - kS1_off * S1b
                   - kC2_on * S1b * E2b
                   + kC2_off * C2;

            // S2b
            k[1] = kS2_on * (N_2 * (Pba + Pbse) / N_Avo 
                   - (S2b + E2b + C1 + C2))* S2 
                   - kS2_off * S2b 
                   - kC1_on * S2b * E1b 
                   + kC1_off*C1;
            
            // E1b
            k[2] = kE1_on * (N_1 * (Pba + Pbse) / N_Avo 
                   - (S1b + E1b + C1 + C2))* E1 
                   - kE1_off * E1b
                   - kC1_on * S2b * E1b
                   + (kC1_off + kC1_cat)*C1
                   + kC2_cat*C2;
            
            // E2b
            k[3] = kE2_on * (N_2 * (Pba + Pbse) / N_Avo 
                   - (S2b + E2b + C1 + C2))* E2
                   - kE2_off * E2b
                   - kC2_on * S1b * E2b
                   + (kC2_off + kC2_cat)*C2
                   + kC1_cat * C1;

            // C1
            k[4] = kC1_on * S2b * E1b 
                   -(kC1_off + kC1_cat)*C1;

            // C2
            k[5] = kC2_on * S1b * E2b 
                   -(kC2_off + kC2_cat)*C2;
        }
};