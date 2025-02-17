// Check if the pressure is negative
if ( min(p).value() < 0 && runTime.time().value() > 1.0 )
{
    
    Info<< "!!! Negative pressure, program terminating !!!" <<endl;
    Info<< "\n!!!!! min(p) = "<< min(p).value() <<" at t = " << 
    runTime.time().value() << " !!!!!"<< endl;
    FatalError<< nl << exit(FatalError);
}


// Check if Theta_T is larger than the maximum allowed threshold
if ( max(Theta_T).value() > pltConst.thresholdThetaMax )
{
    Info<< "!!! Total Platelet Fraction Exceeds "
        << pltConst.thresholdThetaMax <<" !!!" <<endl;
    Info<< "\n!!!!! max(Theta_T) = "<< max(Theta_T).value()
        <<" > " << pltConst.thresholdThetaMax <<" at t = " 
        << runTime.time().value() << " !!!!!"<< endl;
    FatalError<< nl << exit(FatalError);
}

