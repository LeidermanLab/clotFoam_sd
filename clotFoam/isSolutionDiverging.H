// Check if the pressure is negative
if ( min(p).value() < 0 && runTime.time().value() > 1.0 )
{
    
    Info<< "!!! Negative pressure, program terminating !!!" <<endl;
    Info<< "\n!!!!! min(p) = "<< min(p).value() <<" at t = " << 
    runTime.time().value() << " !!!!!"<< endl;
    FatalError<< nl << exit(FatalError);
}

// Check if Theta_T is larger than 1.015
if ( max(Theta_T).value() > 1.015 )
{
    
    Info<< "!!! Total Platelet Fraction Exceeds 1.015 !!!" <<endl;
    Info<< "\n!!!!! max(Theta_T) = "<< max(Theta_T).value()
    <<" > 1.015 at t = " << 
    runTime.time().value() << " !!!!!"<< endl;
    FatalError<< nl << exit(FatalError);
}