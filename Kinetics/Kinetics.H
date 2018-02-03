#ifndef KINETICS_H
#define KINETICS_H

namespace Foam
{
class Kinetics {

public:
    Kinetics
    (
        fvMesh &mesh_, Time &runTime_, volScalarField &alpha2_,
        volScalarField &TS_
    )
    :
    mesh(mesh_),
    runTime(runTime_),
    alpha2(alpha2_),
    TS(TS_)
    {};

    scalar dQdt(scalar Qc, scalar Qo){
        return (mag(Qc - Qo)/max(runTime.deltaTValue(),ROOTVSMALL));
    }

    scalar LliqMax (scalar& tmptr){
        // LliqMax - maximum allowable amount of liquid blowing agent (n-pentane) in liquid
        // @param - double &tmptr, the input temperature.
        // constants:
        const scalar a = 0.0064, h = 0.0551, T0 = 298.0, ww = 17.8;
        scalar lMax;
        if (tmptr > T0)
        {
            scalar tempDummy = Foam::pow((tmptr - T0),2.0);
            lMax = (a + h*Foam::exp((-tempDummy/(2.0*ww*ww))));
        }
        else
        {
            lMax = (a+h);
        }
    return (lMax);
    }

    scalar ddT_LliqMax (double& tmptr){
    // ddT_LliqMax - derivative of LliqMax with respect to temperature
    // @param - double &tmptr, the input temperature.
    // constants:
    const scalar h = 0.0551, T0 = 298.0, ww = 17.8;
    return (
            (-h*(tmptr-T0)*
            Foam::exp((-(Foam::pow((tmptr - T0),2.0))/(2.0*ww*ww))))/(ww*ww)
           );
    }

    scalar xBL (scalar &T, scalar &dxdT){
    // xBL - mole fraction of blowing agent (R-11) in liquid polymer
    // @param - double &T - temperature
    // @param - double &dxdT - constant
    scalar xBL_value;
    xBL_value = dxdT*(T - 300.0) + 0.5;

    if (xBL_value < 0.0)
    {
        xBL_value = scalar(0.0);
    }
    else if (xBL_value > 0.5)
    {
        xBL_value = scalar(0.5);
    }
    else
    {
        xBL_value = xBL_value;
    }
    return (xBL_value);
    }

    scalar wBL_D(scalar &xBL, scalar &M_B, scalar &M_NCO){
    // wBL_D - weight fraction of maximum allowable blowing agent (R-11) in liquid
    // @param - double &xBL - mole fraction of blowing agent (R-11)
    // @param - double &M_B - molecular weight of blowing agent (R-11)
    // @param - double &M_NCO - molecular weight of NCO
        scalar Lm;
        Lm  = (xBL/(1.0 - xBL))*(M_B/M_NCO);

        if (Lm < 1.0e-4)
        {
            Lm = 1.0e-4;
        }
        return (Lm);
    }
    scalar henryCoefficient(scalar &T){
    // henryCoefficient - Henry coefficient for CO2
    // @param - double &T - Temperature
        if (T < 600)
        {
            // constants
            scalar a = 1.771e7;
            scalar b = -1.134e5;
            scalar c = 320.2;
            scalar d = -0.2563;

            return (a + b*T + c*T*T + d*T*T*T);
        }
        else
        {
            return (8.5e6);
        }
    }
    scalar bubbleRadius (const scalar m0, const scalar m1){
    // bubbleRadius - radius of bubbles based on moments
    // @param - const double m0 - moment of order zero
    // @param - const double m1 - moment of order one
        if (m0 != 0.0 && m1 != 0.0)
        {
            scalar R;
            R   = Foam::pow((3.0*m1/(4.0*M_PI*m0)), 1.0/3.0);
            return R;
        }
        else
        {
            return (30e-6);
        }
    }
    scalar partialPressureCO2
        (
         scalar &M_CO2, scalar &M_B, scalar &surfaceTension,
         scalar &wCO2_g, double &wBA_g, double &p, double &R
        ){
    // partialPressureCO2 - partial pressure of CO2
    // @param - double &M_CO2 - molecular weight of CO2
    // @param - double &M_B - molecular weight of blowing agent
    // @param - double &surfaceTension - surface tenstion
    // @param - double &wCO2_g - weight fraction of CO2 in gas
    // @param - double &wBA_g - weight frcation of blowing agent in gas
    // @param - double &p - ambient pressure
    // @param - double &R - bubble radius
        scalar pCO2;

        if (wCO2_g == 0.0)
        {
            pCO2 = ROOTVSMALL;
        }
        else
        {
            pCO2 = ((wCO2_g/M_CO2)/(wBA_g/M_B + wCO2_g/M_CO2))
                  *(p + 2*surfaceTension/Foam::max(R,ROOTVSMALL));
        }
        return (pCO2);
    }
    scalar wCO2Max
           (
            scalar &M_CO2, scalar &M_liq,
            scalar &pCO2, scalar &henryCoeff
           ){
    // wCO2Max - dissolved amount of CO2 in liquid
    // @param - double &M_CO2 - molecular weight of CO2
    // @param - double &M_liq - molecular weight of liquid mixture
    // @param - double &pCO2 - partial pressure of CO2
    // @param - double &henryCoeff - Henry coefficient
        if ((henryCoeff - pCO2 ) > 0.0)
        {
            return ((M_CO2/M_liq)*(pCO2/(henryCoeff - pCO2)));
        }
        else
        {
            Info<< "\nWarning! Invalid wCO2Max value!" << endl;
            Info<< "'wCO2Max' is replaced by a constant." << endl;
            return (4.4e-4);
        }
    }
    scalar creamTemperature(scalar &xBL0, scalar &dxdT)
    {
    // creamTemperature - The temperature that foaming process starts
    // @param - double &xBL0 - initial mole fraction of the blowing agent (R-11)
    // @param - double &dxdT - constant
        return ((xBL0 - 0.5)/dxdT + 300.0);
    }
    scalar arrhenius (scalar& AOH, scalar& EOH, scalar& tempt)
    {
    // arrhenius - Arrhenius function
    // @param - double& AOH - pre-exponential factor
    // @param - double& EOH - activation energy
    // @param - double& tempt - temperature
        const scalar R = 8.3145; // J/mol K
        return (AOH*Foam::exp(-EOH/(R*tempt)));
    }

    scalar BAdilution(scalar& L_l, scalar& rhoPoly, scalar& rhoBL)
    {
    // BAdilution - dilution term for the blowing agent (n-pentane)
    // @param - double& L_l - weight fraction of liquid blowing agent
    // @param - double& rhoPoly - density of liquid mixture (polymer)
    // @param - double& rhoBL - density of blowing agent
        return (1/(1 + L_l*(rhoPoly/(Foam::max(rhoBL,ROOTVSMALL)))));
    }
    scalar exothermicGelling
    (
        scalar& DH_OH,
        scalar& COH_0,
        scalar& dXOHdt,
        scalar& rhoPoly,
        scalar& C_TOT
    )
    {
        return ((-DH_OH*COH_0*dXOHdt)/(rhoPoly*C_TOT));
    }
    scalar exothermicBlowing
    (
        scalar& DH_W,
        scalar& CW_0,
        scalar& dWdt,
        scalar& rhoPoly,
        scalar& C_TOT
    )
    {
        return ((-DH_W*CW_0*dWdt)/(rhoPoly*C_TOT));
    }
    scalar endothermicEvaporation
    (
        scalar& latenth,
        scalar& dLdt,
        scalar& rhoPoly,
        scalar& C_TOT
    )
    {
        return ((latenth*dLdt)/(rhoPoly*C_TOT));
    }
    scalar R11EvaporationRate
    (
        scalar M_B,
        scalar M_NCO,
        scalar xBL_val,
        scalar dxdT,
        scalar latenth
    )
    {
        return ((-(M_B/M_NCO)*(1.0/(Foam::pow((1.0 - xBL_val),2)))*(dxdT))*latenth);
    }
    scalar thermalDiffusivityGas(scalar &T)
    {
        // thermalDiffusivityGas - return thermal diffusivity of gas as a function of T
        // @param = double &T - Temperature
        // constants
        scalar a = 7.860e-10;
        scalar b = 0.000001038;
        scalar c = -0.0001518;

        return (a*T*T + b*T + c);
    }
    scalar npentaneThermalconductivityHighrho(scalar &rho)
    {
        return ((scalar(8.7006e-8)*Foam::pow(rho,2) +
                 scalar(8.4674e-5)*rho + scalar(1.16e-2)));
    }
    scalar npentaneThermalconductivityLowrho(scalar &rho)
    {
        return (scalar(9.3738e-6)*Foam::pow(rho,2) -
                scalar(7.3511e-4)*rho + scalar(2.956e-2));
    }

protected:
    const fvMesh& mesh;
    Time &runTime;
    volScalarField &alpha2;
    volScalarField &TS;

};


}
#endif