/*
 * Copyright (C) 2017 Ruslan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package thomsonsource;

import electronbunch.AbstractElectronBunch;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;
import laserpulse.AbstractLaserPulse;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.RombergIntegrator;
import org.la4j.Vector;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.la4j.Matrix;
import org.la4j.Vectors;
import static thomsonsource.AbstractThomsonSource.MAXIMAL_NUMBER_OF_EVALUATIONS;

/**
 * The main class containing all physics of LEXG in non-linear case
 *
 * @version 1.1
 * @author Ruslan Feshchenko
 */
public final class NonLinearThomsonSource extends AbstractThomsonSource {

    /**
     * Constructor
     *
     * @param l
     * @param b
     * @param n - non-linear order number
     */
    /**
     * The non-linear order number
     */
    private int ordernumber = 1;

    /**
     * The saturating laser intensity
     */
    private double sIntensity = 1;

    /**
     * The list of functional objects for numerical integration
     */
    private List<Function<Double[], Double>> funcarray;

    /**
     * The saturating laser intensity
     *
     * @param l - laser pulse object
     * @param b - electron bunch object
     * @param ornum - the nonlinear order number
     */
    public NonLinearThomsonSource(AbstractLaserPulse l, AbstractElectronBunch b, int ornum) {
        super(l, b);
        this.ordernumber = ornum;
        this.funcarray = generateFunctionList();
        setsIntensity();
        calculateLinearTotalFlux();
        calculateGeometricFactor();
    }

    @Override
    public Object clone() throws CloneNotSupportedException {
        Object ts = super.clone();
        ((NonLinearThomsonSource) ts).funcarray = ((NonLinearThomsonSource) ts).generateFunctionList();
        return ts;
    }

    @Override
    public double directionFrequencyFluxNoSpread(Vector n, Vector v, Vector r, double e) {
        double intensity;
        //If vector r is null then use average intensity
        if (r == null) {
            intensity = lp.getAverageIntensity();
        } else {
            intensity = lp.getIntensity(r);
        }
        double gamma = calculateGamma(n, v, e, intensity);
        double der = calculateGammaDerivative(n, v, e, intensity);
        if (gamma == 0) {
            //returning zero if gamma is zero
            return 0;
        } else {
            return directionFluxBasic(n, v, e, gamma, intensity)
                    * eb.gammaDistribution(gamma) / der;
        }
    }

    @Override
    public double[] directionFrequencyPolarizationNoSpread(Vector n, Vector v, Vector r, double e) {
        double intensity, factor;
        double[] res = new double[]{0, 0, 0, 0};
        //If vector r is null then use average intensity
        if (r == null) {
            intensity = lp.getAverageIntensity();
        } else {
            intensity = lp.getIntensity(r);
        }
        double gamma = calculateGamma(n, v, e, intensity);
        double der = calculateGammaDerivative(n, v, e, intensity);
        if (gamma == 0) {
            //returning zero if gamma is zero
            return res;
        } else {
            res = directionPolarizationBasic(n, v, e, gamma, intensity);
            factor = eb.gammaDistribution(gamma) / der;
            //Multiplying polarization matrix by gamma distribution
            return multiplyPolarizationParameters(res, factor);
        }
    }

    @Override
    public double directionFrequencyPolarizationNoSpread(Vector n, Vector v, Vector r, double e, int index) {
        double intensity;
        //If vector r is null then use average intensity
        if (r == null) {
            intensity = lp.getAverageIntensity();
        } else {
            intensity = lp.getIntensity(r);
        }
        double gamma = calculateGamma(n, v, e, intensity);
        double der = calculateGammaDerivative(n, v, e, intensity);
        if (gamma == 0) {
            //returning zero if gamma is zero
            return 0;
        } else {
            //Multiplying polarization matrix by gamma distribution and returning the result
            return directionPolarizationBasic(n, v, e, gamma, intensity, index) * eb.gammaDistribution(gamma) / der;
        }
    }

    @Override
    public double directionFlux(Vector n, Vector v) {
        return directionFluxBasic(n, v, directionEnergy(n, v), eb.getGamma(), lp.getAverageIntensity());
    }

    @Override
    public double directionEnergy(Vector n, Vector v) {
        return directionEnergyBasic(n, v, eb.getGamma(), lp.getAverageIntensity());
    }

    /**
     * Setting the non-linear order number
     *
     * @return the n
     */
    public int getOrdernumber() {
        return ordernumber;
    }

    /**
     * Getting the non-linear order number
     *
     * @param n the n to set
     */
    public void setOrdernumber(int n) {
        this.ordernumber = n;
    }

    /**
     * @return the current saturating intensity in W/m^2
     */
    public double getsIntensity() {
        return sIntensity;
    }

    /**
     * Sets saturating intensity for a given photon energy, W/m^2
     */
    public void setsIntensity() {
        double k = lp.getPhotonEnergy() / AbstractLaserPulse.HC;
        this.sIntensity = 1e-5 * AbstractLaserPulse.C * k * k * AS * AS / 8 / Math.PI;
    }

    /**
     * Basic method for calculation of X-ray energy
     *
     * @param n
     * @param v
     * @param gamma
     * @param inten
     * @return
     */
    private double directionEnergyBasic(Vector n, Vector v, double gamma, double inten) {
        double mv, M, pr;
        mv = Math.sqrt(1.0 - 1.0 / gamma / gamma);
        pr = n.innerProduct(v);
        M = inten / sIntensity * (1 + pr) * (1 - mv) / 4;
        return ordernumber * (1 + mv) * lp.getPhotonEnergy() / (1 - pr * mv + M);
    }

    /**
     * Basic method for calculation of X-ray flux
     *
     * @param n
     * @param v
     * @param xenergy
     * @param gamma
     * @param inten
     * @return
     */
    private double directionFluxBasic(Vector n, Vector v, double xenergy, double gamma, double inten) {
        //If gamma is zero (negative expression under the root) then return zero
        if (gamma == 0) {
            return 0;
        }
        //If gamma is not zero then proceed
        double mv, M, pr, gamma2, intratio, coef1, coef2, coef3, result;
        double[] K1 = lp.getKA1();
        double[] K2 = lp.getKA2();
        Vector[] A1 = lp.getA1();
        Vector[] A2 = lp.getA2();
        gamma2 = gamma * gamma;
        mv = Math.sqrt(1.0 - 1.0 / gamma2);//Dimesionaless speed
        pr = n.innerProduct(v);
        intratio = inten / sIntensity;
        //Parameter of non-linearity
        M = intratio * (1 + pr) / 4 * (1 - mv);
        //Coefficients
        coef1 = xenergy / lp.getPhotonEnergy()
                * Math.sqrt(intratio) / (1 + mv) / gamma;
        coef2 = xenergy / lp.getPhotonEnergy() * intratio
                * (1 + pr) / Math.pow(gamma * (1 + mv), 2) / 8;
        coef3 = -getLinearTotalFlux() * ordernumber * 3 / 2 / Math.PI
                / Math.pow((1 - pr * mv) * (1 + M), 2) / gamma2;
        //Checking if the radiation is fully polarized and then just calculating intensity
        if (K1[0] == 0 && K2[0] == 0) {
            result = directionFluxBasicAuxiliary(coef1, coef2, intratio, n, new Vector[]{A1[1], A2[1]});
        } else if (K1[1] == 0 && K2[1] == 0) {
            result = directionFluxBasicAuxiliary(coef1, coef2, intratio, n, new Vector[]{A1[0], A2[0]});
        } else {
            //If not fully polarized then everaging over all possible surpepositions of two independant polarizations
            result = (new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                    RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT,
                    RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT)).integrate(MAXIMAL_NUMBER_OF_EVALUATIONS,
                    new UnivariatePolarizationIntegration(coef1, coef2, intratio, n), -Math.PI, Math.PI) / 2 / Math.PI;
        }
        return coef3 * result;
    }

    /**
     * Basic method for calculation of X-ray polarization matrix
     *
     * @param n
     * @param v
     * @param xenergy
     * @param gamma
     * @param inten
     * @return
     */
    private double[] directionPolarizationBasic(Vector n, Vector v, double xenergy, double gamma, double inten) {
        //If gamma is zero (negative expression under the root) then return zero
        double[] res = new double[]{0, 0, 0, 0};
        if (gamma == 0) {
            return res;
        }
        //If gamma is not zero then proceed
        double mv, M, pr, gamma2, intratio, coef1, coef2, coef3;
        double[] result = new double[AbstractThomsonSource.NUMBER_OF_POL_PARAM];
        double[] K1 = lp.getKA1();
        double[] K2 = lp.getKA2();
        Vector[] A1 = lp.getA1();
        Vector[] A2 = lp.getA2();
        gamma2 = gamma * gamma;
        mv = Math.sqrt(1.0 - 1.0 / gamma2);//Dimesionaless speed
        pr = n.innerProduct(v);
        intratio = inten / sIntensity;
        //Parameter of non-linearity
        M = intratio * (1 + pr) / 4 * (1 - mv);
        //Coefficients
        coef1 = xenergy / lp.getPhotonEnergy()
                * Math.sqrt(intratio) / (1 + mv) / gamma;
        coef2 = xenergy / lp.getPhotonEnergy() * intratio
                * (1 + pr) / Math.pow(gamma * (1 + mv), 2) / 8;
        coef3 = -getLinearTotalFlux() * ordernumber * 3 / 2 / Math.PI
                / Math.pow((1 - pr * mv) * (1 + M), 2) / gamma2;
        //Checking if the radiation is fully poolarized and then just calculating intensity
        if (K1[0] == 0 && K2[0] == 0) {
            result = directionPolarizationBasicAuxiliary(coef1, coef2, intratio,
                    n, new Vector[]{A1[1], A2[1]}, gamma, v);
        } else if (K1[1] == 0 && K2[1] == 0) {
            result = directionPolarizationBasicAuxiliary(coef1, coef2, intratio,
                    n, new Vector[]{A1[0], A2[0]}, gamma, v);
        } else {
            //If not fully polarized then everaging over all possible surpepositions of two independant polarizations
            for (int s = 0; s < AbstractThomsonSource.NUMBER_OF_POL_PARAM; s++) {
                result[s] = (new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                        RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT,
                        RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT)).integrate(MAXIMAL_NUMBER_OF_EVALUATIONS,
                        new UnivariatePolarizationIntegrationPolarization(coef1, coef2,
                                intratio, gamma, n, v, s), -Math.PI, Math.PI) / 2 / Math.PI;
            }
        }
        result = multiplyPolarizationParameters(result, coef3);
        return result;
    }

    /**
     * Basic method for calculation of an element of X-ray polarization matrix
     *
     * @param n
     * @param v
     * @param xenergy
     * @param gamma
     * @param inten
     * @param index
     * @return
     */
    private double directionPolarizationBasic(Vector n, Vector v, double xenergy, double gamma, double inten, int index) {
        //If gamma is zero (negative expression under the root) then return zero
        if (gamma == 0) {
            return 0;
        }
        //If gamma is not zero then proceed
        double mv, M, pr, gamma2, intratio, coef1, coef2, coef3;
        double result;
        double[] K1 = lp.getKA1();
        double[] K2 = lp.getKA2();
        Vector[] A1 = lp.getA1();
        Vector[] A2 = lp.getA2();
        gamma2 = gamma * gamma;
        mv = Math.sqrt(1.0 - 1.0 / gamma2);//Dimesionaless speed
        pr = n.innerProduct(v);
        intratio = inten / sIntensity;
        //Parameter of non-linearity
        M = intratio * (1 + pr) / 4 * (1 - mv);
        //Coefficients
        coef1 = xenergy / lp.getPhotonEnergy()
                * Math.sqrt(intratio) / (1 + mv) / gamma;
        coef2 = xenergy / lp.getPhotonEnergy() * intratio
                * (1 + pr) / Math.pow(gamma * (1 + mv), 2) / 8;
        coef3 = -getLinearTotalFlux() * ordernumber * 3 / 2 / Math.PI
                / Math.pow((1 - pr * mv) * (1 + M), 2) / gamma2;
        //Checking if the radiation is fully poolarized and then just calculating intensity
        if (K1[0] == 0 && K2[0] == 0) {
            result = directionPolarizationBasicAuxiliary(coef1, coef2, intratio,
                    n, new Vector[]{A1[1], A2[1]}, gamma, v)[index];
        } else if (K1[1] == 0 && K2[1] == 0) {
            result = directionPolarizationBasicAuxiliary(coef1, coef2, intratio,
                    n, new Vector[]{A1[0], A2[0]}, gamma, v)[index];
        } else {
            //If not fully polarized then averaging over all possible surpepositions of two independant polarizations
            result = coef3 * (new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                    RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT,
                    RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT)).integrate(MAXIMAL_NUMBER_OF_EVALUATIONS,
                    new UnivariatePolarizationIntegrationPolarization(coef1, coef2,
                            intratio, gamma, n, v, index), -Math.PI, Math.PI) / 2 / Math.PI;
        }
        return result;
    }

    /**
     * Basic method for calculation of X-ray flux with a given polarization
     *
     * @param cf1
     * @param cf2
     * @param cf3
     * @param intratio
     * @param n
     * @param B
     * @return
     */
    private double directionFluxBasicAuxiliary(double cf1, double cf2, double intratio,
            Vector n, Vector[] B) {
        double K1, K2, a1, a2, a3, result = 0;
        double[] f = new double[8]; //An array for integrals
        K1 = B[0].innerProduct(B[0]);
        K2 = B[1].innerProduct(B[1]);
        //If the fieled is zero, rerurn zero
        if (K1 != 0 || K2 != 0) {
            a1 = cf1 * n.innerProduct(B[0]);
            a2 = cf1 * n.innerProduct(B[1]);
            a3 = cf2 * (K1 - K2);
            //Calculation of six independent non-linear Fourier integrals necessary for cross-section determination
            for (int t = 2; t < 8; t++) {
                //Creating a Romberg integrator and a UnivariateFunction object and then integrating
                try {
                    f[t] = (new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT,
                            RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT)).integrate(MAXIMAL_NUMBER_OF_EVALUATIONS,
                            new UnivariateFourierHarmonics(a1, a2, a3, t), -Math.PI, Math.PI) / 2 / Math.PI - 1;
                } catch (TooManyEvaluationsException ex) {
                    f[t] = 0;
                }
            }
            //Calculating f_0 using the relation between integrals f_i
            f[0] = -(a1 * f[2] + a2 * f[4] + 2 * a3 * f[6]) / ordernumber;
            f[1] = -(a1 * f[3] + a2 * f[5] + 2 * a3 * f[7]) / ordernumber;
            //Calculating the total flux
            result = ((f[0] * f[0] + f[1] * f[1]) * (1.0 / intratio + (K1 + K2) / 2)
                    - (f[2] * f[2] + f[3] * f[3]) * K1 - (f[4] * f[4] + f[5] * f[5]) * K2
                    + (f[6] * f[0] + f[7] * f[1]) * (K1 - K2) / 2);
        }
        return result;
    }

    /**
     * Basic method for calculation of X-ray polarization matrix with a given
     * polarization
     *
     * @param cf1
     * @param cf2
     * @param intratio
     * @param n
     * @param B
     * @param gamma
     * @param v
     * @return
     */
    private double[] directionPolarizationBasicAuxiliary(double cf1, double cf2, double intratio,
            Vector n, Vector[] B, double gamma, Vector v) {
        double K1, K2, a1, a2, a3, gamma2, mv, pr, sq, sqratio;
        //Auxialiry parameters
        gamma2 = gamma * gamma;
        mv = Math.sqrt(1.0 - 1.0 / gamma2);//Dimesionaless speed
        pr = n.innerProduct(v);
        sq = Math.sqrt(mv * mv - pr * pr);
        sqratio = Math.sqrt(intratio);
        //Arrays for real and imagenary parts of polarization vectors
        double[] pol1 = new double[2];
        double[] pol2 = new double[2];
        //Ortogonal unity vectors
        Vector e1 = crossProduct3D(n, v).divide(sq);
        Vector e2 = crossProduct3D(e1, n);
        double[] result = new double[AbstractThomsonSource.NUMBER_OF_POL_PARAM];
        //A vector for integral coefficients
        double[] f = new double[8]; //An array for integrals
        K1 = B[0].innerProduct(B[0]);
        K2 = B[1].innerProduct(B[1]);
        //If the fieled is zero, rerurn zero
        if (K1 != 0 || K2 != 0) {
            a1 = cf1 * n.innerProduct(B[0]);
            a2 = cf1 * n.innerProduct(B[1]);
            a3 = cf2 * (K1 - K2);
            //Calculation of six independent non-linear Fourier integrals necessary for cross-section determination
            for (int t = 2; t < 8; t++) {
                //Creating a Romberg integrator and a UnivariateFunction object and then integrating
                try {
                    f[t] = (new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT,
                            RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT)).integrate(MAXIMAL_NUMBER_OF_EVALUATIONS,
                            new UnivariateFourierHarmonics(a1, a2, a3, t), -Math.PI, Math.PI) / 2 / Math.PI - 1;
                } catch (TooManyEvaluationsException ex) {
                    f[t] = 0;
                }
            }
            //Calculating f_0 using the relation between integrals f_i
            f[0] = -(a1 * f[2] + a2 * f[4] + 2 * a3 * f[6]) / ordernumber;
            f[1] = -(a1 * f[3] + a2 * f[5] + 2 * a3 * f[7]) / ordernumber;
            //Calculating orthogonal polarization vectors
            pol1[0] = e1.innerProduct(B[0]) * f[2] + e2.innerProduct(B[1]) * f[4];
            pol1[1] = e1.innerProduct(B[0]) * f[3] + e2.innerProduct(B[1]) * f[5];
            pol2[0] = gamma * Math.sqrt(intratio) * (1 - (K1 + K2) / 4 / intratio / gamma2 / mv / (1 + mv)) * sq * f[0]
                    + (n.innerProduct(B[0]) * f[2] + n.innerProduct(B[1]) * f[4]) * pr / sq
                    - (K1 - K2) / 4 / gamma / Math.sqrt(intratio) / mv / (1 + mv) * f[6] * sq;
            pol2[1] = gamma * sqratio * (1 - (K1 + K2) / 4 / intratio / gamma2 / mv / (1 + mv)) * sq * f[1]
                    + (n.innerProduct(B[0]) * f[3] + n.innerProduct(B[1]) * f[5]) * pr / sq
                    - (K1 - K2) / 4 / gamma / sqratio / mv / (1 + mv) * f[7] * sq;
            //Calculating the elements of the polarization matrix
            result[0] = pol1[0] * pol1[0] + pol1[1] * pol1[1] + pol2[0] * pol2[0] + pol2[1] * pol2[1];
            result[1] = pol1[0] * pol2[0] + pol1[1] * pol2[1];
            result[2] = pol1[0] * pol2[1] - pol1[1] * pol2[0];
            result[3] = pol2[0] * pol2[0] + pol2[1] * pol2[1] - pol1[0] * pol1[0] - pol1[1] * pol1[1];
        }
        return result;
    }

    /**
     * A method calculating required gamma factor as function of X-ray energy,
     * laser intensity and other parameters
     *
     * @param n
     * @param v
     * @param e
     * @param inten
     * @return
     */
    private double calculateGamma(Vector n, Vector v, double e, double inten) {
        double rho, pr, fqration, coef;
        pr = n.innerProduct(v);
        rho = inten / sIntensity * (1 + pr) / 4;
        fqration = e / ordernumber / lp.getPhotonEnergy();
        coef = 2 - fqration * (1 - pr);
        if (coef <= 0) {
            //Returning zero if the expression under the root is not positive
            return 0;
        } else {
            return (fqration * (pr + rho) + 1)
                    / Math.sqrt(fqration * (1 + 2 * rho + pr) * coef);
        }
    }

    /**
     * A method calculating derivative of X-ray energy by gamma factor
     * normalized by X-ray energy
     *
     * @param n
     * @param v
     * @param e
     * @param inten
     * @return
     */
    private double calculateGammaDerivative(Vector n, Vector v, double e, double inten) {
        double rho, pr, fqration, coef;
        pr = n.innerProduct(v);
        rho = inten / sIntensity * (1 + pr) / 4;
        fqration = e / ordernumber / lp.getPhotonEnergy();
        coef = 2 - fqration * (1 - pr);
        if (coef <= 0) {
            //Returning unit if the expression under the root is not positive
            return 1;
        } else {
            return Math.sqrt((1 + 2 * rho + pr) * (2 - fqration * (1 - pr)) * fqration)
                    * (2 - fqration * (1 - pr)) / (fqration * (1 + rho) - 1);
        }
    }

    @Override
    public double directionFrequencyBrillianceNoSpread(Vector r0, Vector n, Vector v, double e) throws InterruptedException {
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        //Creating an anonymous class for the integrand
        UnivariateFunction func = new UnivariateFunction() {
            @Override
            public double value(double x) {
                if (n.get(0) + n.get(1) + n.get(2) == 0) {
                    throw new LocalException(x);
                }
                try {
                    return directionFrequencyVolumeFluxNoSpread(r0.add(n.multiply(x)), n, v, e);
                } catch (InterruptedException ex) {
                    return 0;
                }
            }
        };
        //Integrating over a line
        try {
            //If interrupted, throw InterruptedException
            if (Thread.currentThread().isInterrupted()) {
                throw new InterruptedException("directionFrequencyBrillianceNoSpread!");
            }
            return integrator.integrate(AbstractThomsonSource.MAXIMAL_NUMBER_OF_EVALUATIONS, func,
                    r0.fold(Vectors.mkEuclideanNormAccumulator()) - 3 * eb.getLength(), r0.fold(Vectors.mkEuclideanNormAccumulator()) + 3 * eb.getLength());
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
    }

    @Override
    public double directionFrequencyBrillianceSpread(Vector r0, Vector n, Vector v, double e) throws InterruptedException {
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        //Creating an anonymous class for the integrand
        UnivariateFunction func = new UnivariateFunction() {
            @Override
            public double value(double x) {
                if (n.get(0) + n.get(1) + n.get(2) == 0) {
                    throw new LocalException(x);
                }
                try {
                    return directionFrequencyVolumeFluxSpread(r0.add(n.multiply(x)), n, v, e);
                } catch (InterruptedException ex) {
                    Thread.currentThread().interrupt();
                    return 0;
                }
            }
        };
        try {
            //If interrupted, throw InterruptedException
            if (Thread.currentThread().isInterrupted()) {
                throw new InterruptedException("directionFrequencyBrillianceSpread!");
            }
            return integrator.integrate(AbstractThomsonSource.MAXIMAL_NUMBER_OF_EVALUATIONS, func, r0.fold(Vectors.mkEuclideanNormAccumulator()) - 3 * eb.getLength(),
                    r0.fold(Vectors.mkEuclideanNormAccumulator()) + 3 * eb.getLength());
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
    }

    @Override
    public double[] directionFrequencyBrilliancePolarizationNoSpread(Vector r0, Vector n, Vector v, double e) throws InterruptedException {
        double[] array = new double[NUMBER_OF_POL_PARAM];
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateFunction func;
        int[] ia = new int[1];
        for (int i = 0; i < NUMBER_OF_POL_PARAM; i++) {
            //Creating an anonymous class for the integrand
            ia[0] = i;
            func = new UnivariateFunction() {
                @Override
                public double value(double x) {
                    if (n.get(0) + n.get(1) + n.get(2) == 0) {
                        throw new LocalException(x);
                    }
                    try {
                        return directionFrequencyVolumePolarizationNoSpread(r0.add(n.multiply(x)), n, v, e, ia[0]);
                    } catch (InterruptedException ex) {
                        Thread.currentThread().interrupt();
                        return 0;
                    }
                }
            };
            try {
                //If interrupted, throw InterruptedException
                if (Thread.currentThread().isInterrupted()) {
                    throw new InterruptedException("directionFrequencyBrilliancePolarizationNoSpread!");
                }
                array[i] = integrator.integrate(AbstractThomsonSource.MAXIMAL_NUMBER_OF_EVALUATIONS, func,
                        r0.fold(Vectors.mkEuclideanNormAccumulator()) - 3 * eb.getLength(), r0.fold(Vectors.mkEuclideanNormAccumulator()) + 3 * eb.getLength());
            } catch (TooManyEvaluationsException ex) {
                array[i] = 0;
            }
        }
        return array;
    }

    @Override
    public double[] directionFrequencyBrilliancePolarizationSpread(Vector r0, Vector n, Vector v, double e) throws InterruptedException {
        double[] array = new double[NUMBER_OF_POL_PARAM];
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateFunction func;
        int[] ia = new int[1];
        for (int i = 0; i < NUMBER_OF_POL_PARAM; i++) {
            //Creating an anonymous class for the integrand
            ia[0] = i;
            func = new UnivariateFunction() {
                @Override
                public double value(double x) {
                    if (n.get(0) + n.get(1) + n.get(2) == 0) {
                        throw new LocalException(x);
                    }
                    try {
                        return directionFrequencyVolumePolarizationSpread(r0.add(n.multiply(x)), n, v, e, ia[0]);
                    } catch (InterruptedException ex) {
                        Thread.currentThread().interrupt();
                        return 0;
                    }
                }
            };
            try {
                //If interrupted, throw InterruptedException
                if (Thread.currentThread().isInterrupted()) {
                    throw new InterruptedException("directionFrequencyBrilliancePolarizationSpread!");
                }
                array[i] = integrator.integrate(AbstractThomsonSource.MAXIMAL_NUMBER_OF_EVALUATIONS, func,
                        r0.fold(Vectors.mkEuclideanNormAccumulator()) - 3 * eb.getLength(), r0.fold(Vectors.mkEuclideanNormAccumulator()) + 3 * eb.getLength());
            } catch (TooManyEvaluationsException ex) {
                array[i] = 0;
            }
        }
        return array;
    }

    @Override
    public double directionFrequencyBrilliancePolarizationNoSpread(Vector r0, Vector n, Vector v, double e, int index) throws InterruptedException {
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateFunction func;
        //Creating an anonymous class for the integrand
        func = new UnivariateFunction() {
            @Override
            public double value(double x) {
                if (n.get(0) + n.get(1) + n.get(2) == 0) {
                    throw new LocalException(x);
                }
                try {
                    return directionFrequencyVolumePolarizationNoSpread(r0.add(n.multiply(x)), n, v, e, index);
                } catch (InterruptedException ex) {
                    Thread.currentThread().interrupt();
                    return 0;
                }
            }
        };
        try {
            //If interrupted, throw InterruptedException
            if (Thread.currentThread().isInterrupted()) {
                throw new InterruptedException("directionFrequencyBrilliancePolarizationNoSpread!");
            }
            return integrator.integrate(AbstractThomsonSource.MAXIMAL_NUMBER_OF_EVALUATIONS, func,
                    r0.fold(Vectors.mkEuclideanNormAccumulator()) - 3 * eb.getLength(), r0.fold(Vectors.mkEuclideanNormAccumulator()) + 3 * eb.getLength());
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
    }

    @Override
    public double directionFrequencyBrilliancePolarizationSpread(Vector r0, Vector n, Vector v, double e, int index) throws InterruptedException {
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateFunction func;
        //Creating an anonymous class for the integrand
        func = new UnivariateFunction() {
            @Override
            public double value(double x) {
                if (n.get(0) + n.get(1) + n.get(2) == 0) {
                    throw new LocalException(x);
                }
                try {
                    return directionFrequencyVolumePolarizationSpread(r0.add(n.multiply(x)), n, v, e, index);
                } catch (InterruptedException ex) {
                    Thread.currentThread().interrupt();
                    return 0;
                }
            }
        };
        try {
            //If interrupted, throw InterruptedException
            if (Thread.currentThread().isInterrupted()) {
                throw new InterruptedException("directionFrequencyBrilliancePolarizationSpread!");
            }
            return integrator.integrate(AbstractThomsonSource.MAXIMAL_NUMBER_OF_EVALUATIONS, func,
                    r0.fold(Vectors.mkEuclideanNormAccumulator()) - 3 * eb.getLength(), r0.fold(Vectors.mkEuclideanNormAccumulator()) + 3 * eb.getLength());
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
    }

    @Override
    public double directionFrequencyVolumeFluxNoSpread(Vector r, Vector n, Vector v, double e) throws InterruptedException {
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        Vector re = r.copy();
        //Transforming coordinates between laser and electron beam frames
        Matrix T = getTransform(v, lp.getDirection());
        Vector rph = T.multiply(r);
        //Creating an anonymous class for the integrand
        UnivariateFunction func = new UnivariateFunction() {
            @Override
            public double value(double x) {
                //Coordinate transformation between laser and electron beams frames
                re.set(2, re.get(2) - x);
                rph.set(2, rph.get(2) - x);
                return directionFrequencyFluxNoSpread(r, n, v, e) * eb.lSpatialDistribution(re)
                        * lp.lSpatialDistribution(rph);
            }
        };
        try {
            //If interrupted, throw InterruptedException
            if (Thread.currentThread().isInterrupted()) {
                throw new InterruptedException("directionFrequencyVolumeFluxNoSpread!");
            }
            return integrator.integrate(AbstractThomsonSource.MAXIMAL_NUMBER_OF_EVALUATIONS, func, -INT_RANGE * eb.getLength(),
                    INT_RANGE * eb.getLength()) * lp.tSpatialDistribution(rph) * eb.tSpatialDistribution(re);
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
    }

    @Override
    public double[] directionFrequencyVolumePolarizationNoSpread(Vector r, Vector n, Vector v, double e) throws InterruptedException {
        double[] stocks = new double[AbstractThomsonSource.NUMBER_OF_POL_PARAM];
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        Vector re = r.copy();
        //Transforming coordinates between laser and electron beam frames
        Matrix T = getTransform(v, lp.getDirection());
        Vector rph = T.multiply(r);
        int[] ia = new int[1];
        //Intterating over Stocks parameters multiplied by intensity
        for (int i = 0; i < AbstractThomsonSource.NUMBER_OF_POL_PARAM; i++) {
            //Creating an anonymous class for the integrand
            ia[0] = i;
            UnivariateFunction func = new UnivariateFunction() {
                @Override
                public double value(double x) {
                    //Coordinate transformation between laser and electron beams frames
                    re.set(2, re.get(2) - x);
                    rph.set(2, rph.get(2) - x);
                    return directionFrequencyPolarizationNoSpread(r, n, v, e, ia[0]) * eb.lSpatialDistribution(re)
                            * lp.lSpatialDistribution(rph);
                }
            };
            try {
                //If interrupted, throw InterruptedException
                if (Thread.currentThread().isInterrupted()) {
                    throw new InterruptedException("directionFrequencyVolumePolarizationNoSpread!");
                }
                stocks[i] = integrator.integrate(AbstractThomsonSource.MAXIMAL_NUMBER_OF_EVALUATIONS, func, -INT_RANGE * eb.getLength(),
                        INT_RANGE * eb.getLength()) * lp.tSpatialDistribution(rph) * eb.tSpatialDistribution(re);
            } catch (TooManyEvaluationsException ex) {
                stocks[i] = 0;
            }
        }
        return stocks;
    }

    @Override
    public double directionFrequencyVolumePolarizationNoSpread(Vector r, Vector n, Vector v, double e, int index) throws InterruptedException {
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        Vector re = r.copy();
        //Transforming coordinates between laser and electron beam frames
        Matrix T = getTransform(v, lp.getDirection());
        Vector rph = T.multiply(r);
        //Creating an anonymous class for the integrand
        UnivariateFunction func = new UnivariateFunction() {
            @Override
            public double value(double x) {
                //Coordinate transformation between laser and electron beams frames
                re.set(2, re.get(2) - x);
                rph.set(2, rph.get(2) - x);
                return directionFrequencyPolarizationNoSpread(r, n, v, e, index) * eb.lSpatialDistribution(re)
                        * lp.lSpatialDistribution(rph);
            }
        };
        try {
            //If interrupted, throw InterruptedException
            if (Thread.currentThread().isInterrupted()) {
                throw new InterruptedException("directionFrequencyVolumePolarizationNoSpread!");
            }
            return integrator.integrate(AbstractThomsonSource.MAXIMAL_NUMBER_OF_EVALUATIONS, func, -INT_RANGE * eb.getLength(),
                    INT_RANGE * eb.getLength()) * lp.tSpatialDistribution(rph) * eb.tSpatialDistribution(re);
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
    }

    @Override
    public double directionFrequencyVolumeFluxSpread(Vector r, Vector n, Vector v, double e) throws InterruptedException {
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        Vector re = r.copy();
        //Transforming coordinates between laser and electron beam frames
        Matrix T = getTransform(v, lp.getDirection());
        Vector rph = T.multiply(r);
        //Creating an anonymous class for the integrand
        UnivariateFunction func = new UnivariateFunction() {
            @Override
            public double value(double x) {
                //Coordinate transformation between laser and electron beams frames
                re.set(2, re.get(2) - x);
                rph.set(2, rph.get(2) - x);
                try {
                    return directionFrequencyFluxSpread(r, n, v, e) * eb.lSpatialDistribution(re)
                            * lp.lSpatialDistribution(rph);
                } catch (InterruptedException ex) {
                    Thread.currentThread().interrupt();
                    return 0;
                }
            }
        };
        try {
            //If interrupted, throw InterruptedException
            if (Thread.currentThread().isInterrupted()) {
                throw new InterruptedException("Interruption in directionFrequencyFluxSpread!");
            }
            return integrator.integrate(AbstractThomsonSource.MAXIMAL_NUMBER_OF_EVALUATIONS, func, -INT_RANGE * eb.getLength(),
                    INT_RANGE * eb.getLength()) * lp.tSpatialDistribution(rph) * eb.tSpatialDistribution(re);
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
    }

    @Override
    public double[] directionFrequencyVolumePolarizationSpread(Vector r, Vector n, Vector v, double e) throws InterruptedException {
        double[] stocks = new double[AbstractThomsonSource.NUMBER_OF_POL_PARAM];
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        Vector re = r.copy();
        //Transforming coordinates between laser and electron beam frames
        Matrix T = getTransform(v, lp.getDirection());
        Vector rph = T.multiply(r);
        int[] ia = new int[1];
        //Intterating over Stocks parameters multiplied by intensity
        for (int i = 0; i < AbstractThomsonSource.NUMBER_OF_POL_PARAM; i++) {
            //Creating an anonymous class for the integrand
            ia[0] = i;
            UnivariateFunction func = new UnivariateFunction() {
                @Override
                public double value(double x) {
                    //Coordinate transformation between laser and electron beams frames
                    re.set(2, re.get(2) - x);
                    rph.set(2, rph.get(2) - x);
                    try {
                        return directionFrequencyPolarizationSpread(r, n, v, e, ia[0]) * eb.lSpatialDistribution(re)
                                * lp.lSpatialDistribution(rph);
                    } catch (InterruptedException ex) {
                        Thread.currentThread().interrupt();
                        return 0;
                    }
                }
            };
            try {
                stocks[i] = integrator.integrate(AbstractThomsonSource.MAXIMAL_NUMBER_OF_EVALUATIONS, func, -INT_RANGE * eb.getLength(),
                        INT_RANGE * eb.getLength()) * lp.tSpatialDistribution(rph) * eb.tSpatialDistribution(re);
            } catch (TooManyEvaluationsException ex) {
                stocks[i] = 0;
            }
        }
        return stocks;
    }

    @Override
    public double directionFrequencyVolumePolarizationSpread(Vector r, Vector n, Vector v, double e, int index) throws InterruptedException {
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        Vector re = r.copy();
        //Transforming coordinates between laser and electron beam frames
        Matrix T = getTransform(v, lp.getDirection());
        Vector rph = T.multiply(r);
        //Creating an anonymous class for the integrand
        UnivariateFunction func = new UnivariateFunction() {
            @Override
            public double value(double x) {
                //Coordinate transformation between laser and electron beams frames
                re.set(2, re.get(2) - x);
                rph.set(2, rph.get(2) - x);
                try {
                    return directionFrequencyPolarizationSpread(r, n, v, e, index) * eb.lSpatialDistribution(re)
                            * lp.lSpatialDistribution(rph);
                } catch (InterruptedException ex) {
                    Thread.currentThread().interrupt();
                    return 0;
                }
            }
        };
        try {
            return integrator.integrate(AbstractThomsonSource.MAXIMAL_NUMBER_OF_EVALUATIONS, func, -INT_RANGE * eb.getLength(),
                    INT_RANGE * eb.getLength()) * lp.tSpatialDistribution(rph) * eb.tSpatialDistribution(re);
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
    }

    /**
     * An auxiliary class for eight Romberg integrators for Fourier polarization
     * harmonics calculations
     */
    private class UnivariateFourierHarmonics implements UnivariateFunction {

        private final double a1, a2, a3;
        private final Function<Double[], Double> fn;

        public UnivariateFourierHarmonics(double a1, double a2, double a3, int ind) {
            this.a1 = a1;
            this.a2 = a2;
            this.a3 = a3;
            this.fn = funcarray.get(ind);
        }

        @Override
        public double value(double tau) {
            return fn.apply(new Double[]{tau, a1, a2, a3});
        }
    }

    /**
     * An auxiliary class for integration over all possible polarization
     * superpositions for flux calculations
     */
    private class UnivariatePolarizationIntegration implements UnivariateFunction {

        private final double cf1, cf2, intratio;
        private final Vector n;

        public UnivariatePolarizationIntegration(double cf1, double cf2, double intratio, Vector n) {
            this.cf1 = cf1;
            this.cf2 = cf2;
            this.intratio = intratio;
            this.n = n;
        }

        @Override
        public double value(double phase) {
            //Getting phase weighted superposition of polarizations
            Vector[] B = lp.getPhaseWeightedPolarizationVectors2(phase);
            //Returning the calculated intensity
            return directionFluxBasicAuxiliary(cf1, cf2, intratio, n, B);
        }
    }

    /**
     * An auxiliary class for integration over all possible polarization
     * superpositions for polarization matrix calculations
     */
    private class UnivariatePolarizationIntegrationPolarization implements UnivariateFunction {

        private final double cf1, cf2, intratio, gamma;
        private final Vector n, v;
        private final int index;

        public UnivariatePolarizationIntegrationPolarization(double cf1, double cf2,
                double intratio, double gamma, Vector n, Vector v, int index) {
            this.cf1 = cf1;
            this.cf2 = cf2;
            this.intratio = intratio;
            this.n = n;
            this.index = index;
            this.v = v;
            this.gamma = gamma;
        }

        @Override
        public double value(double phase) {
            //Getting phase weighted superposition of polarizations
            Vector[] B = lp.getPhaseWeightedPolarizationVectors2(phase);
            //Returning the calculated intensity
            return directionPolarizationBasicAuxiliary(cf1, cf2, intratio, n, B, gamma, v)[index];
        }
    }

    /**
     * Generates a list of functions for parameters f_i
     *
     * @return
     */
    private List<Function<Double[], Double>> generateFunctionList() {
        List<Function<Double[], Double>> fa = new ArrayList<>();
        //Inoitializeing the array of functions used to calculate non-linear amplitudes
        fa.add(arg -> {
            return 1 + Math.cos(ordernumber * arg[0] + arg[1] * Math.sin(arg[0]) - arg[2] * Math.cos(arg[0]) + arg[3] * Math.sin(2 * arg[0]));
        });
        fa.add(arg -> {
            return 1 + Math.sin(ordernumber * arg[0] + arg[1] * Math.sin(arg[0]) - arg[2] * Math.cos(arg[0]) + arg[3] * Math.sin(2 * arg[0]));
        });
        fa.add(arg -> {
            double cs = Math.cos(arg[0]);
            return 1 + cs * Math.cos(ordernumber * arg[0] + arg[1] * Math.sin(arg[0]) - arg[2] * cs + arg[3] * Math.sin(2 * arg[0]));
        });
        fa.add(arg -> {
            double cs = Math.cos(arg[0]);
            return 1 + cs * Math.sin(ordernumber * arg[0] + arg[1] * Math.sin(arg[0]) - arg[2] * cs + arg[3] * Math.sin(2 * arg[0]));
        });
        fa.add(arg -> {
            double sn = Math.sin(arg[0]);
            return 1 + sn * Math.cos(ordernumber * arg[0] + arg[1] * sn - arg[2] * Math.cos(arg[0]) + arg[3] * Math.sin(2 * arg[0]));
        });
        fa.add(arg -> {
            double sn = Math.sin(arg[0]);
            return 1 + sn * Math.sin(ordernumber * arg[0] + arg[1] * sn - arg[2] * Math.cos(arg[0]) + arg[3] * Math.sin(2 * arg[0]));
        });
        fa.add(arg -> {
            double tau2 = 2 * arg[0];
            return 1 + Math.cos(tau2) * Math.cos(ordernumber * arg[0] + arg[1] * Math.sin(arg[0]) - arg[2] * Math.cos(arg[0]) + arg[3] * Math.sin(2 * arg[0]));
        });
        fa.add(arg -> {
            double tau2 = 2 * arg[0];
            return 1 + Math.cos(tau2) * Math.sin(ordernumber * arg[0] + arg[1] * Math.sin(arg[0]) - arg[2] * Math.cos(arg[0]) + arg[3] * Math.sin(2 * arg[0]));
        });

        return fa;
    }
}
