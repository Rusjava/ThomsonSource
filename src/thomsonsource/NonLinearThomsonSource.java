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
import static thomsonsource.AbstractThomsonSource.MAXIMAL_NUMBER_OF_EVALUATIONS;

/**
 * The main class containing all physics of LEXG in non-linear case
 *
 * @version 1.0
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
    private final List<Function<Double[], Double>> funcarray;

    /**
     * The saturating laser intensity
     *
     * @param l - laser pulse object
     * @param b - electron bunch object
     * @param ordernumber - the nonlinear order number
     */
    public NonLinearThomsonSource(AbstractLaserPulse l, AbstractElectronBunch b, int ordernumber) {
        super(l, b);
        this.ordernumber = ordernumber;
        this.funcarray = new ArrayList<>();
        //Inoitializeing the array of functions used to calculate non-linear amplitudes
        funcarray.add(arg -> {
            return 1 + Math.cos(this.ordernumber * arg[0] + arg[1] * Math.sin(arg[0]) - arg[2] * Math.cos(arg[0]) + arg[3] * Math.sin(2 * arg[0]));
        });
        funcarray.add(arg -> {
            return 1 + Math.sin(this.ordernumber * arg[0] + arg[1] * Math.sin(arg[0]) - arg[2] * Math.cos(arg[0]) + arg[3] * Math.sin(2 * arg[0]));
        });
        funcarray.add(arg -> {
            double cs = Math.cos(arg[0]);
            return 1 + cs * Math.cos(this.ordernumber * arg[0] + arg[1] * Math.sin(arg[0]) - arg[2] * cs + arg[3] * Math.sin(2 * arg[0]));
        });
        funcarray.add(arg -> {
            double cs = Math.cos(arg[0]);
            return 1 + cs * Math.sin(this.ordernumber * arg[0] + arg[1] * Math.sin(arg[0]) - arg[2] * cs + arg[3] * Math.sin(2 * arg[0]));
        });
        funcarray.add(arg -> {
            double sn = Math.sin(arg[0]);
            return 1 + sn * Math.cos(this.ordernumber * arg[0] + arg[1] * sn - arg[2] * Math.cos(arg[0]) + arg[3] * Math.sin(2 * arg[0]));
        });
        funcarray.add(arg -> {
            double sn = Math.sin(arg[0]);
            return 1 + sn * Math.sin(this.ordernumber * arg[0] + arg[1] * sn - arg[2] * Math.cos(arg[0]) + arg[3] * Math.sin(2 * arg[0]));
        });
        funcarray.add(arg -> {
            double tau2 = 2 * arg[0];
            return 1 + Math.cos(tau2) * Math.cos(this.ordernumber * arg[0] + arg[1] * Math.sin(arg[0]) - arg[2] * Math.cos(arg[0]) + arg[3] * Math.sin(2 * arg[0]));
        });
        funcarray.add(arg -> {
            double tau2 = 2 * arg[0];
            return 1 + Math.cos(tau2) * Math.sin(this.ordernumber * arg[0] + arg[1] * Math.sin(arg[0]) - arg[2] * Math.cos(arg[0]) + arg[3] * Math.sin(2 * arg[0]));
        });
        setsIntensity();
        calculateTotalFlux();
        calculateGeometricFactor();
    }

    @Override
    public double directionFrequencyFluxNoSpread(Vector n, Vector v, double e) {
        double K, th;
        th = (1 - n.innerProduct(v)) * 2;
        K = Math.pow((Math.sqrt(e / lp.getPhotonEnergy() / (1 - e * th / lp.getPhotonEnergy() / 4)) - 2 * eb.getGamma()), 2)
                / 4 / Math.pow(eb.getGamma() * eb.getDelgamma(), 2);
        return getTotalFlux() * e * 3.0 / 64 / Math.PI / Math.sqrt(Math.PI) / eb.getDelgamma() / eb.getGamma() / lp.getPhotonEnergy()
                * Math.sqrt(e / lp.getPhotonEnergy()) * (Math.pow((1 - e * th / lp.getPhotonEnergy() / 2), 2) + 1)
                / Math.sqrt(1 - e * th / lp.getPhotonEnergy() / 4) * Math.exp(-K);
    }

    @Override
    public double[] directionFrequencyPolarizationNoSpread(Vector n, Vector v, double e) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public double directionFlux(Vector n, Vector v) {
        double mv, M, pr, gamma2, coef, xenergy, result = 0;
        double[] K1 = lp.getKA1();
        double[] K2 = lp.getKA2();
        Vector[] A1 = lp.getA1();
        Vector[] A2 = lp.getA2();
        double[] f = new double[8]; //An array for integrals
        gamma2 = eb.getGamma() * eb.getGamma();
        mv = Math.sqrt(1.0 - 1.0 / gamma2);//Dimesionaless speed
        pr = n.innerProduct(v);
        xenergy = directionEnergy(n, v);
        coef = xenergy / lp.getPhotonEnergy()
                / Math.sqrt(sIntensity) / (1 + mv) / eb.getGamma();
        //Parameter of non-linearity
        M = lp.getIntensity() / sIntensity * (1 + pr) / 4 / gamma2 / (1 + mv);
        /*
        Cycling over two independent polarizations and adding their intensities
         */
        for (int s = 0; s < 2; s++) {
            if (K1[s] != 0 || K2[s] != 0) {
                double a1 = coef * n.innerProduct(A1[s]);
                double a2 = coef * n.innerProduct(A2[s]);
                double a3 = xenergy / lp.getPhotonEnergy() * (K1[s] - K2[s]) / sIntensity
                        * (1 + pr) / Math.pow(eb.getGamma() * (1 + mv), 2) / 8;
                /*
                Calculation of six independent non-linear Fourier integrals necessary for cross-section determination
                 */
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
                //Returning the total flux
                result += -getTotalFlux() * ordernumber * 3 / 2 / Math.PI
                        / Math.pow((1 - pr * mv) * (1 + M), 2) / gamma2
                        * ((f[0] * f[0] + f[1] * f[1]) * (sIntensity / (K1[s] + K2[s]) + 0.5) - (f[2] * f[2] + f[3] * f[3]) * K1[s] / (K1[s] + K2[s])
                        - (f[4] * f[4] + f[5] * f[5]) * K2[s] / (K1[s] + K2[s]) + (f[6] * f[0] + f[7] * f[1]) * (K1[s] - K2[s]) / (K1[s] + K2[s]) / 2);
            }
        }
        return result;
    }

    @Override
    public double directionEnergy(Vector n, Vector v) {
        return directionEnergyBasic(n, v, eb.getGamma(), lp.getIntensity());
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
     * An auxiliary class for eight Romberg integrators for Fourier harmonics
     * calculations
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
}
