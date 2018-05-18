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
public class NonLinearThomsonSource extends AbstractThomsonSource {

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
     * An array of functional objects
     */
    private final Function<Double[], Double>[] funcarray;

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
        this.funcarray = new Function[8];
        //Inoitializeing the array of functions used to calculate non-linear amplitudes
        funcarray[0] = arg -> {
            return 1 + Math.cos(getOrdernumber() * arg[0] + arg[1] * Math.sin(arg[0]) - arg[2] * Math.cos(arg[0]) + arg[3] * Math.sin(2 * arg[0]));
        };
        funcarray[1] = arg -> {
            return 1 + Math.sin(getOrdernumber() * arg[0] + arg[1] * Math.sin(arg[0]) - arg[2] * Math.cos(arg[0]) + arg[3] * Math.sin(2 * arg[0]));
        };
        funcarray[2] = arg -> {
            double cs = Math.cos(arg[0]);
            return 1 + cs * Math.cos(getOrdernumber() * arg[0] + arg[1] * Math.sin(arg[0]) - arg[2] * cs + arg[3] * Math.sin(2 * arg[0]));
        };
        funcarray[3] = arg -> {
            double cs = Math.cos(arg[0]);
            return 1 + cs * Math.sin(getOrdernumber() * arg[0] + arg[1] * Math.sin(arg[0]) - arg[2] * cs + arg[3] * Math.sin(2 * arg[0]));
        };
        funcarray[4] = arg -> {
            double sn = Math.sin(arg[0]);
            return 1 + sn * Math.cos(getOrdernumber() * arg[0] + arg[1] * sn - arg[2] * Math.cos(arg[0]) + arg[3] * Math.sin(2 * arg[0]));
        };
        funcarray[5] = arg -> {
            double sn = Math.sin(arg[0]);
            return 1 + sn * Math.sin(getOrdernumber() * arg[0] + arg[1] * sn - arg[2] * Math.cos(arg[0]) + arg[3] * Math.sin(2 * arg[0]));
        };
        funcarray[6] = arg -> {
            double tau2 = 2 * arg[0];
            return 1 + Math.cos(tau2) * Math.sin(getOrdernumber() * arg[0] + arg[1] * Math.sin(arg[0]) - arg[2] * Math.cos(arg[0]) + arg[3] * Math.sin(2 * arg[0]));
        };
        funcarray[7] = arg -> {
            double tau2 = 2 * arg[0];
            return 1 + Math.cos(tau2) * Math.sin(getOrdernumber() * arg[0] + arg[1] * Math.sin(arg[0]) - arg[2] * Math.cos(arg[0]) + arg[3] * Math.sin(2 * arg[0]));
        };

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
        double mv, M, pr, gamma2, coef;
        double K1 = lp.getKA1()[0];
        double K2 = lp.getKA2()[0];
        Vector A1 = lp.getA1()[0];
        Vector A2 = lp.getA2()[0];
        double f01, f02, f11, f12, f21, f22, f31, f32;
        gamma2 = eb.getGamma() * eb.getGamma();
        mv = Math.sqrt(1.0 - 1.0 / gamma2);//Dimesionaless speed
        pr = n.innerProduct(v);
        coef = directionEnergy(n, v) / lp.getPhotonEnergy()
                / Math.sqrt(sIntensity) / (1 + mv) / eb.getGamma();
        double a1 = coef * n.innerProduct(A1);
        double a2 = coef * n.innerProduct(A2);
        double a3 = directionEnergy(n, v) / lp.getPhotonEnergy() * (K1 - K2) / sIntensity
                * (1 + pr) / Math.pow(eb.getGamma() * (1 + mv), 2) / 8;

        /*
        Calculating the Fourie harmonics
         */
        // 01 component
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateFunction func = new UnivariateFourierHarmonics01(a1, a2, a3);
        try {
            f01 = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, -Math.PI, Math.PI) / 2 / Math.PI - 1;
        } catch (TooManyEvaluationsException ex) {
            f01 = 0;
        }
        // 02 component
        integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        func = new UnivariateFourierHarmonics02(a1, a2, a3);
        try {
            f02 = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, -Math.PI, Math.PI) / 2 / Math.PI - 1;
        } catch (TooManyEvaluationsException ex) {
            f02 = 0;
        }
        // 11 component
        integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        func = new UnivariateFourierHarmonics11(a1, a2, a3);
        try {
            f11 = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, -Math.PI, Math.PI) / 2 / Math.PI - 1;
        } catch (TooManyEvaluationsException ex) {
            f11 = 0;
        }
        // 12 component
        integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        func = new UnivariateFourierHarmonics12(a1, a2, a3);
        try {
            f12 = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, -Math.PI, Math.PI) / 2 / Math.PI - 1;
        } catch (TooManyEvaluationsException ex) {
            f12 = 0;
        }
        // 21 component
        integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        func = new UnivariateFourierHarmonics21(a1, a2, a3);
        try {
            f21 = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, -Math.PI, Math.PI) / 2 / Math.PI - 1;
        } catch (TooManyEvaluationsException ex) {
            f21 = 0;
        }
        // 22 component
        integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        func = new UnivariateFourierHarmonics22(a1, a2, a3);
        try {
            f22 = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, -Math.PI, Math.PI) / 2 / Math.PI - 1;
        } catch (TooManyEvaluationsException ex) {
            f22 = 0;
        }
        // 31 component
        integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        func = new UnivariateFourierHarmonics31(a1, a2, a3);
        try {
            f31 = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, -Math.PI, Math.PI) / 2 / Math.PI - 1;
        } catch (TooManyEvaluationsException ex) {
            f31 = 0;
        }
        // 32 component
        integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        func = new UnivariateFourierHarmonics32(a1, a2, a3);
        try {
            f32 = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, -Math.PI, Math.PI) / 2 / Math.PI - 1;
        } catch (TooManyEvaluationsException ex) {
            f32 = 0;
        }
        //Parameter of nolinearity
        M = lp.getIntensity() / sIntensity * (1 + pr) / 4 / gamma2 / (1 + mv);
        //Returning the total flux
        return -getTotalFlux() * ordernumber * 3 / 2 / Math.PI
                / Math.pow((1 - pr * mv) * (1 + M), 2) / (K1 + K2) / gamma2
                * ((Math.pow(f01, 2) + Math.pow(f02, 2)) * (sIntensity + (K1 + K2) / 2) - (Math.pow(f11, 2) + Math.pow(f12, 2)) * K1
                - (Math.pow(f21, 2) + Math.pow(f22, 2)) * K2 + (f31 * f01 + f32 * f02) * (K1 - K2) / 2);
    }

    @Override
    public double directionEnergy(Vector n, Vector v) {
        double mv, M, pr, gamma_2;
        gamma_2 = 1.0 / eb.getGamma() / eb.getGamma();
        mv = Math.sqrt(1.0 - gamma_2);
        pr = n.innerProduct(v);
        M = lp.getIntensity() / sIntensity * (1 + pr) * gamma_2 / (1 + mv) / 4;
        return ordernumber * (1 + mv) * lp.getPhotonEnergy() / (1 - pr * mv + M);
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
     * An auxiliary class for Romberg integrator for Fourier harmonics
     * calculations
     */
    private class UnivariateFourierHarmonics implements UnivariateFunction {

        private final double a1, a2, a3;
        private final Function<Double[], Double> fn;

        public UnivariateFourierHarmonics(double a1, double a2, double a3, int ind) {
            this.a1 = a1;
            this.a2 = a2;
            this.a3 = a3;
            this.fn = funcarray[ind];
        }

        @Override
        public double value(double tau) {
            return fn.apply(new Double[]{tau, a1, a2, a3});
        }
    }

    /**
     * An auxiliary class 01 for Romberg integrator for Fourier harmonics
     * calculations
     */
    private class UnivariateFourierHarmonics01 implements UnivariateFunction {

        private final double a1, a2, a3;

        public UnivariateFourierHarmonics01(double a1, double a2, double a3) {
            this.a1 = a1;
            this.a2 = a2;
            this.a3 = a3;
        }

        @Override
        public double value(double tau) {
            return 1 + Math.cos(getOrdernumber() * tau + a1 * Math.sin(tau) - a2 * Math.cos(tau) + a3 * Math.sin(2 * tau));
        }
    }

    /**
     * An auxiliary class 02 for Romberg integrator for Fourier harmonics
     * calculations
     */
    private class UnivariateFourierHarmonics02 implements UnivariateFunction {

        private final double a1, a2, a3;

        public UnivariateFourierHarmonics02(double a1, double a2, double a3) {
            this.a1 = a1;
            this.a2 = a2;
            this.a3 = a3;
        }

        @Override
        public double value(double tau) {
            return 1 + Math.sin(getOrdernumber() * tau + a1 * Math.sin(tau) - a2 * Math.cos(tau) + a3 * Math.sin(2 * tau));
        }
    }

    /**
     * An auxiliary class 11 for Romberg integrator for Fourier harmonics
     * calculations
     */
    private class UnivariateFourierHarmonics11 implements UnivariateFunction {

        private final double a1, a2, a3;

        public UnivariateFourierHarmonics11(double a1, double a2, double a3) {
            this.a1 = a1;
            this.a2 = a2;
            this.a3 = a3;
        }

        @Override
        public double value(double tau) {
            double cs = Math.cos(tau);
            return 1 + cs * Math.cos(getOrdernumber() * tau + a1 * Math.sin(tau) - a2 * cs + a3 * Math.sin(2 * tau));
        }
    }

    /**
     * An auxiliary class 12 for Romberg integrator for Fourier harmonics
     * calculations
     */
    private class UnivariateFourierHarmonics12 implements UnivariateFunction {

        private final double a1, a2, a3;

        public UnivariateFourierHarmonics12(double a1, double a2, double a3) {
            this.a1 = a1;
            this.a2 = a2;
            this.a3 = a3;
        }

        @Override
        public double value(double tau) {
            double cs = Math.cos(tau);
            return 1 + cs * Math.sin(getOrdernumber() * tau + a1 * Math.sin(tau) - a2 * cs + a3 * Math.sin(2 * tau));
        }
    }

    /**
     * An auxiliary class 21 for Romberg integrator for Fourier harmonics
     * calculations
     */
    private class UnivariateFourierHarmonics21 implements UnivariateFunction {

        private final double a1, a2, a3;

        public UnivariateFourierHarmonics21(double a1, double a2, double a3) {
            this.a1 = a1;
            this.a2 = a2;
            this.a3 = a3;
        }

        @Override
        public double value(double tau) {
            double sn = Math.sin(tau);
            return 1 + sn * Math.cos(getOrdernumber() * tau + a1 * sn - a2 * Math.cos(tau) + a3 * Math.sin(2 * tau));
        }
    }

    /**
     * An auxiliary class 22 for Romberg integrator for Fourier harmonics
     * calculations
     */
    private class UnivariateFourierHarmonics22 implements UnivariateFunction {

        private final double a1, a2, a3;

        public UnivariateFourierHarmonics22(double a1, double a2, double a3) {
            this.a1 = a1;
            this.a2 = a2;
            this.a3 = a3;
        }

        @Override
        public double value(double tau) {
            double sn = Math.sin(tau);
            return 1 + sn * Math.sin(getOrdernumber() * tau + a1 * sn - a2 * Math.cos(tau) + a3 * Math.sin(2 * tau));
        }
    }

    /**
     * An auxiliary class 31 for Romberg integrator for Fourier harmonics
     * calculations
     */
    private class UnivariateFourierHarmonics31 implements UnivariateFunction {

        private final double a1, a2, a3;

        public UnivariateFourierHarmonics31(double a1, double a2, double a3) {
            this.a1 = a1;
            this.a2 = a2;
            this.a3 = a3;
        }

        @Override
        public double value(double tau) {
            double tau2 = 2 * tau;
            return 1 + Math.cos(tau2) * Math.cos(getOrdernumber() * tau + a1 * Math.sin(tau) - a2 * Math.cos(tau) + a3 * Math.sin(tau2));
        }
    }

    /**
     * An auxiliary class 32 for Romberg integrator for Fourier harmonics
     * calculations
     */
    private class UnivariateFourierHarmonics32 implements UnivariateFunction {

        private final double a1, a2, a3;

        public UnivariateFourierHarmonics32(double a1, double a2, double a3) {
            this.a1 = a1;
            this.a2 = a2;
            this.a3 = a3;
        }

        @Override
        public double value(double tau) {
            double tau2 = 2 * tau;
            return 1 + Math.cos(tau2) * Math.sin(getOrdernumber() * tau + a1 * Math.sin(tau) - a2 * Math.cos(tau) + a3 * Math.sin(tau2));
        }
    }
}
