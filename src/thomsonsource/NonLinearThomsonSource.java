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

    public NonLinearThomsonSource(AbstractLaserPulse l, AbstractElectronBunch b, int ordernumber) {
        super(l, b);
        this.ordernumber = ordernumber;
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
        double mv, M;
        double K1 = lp.getKA1()[0];
        double K2 = lp.getKA2()[0];
        Vector A1 = lp.getA1()[0];
        Vector A2 = lp.getA2()[0];
        double f01, f02, f11, f12, f21, f22, f31, f32;
        double gamma2 = eb.getGamma() * eb.getGamma();
        mv = Math.sqrt(1.0 - 1.0 / gamma2);//Dimesionaless speed

        double a1 = directionEnergy(n, v) / lp.getPhotonEnergy() * n.innerProduct(A1)
                / Math.sqrt(getsIntensity()) / (1 + mv) / eb.getGamma();
        double a2 = directionEnergy(n, v) / lp.getPhotonEnergy() * n.innerProduct(A2)
                / Math.sqrt(getsIntensity()) / (1 + mv) / eb.getGamma();
        double a3 = directionEnergy(n, v) / lp.getPhotonEnergy() * (K1 - K2) / getsIntensity()
                * (mv + n.innerProduct(v)) / Math.pow(eb.getGamma() * (1 + mv), 2) / mv / 8;
        System.out.println(directionEnergy(n, v) + "\n");
        /*
        Calculating the fourie harmonics
         */
        // 01 component
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateFunction func = new UnivariateFourierHarmonics01(a1, a2, a3);
        try {
            f01 = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, 0, 2 * Math.PI) / 2 / Math.PI;
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
        // 02 component
        integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        func = new UnivariateFourierHarmonics02(a1, a2, a3);
        try {
            f02 = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, 0, 2 * Math.PI) / 2 / Math.PI;
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
        // 11 component
        integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        func = new UnivariateFourierHarmonics11(a1, a2, a3);
        try {
            f11 = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, 0, 2 * Math.PI) / 2 / Math.PI;
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
        // 12 component
        integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        func = new UnivariateFourierHarmonics12(a1, a2, a3);
        try {
            f12 = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, 0, 2 * Math.PI) / 2 / Math.PI;
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
        // 21 component
        integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        func = new UnivariateFourierHarmonics21(a1, a2, a3);
        try {
            f21 = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, 0, 2 * Math.PI) / 2 / Math.PI;
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
        // 22 component
        integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        func = new UnivariateFourierHarmonics22(a1, a2, a3);
        try {
            f22 = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, 0, 2 * Math.PI) / 2 / Math.PI;
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
        // 31 component
        integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        func = new UnivariateFourierHarmonics31(a1, a2, a3);
        try {
            f31 = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, 0, 2 * Math.PI) / 2 / Math.PI;
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
        // 32 component
        integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        func = new UnivariateFourierHarmonics32(a1, a2, a3);
        try {
            f32 = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, 0, 2 * Math.PI);
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }

        //Parameter of nolinearity
        M = lp.getIntensity() / getsIntensity() * (1 + n.innerProduct(v)) / 4 / gamma2 / (1 + mv);
        //Returning the total flux
        return -getTotalFlux() * 3 / 2 / Math.PI * SIGMA_T
                / Math.pow((1 - n.innerProduct(v)) * (1 + M), 2) / (K1 + K2)
                * ((Math.pow(f01, 2) + Math.pow(f02, 2)) * (Math.pow(AS, 2) + (K1 + K2) / 2) - (Math.pow(f11, 2) + Math.pow(f12, 2)) * K1
                - (Math.pow(f21, 2) + Math.pow(f22, 2)) * K2 + (f31 * f01 + f32 * f02) * (K1 - K2) / 2);
    }

    @Override
    public double directionEnergy(Vector n, Vector v) {
        double mv, M, pr, gamma_2;
        gamma_2 = 1.0 / eb.getGamma() / eb.getGamma();
        mv = Math.sqrt(1.0 - gamma_2);
        pr = n.innerProduct(v);
        M = lp.getIntensity() / getsIntensity() * (1 + pr) * gamma_2 / (1 + mv) / 4;
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
        this.sIntensity = 1e-3 * AbstractLaserPulse.C * k * k * AS * AS / 8 / Math.PI;
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
            return Math.cos(getOrdernumber() * tau + a1 * Math.sin(tau) - a2 * Math.cos(tau) + a3 * Math.sin(2 * tau));
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
            return Math.sin(getOrdernumber() * tau + a1 * Math.sin(tau) - a2 * Math.cos(tau) + a3 * Math.sin(2 * tau));
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
            return Math.cos(tau) * Math.cos(getOrdernumber() * tau + a1 * Math.sin(tau) - a2 * Math.cos(tau) + a3 * Math.sin(2 * tau));
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
            return Math.cos(tau) * Math.sin(getOrdernumber() * tau + a1 * Math.sin(tau) - a2 * Math.cos(tau) + a3 * Math.sin(2 * tau));
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
            return Math.sin(tau) * Math.cos(getOrdernumber() * tau + a1 * Math.sin(tau) - a2 * Math.cos(tau) + a3 * Math.sin(2 * tau));
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
            return Math.sin(tau) * Math.sin(getOrdernumber() * tau + a1 * Math.sin(tau) - a2 * Math.cos(tau) + a3 * Math.sin(2 * tau));
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
            return Math.cos(2 * tau) * Math.cos(getOrdernumber() * tau + a1 * Math.sin(tau) - a2 * Math.cos(tau) + a3 * Math.sin(2 * tau));
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
            return Math.cos(2 * tau) * Math.sin(getOrdernumber() * tau + a1 * Math.sin(tau) - a2 * Math.cos(tau) + a3 * Math.sin(2 * tau));
        }
    }
}
