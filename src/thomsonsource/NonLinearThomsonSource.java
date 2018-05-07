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
import org.apache.commons.math3.complex.*;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.la4j.Vectors;
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
        double gamma2;
        double mv, M;
        double K1 = lp.getIntensity() * lp.getPolarization()[0];
        double K2 = lp.getIntensity() * lp.getPolarization()[3];
        double a1 = 1, a2 = 1, a3 = 1;
        double f01, f02, f11, f12, f21, f22, f31, f32;
        Complex f0 = new Complex(1, 0), f1 = new Complex(1, 0), f2 = new Complex(1, 0), f3 = new Complex(1, 0);
        //Calculating the fourie harmonics
        // 01 component
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateFunction func = new UnivariateFourieHarmonics01(a1, a2, a3);
        try {
            f01 = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, 0, 2 * Math.PI);
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
        // 02 component
        integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        func = new UnivariateFourieHarmonics02(a1, a2, a3);
        try {
            f02 = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, 0, 2 * Math.PI);
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
        gamma2 = eb.getGamma() * eb.getGamma();
        mv = Math.sqrt(1.0 - 1.0 / gamma2);
        //Parameter of nolinearity
        M = lp.getIntensity() / getsIntensity() * (1 + n.innerProduct(v)) / 4 / gamma2 / (1 + mv);
        //Returning the total flux
        return -getTotalFlux() * 3 / 2 / Math.PI * SIGMA_T
                / Math.pow((1 - n.innerProduct(v)) * (1 + M), 2) / (K1 + K2)
                * ((Math.pow(f01, 2) + Math.pow(f02, 2)) * (Math.pow(AS, 2) + (K1 + K2) / 2) - Math.pow(f1.abs(), 2) * K1
                - Math.pow(f2.abs(), 2) * K2 + f3.multiply(f0.conjugate()).getReal() * (K1 - K2) / 2);
    }

    @Override
    public double directionEnergy(Vector n, Vector v) {
        double mv, M;
        mv = Math.sqrt(1.0 - 1.0 / eb.getGamma() / eb.getGamma());
        M = lp.getIntensity() / getsIntensity() * (1 + n.innerProduct(v) / mv) / 4 / eb.getGamma() / eb.getGamma() / (1 + mv);
        return ordernumber * (1 + mv) * lp.getPhotonEnergy() / (1 - n.innerProduct(v) * mv + M);
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
        this.sIntensity = 1e-3 * Math.pow(AbstractLaserPulse.C, 3) * k * k * AS * AS / 8 / Math.PI;
    }

    /**
     * An auxiliary class 01 for Romberg integrator for Fourie harmonics
     * calculations
     */
    private class UnivariateFourieHarmonics01 implements UnivariateFunction {

        private final double a1, a2, a3;

        public UnivariateFourieHarmonics01(double a1, double a2, double a3) {
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
     * An auxiliary class 02 for Romberg integrator for Fourie harmonics
     * calculations
     */
    private class UnivariateFourieHarmonics02 implements UnivariateFunction {

        private final double a1, a2, a3;

        public UnivariateFourieHarmonics02(double a1, double a2, double a3) {
            this.a1 = a1;
            this.a2 = a2;
            this.a3 = a3;
        }

        @Override
        public double value(double tau) {
            return Math.sin(getOrdernumber() * tau + a1 * Math.sin(tau) - a2 * Math.cos(tau) + a3 * Math.sin(2 * tau));
        }
    }
}
