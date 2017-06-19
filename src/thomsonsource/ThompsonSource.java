/*
 * Copyright (C) 2015 Ruslan Feshchenko
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

import laserpulse.LaserPulse;
import electronbunch.AbstractElectronBunch;
import electronbunch.ElectronBunch;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.DoubleAdder;
import laserpulse.AbstractLaserPulse;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.la4j.Vector;
import org.apache.commons.math3.analysis.integration.*;
import org.apache.commons.math3.complex.Complex;
import org.la4j.Matrix;
import org.la4j.matrix.dense.Basic1DMatrix;
import org.la4j.Vectors;
import org.la4j.vector.dense.BasicVector;

/**
 * The main class containing all physics of LEXG
 *
 * @author Ruslan Feshchenko
 * @version 3.0
 */
public class ThompsonSource extends AbstractThomsonSource {

    /**
     * Constructor
     *
     * @param l
     * @param b
     */
    public ThompsonSource(LaserPulse l, AbstractElectronBunch b) {
        this.threadNumber = Runtime.getRuntime().availableProcessors();
        this.lp = l;
        this.eb = b;
        this.counter = new AtomicInteger();
        this.partialFlux = new DoubleAdder();
        calculateTotalFlux();
        calculateGeometricFactor();
    }

    /**
     * A method calculating the flux density in a given direction for a given
     * X-ray photon energy without taking into account electron transversal
     * pulse spread
     *
     * @param n direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     */
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

    /**
     * A method calculating the Stocks parameters density in a given direction
     * for a given X-ray photon energy without taking into account electron
     * transversal pulse spread
     *
     * @param n direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     */
    @Override
    public double[] directionFrequencyPolarizationNoSpread(Vector n, Vector v, double e) {
        double[] array = new double[NUMBER_OF_POL_PARAM];
        double K, th, m11, m22, m12, mlt, cs, sn;
        th = (1 - n.innerProduct(v)) * 2;
        mlt = 1 - e * th / lp.getPhotonEnergy() / 2;
        K = Math.pow((Math.sqrt(e / lp.getPhotonEnergy() / (1 - e * th / lp.getPhotonEnergy() / 4)) - 2 * eb.getGamma()), 2)
                / 4 / Math.pow(eb.getGamma() * eb.getDelgamma(), 2);
        m11 = getTotalFlux() * e * 3.0 / 32 / Math.PI / Math.sqrt(Math.PI) / eb.getDelgamma() / eb.getGamma() / lp.getPhotonEnergy()
                * Math.sqrt(e / lp.getPhotonEnergy()) / Math.sqrt(1 - e * th / lp.getPhotonEnergy() / 4) * Math.exp(-K);
        m12 = m11 * mlt;
        m22 = m12 * mlt;
        //Determine the polarization rotation angle
        double vn = v.innerProduct(n);
        double norm = Math.sqrt((1 - vn * vn) * (1 - n.get(0) * n.get(0)));
        if (norm != 0) {
            cs = (v.get(0) - n.get(0) * vn) / norm;
            sn = (n.get(1) * v.get(2) - n.get(2) * v.get(1)) / norm;
        } else {
            cs = 1;
            sn = 0;
        }
        double cs2 = 2 * cs * cs - 1, sn2 = 2 * sn * cs;
        double cs2cs2 = cs2 * cs2, sn2sn2 = sn2 * sn2, cs2sn2 = sn2 * cs2;
        //Calculating Stocks parameters
        array[0] = (m11 + m22 - (cs2 * lp.getPolarization()[0] + sn2 * lp.getPolarization()[1]) * (m11 - m22)) / 2;
        array[1] = (cs2 * (m22 - m11) + lp.getPolarization()[0] * (cs2cs2 * (m11 + m22) + 2 * sn2sn2 * m12)
                + lp.getPolarization()[1] * cs2sn2 * (m11 + m22 - 2 * m12)) / 2;
        array[2] = (sn2 * (m22 - m11) + lp.getPolarization()[1] * (sn2sn2 * (m11 + m22) + 2 * cs2cs2 * m12)
                + lp.getPolarization()[0] * cs2sn2 * (m11 + m22 - 2 * m12)) / 2;
        array[3] = lp.getPolarization()[2] * m12;
        return array;
    }


    /**
     * An auxiliary class for Romberg integrator for flux calculations
     */
    private class UnivariateFrequencyFluxSpreadOuter implements UnivariateFunction {

        private final double e;
        private final Vector n, v0;
        private final BaseAbstractUnivariateIntegrator inergrator;

        public UnivariateFrequencyFluxSpreadOuter(double e, Vector v0, Vector n) {
            this.e = e;
            this.v0 = v0;
            this.n = n;
            this.inergrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                    RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        }

        @Override
        public double value(double phi) {
            UnivariateFunction func
                    = new UnivariateFrequencyFluxSpreadInner(phi, e, v0, n);
            try {
                return inergrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, 0.0, INT_RANGE * eb.getSpread());
            } catch (TooManyEvaluationsException ex) {
                return 0;
            }
        }
    }

    /**
     * An auxiliary class for Romberg integrator for polarization calculations
     */
    private class UnivariateFrequencyPolarizationSpreadOuter implements UnivariateFunction {

        private final double e;
        private final Vector n, v0;
        private final int index;
        private final BaseAbstractUnivariateIntegrator inergrator;

        public UnivariateFrequencyPolarizationSpreadOuter(double e, Vector v0, Vector n, int index) {
            this.e = e;
            this.v0 = v0;
            this.n = n;
            this.index = index;
            this.inergrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                    RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        }

        @Override
        public double value(double phi) {
            if (Thread.currentThread().isInterrupted()) {
                Thread.currentThread().interrupt();
                return 0;
            }
            UnivariateFunction func
                    = new UnivariateFrequencyPolarizationSpreadInner(phi, e, v0, n, index);
            try {
                double rs = inergrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, 0.0, INT_RANGE * eb.getSpread());
                return rs;
            } catch (TooManyEvaluationsException ex) {
                return 0;
            }
        }
    }

    /**
     * An auxiliary class for Romberg integrator for flux calculations
     */
    private class UnivariateFrequencyFluxSpreadInner implements UnivariateFunction {

        private final double snphi, csphi, e;
        private final Vector n, v0;

        public UnivariateFrequencyFluxSpreadInner(double phi, double e, Vector v0, Vector n) {
            this.snphi = Math.sin(phi);
            this.csphi = Math.cos(phi);
            this.e = e;
            this.n = n;
            this.v0 = v0;
        }

        @Override
        public double value(double theta) {
            double u, sn = Math.sin(theta);
            Vector v = new BasicVector(new double[]{sn * csphi, sn * snphi, Math.cos(theta)});
            Vector dv = v.subtract(v0);
            u = theta * directionFrequencyFluxNoSpread(n, v, e) * eb.angleDistribution(dv.get(0), dv.get(1));
            return new Double(u).isNaN() ? 0 : u;
        }
    }

    /**
     * An auxiliary class for Romberg integrator for polarization calculations
     */
    private class UnivariateFrequencyPolarizationSpreadInner implements UnivariateFunction {

        private final double snphi, csphi, e;
        private final int index;
        private final Vector n, v0;

        public UnivariateFrequencyPolarizationSpreadInner(double phi, double e, Vector v0, Vector n, int index) {
            this.snphi = Math.sin(phi);
            this.csphi = Math.cos(phi);
            this.e = e;
            this.n = n;
            this.index = index;
            this.v0 = v0;
        }

        @Override
        public double value(double theta) {
            if (Thread.currentThread().isInterrupted()) {
                Thread.currentThread().interrupt();
                return 0;
            }
            double u, sn = Math.sin(theta);
            Vector v = new BasicVector(new double[]{sn * csphi, sn * snphi, Math.cos(theta)});
            Vector dv = v.subtract(v0);
            u = theta * directionFrequencyPolarizationNoSpread(n, v, e)[index] * eb.angleDistribution(dv.get(0), dv.get(1));
            return new Double(u).isNaN() ? 0 : u;
        }
    }


    /**
     * A method giving the flux density in a given direction
     *
     * @param n direction
     * @param v normalized electron velocity
     * @return
     */
    @Override
    public double directionFlux(Vector n, Vector v) {
        double gamma2, th;
        th = (1 - n.innerProduct(v)) * 2;
        gamma2 = eb.getGamma() * eb.getGamma();
        return getTotalFlux() * 3.0 / 2 / Math.PI * gamma2 * (1 + Math.pow(th * gamma2, 2))
                / Math.pow((1 + gamma2 * th), 4) * getGeometricFactor();
    }

    /**
     * A method calculating X-ray energy in a given direction
     *
     * @param n direction
     * @param v normalized electron velocity
     * @return
     */
    @Override
    public double directionEnergy(Vector n, Vector v) {
        double mv;
        mv = Math.sqrt(1.0 - 1.0 / eb.getGamma() / eb.getGamma());
        return 2 * lp.getPhotonEnergy() / (1 - n.innerProduct(v) * mv);
    }


    /**
     * An auxiliary class for the brilliance calculations
     */
    private class UnivariateVolumeFlux implements UnivariateFunction {

        Vector r0;
        Vector n0;

        public UnivariateVolumeFlux(Vector r0, Vector n0) {
            this.r0 = r0;
            this.n0 = n0;
        }

        @Override
        public double value(double x) {
            Vector r;
            r = r0.add(n0.multiply(x));
            double y = volumeFlux(r);
            if (n0.get(0) + n0.get(1) + n0.get(2) == 0) {
                throw new LocalException(x);
            }
            return y;
        }
    }

    /**
     * A custom exception class
     */
    private static class LocalException extends RuntimeException {

        // The x value that caused the problem.
        private final double x;

        public LocalException(double x) {
            this.x = x;
        }

        public double getX() {
            return x;
        }
    }

}
