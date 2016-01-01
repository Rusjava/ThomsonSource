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

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.DoubleAdder;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.la4j.vector.Vector;
import org.apache.commons.math3.analysis.integration.*;
import org.apache.commons.math3.complex.Complex;
import org.la4j.matrix.Matrix;
import org.la4j.matrix.dense.Basic1DMatrix;
import org.la4j.vector.Vectors;
import org.la4j.vector.dense.BasicVector;

/**
 * The main class containing all physics of LEXG
 *
 * @author Ruslan Feshchenko
 * @version 1.9
 */
public class ThompsonSource implements Cloneable {

    /**
     * Constructor
     *
     * @param l
     * @param b
     */
    public ThompsonSource(LaserPulse l, ElectronBunch b) {
        this.threadNumber = Runtime.getRuntime().availableProcessors();
        this.lp = l;
        this.eb = b;
        this.counter = new AtomicInteger();
        this.partialFlux = new DoubleAdder();
        calculateTotalFlux();
        calculateGeometricFactor();
    }

    /**
     * Multiple for the range
     */
    private final double INT_RANGE = 3;
    /**
     * The number of columns in Shadow files
     */
    public static final int NUMBER_OF_COLUMNS = 18;

    /**
     * Angle range for rays exported for Shadow in the X-direction
     */
    private double rayXAnglerange = 0.0003;

    /**
     * Angle range for rays exported for Shadow in the Y-direction
     */
    private double rayYAnglerange = 0.0003;

    /**
     * Min ray energy
     */
    private double minEnergy = 36 * 1e3 * ElectronBunch.E;

    /**
     * Max ray energy
     */
    private double maxEnergy = 46 * 1e3 * ElectronBunch.E;

    /**
     * Number of points in Monte Carlo calculation of the geometric factor
     */
    private int npGeometricFactor = 5000000;

    /**
     * Maximal number of evaluations in calculations of the brilliance
     */
    public static final int MAXIMAL_NUMBER_OF_EVALUATIONS = 30000;

    /**
     * Precision in calculations of the brilliance
     */
    private double precision = 0.0001;

    /**
     * Normalized total flux from the source
     */
    private double totalFlux;

    /**
     * Geometric factor. Assumes values from 0 to 1
     */
    private double geometricFactor = 1;

    /**
     * Flag - whether or not the electron beam transversal velocity spread is
     * taken into account
     */
    private boolean eSpread = false;

    /**
     * Flux in the phase space volume of ray generation
     */
    private DoubleAdder partialFlux;

    /**
     * Number of used threads
     */
    private int threadNumber;

    /**
     * Counter of ray iterations
     */
    private AtomicInteger counter;

    private LaserPulse lp;
    private ElectronBunch eb;

    private double[] ksi = null;

    @Override
    public Object clone() throws CloneNotSupportedException {
        Object tm = super.clone();
        ((ThompsonSource) tm).eb = (ElectronBunch) this.eb.clone();
        ((ThompsonSource) tm).lp = (LaserPulse) this.lp.clone();
        ((ThompsonSource) tm).counter = new AtomicInteger();
        ((ThompsonSource) tm).partialFlux = new DoubleAdder();
        return tm;
    }

    /**
     * Returns Electron Bunch reference
     *
     * @return
     */
    public ElectronBunch getElectronBunch() {
        return eb;
    }

    /**
     * Returns laser pulse reference
     *
     * @return
     */
    public LaserPulse getLaserPulse() {
        return lp;
    }

    /**
     * The full Thomson cross-section
     */
    public final static double SIGMA_T = 6.65e-29;

    /**
     * A method calculating normalized total flux
     */
    public void calculateTotalFlux() {
        this.totalFlux = SIGMA_T * eb.getNumber() * lp.getPhotonNumber()
                * lp.getFq() / Math.PI / Math.sqrt((lp.getWidth2(0.0) + eb.getxWidth2(0.0))
                        * (lp.getWidth2(0.0) + eb.getyWidth2(0.0)));
    }

    /**
     * A method calculating the geometric factor
     *
     */
    public void calculateGeometricFactor() {
        ExecutorService execs = Executors.newFixedThreadPool(threadNumber);
        // We need to synchronize threads
        CountDownLatch lt = new CountDownLatch(threadNumber);
        // Atomic adder
        DoubleAdder sum = new DoubleAdder();
        double wdx, wdy, len;
        int mult = 2;
        wdx = mult * Math.max(eb.getxWidth(0.0) + Math.abs(eb.getShift().get(0)) / 2, lp.getWidth(0.0) + Math.abs(eb.getShift().get(0)) / 2);
        wdy = mult * Math.max(eb.getyWidth(0.0) + Math.abs(eb.getShift().get(1)) / 2, lp.getWidth(0.0) + Math.abs(eb.getShift().get(1)) / 2);
        len = mult * Math.max(eb.getLength() + Math.abs(eb.getShift().get(2)) / 2, lp.getLength() + Math.abs(eb.getShift().get(2)) / 2);
        final int itNumber = Math.round(getNpGeometricFactor() / threadNumber);
        /*
         Splitting the job into a number of threads
         */
        for (int m = 0; m < threadNumber; m++) {
            execs.execute(() -> {
                double psum = 0;
                Vector iter = new BasicVector(new double[]{0.0, 0.0, 0.0});
                for (int i = 0; i < itNumber; i++) {
                    if (Thread.currentThread().isInterrupted()) {
                        return;
                    }
                    iter.set(0, eb.getShift().get(0) / 2 + wdx * (2 * Math.random() - 1.0));
                    iter.set(1, eb.getShift().get(1) / 2 + wdy * (2 * Math.random() - 1.0));
                    iter.set(2, eb.getShift().get(2) / 2 + len * (2 * Math.random() - 1.0));
                    psum += volumeFlux(iter);
                }
                sum.add(psum);
                lt.countDown();
            });
        }
        try {
            lt.await();
        } catch (InterruptedException ex) {
            Thread.currentThread().interrupt();
        }
        execs.shutdownNow();
        this.geometricFactor = 8 * wdx * wdy * len * sum.sum() / itNumber / threadNumber;    
    }

    /**
     * A method calculating the flux density in a given direction for a given
     * X-ray photon energy
     *
     * @param n direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     */
    public double directionFrequencyFlux(Vector n, Vector v, double e) {
        return iseSpread() ? directionFrequencyFluxSpread(n, v, e) : directionFrequencyFluxNoSpread(n, v, e);
    }

    /**
     * A method calculating the Stocks parameters density in a given direction
     * for a given X-ray photon energy
     *
     * @param n direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     */
    public double[] directionFrequencyPolarization(Vector n, Vector v, double e) {
        return iseSpread() ? directionFrequencyPolarizationSpread(n, v, e) : directionFrequencyPolarizationNoSpread(n, v, e);
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
    public double[] directionFrequencyPolarizationNoSpread(Vector n, Vector v, double e) {
        double[] array = new double[4];
        double K, th, m11, m22, m12, mlt;
        th = (1 - n.innerProduct(v)) * 2;
        mlt = 1 - e * th / lp.getPhotonEnergy() / 2;
        K = Math.exp(-Math.pow((Math.sqrt(e / lp.getPhotonEnergy() / (1 - e * th / lp.getPhotonEnergy() / 4)) - 2 * eb.getGamma()), 2)
                / 4 / Math.pow(eb.getGamma() * eb.getDelgamma(), 2));
        m11 = getTotalFlux() * e * 3.0 / 32 / Math.PI / Math.sqrt(Math.PI) / eb.getDelgamma() / eb.getGamma() / lp.getPhotonEnergy()
                * Math.sqrt(e / lp.getPhotonEnergy()) / Math.sqrt(1 - e * th / lp.getPhotonEnergy() / 4) * K;
        m12 = m11 * mlt;
        m22 = m12 * mlt;
        array[0] = (m22 * (1 + lp.getPolarization()[0]) + m11 * (1 - lp.getPolarization()[0])) / 2;
        array[1] = (m22 * (1 + lp.getPolarization()[0]) - m11 * (1 - lp.getPolarization()[0])) / 2;
        array[2] = m12 * lp.getPolarization()[1];
        array[3] = m12 * lp.getPolarization()[2];
        return array;
    }

    /**
     * A method calculating the flux density in a given direction for a given
     * X-ray photon energy taking into account electron transversal pulse spread
     *
     * @param n direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     */
    public double directionFrequencyFluxSpread(Vector n, Vector v, double e) {
        BaseAbstractUnivariateIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateFunction func
                = new UnivariateFrequencyFluxSpreadOuter(e, v, n);
        try {
            return integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, 0.0, 2 * Math.PI);
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
    }

    /**
     * A method calculating the stocks parameters density in a given direction
     * for a given X-ray photon energy taking into account electron transversal
     * pulse spread
     *
     * @param n direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     */
    public double[] directionFrequencyPolarizationSpread(Vector n, Vector v, double e) {
        BaseAbstractUnivariateIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        double[] array = new double[4];
        UnivariateFunction func;
        for (int i = 0; i < 4; i++) {
            func = new UnivariateFrequencyPolarizationSpreadOuter(e, v, n, i);
            try {
                array[i] = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, 0.0, 2 * Math.PI);
            } catch (TooManyEvaluationsException ex) {
                array[i] = 0;
            }
        }
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
            UnivariateFunction func
                    = new UnivariateFrequencyPolarizationSpreadInner(phi, e, v0, n, index);
            try {
                return inergrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, 0.0, INT_RANGE * eb.getSpread());
            } catch (TooManyEvaluationsException ex) {
                return 0;
            }
        }
    }

    /**
     * An auxiliary class for Romberg integrator for flux calculations
     */
    private class UnivariateFrequencyFluxSpreadInner implements UnivariateFunction {

        private final double phi, e;
        private final Vector n, v0;

        public UnivariateFrequencyFluxSpreadInner(double phi, double e, Vector v0, Vector n) {
            this.phi = phi;
            this.e = e;
            this.n = n;
            this.v0 = v0;
        }

        @Override
        public double value(double theta) {
            double u;
            Vector v = new BasicVector(new double[]{Math.sin(theta) * Math.cos(phi),
                Math.sin(theta) * Math.sin(phi), Math.cos(theta)});
            Vector dv = v.subtract(v0);
            u = theta * directionFrequencyFluxNoSpread(n, v, e)
                    * eb.angleDistribution(dv.get(0), dv.get(1));
            return new Double(u).isNaN() ? 0 : u;
        }
    }

    /**
     * An auxiliary class for Romberg integrator for polarization calculations
     */
    private class UnivariateFrequencyPolarizationSpreadInner implements UnivariateFunction {

        private final double phi, e;
        private final int index;
        private final Vector n, v0;

        public UnivariateFrequencyPolarizationSpreadInner(double phi, double e, Vector v0, Vector n, int index) {
            this.phi = phi;
            this.e = e;
            this.n = n;
            this.index = index;
            this.v0 = v0;
        }

        @Override
        public double value(double theta) {
            double u;
            Vector v = new BasicVector(new double[]{Math.sin(theta) * Math.cos(phi),
                Math.sin(theta) * Math.sin(phi), Math.cos(theta)});
            Vector dv = v.subtract(v0);
            u = theta * directionFrequencyPolarizationNoSpread(n, v, e)[index]
                    * eb.angleDistribution(dv.get(0), dv.get(1));
            return new Double(u).isNaN() ? 0 : u;
        }
    }

    /**
     * A method calculating the flux density in a given direction for a given
     * X-ray photon energy for a given volume element
     *
     * @param r spatial position
     * @param n direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     */
    public double directionFrequencyVolumeFlux(Vector r, Vector n, Vector v, double e) {
        return iseSpread() ? directionFrequencyVolumeFluxSpread(r, n, v, e) : directionFrequencyVolumeFluxNoSpread(r, n, v, e);
    }

    /**
     * A method calculating the Stocks parameters density in a given direction
     * for a given X-ray photon energy for a given volume element
     *
     * @param r spatial position
     * @param n direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     */
    public double[] directionFrequencyVolumePolarization(Vector r, Vector n, Vector v, double e) {
        return iseSpread() ? directionFrequencyVolumePolarizationSpread(r, n, v, e) : directionFrequencyVolumePolarizationNoSpread(r, n, v, e);
    }

    /**
     * A method calculating the flux density in a given direction for a given
     * X-ray photon energy for a given volume element without taking into
     * account electron transversal pulse spread
     *
     * @param r spatial position
     * @param n direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     */
    public double directionFrequencyVolumeFluxNoSpread(Vector r, Vector n, Vector v, double e) {
        return directionFrequencyFluxNoSpread(n, v, e) * volumeFlux(r);
    }

    /**
     * A method calculating the Stocks parameters density in a given direction
     * for a given X-ray photon energy for a given volume element without taking
     * into account electron transversal pulse spread
     *
     * @param r spatial position
     * @param n direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     */
    public double[] directionFrequencyVolumePolarizationNoSpread(Vector r, Vector n, Vector v, double e) {
        double[] stocks = directionFrequencyPolarizationNoSpread(n, v, e);
        double vFlux = volumeFlux(r);
        for (int i = 0; i < 4; i++) {
            stocks[i] *= vFlux;
        }
        return stocks;
    }

    /**
     * A method calculating the flux density in a given direction for a given
     * X-ray photon energy for a given volume element taking into account
     * electron transversal pulse spread
     *
     * @param r spatial position
     * @param n direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     */
    public double directionFrequencyVolumeFluxSpread(Vector r, Vector n, Vector v, double e) {
        return directionFrequencyFluxSpread(n, v, e) * volumeFlux(r);
    }

    /**
     * A method calculating the Stocks parameters density in a given direction
     * for a given X-ray photon energy for a given volume element taking into
     * account electron transversal pulse spread
     *
     * @param r spatial position
     * @param n direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     */
    public double[] directionFrequencyVolumePolarizationSpread(Vector r, Vector n, Vector v, double e) {
        double[] stocks = directionFrequencyPolarizationSpread(n, v, e);
        double vFlux = volumeFlux(r);
        for (int i = 0; i < 4; i++) {
            stocks[i] *= vFlux;
        }
        return stocks;
    }

    /**
     * An auxiliary method calculating volume density of the X-ray source
     *
     * @param r spatial position
     * @return
     */
    public double volumeFlux(Vector r) {
        double u, z0, z, z1, x0, x, x1, y0, y, y1, sn, cs, K, len;
        len = Math.sqrt(lp.getLength() * lp.getLength()
                + eb.getLength() * eb.getLength());
        sn = lp.getDirection().get(1);
        cs = lp.getDirection().get(2);
        x0 = eb.getShift().get(0);
        y0 = eb.getShift().get(1);
        z0 = eb.getShift().get(2);
        x = r.get(0);
        y = r.get(1);
        z = r.get(2);
        x1 = x;
        y1 = -sn * z + cs * y;
        z1 = cs * z + sn * y;
        K = Math.pow((z + z1 - z0 - lp.getDelay()) / len, 2)
                + Math.pow((x - x0), 2) / eb.getxWidth2(z - z0) + Math.pow((y - y0), 2) / eb.getyWidth2(z - z0)
                + (Math.pow(x1, 2) + Math.pow(y1, 2)) / lp.getWidth2(z1);
        u = 2.0 / Math.pow(Math.PI, 1.5) * Math.sqrt((lp.getWidth2(0.0)
                + eb.getxWidth2(0.0)) * (lp.getWidth2(0.0)
                + eb.getyWidth2(0.0))) / len / lp.getWidth2(z1) / eb.getxWidth(z - z0) / eb.getyWidth(z - z0) * Math.exp(-K);
        return new Double(u).isNaN() ? 0 : u;
    }

    /**
     * A method giving the flux density in a given direction
     *
     * @param n direction
     * @param v normalized electron velocity
     * @return
     */
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
    public double directionEnergy(Vector n, Vector v) {
        double mv;
        mv = Math.sqrt(1.0 - 1.0 / eb.getGamma() / eb.getGamma());
        return 2 * lp.getPhotonEnergy() / (1 - n.innerProduct(v) * mv);
    }

    /**
     * A method calculating spectral brilliance in a given direction
     *
     * @param r0 spatial position for brightness
     * @param n direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     */
    public double directionFrequencyBrilliance(Vector r0, Vector n, Vector v, double e) {
        return iseSpread() ? directionFrequencyBrillianceSpread(r0, n, v, e) : directionFrequencyBrillianceNoSpread(r0, n, v, e);
    }

    /**
     * A method calculating spectral brilliance in a given direction with
     * polarization
     *
     * @param r0 spatial position for brightness
     * @param n direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     */
    public double[] directionFrequencyPolarizationBrilliance(Vector r0, Vector n, Vector v, double e) {
        return iseSpread() ? directionFrequencyBrilliancePolarizationSpread(r0, n, v, e) : directionFrequencyBrilliancePolarizationNoSpread(r0, n, v, e);
    }

    /**
     * A method calculating spectral brilliance in a given direction without
     * taking into account electron transversal pulse spread
     *
     * @param r0 spatial position for brightness
     * @param n direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     */
    public double directionFrequencyBrillianceNoSpread(Vector r0, Vector n, Vector v, double e) {
        double u;
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateVolumeFlux func = new UnivariateVolumeFlux(r0, n);
        try {
            u = integrator.integrate(30000, func,
                    r0.fold(Vectors.mkEuclideanNormAccumulator()) - 3 * eb.getLength(),
                    r0.fold(Vectors.mkEuclideanNormAccumulator()) + 3 * eb.getLength());
            return u * directionFrequencyFluxNoSpread(n, v, e);
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
    }

    /**
     * A method calculating spectral brilliance in a given direction without
     * taking into account electron transversal pulse spread but with
     * polarization
     *
     * @param r0 spatial position for brightness
     * @param n direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     */
    public double[] directionFrequencyBrilliancePolarizationNoSpread(Vector r0, Vector n, Vector v, double e) {
        double mlt;
        double[] array = new double[4];
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateVolumeFlux func = new UnivariateVolumeFlux(r0, n);
        try {
            mlt = integrator.integrate(30000, func,
                    r0.fold(Vectors.mkEuclideanNormAccumulator()) - 3 * eb.getLength(),
                    r0.fold(Vectors.mkEuclideanNormAccumulator()) + 3 * eb.getLength());
        } catch (TooManyEvaluationsException ex) {
            mlt = 0;
        }
        for (int i = 0; i < 4; i++) {
            array[i] = mlt * directionFrequencyPolarizationNoSpread(n, v, e)[i];
        }
        return array;
    }

    /**
     * A method calculating spectral brilliance in a given direction taking into
     * account electron transversal pulse spread
     *
     * @param r0 spatial position for brightness
     * @param n direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     */
    public double directionFrequencyBrillianceSpread(Vector r0, Vector n, Vector v, double e) {
        double u;
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateVolumeFlux func = new UnivariateVolumeFlux(r0, n);
        try {
            u = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func,
                    r0.fold(Vectors.mkEuclideanNormAccumulator()) - 3 * eb.getLength(),
                    r0.fold(Vectors.mkEuclideanNormAccumulator()) + 3 * eb.getLength());
            return u * directionFrequencyFluxSpread(n, v, e);
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
    }

    /**
     * A method calculating spectral brilliance in a given direction taking into
     * account electron transversal pulse spread nd with polarization
     *
     * @param r0 spatial position for brightness
     * @param n direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     */
    public double[] directionFrequencyBrilliancePolarizationSpread(Vector r0, Vector n, Vector v, double e) {
        double mlt;
        double[] array = new double[4];
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateVolumeFlux func = new UnivariateVolumeFlux(r0, n);
        try {
            mlt = integrator.integrate(30000, func,
                    r0.fold(Vectors.mkEuclideanNormAccumulator()) - 3 * eb.getLength(),
                    r0.fold(Vectors.mkEuclideanNormAccumulator()) + 3 * eb.getLength());
        } catch (TooManyEvaluationsException ex) {
            mlt = 0;
        }
        for (int i = 0; i < 4; i++) {
            array[i] = mlt * directionFrequencyPolarizationSpread(n, v, e)[i];
        }
        return array;
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

    /**
     * Returning a random ray
     *
     * @return an array with ray parameters
     * @throws java.lang.InterruptedException
     */
    public double[] getRay() throws InterruptedException {
        double[] ray = new double[NUMBER_OF_COLUMNS];
        Matrix T;
        Vector n = new BasicVector(new double[]{0.0, 0.0, 1.0});
        Vector r = new BasicVector(new double[]{0.0, 0.0, 0.0});
        Vector n0 = new BasicVector(new double[]{0.0, 1.0, 0.0}), As;
        double prob0, prob, EMax, mult = 2, factor, sum = 0;
        double[] pol, polParam;
        EMax = directionEnergy(n, n);
        factor = 64 * Math.max(eb.getxWidth(0.0), lp.getWidth(0.0)) * Math.max(eb.getyWidth(0.0), lp.getWidth(0.0))
                * Math.max(eb.getLength(), lp.getLength()) * 4 * rayXAnglerange * rayYAnglerange
                * (maxEnergy - minEnergy);
        prob0 = directionFrequencyVolumePolarizationNoSpread(r, n, new BasicVector(new double[]{0.0, 0.0, 1.0}), EMax)[0];
        if (iseSpread()) {
            prob0 *= eb.angleDistribution(0, 0);
            factor *= 4 * mult * mult * eb.getXSpread() * eb.getYSpread();
        }
        do {
            if (Thread.currentThread().isInterrupted()) {
                throw new InterruptedException();
            }
            ray[0] = 2 * (2 * Math.random() - 1.0) * Math.max(eb.getxWidth(0.0), lp.getWidth(0.0));
            ray[2] = 2 * (2 * Math.random() - 1.0) * Math.max(eb.getyWidth(0.0), lp.getWidth(0.0));
            ray[1] = 2 * (2 * Math.random() - 1.0) * Math.max(eb.getLength(), lp.getLength());
            r.set(0, ray[0]);
            r.set(1, ray[2]);
            r.set(2, ray[1]);
            ray[3] = rayXAnglerange * (2 * Math.random() - 1.0);
            ray[5] = rayYAnglerange * (2 * Math.random() - 1.0);
            n.set(0, ray[3]);
            n.set(1, ray[5]);
            n.set(2, 1.0);
            n = n.divide(n.fold(Vectors.mkEuclideanNormAccumulator()));
            ray[3] = n.get(0);
            ray[5] = n.get(1);
            ray[4] = n.get(2);
            ray[10] = Math.random() * (maxEnergy - minEnergy) + minEnergy;
            if (iseSpread()) {
                double thetax = mult * eb.getXSpread() * (2 * Math.random() - 1);
                double thetay = mult * eb.getYSpread() * (2 * Math.random() - 1);
                Vector v = new BasicVector(new double[]{thetax, thetay, Math.sqrt(1 - thetax * thetax - thetay * thetay)});
                polParam = directionFrequencyVolumePolarizationNoSpread(r, n, v, ray[10]);
                prob = polParam[0] * eb.angleDistribution(thetax, thetay);
            } else {
                polParam = directionFrequencyVolumePolarizationNoSpread(r, n, new BasicVector(new double[]{0.0, 0.0, 1.0}), ray[10]);
                prob = polParam[0];
            }
            if (!new Double(prob).isNaN()) {
                sum += prob / ray[10];
            }
            counter.incrementAndGet();
        } while (prob / prob0 < Math.random() || (new Double(prob)).isNaN());
        // Calculation of the rotated polarization vector and getting the full polarizaation state
        n = new BasicVector(new double[]{ray[3], ray[4], ray[5]});
        T = getTransform(n, n0);
        //Checking if polarization is pre-specified
        pol = (ksi != null) ? getPolarization(ksi)
                : getPolarization(new double[]{polParam[1] / polParam[0], polParam[2] / polParam[0], polParam[3] / polParam[0]});
        //Rotating the ray electrical vectors
        As = T.multiply(new BasicVector(new double[]{1.0, 0.0, 0.0})).multiply(pol[0]);
        ray[6] = As.get(0);
        ray[7] = As.get(1);
        ray[8] = As.get(2);
        As = T.multiply(new BasicVector(new double[]{0.0, 0.0, 1.0})).multiply(pol[1]);
        ray[15] = As.get(0);
        ray[16] = As.get(1);
        ray[17] = As.get(2);
        //Setting other columns
        ray[9] = 1.0;
        ray[13] = pol[2];
        ray[14] = pol[3];
        partialFlux.add(sum * factor);
        return ray;
    }

    /**
     * Setting ranges for the Shadow ray generation
     *
     * @param xangle angle in the X direction
     * @param yangle angle in the Y direction
     * @param minEn minimum ray energy
     * @param maxEn maximum ray energy
     */
    public void setRayRanges(double xangle, double yangle, double minEn, double maxEn) {
        rayXAnglerange = xangle;
        rayYAnglerange = yangle;
        minEnergy = minEn;
        maxEnergy = maxEn;
    }

    /**
     * Returning the matrix of 3D rotation based on two unity vectors
     *
     * @param n
     * @param n0
     * @return transformation matrix
     */
    protected Matrix getTransform(Vector n, Vector n0) {
        Matrix D, A, I = new Basic1DMatrix(3, 3);
        I.set(0, 0, 1.0);
        I.set(1, 1, 1.0);
        I.set(2, 2, 1.0);
        double innerProduct;
        innerProduct = n.innerProduct(n0);
        D = n.outerProduct(n0).add(n0.outerProduct(n)).multiply(innerProduct).subtract(n.outerProduct(n)
                .add(n0.outerProduct(n0))).divide(innerProduct * innerProduct - 1.0);
        A = n.outerProduct(n0).subtract(n0.outerProduct(n)).add(I.multiply(innerProduct));
        return I.subtract(D).multiply(1 - innerProduct).add(A);
    }

    /**
     * Setting Stocks parameters
     *
     * @param ksi
     */
    public void setPolarization(double[] ksi) {
        this.ksi = ksi;
    }

    /**
     * Number of points in Monte Carlo calculation of the geometric factor
     *
     * @return the npGeometricFactor
     */
    public int getNpGeometricFactor() {
        return npGeometricFactor;
    }

    /**
     * Number of points in Monte Carlo calculation of the geometric factor
     *
     * @param npGeometricFactor the npGeometricFactor to set
     */
    public void setNpGeometricFactor(int npGeometricFactor) {
        this.npGeometricFactor = npGeometricFactor;
    }

    /**
     * Precision in calculations of the brilliance
     *
     * @return the precision
     */
    public double getPrecision() {
        return precision;
    }

    /**
     * Precision in calculations of the brilliance
     *
     * @param precision the precision to set
     */
    public void setPrecision(double precision) {
        this.precision = precision;
    }

    /**
     * Normalized total flux from the source
     *
     * @return the totalFlux
     */
    public double getTotalFlux() {
        return totalFlux;
    }

    /**
     * Geometric factor. Assumes values from 0 to 1
     *
     * @return the geometricFactor
     */
    public double getGeometricFactor() {
        return geometricFactor;
    }

    /**
     * Flag - whether or not the electron beam transversal velocity spread is
     * taken into account
     *
     * @return the eSpread
     */
    public boolean iseSpread() {
        return eSpread;
    }

    /**
     * Flag - whether or not the electron beam transversal velocity spread is
     * taken into account
     *
     * @param eSpread the eSpread to set
     */
    public void seteSpread(boolean eSpread) {
        this.eSpread = eSpread;
    }

    /**
     * Flux in the phase space volume of ray generation
     *
     * @return the partialFlux
     */
    public double getPartialFlux() {
        return partialFlux.sum();
    }

    /**
     * Counter of ray iterations
     *
     * @return the counter
     */
    public int getCounter() {
        return counter.get();
    }

    /**
     * Setting the number of used threads
     *
     * @param threadNumber the threadNumber to set
     */
    public void setThreadNumber(int threadNumber) {
        this.threadNumber = threadNumber;
    }

    /**
     * Getting the number of used threads
     *
     * @return
     */
    public int getThreadNumber() {
        return threadNumber;
    }

    /**
     * Calculation of random amplitudes and phases for an arbitrary state of
     * polarization
     *
     * @return
     */
    private double[] getPolarization(double[] ksiVector) {
        double[] pol = new double[4];
        Complex phase1 = Complex.I.multiply(Math.random() * 2 * Math.PI).exp();
        Complex phase2 = Complex.I.multiply(Math.random() * 2 * Math.PI).exp();
        double p = Math.sqrt(ksiVector[0] * ksiVector[0] + ksiVector[1] * ksiVector[1] + ksiVector[2] * ksiVector[2]);
        //If p > 1 reducing it yo 1
        p = p > 1 ? 1 : p;
        //Special case when the denomonator is zero
        if (ksiVector[0] == -p) {
            pol[0] = Math.sqrt((1 - ksiVector[0]) / 2);
            pol[1] = Math.sqrt((1 + ksiVector[0]) / 2);
            pol[2] = Math.random() * 2 * Math.PI;
            pol[3] = Math.random() * 2 * Math.PI;
        } else {
            //General case
            double k1 = Math.sqrt(1 - p);
            double k2 = Math.sqrt(1 + p);
            double coef = Math.sqrt((p + ksiVector[0]) / p) / 2;
            Complex ksiD = new Complex(ksiVector[1], ksiVector[2]);
            Complex e1 = phase1.multiply(k1).add(phase2.multiply(ksiD).multiply(k2).divide(p + ksiVector[0])).multiply(coef);
            Complex e2 = phase2.multiply(k2).subtract(phase1.multiply(ksiD.conjugate()).multiply(k1).divide(p + ksiVector[0])).multiply(coef);
            pol[0] = e1.abs();
            pol[1] = e2.abs();
            pol[2] = e1.getArgument();
            pol[3] = e2.getArgument();
        }
        return pol;
    }
}
