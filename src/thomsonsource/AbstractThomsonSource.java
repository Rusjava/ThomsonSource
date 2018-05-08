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
import electronbunch.GaussianElectronBunch;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.DoubleAdder;
import laserpulse.AbstractLaserPulse;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.BaseAbstractUnivariateIntegrator;
import org.apache.commons.math3.analysis.integration.RombergIntegrator;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.la4j.Matrix;
import org.la4j.Vector;
import org.la4j.Vectors;
import org.la4j.matrix.dense.Basic1DMatrix;
import org.la4j.vector.dense.BasicVector;

/**
 * An abstract class for Thomson source. Methods that calculated scattering by
 * one electron need to be defined.
 *
 * @version 1.0
 * @author Ruslan Feshchenko
 */
public abstract class AbstractThomsonSource implements Cloneable {

    /**
     * Constructor
     *
     * @param l
     * @param b
     */
    public AbstractThomsonSource(AbstractLaserPulse l, AbstractElectronBunch b) {
        this.threadNumber = Runtime.getRuntime().availableProcessors();
        this.lp = l;
        this.eb = b;
        this.counter = new AtomicInteger();
        this.partialFlux = new DoubleAdder();
    }

    /**
     * Multiple for the range
     */
    private static final double INT_RANGE = 3;
    /**
     * The number of columns in Shadow files
     */
    public static final int NUMBER_OF_COLUMNS = 18;
    /**
     * The number of polarization parameters
     */
    public static final int NUMBER_OF_POL_PARAM = 4;
    /**
     * Maximal number of evaluations in calculations of the brilliance and
     * polarization
     */
    public static final int MAXIMAL_NUMBER_OF_EVALUATIONS = 1000000;
    /**
     * The full Thomson cross-section
     */
    public static final double SIGMA_T = 6.65e-29;
    /**
     * The saturating potential in CGS system
     */
    public static final double AS = 1.704509e3;
        /**
     * Angle range for rays exported for Shadow in the X-direction
     */
    protected double rayXAnglerange = 0.0003;
    /**
     * Angle range for rays exported for Shadow in the Y-direction
     */
    protected double rayYAnglerange = 0.0003;
    /**
     * Min ray energy
     */
    protected double minEnergy = 36 * 1e3 * GaussianElectronBunch.E;
    /**
     * Max ray energy
     */
    protected double maxEnergy = 46 * 1e3 * GaussianElectronBunch.E;
    /**
     * Number of points in Monte Carlo calculation of the geometric factor
     */
    protected int npGeometricFactor = 50000;
    /**
     * Precision in calculations of the brilliance
     */
    protected double precision = 0.0001;
    /**
     * Normalized total flux from the source
     */
    protected double totalFlux;
    /**
     * Geometric factor. Assumes values from 0 to 1
     */
    protected double geometricFactor = 1;
    /**
     * Flag - whether or not the electron beam transversal velocity spread is
     * taken into account
     */
    protected boolean eSpread = false;
    /**
     * Flux in the phase space volume of ray generation
     */
    protected DoubleAdder partialFlux;
    /**
     * Number of used threads
     */
    protected int threadNumber;
    /**
     * Counter of ray iterations
     */
    protected AtomicInteger counter;

    /**
     * A pointer for laser pulse
     */
    protected AbstractLaserPulse lp;

    /**
     * A pointer for electron bunch
     */
    protected AbstractElectronBunch eb;

    /**
     * An array for laser pulse polarization state
     */
    protected double[] ksi = null;

    @Override
    public Object clone() throws CloneNotSupportedException {
        Object tm = super.clone();
        ((AbstractThomsonSource) tm).eb = (AbstractElectronBunch) this.eb.clone();
        ((AbstractThomsonSource) tm).lp = (AbstractLaserPulse) this.lp.clone();
        ((AbstractThomsonSource) tm).counter = new AtomicInteger();
        ((AbstractThomsonSource) tm).partialFlux = new DoubleAdder();
        ((AbstractThomsonSource) tm).ksi = (double[]) ksi.clone();
        return tm;
    }

    /**
     * Returns Electron Bunch reference
     *
     * @return
     */
    public AbstractElectronBunch getElectronBunch() {
        return eb;
    }

    /**
     * Returns laser pulse reference
     *
     * @return
     */
    public AbstractLaserPulse getLaserPulse() {
        return lp;
    }

    /**
     * A method calculating normalized total flux
     */
    public final void calculateTotalFlux() {
        this.totalFlux = SIGMA_T * eb.getNumber() * lp.getPhotonNumber() * lp.getFq() / Math.PI / Math.sqrt((lp.getWidth2(0.0) + eb.getxWidth2(0.0)) * (lp.getWidth2(0.0) + eb.getyWidth2(0.0)));
    }

    /**
     * A method calculating the geometric factor
     *
     */
    public final void calculateGeometricFactor() {
        ExecutorService execs = Executors.newFixedThreadPool(threadNumber);
        // We need to synchronize threads
        CountDownLatch lt = new CountDownLatch(threadNumber);
        // Atomic adder
        DoubleAdder sum = new DoubleAdder();
        double wdx;
        double wdy;
        double len;
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
     * @throws java.lang.InterruptedException
     */
    public double[] directionFrequencyPolarization(Vector n, Vector v, double e) throws InterruptedException {
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
    public abstract double directionFrequencyFluxNoSpread(Vector n, Vector v, double e);

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
    public abstract double[] directionFrequencyPolarizationNoSpread(Vector n, Vector v, double e);

    /**
     * A method calculating the flux density in a given direction for a given
     * X-ray photon energy taking into account electron transversal pulse spread
     *
     * @param n direction
     * @param v0 normalized electron velocity
     * @param e X-ray energy
     * @return
     */
    public double directionFrequencyFluxSpread(Vector n, Vector v0, double e) {
        BaseAbstractUnivariateIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateFunction func = new UnivariateFrequencyFluxSpreadOuter(e, v0, n);
        try {
            return integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, 0.0, 2 * Math.PI);
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
    }

    /**
     * A multi-threaded method calculating the full polarization tensor density
     * in a given direction for a given X-ray photon energy taking into account
     * electron transversal pulse spread
     *
     * @param n direction
     * @param v0 normalized electron velocity
     * @param e X-ray energy
     * @return
     * @throws java.lang.InterruptedException
     */
    public double[] directionFrequencyPolarizationSpread(final Vector n, final Vector v0, final double e) throws InterruptedException {
        //An array for results
        double[] array = new double[NUMBER_OF_POL_PARAM];
        //Creating a latch for threads
        CountDownLatch lt = new CountDownLatch(NUMBER_OF_POL_PARAM);
        //Creating a pool of threads for calculations
        ExecutorService execs = Executors.newFixedThreadPool(threadNumber);
        //Calculating the polarization tensor elements
        for (int i = 0; i < NUMBER_OF_POL_PARAM; i++) {
            int[] ia = new int[]{i};
            execs.execute(() -> {
                UnivariateFunction func = new UnivariateFrequencyPolarizationSpreadOuter(e, v0, n, ia[0]);
                try {
                    //Creating a separate inegrator for each thread
                    array[ia[0]] = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT).integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, 0.0, 2 * Math.PI);
                } catch (TooManyEvaluationsException ex) {
                    array[ia[0]] = 0;
                }
                lt.countDown();
            });
        }
        //Waiting for an interruption and shuting down threads if interrupted
        try {
            lt.await();
        } catch (InterruptedException ex) {
            execs.shutdownNow();
            throw ex;
        }
        execs.shutdownNow();
        return array;
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
     * @throws java.lang.InterruptedException
     */
    public double[] directionFrequencyVolumePolarization(Vector r, Vector n, Vector v, double e) throws InterruptedException {
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
        for (int i = 0; i < NUMBER_OF_POL_PARAM; i++) {
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
     * @throws java.lang.InterruptedException
     */
    public double[] directionFrequencyVolumePolarizationSpread(Vector r, Vector n, Vector v, double e) throws InterruptedException {
        double[] stocks = directionFrequencyPolarizationSpread(n, v, e);
        double vFlux = volumeFlux(r);
        for (int i = 0; i < NUMBER_OF_POL_PARAM; i++) {
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
        double u;
        double z;
        double z1;
        double x;
        double x1;
        double y;
        double y1;
        double sn;
        double cs;
        double K;
        double len;
        len = Math.sqrt(lp.getLength() * lp.getLength() + eb.getLength() * eb.getLength());
        sn = lp.getDirection().get(1);
        cs = lp.getDirection().get(2);
        x = r.get(0);
        y = r.get(1);
        z = r.get(2);
        x1 = x;
        y1 = -sn * z + cs * y;
        z1 = cs * z + sn * y;
        K = Math.pow((z + z1 - eb.getShift().get(2) - lp.getDelay()) / len, 2);
        u = 2.0 * Math.sqrt(Math.PI) * Math.sqrt((lp.getWidth2(0.0) + eb.getxWidth2(0.0)) * (lp.getWidth2(0.0) + eb.getyWidth2(0.0))) / len * Math.exp(-K) * eb.tSpatialDistribution(r) * lp.tSpatialDistribution(new BasicVector(new double[]{x1, y1, z1}));
        return new Double(u).isNaN() ? 0 : u;
    }

    /**
     * A method giving the flux density in a given direction
     *
     * @param n direction
     * @param v normalized electron velocity
     * @return
     */
    public abstract double directionFlux(Vector n, Vector v);

    /**
     * A method calculating X-ray energy in a given direction
     *
     * @param n direction
     * @param v normalized electron velocity
     * @return
     */
    public abstract double directionEnergy(Vector n, Vector v);

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
     * @throws java.lang.InterruptedException
     */
    public double[] directionFrequencyPolarizationBrilliance(Vector r0, Vector n, Vector v, double e) throws InterruptedException {
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
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateVolumeFlux func = new UnivariateVolumeFlux(r0, n);
        try {
            u = integrator.integrate(30000, func, r0.fold(Vectors.mkEuclideanNormAccumulator()) - 3 * eb.getLength(), r0.fold(Vectors.mkEuclideanNormAccumulator()) + 3 * eb.getLength());
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
        double[] array = new double[NUMBER_OF_POL_PARAM];
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateVolumeFlux func = new UnivariateVolumeFlux(r0, n);
        try {
            mlt = integrator.integrate(30000, func, r0.fold(Vectors.mkEuclideanNormAccumulator()) - 3 * eb.getLength(), r0.fold(Vectors.mkEuclideanNormAccumulator()) + 3 * eb.getLength());
        } catch (TooManyEvaluationsException ex) {
            mlt = 0;
        }
        for (int i = 0; i < NUMBER_OF_POL_PARAM; i++) {
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
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateVolumeFlux func = new UnivariateVolumeFlux(r0, n);
        try {
            u = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, r0.fold(Vectors.mkEuclideanNormAccumulator()) - 3 * eb.getLength(), r0.fold(Vectors.mkEuclideanNormAccumulator()) + 3 * eb.getLength());
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
     * @throws java.lang.InterruptedException
     */
    public double[] directionFrequencyBrilliancePolarizationSpread(Vector r0, Vector n, Vector v, double e) throws InterruptedException {
        double mlt;
        double[] array = new double[NUMBER_OF_POL_PARAM];
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateVolumeFlux func = new UnivariateVolumeFlux(r0, n);
        try {
            mlt = integrator.integrate(30000, func, r0.fold(Vectors.mkEuclideanNormAccumulator()) - 3 * eb.getLength(), r0.fold(Vectors.mkEuclideanNormAccumulator()) + 3 * eb.getLength());
        } catch (TooManyEvaluationsException ex) {
            mlt = 0;
        }
        for (int i = 0; i < NUMBER_OF_POL_PARAM; i++) {
            array[i] = mlt * directionFrequencyPolarizationSpread(n, v, e)[i];
        }
        return array;
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
        Vector n0 = new BasicVector(new double[]{0.0, 1.0, 0.0});
        Vector As;
        double prob0;
        double prob;
        double EMax;
        double mult = 2;
        double factor;
        double sum = 0;
        double[] pol;
        double[] polParam;
        EMax = directionEnergy(n, n);
        factor = 64 * Math.max(eb.getxWidth(0.0), lp.getWidth(0.0)) * Math.max(eb.getyWidth(0.0), lp.getWidth(0.0)) * Math.max(eb.getLength(), lp.getLength()) * 4 * rayXAnglerange * rayYAnglerange * (maxEnergy - minEnergy);
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
                if (ksi == null) {
                    polParam = directionFrequencyVolumePolarizationNoSpread(r, n, v, ray[10]);
                } else {
                    polParam = new double[]{directionFrequencyVolumeFluxNoSpread(r, n, v, ray[10]), ksi[0], ksi[1], ksi[2]};
                }
                prob = polParam[0] * eb.angleDistribution(thetax, thetay);
            } else {
                if (ksi == null) {
                    polParam = directionFrequencyVolumePolarizationNoSpread(r, n, new BasicVector(new double[]{0.0, 0.0, 1.0}), ray[10]);
                } else {
                    polParam = new double[]{directionFrequencyVolumeFluxNoSpread(r, n, new BasicVector(new double[]{0.0, 0.0, 1.0}), ray[10]), ksi[0], ksi[1], ksi[2]};
                }
                prob = polParam[0];
            }
            if (!new Double(prob).isNaN()) {
                sum += prob / ray[10];
            }
            counter.incrementAndGet();
        } while (prob / prob0 < Math.random() || (new Double(prob)).isNaN());
        // Calculating the rotated polarization vector and getting the full polarizaation state
        n = new BasicVector(new double[]{ray[3], ray[4], ray[5]});
        T = getTransform(n, n0);
        //Checking if polarization is pre-specified
        pol = (ksi != null) ? getPolarization(ksi) : getPolarization(new double[]{polParam[1] / polParam[0], polParam[2] / polParam[0], polParam[3] / polParam[0]});
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
        Matrix D;
        Matrix A;
        Matrix I = new Basic1DMatrix(3, 3);
        I.set(0, 0, 1.0);
        I.set(1, 1, 1.0);
        I.set(2, 2, 1.0);
        double innerProduct;
        innerProduct = n.innerProduct(n0);
        D = n.outerProduct(n0).add(n0.outerProduct(n)).multiply(innerProduct).subtract(n.outerProduct(n).add(n0.outerProduct(n0))).divide(innerProduct * innerProduct - 1.0);
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
     * Approximate geometric factor calculated neglecting the "hourglass"
     * effect. It assumes values from 0 to 1.
     *
     * @return the approximate geometricFactor
     */
    public double getApproxGeometricFactor() {
        double cs2 = (1 + lp.getDirection().innerProduct(new BasicVector(new double[]{0, 0, 1}))) / 2;
        double sn2 = (1 - lp.getDirection().innerProduct(new BasicVector(new double[]{0, 0, 1}))) / 2;
        double w2 = lp.getWidth2(0) + eb.getWidth2(0);
        double l2 = lp.getLength() * lp.getLength() + eb.getLength() * eb.getLength();
        return Math.sqrt(w2 / (l2 * sn2 + w2 * cs2) / cs2);
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
     * @param ksiVector
     * @return
     */
    protected double[] getPolarization(double[] ksiVector) {
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
    public static class LocalException extends RuntimeException {

        // The x value that caused the problem.
        private final double x;

        /**
         * A class for exceptions in custom functions
         *
         * @param x
         */
        public LocalException(double x) {
            this.x = x;
        }

        public double getX() {
            return x;
        }
    }
}
