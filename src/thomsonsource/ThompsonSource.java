/*
 * Copyright (C) 2023 Ruslan Feshchenko
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
import java.util.concurrent.atomic.*;
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
 * @version 2.81
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
        this.montecarlocounter = new AtomicLong();
        this.partialFlux = new DoubleAdder();
        calculateTotalFlux();
        calculateGeometricFactor();
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
     * A numerical shift to improve convergence of polarization angular
     * integrals
     */
    public static final double SHIFT = 1e10;

    /**
     * Angle range for rays exported for Shadow in the X-direction
     */
    private double rayXAnglerange = 0.0005;

    /**
     * Angle range for rays exported for Shadow in the Y-direction
     */
    private double rayYAnglerange = 0.0005;

    /**
     * Min ray energy
     */
    private double minEnergy = 25 * 1e3 * ElectronBunch.E;

    /**
     * Max ray energy
     */
    private double maxEnergy = 35 * 1e3 * ElectronBunch.E;

    /**
     * Number of points in Monte Carlo calculation of the geometric factor
     */
    private int npGeometricFactor = 50000;

    /**
     * Number of points in Monte Carlo calculation of the emittance averaging
     */
    private int npEmittance = 30000;

    /**
     * Maximal number of evaluations in calculations of the brilliance and
     * polarization
     */
    public static final int MAXIMAL_NUMBER_OF_EVALUATIONS = 1000000;

    /**
     * A shift factor to improve numerical integral convergence in polarization
     * calculations
     */
    private double shiftfactor = 1;
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
     * Flag - whether or not the Monte-Carlo method is used to calculate the
     * directional integral
     */
    private boolean IsMonteCarlo = true;

    /**
     * Flag - whether or not the Compton cross-sections are used
     */
    private boolean IsCompton = true;

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
    private AtomicLong montecarlocounter;

    private LaserPulse lp;
    private ElectronBunch eb;

    private double[] ksi = null;

    @Override
    public Object clone() throws CloneNotSupportedException {
        Object tm = super.clone();
        ((ThompsonSource) tm).eb = (ElectronBunch) this.eb.clone();
        ((ThompsonSource) tm).lp = (LaserPulse) this.lp.clone();
        ((ThompsonSource) tm).montecarlocounter = new AtomicLong();
        ((ThompsonSource) tm).partialFlux = new DoubleAdder();
        if (ksi != null) {
            ((ThompsonSource) tm).ksi = (double[]) ksi.clone();
        }
        return tm;
    }

    /**
     * Returns electron bunch reference
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
    public final void calculateTotalFlux() {
        this.totalFlux = SIGMA_T * eb.getNumber() * lp.getPhotonNumber()
                * lp.getFq() / Math.PI / Math.sqrt((lp.getWidth2(0.0) + eb.getxWidth2(0.0))
                * (lp.getWidth2(0.0) + eb.getyWidth2(0.0)));
    }

    /**
     * A method calculating normalized total flux within a certain angle
     *
     * @param maxAngle
     * @return
     */
    public final double calculateAngleTotalFlux(double maxAngle) {
        double gamma2 = eb.getGamma() * eb.getGamma();
        double v = Math.sqrt(1 - 1 / gamma2);
        double v2 = v * v;
        double cs = Math.cos(maxAngle);
        if (!IsCompton) {
            return 3.0 / 4 / gamma2 * ((1 - cs) / (1 - v * cs) / (1 - v)
                    * (5.0 / 6 + 1.0 / 6 / v2 - 1.0 / 6 / gamma2 / v2 * (1 - v2 * cs) / (1 - v) / (1 - v * cs))
                    + 1.0 / 6 / gamma2 / v2 * (1 - cs * cs) / Math.pow((1 - v * cs), 3))
                    * getTotalFlux() * getGeometricFactor();
        } else {
            return 3.0 / 4 / gamma2 * ((1 - cs) / (1 - v * cs) / (1 - v)
                    * (5.0 / 6 + 1.0 / 6 / v2 - 1.0 / 6 / gamma2 / v2 * (1 - v2 * cs) / (1 - v) / (1 - v * cs))
                    + 1.0 / 6 / gamma2 / v2 * (1 - cs * cs) / Math.pow((1 - v * cs), 3))
                    * getTotalFlux() * getGeometricFactor();
        }
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
        double wdx, wdy, len;
        int mult = 3;
        wdx = mult * (eb.getxWidth(0.0) * lp.getWidth(0.0) / Math.sqrt(eb.getxWidth2(0.0) + lp.getWidth2(0.0)));
        wdy = mult * (eb.getyWidth(0.0) * lp.getWidth(0.0) / Math.sqrt(eb.getyWidth2(0.0) + lp.getWidth2(0.0)));
        len = mult * eb.getLength() * lp.getLength() / Math.sqrt(eb.getLength() * eb.getLength() + lp.getLength() * lp.getLength());
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
        return iseSpread() ? directionFrequencyFluxSpreadIntegral(n, v, e) : directionFrequencyFluxNoSpread(n, v, e);
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
        return iseSpread() ? (isIsMonteCarlo() ? directionFrequencyPolarizationSpreadMonteCarlo(n, v, e) : directionFrequencyPolarizationSpreadIntegral(n, v, e)) : directionFrequencyPolarizationNoSpread(n, v, e);
    }

    /**
     * A method calculating the flux density in a given direction for a given
     * X-ray photon energy without taking into account electron transversal
     * momentum spread
     *
     * @param n observation direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     */
    public double directionFrequencyFluxNoSpread(Vector n, Vector v, double e) {
        double K, th, res, ac, koef, gamma, gamma2;
        th = (1 - n.innerProduct(v)) * 2;
        if (!IsCompton) {
            K = Math.pow((Math.sqrt(e / lp.getPhotonEnergy() / (1 - e * th / lp.getPhotonEnergy() / 4)) - 2 * eb.getGamma()), 2)
                    / 4 / Math.pow(eb.getGamma() * eb.getDelgamma(), 2);
            res = getTotalFlux() * e * 3.0 / 64 / Math.PI / Math.sqrt(Math.PI) / eb.getDelgamma() / eb.getGamma() / lp.getPhotonEnergy()
                    * Math.sqrt(e / lp.getPhotonEnergy()) * (Math.pow((1 - e * th / lp.getPhotonEnergy() / 2), 2) + 1)
                    / Math.sqrt(1 - e * th / lp.getPhotonEnergy() / 4) * Math.exp(-K);
        } else {
            ac = lp.getPhotonEnergy() / (ElectronBunch.MC2 * ElectronBunch.E * 1e6);
            koef = 4 * lp.getPhotonEnergy() / e - th;
            gamma = (2 * ac + Math.sqrt(4 * ac * ac + koef)) / koef;
            gamma2 = gamma * gamma;
            res = getTotalFlux() * e * 1.5 / Math.pow(Math.PI, 1.5) / eb.getDelgamma() / eb.getGamma() * lp.getPhotonEnergy() / Math.pow(e, 2)
                    * Math.pow(eb.getGamma(), 5) / Math.pow((1 + gamma2 * th + 4 * ac * eb.getGamma()), 2) / (1 + 2 * gamma * ac)
                    * (1 + Math.pow((1 - gamma2 * th) / (1 + gamma2 * th), 2)) * Math.exp(-Math.pow((gamma - eb.getGamma()) / eb.getDelgamma() / eb.getGamma(), 2));
        }
        return new Double(res).isNaN() ? 0 : res;
    }

    /**
     * A method calculating the Stocks parameters density in a given direction
     * for a given X-ray photon energy without taking into account electron
     * transversal momentum spread
     *
     * @param n observation direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     */
    public double[] directionFrequencyPolarizationNoSpread(Vector n, Vector v, double e) {
        double[] res = new double[NUMBER_OF_POL_PARAM];
        double K, th, m11, m22, m12, mlt, cs, sn, ac, koef, gamma, gamma2;
        th = (1 - n.innerProduct(v)) * 2;
        if (!IsCompton) {
            K = Math.pow((Math.sqrt(e / lp.getPhotonEnergy() / (1 - e * th / lp.getPhotonEnergy() / 4)) - 2 * eb.getGamma()), 2)
                    / 4 / Math.pow(eb.getGamma() * eb.getDelgamma(), 2);
            m11 = getTotalFlux() * e * 3.0 / 32 / Math.PI / Math.sqrt(Math.PI) / eb.getDelgamma() / eb.getGamma() / lp.getPhotonEnergy()
                    * Math.sqrt(e / lp.getPhotonEnergy()) / Math.sqrt(1 - e * th / lp.getPhotonEnergy() / 4) * Math.exp(-K);
            mlt = 1 - e * th / lp.getPhotonEnergy() / 2;
        } else {
            ac = lp.getPhotonEnergy() / (ElectronBunch.MC2 * ElectronBunch.E * 1e6);
            koef = 4 * lp.getPhotonEnergy() / e - th;
            gamma = (2 * ac + Math.sqrt(4 * ac * ac + koef)) / koef;
            gamma2 = gamma * gamma;
            m11 = getTotalFlux() * e * 3.0 / Math.pow(Math.PI, 1.5) / eb.getDelgamma() / eb.getGamma() * lp.getPhotonEnergy() / Math.pow(e, 2)
                    * Math.pow(eb.getGamma(), 5) / Math.pow((1 + gamma2 * th + 4 * ac * eb.getGamma()), 2) / (1 + 2 * gamma * ac)
                    * Math.exp(-Math.pow((gamma - eb.getGamma()) / eb.getDelgamma() / eb.getGamma(), 2));
            mlt = (1 - gamma2 * th) / (1 + gamma2 * th);
        }
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
        res[0] = (m11 + m22 - (cs2 * lp.getPolarization()[2] + sn2 * lp.getPolarization()[0]) * (m11 - m22)) / 2;
        res[3] = (cs2 * (m22 - m11) + lp.getPolarization()[2] * (cs2cs2 * (m11 + m22) + 2 * sn2sn2 * m12)
                + lp.getPolarization()[0] * cs2sn2 * (m11 + m22 - 2 * m12)) / 2;
        res[1] = (sn2 * (m22 - m11) + lp.getPolarization()[0] * (sn2sn2 * (m11 + m22) + 2 * cs2cs2 * m12)
                + lp.getPolarization()[2] * cs2sn2 * (m11 + m22 - 2 * m12)) / 2;
        res[2] = lp.getPolarization()[1] * m12;

        // If the intensity is NaN, zero or less than zero then set all Stocks intensities to zero
        // If a Stocks intensity is NaN then set it to zero
        for (int i = 0; i < NUMBER_OF_POL_PARAM; i++) {
            if (new Double(res[i]).isNaN() || new Double(res[0]).isNaN() || res[0] <= 0) {
                res[i] = 0;
            }
        }
        return res;
    }

    /**
     * A method calculating the flux density in a given direction for a given
     * X-ray photon energy taking into account electron transversal momentum
     * spread and using the conventional integration
     *
     * @param n observation direction
     * @param v0 normalized mean electron velocity
     * @param e X-ray energy
     * @return
     */
    public double directionFrequencyFluxSpreadIntegral(Vector n, Vector v0, double e) {
        BaseAbstractUnivariateIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        double res;
        UnivariateFunction func
                = new UnivariateFrequencyFluxSpreadOuter(e, v0, n);
        try {
            res = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, -INT_RANGE * eb.getXSpread(), INT_RANGE * eb.getXSpread())
                    - getShiftfactor() * SHIFT * 4 * INT_RANGE * INT_RANGE * eb.getXSpread() * eb.getYSpread();
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
        return new Double(res).isNaN() ? 0 : res;
    }

    /**
     * A method calculating the flux density in a given direction for a given
     * X-ray photon energy taking into account electron transversal momentum
     * spread and using the Monte-Carlo method
     *
     * @param n observation direction
     * @param v0 normalized mean electron velocity
     * @param e X-ray energy
     * @return
     */
    public double directionFrequencyFluxSpreadMonteCarlo(Vector n, Vector v0, double e) {
        ExecutorService execs = Executors.newFixedThreadPool(threadNumber);
        // We need to synchronize threads
        CountDownLatch lt = new CountDownLatch(threadNumber);
        // Atomic adder
        DoubleAdder sum = new DoubleAdder();
        double res;
        final int itNumber = Math.round(getNpEmittance() / threadNumber);

        // Splitting the job into a number of threads
        for (int m = 0; m < threadNumber; m++) {
            execs.execute(() -> {
                double rx, ry, tm, psum = 0;
                Vector dv, v = new BasicVector(new double[]{0.0, 0.0, 0.0});
                //Calculating a partial sum
                for (int i = 0; i < itNumber; i++) {
                    if (Thread.currentThread().isInterrupted()) {
                        return;
                    }
                    rx = (2 * Math.random() - 1) * INT_RANGE * eb.getXSpread();
                    ry = (2 * Math.random() - 1) * INT_RANGE * eb.getYSpread();
                    v.set(0, rx);
                    v.set(1, ry);
                    v.set(2, Math.sqrt(1 - rx * rx - ry * ry));
                    dv = v.subtract(v0);
                    tm = directionFrequencyFluxNoSpread(n, v, e) * eb.angleDistribution(dv.get(0), dv.get(1));
                    psum += new Double(tm).isNaN() ? 0 : tm;
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
        // Outputting the final result
        res = 4 * INT_RANGE * INT_RANGE * eb.getXSpread() * eb.getYSpread() * sum.sum() / itNumber / threadNumber;
        return new Double(res).isNaN() ? 0 : res;
    }

    /**
     * A multi-threaded method calculating the full polarization tensor density
     * in a given direction for a given X-ray photon energy taking into account
     * electron transversal momentum spread using the conventional integration
     *
     * @param n observation direction
     * @param v0 normalized mean electron velocity
     * @param e X-ray energy
     * @return
     * @throws java.lang.InterruptedException
     */
    public double[] directionFrequencyPolarizationSpreadIntegral(final Vector n, final Vector v0, final double e) throws InterruptedException {
        //An res for results
        double[] res = new double[NUMBER_OF_POL_PARAM];
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
                    res[ia[0]] = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                            RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT).
                            integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, -INT_RANGE * eb.getXSpread(), INT_RANGE * eb.getXSpread())
                            - getShiftfactor() * SHIFT * 4 * INT_RANGE * INT_RANGE * eb.getXSpread() * eb.getYSpread();
                } catch (TooManyEvaluationsException ex) {
                    res[ia[0]] = 0;
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

        // If the intensity is NaN, zero or less than zero then all Stocks intensities to zero
        // If a Stocks intensity is NaN then set it to zero
        for (int i = 0; i < NUMBER_OF_POL_PARAM; i++) {
            if (new Double(res[i]).isNaN() || new Double(res[0]).isNaN() || res[0] <= 0) {
                res[i] = 0;
            }
        }
        return res;
    }

    /**
     * A multi-threaded method calculating the full polarization tensor density
     * in a given direction for a given X-ray photon energy taking into account
     * electron transversal momentum spread using the Monte-Carlo method
     *
     * @param n observation direction
     * @param v0 normalized mean electron velocity
     * @param e X-ray energy
     * @return
     * @throws java.lang.InterruptedException
     */
    public double[] directionFrequencyPolarizationSpreadMonteCarlo(final Vector n, final Vector v0, final double e) throws InterruptedException {
        //An res for results
        double[] res = new double[NUMBER_OF_POL_PARAM];
        //Creating a pool of threads for calculations
        ExecutorService execs = Executors.newFixedThreadPool(threadNumber);
        //The number of threads used to calculate Stocks parameters
        final int itNumber = Math.round(getNpEmittance() / threadNumber);
        //Calculating the polarization tensor elements
        for (int i = 0; i < NUMBER_OF_POL_PARAM; i++) {
            int[] ia = new int[]{i};
            // An atomic adder
            DoubleAdder sum = new DoubleAdder();
            //Creating a latch for the threads
            CountDownLatch lt = new CountDownLatch(threadNumber);
            for (int m = 0; m < threadNumber; m++) {
                execs.execute(() -> {
                    double rx, ry, rth, tm, psum = 0;
                    Vector dv, v = new BasicVector(new double[]{0.0, 0.0, 0.0});
                    //Calculating a partial sum
                    for (int p = 0; p < itNumber; p++) {
                        if (Thread.currentThread().isInterrupted()) {
                            return;
                        }
                        rx = (2 * Math.random() - 1) * INT_RANGE * eb.getXSpread();
                        ry = (2 * Math.random() - 1) * INT_RANGE * eb.getYSpread();
                        v.set(0, rx);
                        v.set(1, ry);
                        v.set(2, Math.sqrt(1 - rx * rx - ry * ry));
                        dv = v.subtract(v0);
                        tm = directionFrequencyPolarizationNoSpread(n, v, e)[ia[0]] * eb.angleDistribution(dv.get(0), dv.get(1));
                        psum += new Double(tm).isNaN() ? 0 : tm;
                    }
                    //Adding to the full sum
                    sum.add(psum);
                    //Counting down the latch
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
            //Outputting the result for the i-th Stocks intensity
            res[i] = 4 * INT_RANGE * INT_RANGE * eb.getXSpread() * eb.getYSpread() * sum.sum() / itNumber / threadNumber;
        }
        //Shutting down the execution services
        execs.shutdownNow();

        // If the intensity is NaN, zero or less than zero then all Stocks intensities to zero
        // If a Stocks intensity is NaN then set it to zero
        for (int i = 0; i < NUMBER_OF_POL_PARAM; i++) {
            if (new Double(res[i]).isNaN() || new Double(res[0]).isNaN() || res[0] <= 0) {
                res[i] = 0;
            }
        }
        return res;
    }

    /**
     * Getting a numerical factor to improve integral convergence
     *
     * @return the factor
     */
    public double getShiftfactor() {
        return shiftfactor;
    }

    /**
     * Setting a numerical factor to improve integral convergence
     *
     * @param shift the factor to set
     */
    public void setShiftfactor(double shift) {
        this.shiftfactor = shift;
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
        public double value(double x) {
            double tmp;
            if (Thread.currentThread().isInterrupted()) {
                Thread.currentThread().interrupt();
                return 0;
            }
            UnivariateFunction func
                    = new UnivariateFrequencyFluxSpreadInner(x, e, v0, n);
            try {
                tmp = inergrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, -INT_RANGE * eb.getYSpread(), INT_RANGE * eb.getYSpread());
            } catch (TooManyEvaluationsException ex) {
                return 0;
            }
            //If the integral is NaN then set it as zero
            return new Double(tmp).isNaN() ? 0 : tmp;
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
        public double value(double x) {
            double tmp;
            if (Thread.currentThread().isInterrupted()) {
                Thread.currentThread().interrupt();
                return 0;
            }
            UnivariateFunction func
                    = new UnivariateFrequencyPolarizationSpreadInner(x, e, v0, n, index);
            try {
                tmp = inergrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, -INT_RANGE * eb.getYSpread(), INT_RANGE * eb.getYSpread());
            } catch (TooManyEvaluationsException ex) {
                return 0;
            }
            //If the integral is NaN then set it as zero
            return new Double(tmp).isNaN() ? 0 : tmp;
        }
    }

    /**
     * An auxiliary class for Romberg integrator for flux calculations
     */
    private class UnivariateFrequencyFluxSpreadInner implements UnivariateFunction {

        private final double x, e;
        private final Vector n, v0;

        public UnivariateFrequencyFluxSpreadInner(double x, double e, Vector v0, Vector n) {
            this.x = x;
            this.e = e;
            this.n = n;
            this.v0 = v0;
        }

        @Override
        public double value(double y) {
            if (Thread.currentThread().isInterrupted()) {
                Thread.currentThread().interrupt();
                return 0;
            }
            double u;
            Vector v = new BasicVector(new double[]{x, y, Math.sqrt(1 - y * y - x * x)});
            Vector dv = v.subtract(v0);
            u = directionFrequencyFluxNoSpread(n, v, e) * eb.angleDistribution(dv.get(0), dv.get(1))
                    + getShiftfactor() * SHIFT;
            return new Double(u).isNaN() ? 0 : u;
        }
    }

    /**
     * An auxiliary class for Romberg integrator for polarization calculations
     */
    private class UnivariateFrequencyPolarizationSpreadInner implements UnivariateFunction {

        private final double x, e;
        private final int index;
        private final Vector n, v0;

        public UnivariateFrequencyPolarizationSpreadInner(double x, double e, Vector v0, Vector n, int index) {
            this.x = x;
            this.e = e;
            this.n = n;
            this.index = index;
            this.v0 = v0;
        }

        @Override
        public double value(double y) {
            if (Thread.currentThread().isInterrupted()) {
                Thread.currentThread().interrupt();
                return 0;
            }
            double u;
            Vector v = new BasicVector(new double[]{x, y, Math.sqrt(1 - y * y - x * x)});
            Vector dv = v.subtract(v0);
            // Normalization to the peak value
            u = directionFrequencyPolarizationNoSpread(n, v, e)[index] * eb.angleDistribution(dv.get(0), dv.get(1))
                    + getShiftfactor() * SHIFT;
            return new Double(u).isNaN() ? 0 : u;
        }
    }

    /**
     * A method calculating the flux density in a given direction for a given
     * X-ray photon energy for a given volume element
     *
     * @param r spatial position
     * @param n observation direction
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
     * @param n observation direction
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
     * account electron transversal momentum spread
     *
     * @param r spatial position
     * @param n observation direction
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
     * into account electron transversal momentum spread
     *
     * @param r spatial position
     * @param n observation direction
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
     * @param n observation direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     */
    public double directionFrequencyVolumeFluxSpread(Vector r, Vector n, Vector v, double e) {
        return isIsMonteCarlo() ? directionFrequencyFluxSpreadMonteCarlo(n, v, e) * volumeFlux(r) : directionFrequencyFluxSpreadIntegral(n, v, e) * volumeFlux(r);
    }

    /**
     * A method calculating the Stocks parameters density in a given direction
     * for a given X-ray photon energy for a given volume element taking into
     * account electron transversal pulse spread
     *
     * @param r spatial position
     * @param n observation direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     * @throws java.lang.InterruptedException
     */
    public double[] directionFrequencyVolumePolarizationSpread(Vector r, Vector n, Vector v, double e) throws InterruptedException {
        double[] stocks = isIsMonteCarlo() ? directionFrequencyPolarizationSpreadMonteCarlo(n, v, e) : directionFrequencyPolarizationSpreadIntegral(n, v, e);
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
     * @param n observation direction
     * @param v normalized electron velocity
     * @return
     */
    public double directionFlux(Vector n, Vector v) {
        double gamma2, th;
        th = (1 - n.innerProduct(v)) * 2;
        gamma2 = eb.getGamma() * eb.getGamma();
        if (!IsCompton) {
            return getTotalFlux() * 3.0 / 2 / Math.PI * gamma2 * (1 + Math.pow(th * gamma2, 2))
                    / Math.pow((1 + gamma2 * th), 4) * getGeometricFactor();
        } else {
            return getTotalFlux() * 3.0 / 2 / Math.PI * gamma2 * (1 + Math.pow(th * gamma2, 2))
                    / Math.pow((1 + gamma2 * th), 2)
                    / Math.pow((1 + gamma2 * th + 4 * lp.getPhotonEnergy() / (ElectronBunch.MC2 * ElectronBunch.E * 1e6) * eb.getGamma()), 2) * getGeometricFactor();
        }
    }

    /**
     * A method calculating X-ray energy in a given direction
     *
     * @param n observation direction
     * @param v normalized electron velocity
     * @return
     */
    public double directionEnergy(Vector n, Vector v) {
        double mv, cs;
        mv = Math.sqrt(1.0 - 1.0 / eb.getGamma() / eb.getGamma());
        cs = n.innerProduct(v);
        if (!IsCompton) {
            return (1 + mv) * lp.getPhotonEnergy() / (1 - cs * mv);
        } else {
            return (1 + mv) * lp.getPhotonEnergy() / (1 - cs * mv + lp.getPhotonEnergy() / (ElectronBunch.MC2 * ElectronBunch.E * 1e6) / eb.getGamma() * (1 + cs));
        }
    }

    /**
     * A method calculating spectral brilliance in a given direction
     *
     * @param r0 spatial position for brightness
     * @param n observation direction
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
     * @param n observation direction
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
     * taking into account electron transversal momentum spread
     *
     * @param r0 spatial position for brightness
     * @param n observation direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     */
    public double directionFrequencyBrillianceNoSpread(Vector r0, Vector n, Vector v, double e) {
        double u;
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateVolumeFlux func = new UnivariateVolumeFlux(r0, n);
        //Semiwidth of the integration interval
        double semiwidth = INT_RANGE * lp.getLength() * lp.getWidth(0)
                / Math.sqrt(Math.pow(lp.getLength() * n.get(0), 2) + lp.getWidth2(0) * Math.pow(n.get(2), 2)) / 2;

        try {
            u = integrator.integrate(30000, func,
                    r0.fold(Vectors.mkEuclideanNormAccumulator()) - semiwidth,
                    r0.fold(Vectors.mkEuclideanNormAccumulator()) + semiwidth);

            return new Double(u).isNaN() ? 0 : u * directionFrequencyFluxNoSpread(n, v, e);
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
    }

    /**
     * A method calculating spectral brilliance in a given direction without
     * taking into account electron transversal momentum spread but with
     * polarization
     *
     * @param r0 spatial position for brightness
     * @param n observation direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     */
    public double[] directionFrequencyBrilliancePolarizationNoSpread(Vector r0, Vector n, Vector v, double e) {
        double mlt;
        double[] array = new double[NUMBER_OF_POL_PARAM];
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateVolumeFlux func = new UnivariateVolumeFlux(r0, n);
        //Semiwidth of the integration interval
        double semiwidth = INT_RANGE * lp.getLength() * lp.getWidth(0)
                / Math.sqrt(Math.pow(lp.getLength() * n.get(0), 2) + lp.getWidth2(0) * Math.pow(n.get(2), 2)) / 2;

        try {
            mlt = integrator.integrate(30000, func,
                    r0.fold(Vectors.mkEuclideanNormAccumulator()) - semiwidth,
                    r0.fold(Vectors.mkEuclideanNormAccumulator()) + semiwidth);
            mlt = new Double(mlt).isNaN() ? 0 : mlt;
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
     * account electron transversal momentum spread
     *
     * @param r0 spatial position for brightness
     * @param n observation direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     */
    public double directionFrequencyBrillianceSpread(Vector r0, Vector n, Vector v, double e) {
        double u;
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateVolumeFlux func = new UnivariateVolumeFlux(r0, n);
        //Semiwidth of the integration interval
        double semiwidth = INT_RANGE * lp.getLength() * lp.getWidth(0)
                / Math.sqrt(Math.pow(lp.getLength() * n.get(0), 2) + lp.getWidth2(0) * Math.pow(n.get(2), 2)) / 2;

        try {
            u = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func,
                    r0.fold(Vectors.mkEuclideanNormAccumulator()) - semiwidth,
                    r0.fold(Vectors.mkEuclideanNormAccumulator()) + semiwidth);

            return new Double(u).isNaN() ? 0 : u * (isIsMonteCarlo() ? directionFrequencyFluxSpreadMonteCarlo(n, v, e) : directionFrequencyFluxSpreadIntegral(n, v, e));
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
    }

    /**
     * A method calculating spectral brilliance in a given direction taking into
     * account the electron transversal momentum spread and with polarization
     *
     * @param r0 spatial position for brightness
     * @param n observation direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     * @throws java.lang.InterruptedException
     */
    public double[] directionFrequencyBrilliancePolarizationSpread(Vector r0, Vector n, Vector v, double e) throws InterruptedException {
        double mlt;
        double[] array = new double[NUMBER_OF_POL_PARAM];
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateVolumeFlux func = new UnivariateVolumeFlux(r0, n);
        //Semiwidth of the integration interval
        double semiwidth = INT_RANGE * lp.getLength() * lp.getWidth(0)
                / Math.sqrt(Math.pow(lp.getLength() * n.get(0), 2) + lp.getWidth2(0) * Math.pow(n.get(2), 2)) / 2;

        try {
            mlt = integrator.integrate(30000, func,
                    r0.fold(Vectors.mkEuclideanNormAccumulator()) - semiwidth,
                    r0.fold(Vectors.mkEuclideanNormAccumulator()) + semiwidth);
            mlt = new Double(mlt).isNaN() ? 0 : mlt;
        } catch (TooManyEvaluationsException ex) {
            mlt = 0;
        }
        for (int i = 0; i < NUMBER_OF_POL_PARAM; i++) {
            array[i] = mlt * (isIsMonteCarlo() ? directionFrequencyPolarizationSpreadMonteCarlo(n, v, e)[i] : directionFrequencyPolarizationSpreadIntegral(n, v, e)[i]);
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
            if (Thread.currentThread().isInterrupted()) {
                Thread.currentThread().interrupt();
                return 0;
            }
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
     * @return an res with ray parameters
     * @throws java.lang.InterruptedException
     */
    public double[] getRay() throws InterruptedException {
        double[] ray = new double[NUMBER_OF_COLUMNS];
        Matrix T;
        Vector n = new BasicVector(new double[]{0.0, 0.0, 1.0});
        Vector r = new BasicVector(new double[]{0.0, 0.0, 0.0});
        Vector n0 = new BasicVector(new double[]{0.0, 1.0, 0.0}), As;
        double prob0, prob, EMax, MULT = 2, factor, sum = 0;
        int lcn = 0;
        double[] pol, polParam;
        EMax = directionEnergy(n, n);
        factor = 32 * MULT * MULT * MULT * Math.max(eb.getxWidth(0.0), lp.getWidth(0.0)) * Math.max(eb.getyWidth(0.0), lp.getWidth(0.0))
                * Math.max(eb.getLength(), lp.getLength()) * rayXAnglerange * rayYAnglerange
                * (maxEnergy - minEnergy);
        prob0 = directionFrequencyVolumePolarizationNoSpread(r, n, new BasicVector(new double[]{0.0, 0.0, 1.0}), EMax)[0] / EMax;
        if (iseSpread()) {
            prob0 *= eb.angleDistribution(0, 0);
            factor *= 4 * MULT * MULT * eb.getXSpread() * eb.getYSpread();
        }
        do {
            if (Thread.currentThread().isInterrupted()) {
                throw new InterruptedException();
            }
            ray[0] = MULT * (2 * Math.random() - 1.0) * Math.max(eb.getxWidth(0.0), lp.getWidth(0.0));
            ray[2] = MULT * (2 * Math.random() - 1.0) * Math.max(eb.getyWidth(0.0), lp.getWidth(0.0));
            ray[1] = MULT * (2 * Math.random() - 1.0) * Math.max(eb.getLength(), lp.getLength());
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
                double thetax = MULT * eb.getXSpread() * (2 * Math.random() - 1);
                double thetay = MULT * eb.getYSpread() * (2 * Math.random() - 1);
                Vector v = new BasicVector(new double[]{thetax, thetay, Math.sqrt(1 - thetax * thetax - thetay * thetay)});
                if (ksi == null) {
                    polParam = directionFrequencyVolumePolarizationNoSpread(r, n, v, ray[10]);
                } else {
                    polParam = new double[]{directionFrequencyVolumeFluxNoSpread(r, n, v, ray[10]), ksi[0], ksi[1], ksi[2]};
                }
                prob = polParam[0] * eb.angleDistribution(thetax, thetay) / ray[10];
            } else {
                if (ksi == null) {
                    polParam = directionFrequencyVolumePolarizationNoSpread(r, n, new BasicVector(new double[]{0.0, 0.0, 1.0}), ray[10]);
                } else {
                    polParam = new double[]{directionFrequencyVolumeFluxNoSpread(r, n, new BasicVector(new double[]{0.0, 0.0, 1.0}), ray[10]), ksi[0], ksi[1], ksi[2]};
                }
                prob = polParam[0] / ray[10];
            }
            if (!new Double(prob).isNaN()) {
                sum += prob;
            }
            lcn++;
        } while (prob / prob0 < Math.random() || (new Double(prob)).isNaN());

        // Calculating the rotated polarization vector and getting the full polarizaation state
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
        //Incrementing sum and inegration point counter
        partialFlux.add(sum * factor);
        montecarlocounter.addAndGet(lcn);
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
     * Counter of integral ray iterations
     *
     * @return the montecarlocounter
     */
    public long getMonteCarloCounter() {
        return montecarlocounter.get() == 0 ? 1 : montecarlocounter.get();
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
        //If p > 1 reducing it to 1
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
     * Getting the number of points in the emittance Monte-Carlo integration
     *
     * @return the npEmittance
     */
    public int getNpEmittance() {
        return npEmittance;
    }

    /**
     * Setting the number of points in the emittance Monte-Carlo integration
     *
     * @param npEmittance the npEmittance to set
     */
    public void setNpEmittance(int npEmittance) {
        this.npEmittance = npEmittance;
    }

    /**
     * Getting the IsMonteCarlo flag
     *
     * @return the IsMonteCarlo
     */
    public boolean isIsMonteCarlo() {
        return IsMonteCarlo;
    }

    /**
     * Setting the IsMonteCarlo flag
     *
     * @param IsMonteCarlo the IsMonteCarlo to set
     */
    public void setIsMonteCarlo(boolean IsMonteCarlo) {
        this.IsMonteCarlo = IsMonteCarlo;
    }

    /**
     * @return the IsCompton
     */
    public boolean isIsCompton() {
        return IsCompton;
    }

    /**
     * @param IsCompton the IsCompton to set
     */
    public void setIsCompton(boolean IsCompton) {
        this.IsCompton = IsCompton;
    }
}
