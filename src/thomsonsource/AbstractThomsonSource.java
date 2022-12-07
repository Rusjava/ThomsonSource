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
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Properties;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.DoubleAdder;
import java.util.function.Consumer;
import laserpulse.AbstractLaserPulse;
import laserpulse.GaussianLaserPulse;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.BaseAbstractUnivariateIntegrator;
import org.apache.commons.math3.analysis.integration.RombergIntegrator;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.la4j.Matrix;
import org.la4j.Vector;
import org.la4j.Vectors;
import org.la4j.vector.dense.BasicVector;
import shadowfileconverter.ShadowFiles;

/**
 * An abstract class for Thomson source. Methods that calculated scattering by
 * one electron need to be defined.
 *
 * @version 1.32
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
        this.montecarlocounter = new AtomicInteger();
        this.rayCounter = new AtomicInteger();
        this.partialFlux = new DoubleAdder();
        this.ksi = new double[]{0, 0, -1};
        this.paramNames = new String[]{"Electron_energy_MeV", "Electron_bunch_charge_nQ",
            "Electron_bunch_relative_energy_spread", "Electron_bunch_length_ps",
            "X-emittance_mm*mrad", "Y-emittance_mm*mrad", "Beta-x_function_mm", "Beta-y_function_mm", "Photon_energy_eV",
            "Pulse_energy_mJ", "Laser_pulse_length_ps", "Rayleigh_length_mm",
            "Pulse_frequency_MHz", "Delay_ps", "X-shift_mm",
            "Y-shift_mm", "Z-shift_mm", "Laser-electron_direction_x",
            "Laser-electron_direction_y", "Laser-electron_direction_z", "Laser-electron_angle_mrad"};
    }
    /**
     * Names of Thomson source parameters
     */
    private final String[] paramNames;
    /**
     * Multiple for the range
     */
    protected static final double INT_RANGE = 3;
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
     * A numerical shift to improve convergence of integrals
     */
    public static final double SHIFT = 1e10;
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
    protected double minEnergy = 25 * 1e3 * GaussianElectronBunch.E;
    /**
     * Max ray energy
     */
    protected double maxEnergy = 35 * 1e3 * GaussianElectronBunch.E;
    /**
     * Number of points in Monte Carlo calculation of the geometric factor
     */
    protected int npGeometricFactor = 50000;
    /**
     * A shift factor to improve numerical integral convergence in polarization
     * calculations
     */
    private double shiftfactor = 1;
    /**
     * Precision in calculations of the brilliance
     */
    protected double precision = 0.001;
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
     * Counter of Monte-Carlo iterations
     */
    protected AtomicInteger montecarlocounter;
    /**
     * Counter of rays
     */
    private AtomicInteger rayCounter;

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
        ((AbstractThomsonSource) tm).montecarlocounter = new AtomicInteger();
        ((AbstractThomsonSource) tm).rayCounter = new AtomicInteger();
        ((AbstractThomsonSource) tm).partialFlux = new DoubleAdder();
        if (ksi != null) {
            ((AbstractThomsonSource) tm).ksi = (double[]) ksi.clone();
        }
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
     * A method calculating the linear total flux
     */
    public final void calculateLinearTotalFlux() {
        this.totalFlux = SIGMA_T * eb.getNumber() * lp.getPhotonNumber() * lp.getFq()
                / Math.PI / Math.sqrt((lp.getWidth2(0.0) + eb.getxWidth2(0.0)) * (lp.getWidth2(0.0) + eb.getyWidth2(0.0)));
    }

    /**
     * A method calculating normalized total flux within a certain angle
     *
     * @param maxAngle
     * @return
     */
    public final double calculateAngleLinearTotalFlux(double maxAngle) {
        double gamma2 = eb.getGamma() * eb.getGamma();
        double v = Math.sqrt(1 - 1 / gamma2);
        double cs = Math.cos(maxAngle);
        return 3.0 / 4 / gamma2 * ((1 - cs) / (1 - v * cs) / (1 - v)
                * (5.0 / 6 + 1.0 / 6 / v / v - 1.0 / 6 / gamma2 / v / v * (1 - v * v * cs) / (1 - v) / (1 - v * cs))
                + 1.0 / 6 / gamma2 / v * (1 - cs * cs) / Math.pow((1 - v * cs), 3))
                * getLinearTotalFlux() * getGeometricFactor();
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
     * @param n viewing direction
     * @param v normalized electron velocity
     * @param r spatial position
     * @param e X-ray energy
     * @return
     * @throws java.lang.InterruptedException
     */
    public double directionFrequencyFlux(Vector n, Vector v, Vector r, double e) throws InterruptedException {
        return iseSpread() ? directionFrequencyFluxSpread(n, v, r, e) : directionFrequencyFluxNoSpread(n, v, r, e);
    }

    /**
     * A method calculating the Stocks parameters density in a given direction
     * for a given X-ray photon energy
     *
     * @param n viewing direction
     * @param v normalized electron velocity
     * @param r spatial position
     * @param e X-ray energy
     * @return
     * @throws java.lang.InterruptedException
     */
    public double[] directionFrequencyPolarization(Vector n, Vector v, Vector r, double e) throws InterruptedException {
        return iseSpread() ? directionFrequencyPolarizationSpread(n, v, r, e) : directionFrequencyPolarizationNoSpread(n, v, r, e);
    }

    /**
     * A method calculating a Stocks parameter density in a given direction for
     * a given X-ray photon energy
     *
     * @param n viewing direction
     * @param v normalized electron velocity
     * @param r spatial position
     * @param e X-ray energy
     * @param index
     * @return
     * @throws java.lang.InterruptedException
     */
    public double directionFrequencyPolarization(Vector n, Vector v, Vector r, double e, int index) throws InterruptedException {
        return iseSpread() ? directionFrequencyPolarizationSpread(n, v, r, e, index) : directionFrequencyPolarizationNoSpread(n, v, r, e, index);
    }

    /**
     * A method calculating the flux density in a given direction for a given
     * X-ray photon energy without taking into account electron transversal
     * pulse spread
     *
     * @param n viewing direction
     * @param v normalized electron velocity
     * @param r spatial position in the laser coordinates
     * @param e X-ray energy
     * @return
     */
    public abstract double directionFrequencyFluxNoSpread(Vector n, Vector v, Vector r, double e);

    /**
     * A method calculating the Stocks parameters density in a given direction
     * for a given X-ray photon energy without taking into account electron
     * transversal pulse spread
     *
     * @param n viewing direction
     * @param v normalized electron velocity
     * @param r spatial position in the laser coordinates
     * @param e X-ray energy
     * @return
     */
    public abstract double[] directionFrequencyPolarizationNoSpread(Vector n, Vector v, Vector r, double e);

    /**
     * A method calculating a Stocks parameter density in a given direction for
     * a given X-ray photon energy without taking into account electron
     * transversal pulse spread
     *
     * @param n viewing direction
     * @param v normalized electron velocity
     * @param r spatial position in the laser coordinates
     * @param e X-ray energy
     * @param index polarization matrix element
     * @return
     */
    public abstract double directionFrequencyPolarizationNoSpread(Vector n, Vector v, Vector r, double e, int index);

    /**
     * A method calculating the flux density in a given direction for a given
     * X-ray photon energy taking into account the electron transversal pulse
     * spread
     *
     * @param n viewing direction
     * @param v0 normalized electron velocity
     * @param r spatial position
     * @param e X-ray energy
     * @return
     * @throws java.lang.InterruptedException
     */
    public double directionFrequencyFluxSpread(Vector n, Vector v0, Vector r, double e) throws InterruptedException {
        BaseAbstractUnivariateIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateFunction func = new UnivariateFrequencyFluxSpreadOuter(e, v0, r, n);
        double tmp;
        try {
            //If interrupted, throw InterruptedException
            if (Thread.currentThread().isInterrupted()) {
                throw new InterruptedException("Interruption in directionFrequencyFluxSpread!");
            }
            tmp = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, 0.0, 2 * Math.PI)
                    - getShiftfactor() * SHIFT * 2 * Math.PI * (1 - Math.cos(INT_RANGE * eb.getSpread()));
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
        return new Double(tmp).isNaN() ? 0 : tmp;
    }

    /**
     * A method calculating the full polarization tensor density in a given
     * direction for a given X-ray photon energy taking into account the
     * electron transversal pulse spread
     *
     * @param n viewing direction
     * @param v0 normalized electron velocity
     * @param r spatial position
     * @param e X-ray energy
     * @return
     * @throws java.lang.InterruptedException
     */
    public double[] directionFrequencyPolarizationSpread(final Vector n, final Vector v0, final Vector r, final double e) throws InterruptedException {
        //An array for results
        double[] array = new double[NUMBER_OF_POL_PARAM];
        //Calculating the polarization tensor elements
        for (int i = 0; i < NUMBER_OF_POL_PARAM; i++) {
            try {
                //Calculating elements of polarization matrix
                array[i] = directionFrequencyPolarizationSpread(n, v0, r, e, i);
            } catch (TooManyEvaluationsException ex) {
                array[i] = 0;
            }
        }

        //If the intensity is NaN, zero or less than zero then all Stocks intensities to zero
        // If a Stocks intensity is NaN then set it to zero
        for (int i = 1; i < NUMBER_OF_POL_PARAM; i++) {
            if (new Double(array[i]).isNaN() || new Double(array[0]).isNaN() || array[0] <= 0) {
                array[i] = 0;
            }
        }
        return array;
    }

    /**
     * A multi-threaded method calculating a Stocks parameter density in a given
     * direction for a given X-ray photon energy taking into account the
     * electron transversal pulse spread
     *
     * @param n viewing direction
     * @param v0 normalized electron velocity
     * @param r spatial position
     * @param e X-ray energy
     * @param index polarization matrix element
     * @return
     * @throws java.lang.InterruptedException
     */
    public double directionFrequencyPolarizationSpread(final Vector n, final Vector v0, final Vector r,
            final double e, int index) throws InterruptedException {
        //The function to integrate
        UnivariateFunction func = new UnivariateFrequencyPolarizationSpreadOuter(e, v0, n, r, index);
        double res;
        try {
            //If interrupted, throw InterruptedException
            if (Thread.currentThread().isInterrupted()) {
                throw new InterruptedException("Interruption in directionFrequencyPolarizationSpread!");
            }
            //Creating an inegrator and integrating
            res = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                    RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT)
                    .integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, 0.0, 2 * Math.PI)
                    - getShiftfactor() * SHIFT * 2 * Math.PI * (1 - Math.cos(INT_RANGE * eb.getSpread()));
        } catch (TooManyEvaluationsException ex) {
            res = 0;
        }
        return new Double(res).isNaN() ? 0 : res;
    }

    /**
     * A method calculating the flux density in a given direction for a given
     * X-ray photon energy for a given volume element
     *
     * @param r spatial position
     * @param n viewing direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     * @throws java.lang.InterruptedException
     */
    public double directionFrequencyVolumeFlux(Vector r, Vector n, Vector v, double e) throws InterruptedException {
        return iseSpread() ? directionFrequencyVolumeFluxSpread(r, n, v, e) : directionFrequencyVolumeFluxNoSpread(r, n, v, e);
    }

    /**
     * A method calculating the Stocks parameters density in a given direction
     * for a given X-ray photon energy for a given volume element
     *
     * @param r spatial position
     * @param n viewing direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     * @throws java.lang.InterruptedException
     */
    public double[] directionFrequencyVolumePolarization(Vector r, Vector n, Vector v, double e) throws InterruptedException {
        return iseSpread() ? directionFrequencyVolumePolarizationSpread(r, n, v, e) : directionFrequencyVolumePolarizationNoSpread(r, n, v, e);
    }

    /**
     * A method calculating a Stocks parameter density in a given direction for
     * a given X-ray photon energy for a given volume element
     *
     * @param r spatial position
     * @param n viewing direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @param index polarization matrix element
     * @return
     * @throws java.lang.InterruptedException
     */
    public double directionFrequencyVolumePolarization(Vector r, Vector n, Vector v, double e, int index) throws InterruptedException {
        return iseSpread() ? directionFrequencyVolumePolarizationSpread(r, n, v, e, index) : directionFrequencyVolumePolarizationNoSpread(r, n, v, e, index);
    }

    /**
     * A method calculating the flux density in a given direction for a given
     * X-ray photon energy for a given volume element
     *
     * @param r spatial position
     * @param n viewing direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     * @throws java.lang.InterruptedException
     */
    abstract public double directionFrequencyVolumeFluxNoSpread(Vector r, Vector n, Vector v, double e) throws InterruptedException;

    /**
     * A method calculating the Stocks parameters density in a given direction
     * for a given X-ray photon energy for a given volume element without taking
     * into account the electron transversal pulse spread
     *
     * @param r spatial position
     * @param n viewing direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     * @throws java.lang.InterruptedException
     */
    abstract public double[] directionFrequencyVolumePolarizationNoSpread(Vector r, Vector n, Vector v, double e) throws InterruptedException;

    /**
     * A method calculating a Stocks parameter density in a given direction for
     * a given X-ray photon energy for a given volume element without taking
     * into account the electron transversal pulse spread
     *
     * @param r spatial position
     * @param n viewing direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @param index polarization matrix element
     * @return
     * @throws java.lang.InterruptedException
     */
    abstract public double directionFrequencyVolumePolarizationNoSpread(Vector r, Vector n, Vector v, double e, int index) throws InterruptedException;

    /**
     * A method calculating the flux density in a given direction for a given
     * X-ray photon energy for a given volume element taking into account the
     * electron transversal pulse spread
     *
     * @param r spatial position
     * @param n viewing direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     * @throws java.lang.InterruptedException
     */
    abstract public double directionFrequencyVolumeFluxSpread(Vector r, Vector n, Vector v, double e) throws InterruptedException;

    /**
     * A method calculating the Stocks parameters density in a given direction
     * for a given X-ray photon energy for a given volume element taking into
     * account the electron transversal pulse spread
     *
     * @param r spatial position
     * @param n viewing direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     * @throws java.lang.InterruptedException
     */
    abstract public double[] directionFrequencyVolumePolarizationSpread(Vector r, Vector n, Vector v, double e) throws InterruptedException;

    /**
     * A method calculating a Stocks parameter density in a given direction for
     * a given X-ray photon energy for a given volume element taking into
     * account the electron transversal pulse spread
     *
     * @param r spatial position
     * @param n viewing direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @param index polarization matrix element
     * @return
     * @throws java.lang.InterruptedException
     */
    abstract public double directionFrequencyVolumePolarizationSpread(Vector r, Vector n, Vector v, double e, int index) throws InterruptedException;

    /**
     * An auxiliary method calculating volume density of the X-ray source
     *
     * @param r spatial position
     * @return
     */
    public double volumeFlux(Vector r) {
        //If r is null then use zero coordinates
        if (r == null) {
            r = new BasicVector(new double[]{0, 0, 0});
        }
        double u;
        double len = Math.sqrt(lp.getLength() * lp.getLength() + eb.getLength() * eb.getLength());
        Vector r1 = lp.getTransformedCoordinates(r);
        double K = Math.pow((r.get(2) - r1.get(2) - eb.getShift().get(2) + lp.getDelay()) / len, 2);
        u = 2.0 * Math.sqrt(Math.PI) * Math.sqrt((lp.getWidth2(0.0) + eb.getxWidth2(0.0))
                * (lp.getWidth2(0.0) + eb.getyWidth2(0.0))) / len * Math.exp(-K) * eb.tSpatialDistribution(r)
                * lp.tSpatialDistribution(r1);
        return new Double(u).isNaN() ? 0 : u;
    }

    /**
     * A method giving the flux density in a given direction
     *
     * @param n viewing direction
     * @param v normalized electron velocity
     * @return
     */
    public abstract double directionFlux(Vector n, Vector v);

    /**
     * A method calculating X-ray energy in a given direction
     *
     * @param n viewing direction
     * @param v normalized (to unity) electron velocity
     * @return
     */
    public abstract double directionEnergy(Vector n, Vector v);

    /**
     * A method calculating spectral brilliance in a given direction
     *
     * @param r0 spatial position for brightness
     * @param n viewing direction
     * @param v normalized (to unity) electron velocity
     * @param e X-ray energy
     * @return
     * @throws java.lang.InterruptedException
     */
    public double directionFrequencyBrilliance(Vector r0, Vector n, Vector v, double e) throws InterruptedException {
        return iseSpread() ? directionFrequencyBrillianceSpread(r0, n, v, e) : directionFrequencyBrillianceNoSpread(r0, n, v, e);
    }

    /**
     * A method calculating spectral brilliance in a given direction with
     * polarization
     *
     * @param r0 spatial position for brightness
     * @param n viewing direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     * @throws java.lang.InterruptedException
     */
    public double[] directionFrequencyPolarizationBrilliance(Vector r0, Vector n, Vector v, double e) throws InterruptedException {
        return iseSpread() ? directionFrequencyBrilliancePolarizationSpread(r0, n, v, e) : directionFrequencyBrilliancePolarizationNoSpread(r0, n, v, e);
    }

    /**
     * A method calculating spectral brilliance in a given direction with
     * polarization
     *
     * @param r0 spatial position for brightness
     * @param n viewing direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @param index polarization matrix element
     * @return
     * @throws java.lang.InterruptedException
     */
    public double directionFrequencyPolarizationBrilliance(Vector r0, Vector n, Vector v, double e, int index) throws InterruptedException {
        return iseSpread() ? directionFrequencyBrilliancePolarizationSpread(r0, n, v, e, index) : directionFrequencyBrilliancePolarizationNoSpread(r0, n, v, e, index);
    }

    /**
     * A method calculating spectral brilliance in a given direction without
     * taking into account electron transversal pulse spread
     *
     * @param r0 spatial position for brightness
     * @param n viewing direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     * @throws java.lang.InterruptedException
     */
    abstract public double directionFrequencyBrillianceNoSpread(Vector r0, Vector n, Vector v, double e) throws InterruptedException;

    /**
     * A method calculating spectral brilliance in a given direction without
     * taking into account electron transversal pulse spread but with
     * polarization
     *
     * @param r0 spatial position for brightness
     * @param n viewing direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     * @throws java.lang.InterruptedException
     */
    abstract public double[] directionFrequencyBrilliancePolarizationNoSpread(Vector r0, Vector n, Vector v, double e) throws InterruptedException;

    /**
     * A method calculating spectral brilliance in a given direction without
     * taking into account electron transversal pulse spread but with
     * polarization
     *
     * @param r0 spatial position for brightness
     * @param n viewing direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @param index polarization matrix element
     * @return
     * @throws java.lang.InterruptedException
     */
    abstract public double directionFrequencyBrilliancePolarizationNoSpread(Vector r0, Vector n, Vector v, double e, int index) throws InterruptedException;

    /**
     * A method calculating spectral brilliance in a given direction taking into
     * account electron transversal pulse spread
     *
     * @param r0 spatial position for brightness
     * @param n viewing direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     * @throws java.lang.InterruptedException
     */
    abstract public double directionFrequencyBrillianceSpread(Vector r0, Vector n, Vector v, double e) throws InterruptedException;

    /**
     * A method calculating spectral brilliance in a given direction taking into
     * account electron transversal pulse spread nd with polarization
     *
     * @param r0 spatial position for brightness
     * @param n viewing direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @return
     * @throws java.lang.InterruptedException
     */
    abstract public double[] directionFrequencyBrilliancePolarizationSpread(Vector r0, Vector n, Vector v, double e) throws InterruptedException;

    /**
     * A method calculating spectral brilliance in a given direction taking into
     * account electron transversal pulse spread nd with polarization
     *
     * @param r0 spatial position for brilliance calculations
     * @param n viewing direction
     * @param v normalized electron velocity
     * @param e X-ray energy
     * @param index polarization matrix element
     * @return
     * @throws java.lang.InterruptedException
     */
    abstract public double directionFrequencyBrilliancePolarizationSpread(Vector r0, Vector n, Vector v, double e, int index) throws InterruptedException;

    /**
     * An auxiliary function calculating the brilliance integral along the given
     * direction
     *
     * @param r0 spatial position for brightness
     * @param n viewing direction
     * @param func the function to integrate over
     * @param index the order number
     * @return
     * @throws InterruptedException
     */
    protected double directionIntegralBasic(Vector r0, Vector n, UnivariateFunction func, int index) throws InterruptedException {
        //return directionFrequencyVolumeFluxNoSpread(r0, n, v, e);
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);

        //Semiwidth of the integration interval
        double semiwidth = INT_RANGE * lp.getLength() * lp.getWidth(0)
                / Math.sqrt(Math.pow(lp.getLength() * n.get(0), 2) + lp.getWidth2(0) * Math.pow(n.get(2), 2)) / 2 / Math.sqrt(index);
        //Integrating over a line to calculate spectral brilliance
        try {
            //If interrupted, throw InterruptedException
            if (Thread.currentThread().isInterrupted()) {
                throw new InterruptedException("directionIntegralBasic!");
            }
            double u = integrator.integrate(AbstractThomsonSource.MAXIMAL_NUMBER_OF_EVALUATIONS, func,
                    r0.fold(Vectors.mkEuclideanNormAccumulator()) - semiwidth, r0.fold(Vectors.mkEuclideanNormAccumulator()) + semiwidth)
                    - 2 * semiwidth * SHIFT * getShiftfactor();
            return u;
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
    }

    /**
     * Calculates the integral over time for the laser-electron interaction
     *
     * @param r spatial position
     * @param n viewing direction
     * @param func the function to integrate
     * @return
     * @throws java.lang.InterruptedException
     */
    protected double timeIntegralBasic(Vector r, Vector n, UnivariateFunction func) throws InterruptedException {
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        double tmp;
        //Transforming coordinates between laser and electron beam frames
        Vector rph = lp.getTransformedCoordinates(r);

        //Defining the upper nad lower integration limits
        double semilength = INT_RANGE * eb.getLength() * lp.getLength() / Math.sqrt(eb.getLength() * eb.getLength() + lp.getLength() * lp.getLength());
        double shft = ((rph.get(2) - lp.getDelay()) * eb.getLength() * eb.getLength() + (r.get(2) - eb.getShift().get(2)) * lp.getLength() * lp.getLength())
                / (eb.getLength() * eb.getLength() + lp.getLength() * lp.getLength());

        double zmin = -semilength + shft;
        double zmax = semilength + shft;
        //Integrating by time
        try {
            //If interrupted, throw InterruptedException
            if (Thread.currentThread().isInterrupted()) {
                throw new InterruptedException("timeIntegralBasic!");
            }
            tmp = 2.0 * Math.PI * Math.sqrt((lp.getWidth2(0.0) + eb.getxWidth2(0.0))
                    * (lp.getWidth2(0.0) + eb.getyWidth2(0.0))) * integrator.integrate(AbstractThomsonSource.MAXIMAL_NUMBER_OF_EVALUATIONS, func, zmin, zmax)
                    * lp.tSpatialDistribution(rph) * eb.tSpatialDistribution(r);

            //Testing if NaN, then return zero
            return new Double(tmp).isNaN() ? 0 : tmp;
        } catch (TooManyEvaluationsException ex) {
            return 0;
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
        Vector n0 = new BasicVector(new double[]{0.0, 1.0, 0.0});
        Vector As;
        double prob0;
        double prob;
        double EMax;
        double MULT = 2;
        double factor;
        double sum = 0;
        double[] pol;
        double[] polParam;
        EMax = directionEnergy(n, n);
        factor = 32 * MULT * MULT * MULT * Math.max(eb.getxWidth(0.0), lp.getWidth(0.0)) * Math.max(eb.getyWidth(0.0), lp.getWidth(0.0)) * Math.max(eb.getLength(), lp.getLength())
                * rayXAnglerange * rayYAnglerange * (maxEnergy - minEnergy);
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
            //Incrementing the ray montecarlocounter
            montecarlocounter.incrementAndGet();
        } while (prob / prob0 < Math.random() || (new Double(prob)).isNaN());
        // Calculating the rotated polarization vector and getting the full polarizaation state
        //n is defined here with y in longitudinal direction
        n = new BasicVector(new double[]{ray[3], ray[4], ray[5]});
        T = get3DTransform(n, n0);
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
        //Adding to the Monte-Carlo sum and counter
        partialFlux.add(sum * factor);
        rayCounter.incrementAndGet();
        return ray;
    }

    /**
     * Writing a specific number of rays into a file
     *
     * @param shadowFile The handle of the Shadow binary ray file
     * @param numberOfRays The number of rays to generate
     * @param con Consumer to update the progress bar
     * @return
     * @throws java.lang.InterruptedException
     */
    public boolean writeRays(ShadowFiles shadowFile, int numberOfRays, Consumer<Integer> con) throws InterruptedException {
        // Creating a pool of threads, a lock, an atomic interger and a latch
        ExecutorService excs = Executors.newFixedThreadPool(getThreadNumber());
        CountDownLatch lt = new CountDownLatch(getThreadNumber());
        int rayNumber = getThreadNumber() * (numberOfRays / getThreadNumber());
        //Creating the number of tasks equal to the number of threads
        for (int th = 0; th < getThreadNumber(); th++) {
            if (Thread.currentThread().isInterrupted()) {
                excs.shutdownNow();
                throw new InterruptedException("writeRays method interrupted!");
            }
            //Creating multiple threads to accelerate calculations
            excs.execute(() -> {
                for (int i = 0; i < rayNumber / getThreadNumber(); i++) {
                    try {
                        //Getting a ray
                        double[] ray = getRay();
                        //Units conversions
                        ray[0] *= 1e2;
                        ray[1] *= 1e2;
                        ray[2] *= 1e2;
                        ray[10] *= 1e-2 / GaussianLaserPulse.HC;
                        ray[11] = i;
                        shadowFile.write(ray);
                        //Updating the progress bar
                        con.accept((int) 100 * (rayCounter.get() + 1) / rayNumber);
                    } catch (IOException | InterruptedException ex) {
                        break;
                    }
                }
                lt.countDown();
            });
        }
        //If interrupted shutdown threads and rethrow the exception
        try {
            lt.await();
        } catch (InterruptedException ex) {
            excs.shutdownNow();
            throw ex;
        }
        excs.shutdownNow();
        return true;
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
     * @param n vector 1
     * @param n0 vector 2
     * @return transformation matrix
     */
    protected Matrix get3DTransform(Vector n, Vector n0) {
        Matrix D;
        Matrix A;
        Matrix I = Matrix.identity(3);
        double innerProduct = n.innerProduct(n0);
        D = n.outerProduct(n0).add(n0.outerProduct(n)).multiply(innerProduct).subtract(n.outerProduct(n).add(n0.outerProduct(n0))).divide(innerProduct + 1.0);
        A = n.outerProduct(n0).subtract(n0.outerProduct(n)).add(I.multiply(innerProduct));
        return I.multiply(1 - innerProduct).add(D).add(A);
    }

    /**
     * Returning the matrix of 2D rotation based on two unity vectors
     *
     * @param n vector 1
     * @param n0 vector 2
     * @return transformation matrix
     */
    protected Matrix get2DTransform(Vector n, Vector n0) {
        Matrix I = Matrix.identity(2);
        double cs = n.innerProduct(n0);
        double sn = Math.sqrt(1 - cs * cs);
        I = I.multiply(cs);
        I.set(0, 1, sn);
        I.set(1, 0, -sn);
        return I;
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
    public double getLinearTotalFlux() {
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
        return partialFlux.sum() / (montecarlocounter.get() == 0 ? 1 : montecarlocounter.get());
    }

    /**
     * Counter of ray iterations
     *
     * @return the montecarlocounter
     */
    public int getCounter() {
        return montecarlocounter.get();
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
     * Save Thomson source properties into a file
     *
     * @param file file to save parameters to
     * @throws java.io.IOException
     */
    public void savePropperties(File file) throws IOException {
        //Creating Properties object to store program parameters
        Properties prop = new Properties();
        try (FileWriter fw = new FileWriter(file, false)) {
            prop.setProperty(paramNames[0], Double.toString(eb.getGamma() * 0.512));
            prop.setProperty(paramNames[1], Double.toString(eb.getNumber() * GaussianElectronBunch.E * 1e9));
            prop.setProperty(paramNames[2], Double.toString(eb.getDelGamma()));
            prop.setProperty(paramNames[3], Double.toString(eb.getLength() * 2 / 3e-4));
            prop.setProperty(paramNames[4], Double.toString(eb.getEpsx() * 1e6));
            prop.setProperty(paramNames[5], Double.toString(eb.getEpsy() * 1e6));
            prop.setProperty(paramNames[6], Double.toString(eb.getBetax() * 1e3));
            prop.setProperty(paramNames[7], Double.toString(eb.getBetay() * 1e3));
            prop.setProperty(paramNames[8], Double.toString(lp.getPhotonEnergy() / GaussianElectronBunch.E));
            prop.setProperty(paramNames[9], Double.toString(lp.getPulseEnergy() * 1e3));
            prop.setProperty(paramNames[10], Double.toString(lp.getLength() * 2 / 3e-4));
            prop.setProperty(paramNames[11], Double.toString(lp.getRlength() * 1e3));
            prop.setProperty(paramNames[12], Double.toString(lp.getFq() * 1e-6));
            prop.setProperty(paramNames[13], Double.toString(lp.getDelay() / 3e-4));
            prop.setProperty(paramNames[14], Double.toString(eb.getShift().get(0) * 1e3));
            prop.setProperty(paramNames[15], Double.toString(eb.getShift().get(1) * 1e3));
            prop.setProperty(paramNames[16], Double.toString(eb.getShift().get(2) * 1e3));
            prop.setProperty(paramNames[17], Double.toString(lp.getDirection().get(0)));
            prop.setProperty(paramNames[18], Double.toString(lp.getDirection().get(1)));
            prop.setProperty(paramNames[19], Double.toString(lp.getDirection().get(2)));
            prop.store(fw, "Thomson source parameters");
        } catch (IOException e) {
            throw e;
        }
    }

    /**
     * Reads the file and loads the parameters into object field
     *
     * @param file file date to be read from
     * @return
     * @throws IOException
     * @throws NumberFormatException
     */
    public Properties loadProperties(File file) throws IOException, NumberFormatException {
        Properties prop = new Properties();
        try (FileReader fr = new FileReader(file)) {
            prop.load(fr);
        } catch (IOException e) {
            throw e;
        }
        try {
            eb.setGamma(Float.parseFloat(prop.getProperty(paramNames[0], "0")) / 0.512);
            eb.setNumber(Float.parseFloat(prop.getProperty(paramNames[1], "0")) / GaussianElectronBunch.E * 1e-9);
            eb.setDelgamma(Float.parseFloat(prop.getProperty(paramNames[2], "0")));
            eb.setLength(Float.parseFloat(prop.getProperty(paramNames[3], "0")) / 2 * 3e-4);
            eb.setEpsx(Float.parseFloat(prop.getProperty(paramNames[4], "0")) / 1e6);
            eb.setEpsy(Float.parseFloat(prop.getProperty(paramNames[5], "0")) / 1e6);
            eb.setBetax(Float.parseFloat(prop.getProperty(paramNames[6], "0")) * 1e-3);
            eb.setBetay(Float.parseFloat(prop.getProperty(paramNames[7], "0")) * 1e-3);
            lp.setPhotonEnergy(Float.parseFloat(prop.getProperty(paramNames[8], "0")) * GaussianElectronBunch.E);
            lp.setPulseEnergy(Float.parseFloat(prop.getProperty(paramNames[9], "0")) * 1e-3);
            lp.setLength(Float.parseFloat(prop.getProperty(paramNames[10], "0")) / 2 * 3e-4);
            lp.setRlength(Float.parseFloat(prop.getProperty(paramNames[11], "0")) * 1e-3);
            lp.setFq(Float.parseFloat(prop.getProperty(paramNames[12], "0")) * 1e6);
            lp.setDelay(Float.parseFloat(prop.getProperty(paramNames[13], "0")) * 3e-4);
            Vector sht = new BasicVector(new double[3]);
            sht.set(0, Float.parseFloat(prop.getProperty(paramNames[14], "0")) * 1e-3);
            sht.set(1, Float.parseFloat(prop.getProperty(paramNames[15], "0")) * 1e-3);
            sht.set(2, Float.parseFloat(prop.getProperty(paramNames[16], "0")) * 1e-3);
            eb.setShift(sht);
            Vector dir = new BasicVector(new double[3]);
            dir.set(0, Float.parseFloat(prop.getProperty(paramNames[17], "0")) * 1e-3);
            dir.set(1, Float.parseFloat(prop.getProperty(paramNames[18], "0")) * 1e-3);
            dir.set(2, Float.parseFloat(prop.getProperty(paramNames[19], "0")) * 1e-3);
            dir.set(1, Math.sin(Float.parseFloat(prop.getProperty(paramNames[20], "0")) * 1e-3));
            dir.set(2, Math.cos(Float.parseFloat(prop.getProperty(paramNames[20], "0")) * 1e-3));
            lp.setDirection(dir);
        } catch (NumberFormatException e) {
            throw e;
        }
        return prop;
    }

    /**
     * @return the paramNames
     */
    public String[] getParamNames() {
        return paramNames;
    }

    /**
     * @return the rayCounter
     */
    public AtomicInteger getRayCounter() {
        return rayCounter;
    }

    /**
     *
     * @param polParam polarization parameter vector to multiply
     * @param factor factor to multiply by
     * @return
     */
    protected double[] multiplyPolarizationParameters(double[] polParam, double factor) {
        for (int s = 0; s < AbstractThomsonSource.NUMBER_OF_POL_PARAM; s++) {
            polParam[s] *= factor;
        }
        return polParam;
    }

    /**
     * Calculates vector product of two vectors in 3D space
     *
     * @param a
     * @param b
     * @return
     */
    protected Vector crossProduct3D(Vector a, Vector b) {
        return new BasicVector(new double[]{a.get(1) * b.get(2) - a.get(2) * b.get(1),
            a.get(2) * b.get(0) - a.get(0) * b.get(2), a.get(0) * b.get(1) - a.get(1) * b.get(0)});
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
        private final Vector n, r, v0;
        private final BaseAbstractUnivariateIntegrator inergrator;

        public UnivariateFrequencyFluxSpreadOuter(double e, Vector v0, Vector r, Vector n) {
            this.e = e;
            this.v0 = v0;
            this.n = n;
            this.r = r;
            this.inergrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                    RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        }

        @Override
        public double value(double phi) {
            double tmp;
            if (Thread.currentThread().isInterrupted()) {
                Thread.currentThread().interrupt();
                return 0;
            }
            try {
                tmp = inergrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, new UnivariateFrequencyFluxSpreadInner(phi, e, v0, r, n),
                        0.0, INT_RANGE * eb.getSpread());
            } catch (TooManyEvaluationsException ex) {
                return 0;
            }
            return new Double(tmp).isNaN() ? 0 : tmp;
        }
    }

    /**
     * An auxiliary class for Romberg integrator for flux calculations
     */
    private class UnivariateFrequencyFluxSpreadInner implements UnivariateFunction {

        private final double snphi, csphi, e;
        private final Vector n, r, v0;

        public UnivariateFrequencyFluxSpreadInner(double phi, double e, Vector v0, Vector r, Vector n) {
            this.snphi = Math.sin(phi);
            this.csphi = Math.cos(phi);
            this.e = e;
            this.n = n;
            this.r = r;
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
            u = sn * directionFrequencyFluxNoSpread(n, v, r, e) * eb.angleDistribution(dv.get(0), dv.get(1))
                    + getShiftfactor() * SHIFT;
            return new Double(u).isNaN() ? 0 : u;
        }
    }

    /**
     * An auxiliary class for Romberg integrator for polarization calculations
     */
    private class UnivariateFrequencyPolarizationSpreadOuter implements UnivariateFunction {

        private final double e;
        private final Vector n, v0, r;
        private final int index;
        private final BaseAbstractUnivariateIntegrator inergrator;

        public UnivariateFrequencyPolarizationSpreadOuter(double e, Vector v0, Vector n, Vector r, int index) {
            this.e = e;
            this.v0 = v0;
            this.n = n;
            this.r = r;
            this.index = index;
            this.inergrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                    RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        }

        @Override
        public double value(double phi) {
            double tmp;
            if (Thread.currentThread().isInterrupted()) {
                Thread.currentThread().interrupt();
                return 0;
            }
            try {
                tmp = inergrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS,
                        new UnivariateFrequencyPolarizationSpreadInner(phi, e, v0, n, r, index), 0, INT_RANGE * eb.getSpread());
            } catch (TooManyEvaluationsException ex) {
                return 0;
            }
            return new Double(tmp).isNaN() ? 0 : tmp;
        }
    }

    /**
     * An auxiliary class for Romberg integrator for polarization calculations
     */
    private class UnivariateFrequencyPolarizationSpreadInner implements UnivariateFunction {

        private final double snphi, csphi, e;
        private final int index;
        private final Vector n, v0, r;

        public UnivariateFrequencyPolarizationSpreadInner(double phi, double e, Vector v0, Vector n, Vector r, int index) {
            this.snphi = Math.sin(phi);
            this.csphi = Math.cos(phi);
            this.e = e;
            this.n = n;
            this.r = r;
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
            u = sn * directionFrequencyPolarizationNoSpread(n, v, r, e)[index] * eb.angleDistribution(dv.get(0), dv.get(1))
                    + getShiftfactor() * SHIFT;
            return new Double(u).isNaN() ? 0 : u;
        }
    }

    /**
     * An auxiliary class for the brilliance calculations
     */
    protected class UnivariateVolumeFlux implements UnivariateFunction {

        Vector r0;
        Vector n0;

        public UnivariateVolumeFlux(Vector r0, Vector n0) {
            this.r0 = r0;
            this.n0 = n0;
        }

        @Override
        public double value(double x) {
            if (n0.get(0) + n0.get(1) + n0.get(2) == 0) {
                throw new LocalException(x);
            }
            return volumeFlux(r0.add(n0.multiply(x)));
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
