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

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.la4j.vector.Vector;
import org.apache.commons.math3.analysis.integration.*;
import org.la4j.matrix.Matrix;
import org.la4j.matrix.dense.Basic1DMatrix;
import org.la4j.vector.Vectors;
import org.la4j.vector.dense.BasicVector;

/*
 * The main class contating all physics of LEXG
 *
 * @author Ruslan Feshchenko
 * @version 1.2
 */
public class ThompsonSource implements Cloneable {

    ThompsonSource(LaserPulse l, ElectronBunch b) {
        this.lp = l;
        this.eb = b;
        calculateTotalFlux();
        calculateGeometricFactor();
    }

    /**
     * The number of columns in Shadow files
     */
    public static final int NUMBER_OF_COLUMNS = 18;

    /**
     * Angle range for rays exported for Shadow in the X-direction
     */
    private double rayXAnglerange = 0.05;

    /**
     * Angle range for rays exported for Shadow in the Y-direction
     */
    private double rayYAnglerange = 0.05;

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
    public int npGeometricFactor = 5000000;

    /**
     * Maximal number of evaluations in calculations of the brilliance
     */
    static final public int MAXIMAL_NUMBER_OF_EVALUATIONS = 30000;

    /**
     * Precision in calculations of the brilliance
     */
    public double precision = 0.0001;

    /**
     * Normalized total flux from the source
     */
    public double totalFlux;

    /**
     * Geometric factor. Assumes values from 0 to 1
     */
    public double geometricFactor = 1;

    /**
     * Flag - whether or not the electron beam transversal velocity spread is
     * taken into account
     */
    public boolean eSpread = false;
    
    /**
     * Flux in the phase space volume of ray generation
     */
    public double partialFlux;
    
    /**
     * Counter of ray iterations
     */
    public int counter;

    private LaserPulse lp;
    private ElectronBunch eb;

    @Override
    public Object clone() throws CloneNotSupportedException {
        Object tm = super.clone();
        ((ThompsonSource) tm).eb = (ElectronBunch) this.eb.clone();
        ((ThompsonSource) tm).lp = (LaserPulse) this.lp.clone();
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
     * Returns Laser Pulse reference
     *
     * @return
     */
    public LaserPulse getLaserPulse() {
        return lp;
    }

    /**
     * The full Thompson cross-section
     */
    public final static double SIGMA_T = 6.65e-29; /* Thompson cross-section, m2 */


    /**
     * A method calculating normalized total flux
     */
    public void calculateTotalFlux() {
        this.totalFlux = SIGMA_T * eb.number * lp.getPhotonNumber()
                * lp.fq / Math.PI / Math.sqrt((lp.getWidth2(0.0) + eb.getxWidth2(0.0))
                        * (lp.getWidth2(0.0) + eb.getyWidth2(0.0)));
    }

    /**
     * A method calculating the geometric factor
     */
    public void calculateGeometricFactor() {
        Vector iter = new BasicVector(new double[]{0.0, 0.0, 0.0});
        double sum = 0, wdx, wdy, len;
        int mult = 2;
        wdx = mult * Math.max(eb.getxWidth(0.0) + Math.abs(eb.shift.get(0)) / 2, lp.getWidth(0.0) + Math.abs(eb.shift.get(0)) / 2);
        wdy = mult * Math.max(eb.getyWidth(0.0) + Math.abs(eb.shift.get(1)) / 2, lp.getWidth(0.0) + Math.abs(eb.shift.get(1)) / 2);
        len = mult * Math.max(eb.length + Math.abs(eb.shift.get(2)) / 2, lp.length + Math.abs(eb.shift.get(2)) / 2);

        for (int i = 0; i < npGeometricFactor; i++) {
            iter.set(0, eb.shift.get(0) / 2 + wdx * (2 * Math.random() - 1.0));
            iter.set(1, eb.shift.get(1) / 2 + wdy * (2 * Math.random() - 1.0));
            iter.set(2, eb.shift.get(2) / 2 + len * (2 * Math.random() - 1.0));
            sum += volumeFlux(iter);
        }
        this.geometricFactor = 8 * wdx * wdy * len * sum / npGeometricFactor;
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
        if (eSpread) {
            return directionFrequencyFluxSpread(n, v, e);
        } else {
            return directionFrequencyFluxNoSpread(n, v, e);
        }
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
        double u, K, th;
        th = (1 - n.innerProduct(v)) * 2;
        K = Math.pow((Math.sqrt(e / lp.getPhotonEnergy() / (1 - e * th / lp.getPhotonEnergy() / 4)) - 2 * eb.getGamma()), 2)
                / 4 / Math.pow(eb.getGamma() * eb.delgamma, 2);
        u = totalFlux * e * 3.0 / 64 / Math.PI / Math.sqrt(Math.PI) / eb.delgamma / eb.getGamma() / lp.getPhotonEnergy()
                * Math.sqrt(e / lp.getPhotonEnergy()) * (Math.pow((1 - e * th / lp.getPhotonEnergy() / 2), 2) + 1)
                / Math.sqrt(1 - e * th / lp.getPhotonEnergy() / 4) * Math.exp(-K);
        return u;
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
        double u;
        BaseAbstractUnivariateIntegrator integrator = new RombergIntegrator(precision, RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateFunction func
                = new UnivariateFrequencyFluxSpreadOuter(e, v, n);
        try {
            u = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, 0.0, 2 * Math.PI);
            return u;
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
    }

    class UnivariateFrequencyFluxSpreadOuter implements UnivariateFunction {

        private final double e;
        private final Vector n, v0;
        private final BaseAbstractUnivariateIntegrator inergrator;

        public UnivariateFrequencyFluxSpreadOuter(double e, Vector v0, Vector n) {
            this.e = e;
            this.v0 = v0;
            this.n = n;
            this.inergrator = new RombergIntegrator(precision, RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                    RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        }

        @Override
        public double value(double phi) {
            double u;
            UnivariateFunction func
                    = new UnivariateFrequencyFluxSpreadInner(phi, e, v0, n);
            try {
                u = inergrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, 0.0, 3 * eb.getSpread());
                return u;
            } catch (TooManyEvaluationsException ex) {
                return 0;
            }
        }
    }

    class UnivariateFrequencyFluxSpreadInner implements UnivariateFunction {

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
            Vector dv;
            Vector v = new BasicVector(new double[]{Math.sin(theta) * Math.cos(phi),
                Math.sin(theta) * Math.sin(phi), Math.cos(theta)});
            dv = v.subtract(v0);
            u = theta * directionFrequencyFluxNoSpread(n, v, e)
                    * eb.angleDistribution(dv.get(0), dv.get(1));
            if ((new Double(u)).isNaN()) {
                return 0;
            }
            return u;
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
        if (eSpread) {
            return directionFrequencyVolumeFluxSpread(r, n, v, e);
        } else {
            return directionFrequencyVolumeFluxNoSpread(r, n, v, e);
        }
    }

    /**
     * A method calculating the flux density in a given direction for a given
     * X-ray photon energy for a given volume element without taking into
     * account electron transversial pulse spread
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
     * An auxiliary method calculating volume density of the X-ray source
     *
     * @param r spatial position
     * @return
     */
    public double volumeFlux(Vector r) {
        double u, z0, z, z1, x0, x, x1, y0, y, y1, sn, cs, K, len;
        len = Math.sqrt(lp.length * lp.length
                + eb.length * eb.length);
        sn = lp.direction.get(1);
        cs = lp.direction.get(2);
        x0 = eb.shift.get(0);
        y0 = eb.shift.get(1);
        z0 = eb.shift.get(2);
        x = r.get(0);
        y = r.get(1);
        z = r.get(2);
        x1 = x;
        y1 = -sn * z + cs * y;
        z1 = cs * z + sn * y;
        K = Math.pow((z + z1 - z0 - lp.delay) / len, 2)
                + Math.pow((x - x0), 2) / eb.getxWidth2(z - z0) + Math.pow((y - y0), 2) / eb.getyWidth2(z - z0)
                + (Math.pow(x1, 2) + Math.pow(y1, 2)) / lp.getWidth2(z1);
        u = 2.0 / Math.pow(Math.PI, 1.5) * Math.sqrt((lp.getWidth2(0.0)
                + eb.getxWidth2(0.0)) * (lp.getWidth2(0.0)
                + eb.getyWidth2(0.0))) / len / lp.getWidth2(z1) / eb.getxWidth(z - z0) / eb.getyWidth(z - z0) * Math.exp(-K);
        if ((new Double(u)).isNaN()) {
            u = 0;
        }
        return u;
    }

    /**
     * A method giving the flux density in a given direction
     *
     * @param n direction
     * @param v normalized electron velocity
     * @return
     */
    public double directionFlux(Vector n, Vector v) {
        double u, gamma2, th;
        th = (1 - n.innerProduct(v)) * 2;
        gamma2 = eb.getGamma() * eb.getGamma();
        u = totalFlux * 3.0 / 2 / Math.PI * gamma2 * (1 + Math.pow(th * gamma2, 2))
                / Math.pow((1 + gamma2 * th), 4) * geometricFactor;
        return u;
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
        if (eSpread) {
            return directionFrequencyBrillianceSpread(r0, n, v, e);
        } else {
            return directionFrequencyBrillianceNoSpread(r0, n, v, e);
        }
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
        RombergIntegrator integrator = new RombergIntegrator(precision, RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateVolumeFlux func = new UnivariateVolumeFlux(r0, n);
        try {
            u = integrator.integrate(30000, func,
                    r0.fold(Vectors.mkEuclideanNormAccumulator()) - 3 * eb.length,
                    r0.fold(Vectors.mkEuclideanNormAccumulator()) + 3 * eb.length);
            u = u * directionFrequencyFluxNoSpread(n, v, e);
            return u;
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
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
        RombergIntegrator integrator = new RombergIntegrator(precision, RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateVolumeFlux func = new UnivariateVolumeFlux(r0, n);
        try {
            u = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func,
                    r0.fold(Vectors.mkEuclideanNormAccumulator()) - 3 * eb.length,
                    r0.fold(Vectors.mkEuclideanNormAccumulator()) + 3 * eb.length);
            u = u * directionFrequencyFluxSpread(n, v, e);
            return u;
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
    }

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
            if (n0.get(0) + n0.get(1) + n0.get(2) == 0f) {
                throw new LocalException(x);
            }
            return y;
        }
    }

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
     */
    public double[] getRay() {
        double[] ray = new double[NUMBER_OF_COLUMNS];
        Matrix M, D, T, A, I = new Basic1DMatrix(3, 3);
        I.set(0, 0, 1.0);
        I.set(1, 1, 1.0);
        I.set(2, 2, 1.0);
        Vector n = new BasicVector(new double[]{0.0, 0.0, 1.0});
        Vector r = new BasicVector(new double[]{0.0, 0.0, 0.0});
        Vector n0 = new BasicVector(new double[]{0.0, 1.0, 0.0}), As;
        double prob0, prob, EMax, mult = 2, innerProduct, factor, sum=0;
        EMax = directionEnergy(n, n);
        factor = 64 * Math.max(eb.getxWidth(0.0), lp.getWidth(0.0)) * Math.max(eb.getyWidth(0.0), lp.getWidth(0.0))
                * Math.max(eb.length, lp.length) * 4 * rayXAnglerange * rayYAnglerange
                * (maxEnergy / minEnergy - 1);
        prob0 = directionFrequencyVolumeFluxNoSpread(r, n, new BasicVector(new double[]{0.0, 0.0, 1.0}), EMax);
        if (eSpread) {
            prob0 *= eb.angleDistribution(0, 0);
            factor *= 4 * mult * mult * eb.getXSpread() * eb.getYSpread();
        }
        do {
            ray[0] = 2 * (2 * Math.random() - 1.0) * Math.max(eb.getxWidth(0.0), lp.getWidth(0.0));
            ray[2] = 2 * (2 * Math.random() - 1.0) * Math.max(eb.getyWidth(0.0), lp.getWidth(0.0));
            ray[1] = 2 * (2 * Math.random() - 1.0) * Math.max(eb.length, lp.length);
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
            if (eSpread) {
                double thetax = mult * eb.getXSpread() * (2 * Math.random() - 1);
                double thetay = mult * eb.getYSpread() * (2 * Math.random() - 1);
                Vector v = new BasicVector(new double[]{thetax, thetay, Math.sqrt(1 - thetax * thetax - thetay * thetay)});
                prob = directionFrequencyVolumeFluxNoSpread(r, n, v, ray[10]) * eb.angleDistribution(thetax, thetay);
            } else {
                prob = directionFrequencyVolumeFluxNoSpread(r, n, new BasicVector(new double[]{0.0, 0.0, 1.0}), ray[10]);
            }
            sum += prob;
            counter++;
        } while (prob / prob0 < Math.random() || (new Double(prob)).isNaN());
        // Calculation of the rotated polarization vector
        n = new BasicVector(new double[]{ray[3], ray[4], ray[5]});
        innerProduct = n.innerProduct(n0);
        D = n.outerProduct(n0).add(n0.outerProduct(n)).multiply(innerProduct).subtract(n.outerProduct(n).
                add(n0.outerProduct(n0)).divide(innerProduct * innerProduct - 1.0));
        A = n.outerProduct(n0).subtract(n0.outerProduct(n)).add(I.multiply(innerProduct));
        T = I.subtract(D).multiply(1 - innerProduct).add(A);
        As = T.multiply(new BasicVector(new double[]{1.0, 0.0, 0.0}));
        ray[6] = As.get(0);
        ray[7] = As.get(1);
        ray[8] = As.get(2);
        //Setting other columns
        ray[9] = 1.0;
        ray[13] = Math.random() * 2 * Math.PI;
        partialFlux += sum * factor;
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
}
