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

import electronbunch.AbstractElectronBunch;
import laserpulse.AbstractLaserPulse;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.RombergIntegrator;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.la4j.Vector;
import org.la4j.Vectors;

/**
 * The main class containing all physics of LEXG in linear approximation
 *
 * @author Ruslan Feshchenko
 * @version 3.0
 */
public class LinearThomsonSource extends AbstractThomsonSource {

    /**
     * Constructor
     *
     * @param l
     * @param b
     */
    public LinearThomsonSource(AbstractLaserPulse l, AbstractElectronBunch b) {
        super(l,b);
        calculateLinearTotalFlux();
        calculateGeometricFactor();
    }

    @Override
    public double directionFrequencyFluxNoSpread(Vector n, Vector v, Vector r, double e) {
        double K, th;
        th = (1 - n.innerProduct(v)) * 2;
        K = Math.pow((Math.sqrt(e / lp.getPhotonEnergy() / (1 - e * th / lp.getPhotonEnergy() / 4)) - 2 * eb.getGamma()), 2)
                / 4 / Math.pow(eb.getGamma() * eb.getDelGamma(), 2);
        return getLinearTotalFlux() * e * 3.0 / 64 / Math.PI / Math.sqrt(Math.PI) / eb.getDelGamma() / eb.getGamma() / lp.getPhotonEnergy()
                * Math.sqrt(e / lp.getPhotonEnergy()) * (Math.pow((1 - e * th / lp.getPhotonEnergy() / 2), 2) + 1)
                / Math.sqrt(1 - e * th / lp.getPhotonEnergy() / 4) * Math.exp(-K);
    }

    @Override
    public double[] directionFrequencyPolarizationNoSpread(Vector n, Vector v, Vector r, double e) {
        double[] array = new double[NUMBER_OF_POL_PARAM];
        double K, th, m11, m22, m12, mlt, cs, sn;
        th = (1 - n.innerProduct(v)) * 2;
        mlt = 1 - e * th / lp.getPhotonEnergy() / 2;
        K = Math.pow((Math.sqrt(e / lp.getPhotonEnergy() / (1 - e * th / lp.getPhotonEnergy() / 4)) - 2 * eb.getGamma()), 2)
                / 4 / Math.pow(eb.getGamma() * eb.getDelGamma(), 2);
        m11 = getLinearTotalFlux() * e * 3.0 / 32 / Math.PI / Math.sqrt(Math.PI) / eb.getDelGamma() / eb.getGamma() / lp.getPhotonEnergy()
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

    @Override
    public double directionFlux(Vector n, Vector v) {
        double gamma2, th;
        th = (1 - n.innerProduct(v)) * 2;
        gamma2 = eb.getGamma() * eb.getGamma();
        return getLinearTotalFlux() * 3.0 / 2 / Math.PI * gamma2 * (1 + Math.pow(th * gamma2, 2))
                / Math.pow((1 + gamma2 * th), 4) * getGeometricFactor();
    }

    @Override
    public double directionEnergy(Vector n, Vector v) {
        double mv;
        mv = Math.sqrt(1.0 - 1.0 / eb.getGamma() / eb.getGamma());
        return 2 * lp.getPhotonEnergy() / (1 - n.innerProduct(v) * mv);
    }

    @Override
    public double directionFrequencyBrillianceNoSpread(Vector r0, Vector n, Vector v, double e) {
        double u;
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateFunction func = new UnivariateVolumeFlux(r0, n);
        try {
            u = integrator.integrate(30000, func, r0.fold(Vectors.mkEuclideanNormAccumulator()) - 3 * eb.getLength(), r0.fold(Vectors.mkEuclideanNormAccumulator()) + 3 * eb.getLength());
            return u * directionFrequencyFluxNoSpread(n, v, null, e);
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
    }

    @Override
    public double directionFrequencyBrillianceSpread(Vector r0, Vector n, Vector v, double e) {
        double u;
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateFunction func = new UnivariateVolumeFlux(r0, n);
        try {
            u = integrator.integrate(MAXIMAL_NUMBER_OF_EVALUATIONS, func, r0.fold(Vectors.mkEuclideanNormAccumulator()) - 3 * eb.getLength(), r0.fold(Vectors.mkEuclideanNormAccumulator()) + 3 * eb.getLength());
            return u * directionFrequencyFluxSpread(n, v, null, e);
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }
    }

    @Override
    public double[] directionFrequencyBrilliancePolarizationNoSpread(Vector r0, Vector n, Vector v, double e) {
        double mlt;
        double[] array = new double[AbstractThomsonSource.NUMBER_OF_POL_PARAM];
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateFunction func = new UnivariateVolumeFlux(r0, n);
        try {
            mlt = integrator.integrate(30000, func, r0.fold(Vectors.mkEuclideanNormAccumulator()) - 3 * eb.getLength(), r0.fold(Vectors.mkEuclideanNormAccumulator()) + 3 * eb.getLength());
        } catch (TooManyEvaluationsException ex) {
            mlt = 0;
        }
        for (int i = 0; i < AbstractThomsonSource.NUMBER_OF_POL_PARAM; i++) {
            array[i] = mlt * directionFrequencyPolarizationNoSpread(n, v, null, e)[i];
        }
        return array;
    }

    @Override
    public double[] directionFrequencyBrilliancePolarizationSpread(Vector r0, Vector n, Vector v, double e) throws InterruptedException {
        double mlt;
        double[] array = new double[AbstractThomsonSource.NUMBER_OF_POL_PARAM];
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateFunction func = new UnivariateVolumeFlux(r0, n);
        try {
            mlt = integrator.integrate(30000, func, r0.fold(Vectors.mkEuclideanNormAccumulator()) - 3 * eb.getLength(), r0.fold(Vectors.mkEuclideanNormAccumulator()) + 3 * eb.getLength());
        } catch (TooManyEvaluationsException ex) {
            mlt = 0;
        }
        for (int i = 0; i < AbstractThomsonSource.NUMBER_OF_POL_PARAM; i++) {
            array[i] = mlt * directionFrequencyPolarizationSpread(n, v, null, e)[i];
        }
        return array;
    }

    @Override
    public double directionFrequencyPolarizationNoSpread(Vector n, Vector v, Vector r, double e, int index) {
        return directionFrequencyPolarizationNoSpread(n, v, r, e)[index];
    }

    @Override
    public double directionFrequencyBrilliancePolarizationNoSpread(Vector r0, Vector n, Vector v, double e, int index) {
        double mlt;
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateFunction func = new UnivariateVolumeFlux(r0, n);
        try {
            mlt = integrator.integrate(30000, func, r0.fold(Vectors.mkEuclideanNormAccumulator()) - 3 * eb.getLength(), r0.fold(Vectors.mkEuclideanNormAccumulator()) + 3 * eb.getLength());
        } catch (TooManyEvaluationsException ex) {
            mlt = 0;
        }
        return mlt * directionFrequencyPolarizationNoSpread(n, v, null, e)[index];
    }

    @Override
    public double directionFrequencyBrilliancePolarizationSpread(Vector r0, Vector n, Vector v, double e, int index) throws InterruptedException {
        double mlt;
        RombergIntegrator integrator = new RombergIntegrator(getPrecision(), RombergIntegrator.DEFAULT_ABSOLUTE_ACCURACY, RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT, RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
        UnivariateFunction func = new UnivariateVolumeFlux(r0, n);
        try {
            mlt = integrator.integrate(30000, func, r0.fold(Vectors.mkEuclideanNormAccumulator()) - 3 * eb.getLength(), r0.fold(Vectors.mkEuclideanNormAccumulator()) + 3 * eb.getLength());
        } catch (TooManyEvaluationsException ex) {
            mlt = 0;
        }
        return mlt * directionFrequencyPolarizationSpread(n, v, null, e)[index];
    }

    /**
     * An auxiliary class for the brilliance calculations
     */
    private class UnivariateVolumeFlux implements UnivariateFunction {

        Vector r0;
        Vector n0;

        public UnivariateVolumeFlux(Vector r0, Vector n0) {
            super();
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
}
