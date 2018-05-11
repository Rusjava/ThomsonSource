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
package laserpulse;

import electronbunch.GaussianElectronBunch;
import org.la4j.Vector;
import org.la4j.Vectors;
import org.la4j.vector.dense.BasicVector;

/**
 *
 * @author Ruslan
 * @version 2.0
 */
public abstract class AbstractLaserPulse implements Cloneable {

    /**
     * hc constant
     */
    public static final double HC = 3.1614e-26;
    /**
     * The speed of light, m/c
     */
    public static final double C = 299792458;
    /**
     * Photon energy, J
     */
    private double photonenergy;
    /**
     * Number of photon in laser pulse
     */
    private double number;
    /**
     * Laser pulse semi-length, m
     */
    private double length = 4.5e-3;
    /**
     * Mean direction of the laser pulse
     */
    private Vector direction;
    /**
     * Laser pulse Rayleigh length, m
     */
    private double rlength = 5.4e-3;
    /**
     * Pulse frequency, 1/s
     */
    private double fq = 7.9e7;
    /**
     * Laser pulse delay, m
     */
    private double delay = 0;
    /**
     * Stocks parameters
     */
    private double ksi1 = 0;
    private double ksi2 = 0;
    private double ksi3 = 0;

    private double[] KA1;
    private double[] KA2;
    private double p;
    private Vector[] A1;
    private Vector[] A2;

    protected double rk;

    public AbstractLaserPulse() {
        this.photonenergy = 1.1 * GaussianElectronBunch.E;
        this.setPulseEnergy(2.0e-2);
        this.rk = HC / this.photonenergy;
        this.direction = new BasicVector(new double[]{0.0, 0.0, 1.0});
        this.KA1 = new double[2];
        this.KA2 = new double[2];
        this.A1 = new Vector[]{new BasicVector(new double[]{1.0, 0.0}), new BasicVector(new double[]{1.0, 0.0})};
        this.A2 = new Vector[]{new BasicVector(new double[]{1.0, 0.0}), new BasicVector(new double[]{1.0, 0.0})};
    }

    @Override
    public Object clone() throws CloneNotSupportedException {
        Object tm = super.clone();
        ((AbstractLaserPulse) tm).setDirection(new BasicVector(new double[]{0.0, 0.0, 0.0}));
        ((AbstractLaserPulse) tm).getDirection().set(0, this.getDirection().get(0));
        ((AbstractLaserPulse) tm).getDirection().set(1, this.getDirection().get(1));
        ((AbstractLaserPulse) tm).getDirection().set(2, this.getDirection().get(2));
        return tm;
    }

    /**
     * Returns the laser bunch width at position z, m
     *
     * @param z coordinate, m
     * @return
     */
    public abstract double getWidth(double z);

    /**
     * Sets the laser bunch width
     *
     * @param w width, m
     */
    public void setWidth(double w) {
        setRlength(2 * w * w / rk);
    }

    /**
     * Returns the laser bunch width at position z, m
     *
     * @param z coordinate z, m
     * @return
     */
    public abstract double getWidth2(double z);

    /**
     * Sets the energy of laser photons
     *
     * @param e photon energy, J
     */
    public void setPhotonEnergy(double e) {
        this.photonenergy = e;
        this.rk = HC / this.photonenergy;
    }

    /**
     * Return the energy of laser photons, J
     *
     * @return
     */
    public double getPhotonEnergy() {
        return this.photonenergy;
    }

    /**
     * Sets the number of photons in laser pulse
     *
     * @param n number of photons
     */
    public void setPhotonNumber(double n) {
        this.number = n;
    }

    /**
     * Returns the number of photons in laser pulse
     *
     * @return
     */
    public double getPhotonNumber() {
        return this.number;
    }

    /**
     * Sets the energy of laser pulse (by setting the photon number)
     *
     * @param e energy of laser pulse, J
     */
    public final void setPulseEnergy(double e) {
        this.number = e / this.photonenergy;
    }

    /**
     * Returns the energy of laser pulse
     *
     * @return
     */
    public double getPulseEnergy() {
        return this.number * this.photonenergy;
    }

    /**
     * Laser pulse semi-length, m
     *
     * @return the length
     */
    public double getLength() {
        return length;
    }

    /**
     * Laser pulse semi-length, m
     *
     * @param length the length to set
     */
    public void setLength(double length) {
        this.length = length;
    }

    /**
     * Mean direction of the laser pulse
     *
     * @return the direction
     */
    public Vector getDirection() {
        return direction;
    }

    /**
     * Mean direction of the laser pulse
     *
     * @param direction the direction to set
     */
    public void setDirection(Vector direction) {
        this.direction = direction;
    }

    /**
     * Laser pulse Rayleigh length, m
     *
     * @return the rlength
     */
    public double getRlength() {
        return rlength;
    }

    /**
     * Laser pulse Rayleigh length, m
     *
     * @param rlength the rlength to set
     */
    public void setRlength(double rlength) {
        this.rlength = rlength;
    }

    /**
     * Pulse frequency, 1/s
     *
     * @return the fq
     */
    public double getFq() {
        return fq;
    }

    /**
     * Pulse frequency, 1/s
     *
     * @param fq the fq to set
     */
    public void setFq(double fq) {
        this.fq = fq;
    }

    /**
     * Laser pulse delay, m
     *
     * @return the delay
     */
    public double getDelay() {
        return delay;
    }

    /**
     * Laser pulse delay, m
     *
     * @param delay the delay to set
     */
    public void setDelay(double delay) {
        this.delay = delay;
    }

    /**
     * Setting laser polarization state
     *
     * @param ksi1
     * @param ksi2
     * @param ksi3
     */
    public void setPolarization(double ksi1, double ksi2, double ksi3) {
        this.ksi1 = ksi1;
        this.ksi2 = ksi2;
        this.ksi3 = ksi3;
        this.p = Math.sqrt(ksi1 * ksi1 + ksi2 * ksi2 + ksi2 * ksi3);
        double[] t = new double[]{p + ksi3, p - ksi3};
        double[] kappa = new double[]{1 + p, 1 - p};
        Vector[] AA1 = new Vector[2];
        Vector[] AA2 = new Vector[2];
        double a1, a2, M, c1, c2, h, coef;

        //Auxialiry paremeters
        for (int s = 0; s < 2; s++) {
            h = Math.sqrt(ksi1 * ksi1 + ksi2 * ksi2 + t[s] * t[s]);
            AA1[s] = new BasicVector(new double[]{ksi1, t[s]});
            AA2[s] = new BasicVector(new double[]{-ksi2, 0});
            AA1[s] = AA1[s].divide(h);
            AA2[s] = AA2[s].divide(h);
        }
        //Orthogonal polarization vectors
        for (int s = 0; s < 2; s++) {
            if (AA1[s].innerProduct(AA2[s]) == 0) {
                A1[s] = AA1[s].copy();
                A2[s] = AA2[s].copy();
            } else {
                a1 = AA1[s].fold(Vectors.mkEuclideanNormAccumulator());
                a2 = AA2[s].fold(Vectors.mkEuclideanNormAccumulator());
                M = (a1 - a2)
                        / Math.sqrt(Math.pow(a1, 2) + Math.pow(a2, 2) - 2 * a1 * a2
                                + 4 * Math.pow(AA1[s].innerProduct(AA2[s]), 2));
                c1 = (1 + M) / 2;
                c2 = (1 - M) / 2;
                this.A1[s].set(0, c1 * AA1[s].get(0) - c2 * AA2[s].get(0));
                this.A1[s].set(1, c1 * AA1[s].get(1) - c2 * AA2[s].get(1));
                this.A2[s].set(0, c1 * AA2[s].get(0) + c2 * AA1[s].get(0));
                this.A2[s].set(1, c1 * AA2[s].get(1) + c2 * AA1[s].get(1));
                coef=Math.sqrt(4*Math.PI*getIntensity()*kappa[s]/this.getPhotonEnergy()/C*HC);
                this.A1[s].divide(coef);
                this.A2[s].divide(coef);
            }
            //Intensities of orthogonal polarizations
            this.KA1[s] = A1[s].fold(Vectors.mkEuclideanNormAccumulator());
            this.KA2[s] = A2[s].fold(Vectors.mkEuclideanNormAccumulator());
        }
    }

    /**
     * Getting laser polarization state vector
     *
     * @return
     */
    public double[] getPolarization() {
        return new double[]{ksi1, ksi2, ksi3};
    }

    /**
     * The transversal spatial distribution of photons in the pulse
     *
     * @param r
     * @return
     */
    public abstract double tSpatialDistribution(Vector r);

    /**
     * The longitudinal spatial distribution of photons in the pulse
     *
     * @param r
     * @return
     */
    public abstract double lSpatialDistribution(Vector r);

    /**
     * Getting an average pulse intensity
     *
     * @return the intensity, W/m^2
     */
    public double getIntensity() {
        return getPulseEnergy() / getLength() / Math.PI / getWidth2(0);
    }

    /**
     * @return the KA1
     */
    public double[] getKA1() {
        return KA1;
    }

    /**
     * @return the KA2
     */
    public double[] getKA2() {
        return KA2;
    }

    /**
     * @return the p
     */
    public double getP() {
        return p;
    }

    /**
     * @return the A1
     */
    public Vector[] getA1() {
        return A1;
    }

    /**
     * @return the A2
     */
    public Vector[] getA2() {
        return A2;
    }
}
