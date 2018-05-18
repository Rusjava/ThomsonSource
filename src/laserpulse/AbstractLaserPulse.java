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
    private double delay;
    /**
     * Stocks parameters
     */
    private double ksi1;
    private double ksi2;
    private double ksi3;

    private double[] KA1;
    private double[] KA2;
    private double p;
    private Vector[] A1;
    private Vector[] A2;
    private double intensity = 1; //Average laser pulse intensity
    protected double rk;

    public AbstractLaserPulse() {
        this.photonenergy = 1.1 * GaussianElectronBunch.E;
        this.number = 2.0e-2 / photonenergy;
        this.rk = HC / photonenergy;
        this.direction = new BasicVector(new double[]{0.0, 0.0, 1.0});
        this.KA1 = new double[2];
        this.KA2 = new double[2];
        this.A1 = new Vector[]{new BasicVector(new double[]{0, 0, 0}), new BasicVector(new double[]{0, 0, 0})};
        this.A2 = new Vector[]{new BasicVector(new double[]{0, 0, 0}), new BasicVector(new double[]{0, 0, 0})};
        setIntensity();
        setPolarization(ksi1, ksi2, ksi3);
    }

    @Override
    public Object clone() throws CloneNotSupportedException {
        Object tm = super.clone();
        ((AbstractLaserPulse) tm).setDirection(new BasicVector(new double[]{0.0, 0.0, 0.0}));
        ((AbstractLaserPulse) tm).getDirection().set(0, this.getDirection().get(0));
        ((AbstractLaserPulse) tm).getDirection().set(1, this.getDirection().get(1));
        ((AbstractLaserPulse) tm).getDirection().set(2, this.getDirection().get(2));
        ((AbstractLaserPulse) tm).KA1 = new double[]{KA1[0], KA1[1]};
        ((AbstractLaserPulse) tm).KA2 = new double[]{KA2[0], KA2[1]};
        ((AbstractLaserPulse) tm).A1 = new Vector[]{A1[0].copy(), A1[1].copy(), A1[2].copy()};
        ((AbstractLaserPulse) tm).A2 = new Vector[]{A2[0].copy(), A2[1].copy(), A2[2].copy()};
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
        setIntensity();
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
        setIntensity();
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
        setIntensity();
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
        setIntensity();
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
    public final void setPolarization(double ksi1, double ksi2, double ksi3) {
        this.ksi1 = ksi1;
        this.ksi2 = ksi2;
        this.ksi3 = ksi3;
        this.p = Math.sqrt(ksi1 * ksi1 + ksi2 * ksi2 + ksi2 * ksi3);
        double[] t = new double[]{ksi3 + p, ksi3 - p};
        double[] kappa = new double[]{1 + p, 1 - p};
        Vector[] AA1 = new Vector[2];
        Vector[] AA2 = new Vector[2];
        double a1, a2, M, c1, c2, h, coef;

        //Auxialiry parameters
        for (int s = 0; s < 2; s++) {
            if (p == 0) {
                AA1[s] = (new BasicVector(new double[]{1, 1 - 2 * s})).divide(Math.sqrt(2));
                AA2[s] = new BasicVector(new double[]{0, 0});
            } else {
                h = Math.sqrt(ksi1 * ksi1 + ksi2 * ksi2 + t[s] * t[s]);
                AA1[s] = new BasicVector(new double[]{ksi1, t[s]});
                AA2[s] = new BasicVector(new double[]{-ksi2, 0});
                AA1[s] = AA1[s].divide(h);
                AA2[s] = AA2[s].divide(h);
            }
        }
        //Orthogonal polarization vectors
        for (int s = 0; s < 2; s++) {
            if (AA1[s].innerProduct(AA2[s]) == 0) {
                this.A1[s].set(0, AA1[s].get(0));
                this.A1[s].set(1, AA1[s].get(1));
                this.A2[s].set(0, AA2[s].get(0));
                this.A2[s].set(1, AA2[s].get(1));
            } else {
                a1 = AA1[s].innerProduct(AA1[s]);
                a2 = AA2[s].innerProduct(AA2[s]);
                M = Math.abs(a1 - a2) / Math.sqrt(a1 * a1 + a2 * a2 - 2 * a1 * a2
                        + 4 * Math.pow(AA1[s].innerProduct(AA2[s]), 2));
                c1 = (1 + M) / 2;
                c2 = (1 - M) / 2;
                this.A1[s].set(0, c1 * AA1[s].get(0) - c2 * AA2[s].get(0));
                this.A1[s].set(1, c1 * AA1[s].get(1) - c2 * AA2[s].get(1));
                this.A2[s].set(0, c1 * AA2[s].get(0) + c2 * AA1[s].get(0));
                this.A2[s].set(1, c1 * AA2[s].get(1) + c2 * AA1[s].get(1));
            }
            coef = Math.sqrt(intensity * kappa[s] / 2);
            this.A1[s] = A1[s].multiply(coef);
            this.A2[s] = A2[s].multiply(coef);
            //Intensities of orthogonal polarizations
            this.KA1[s] = A1[s].innerProduct(A1[s]);
            this.KA2[s] = A2[s].innerProduct(A2[s]);
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
        return intensity;
    }

    /**
     * Setting an average pulse intensity and updating polarization vectors
     *
     */
    public final void setIntensity() {
        double newIntensity = getPulseEnergy() / getLength() / Math.PI / getWidth2(0) * AbstractLaserPulse.C;
        for (int s = 0; s < 2; s++) {
            this.A1[s] = A1[s].multiply(Math.sqrt(newIntensity / intensity));
            this.A2[s] = A2[s].multiply(Math.sqrt(newIntensity / intensity));
            this.KA1[s] = A1[s].innerProduct(A1[s]);
            this.KA2[s] = A2[s].innerProduct(A2[s]);
        }
        this.intensity = newIntensity;
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
