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

import electronbunch.ElectronBunch;
import org.la4j.Vector;
import org.la4j.vector.dense.BasicVector;

/**
 *
 * @author Ruslan
 * @version 1.0
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
    
    /**
     * Laser beam intensity, W/m^2
     */
    private double intensity = 0;

    protected double rk;

    public AbstractLaserPulse() {
        this.photonenergy = 1.1 * ElectronBunch.E;
        this.setPulseEnergy(2.0e-2);
        this.rk = HC / this.photonenergy;
        this.direction = new BasicVector(new double[]{0.0, 0.0, 1.0});
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
     * @return the intensity, W/m^2
     */
    public double getIntensity() {
        return intensity;
    }

    /**
     * @param intensity the intensity to set, W/m^2
     */
    public void setIntensity(double intensity) {
        this.intensity = intensity;
    }
}
