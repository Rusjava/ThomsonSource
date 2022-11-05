/*
 * Copyright (BB) 2017 Ruslan
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
 * @version 2.01
 */
public abstract class AbstractLaserPulse implements Cloneable {

    /**
     * hc constant, kg*m^2/s
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
    private double length = 1.5e-3;
    /**
     * Mean direction of the laser pulse
     */
    private Vector direction;
    /**
     * Laser pulse Rayleigh length, m
     */
    private double rlength = 0.35e-3;
    /**
     * Pulse frequency, 1/s
     */
    private double fq = 1e3;
    /**
     * Laser pulse delay, m
     */
    private double delay;
    /**
     * Stocks parameters
     */
    private double ksi1 = 1;
    private double ksi2 = 0;
    private double ksi3 = 0;

    private double[] KA1;
    private double[] KA2;
    private double p;
    private Vector[] A1;
    private Vector[] A2;
    private double intensity = 1; //Average laser pulse intensity
    protected double rk;

    public AbstractLaserPulse() {
        this.photonenergy = 1.204 * GaussianElectronBunch.E;
        setPulseEnergy(0.1);
        this.rk = HC / photonenergy;
        this.direction = new BasicVector(new double[]{0.0, Math.sin(0.052), Math.cos(0.052)});
        this.KA1 = new double[2];
        this.KA2 = new double[2];
        this.A1 = new Vector[]{new BasicVector(new double[]{0, 0, 0}), new BasicVector(new double[]{0, 0, 0})};
        this.A2 = new Vector[]{new BasicVector(new double[]{0, 0, 0}), new BasicVector(new double[]{0, 0, 0})};
        setAverageIntensity();
        setPolarization(ksi1, ksi2, ksi3);
    }

    @Override
    public Object clone() throws CloneNotSupportedException {
        Object tm = super.clone();
        ((AbstractLaserPulse) tm).setDirection(new BasicVector(new double[]{this.getDirection().get(0),
            this.getDirection().get(1), this.getDirection().get(2)}));
        ((AbstractLaserPulse) tm).KA1 = new double[]{KA1[0], KA1[1]};
        ((AbstractLaserPulse) tm).KA2 = new double[]{KA2[0], KA2[1]};
        ((AbstractLaserPulse) tm).A1 = new Vector[]{A1[0].copy(), A1[1].copy()};
        ((AbstractLaserPulse) tm).A2 = new Vector[]{A2[0].copy(), A2[1].copy()};
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
        setAverageIntensity();
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
        setAverageIntensity();
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
        setAverageIntensity();
    }

    /**
     * Mean direction of the laser pulse
     *
     * @return the direction
     */
    public final Vector getDirection() {
        return direction.copy();
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
        setAverageIntensity();
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
        this.p = Math.sqrt(ksi1 * ksi1 + ksi2 * ksi2 + ksi3 * ksi3);
        double[] t = new double[]{ksi3 + p, ksi3 - p};
        //An array for polarization intensities
        double[] coef = new double[]{Math.sqrt((1 + p) / 2), Math.sqrt((1 - p) / 2)};
        //Arrays for real and imagenery parts of polarization vectors
        Vector AA1, AA2;
        Vector[] B;
        double h;
        //Auxialiry parameters
        for (int s = 0; s < 2; s++) {
            h = Math.sqrt(ksi1 * ksi1 + ksi2 * ksi2 + t[s] * t[s]);
            if (p == 0) {
                //If unpolarized then special treatment
                AA1 = (new BasicVector(new double[]{1, 1 - 2 * s, 0})).divide(Math.sqrt(2));
                AA2 = new BasicVector(new double[]{0, 0, 0});
            } else {
                if (h == 0) {
                    AA1 = new BasicVector(new double[]{1, 0, 0});
                    AA2 = new BasicVector(new double[]{0, 0, 0});
                } else {
                    AA1 = (new BasicVector(new double[]{ksi1, t[s], 0})).divide(h);
                    AA2 = (new BasicVector(new double[]{-ksi2, 0, 0})).divide(h);
                }
            }
            //Orthogonalizing real and imagenery parts of polarization vectors
            B = getOrthogonalPolarizationVectors(AA1, AA2);
            this.A1[s] = B[0];
            this.A2[s] = B[1];
            //Normalizing to the intensity
            this.A1[s] = A1[s].multiply(coef[s]);
            this.A2[s] = A2[s].multiply(coef[s]);
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
    public double getAverageIntensity() {
        return intensity;
    }

    /**
     * Setting an average pulse intensity and updating polarization vectors
     *
     */
    public final void setAverageIntensity() {
        this.intensity = getPulseEnergy() / getLength() / Math.PI / getWidth2(0) * AbstractLaserPulse.C;
    }

    /**
     * Getting coordinate dependant pulse intensity
     *
     * @param r
     * @return
     */
    public double getIntensity(Vector r) {
        return getPulseEnergy() * tSpatialDistribution(r) * lSpatialDistribution(r) * AbstractLaserPulse.C;
    }

    /**
     * @return the KA1
     */
    public double[] getKA1() {
        return new double[]{KA1[0], KA1[1]};
    }

    /**
     * @return the KA2
     */
    public double[] getKA2() {
        return new double[]{KA2[0], KA2[1]};
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
        return new Vector[]{A1[0].copy(), A1[1].copy()};
    }

    /**
     * @return the A2
     */
    public Vector[] getA2() {
        return new Vector[]{A2[0].copy(), A2[1].copy()};
    }

    /**
     * Orthogonalization of the real and imaginary parts of a complex vector
     *
     * @param AA1 real part
     * @param AA2 imaginary part
     * @return array of orthogonal vectors
     */
    public Vector[] getOrthogonalPolarizationVectors(Vector AA1, Vector AA2) {
        Vector B1 = new BasicVector(new double[]{0, 0, 0});
        Vector B2 = new BasicVector(new double[]{0, 0, 0});
        double c1, c2, M, a1, a2;
        if (AA1.innerProduct(AA2) == 0) {
            //If already orthogonal, return initial vectors
            return new Vector[]{AA1, AA2};
        } else {
            //If not orthogonal, then making them orthogonal
            a1 = AA1.innerProduct(AA1);
            a2 = AA2.innerProduct(AA2);
            M = Math.abs(a1 - a2) / Math.sqrt(a1 * a1 + a2 * a2 - 2 * a1 * a2
                    + 4 * Math.pow(AA1.innerProduct(AA2), 2));
            c1 = Math.sqrt((1 + M) / 2);
            c2 = Math.sqrt((1 - M) / 2);
            B1.set(0, c1 * AA1.get(0) - c2 * AA2.get(0));
            B1.set(1, c1 * AA1.get(1) - c2 * AA2.get(1));
            B2.set(0, c1 * AA2.get(0) + c2 * AA1.get(0));
            B2.set(1, c1 * AA2.get(1) + c2 * AA1.get(1));
            return new Vector[]{B1, B2};
        }
    }

    /**
     * Returning orthogonalized real and imaginary parts of a linear
     * superposition of two independent polarization states
     *
     * @param phase mutual phase of two polarizations
     * @return an array of polarization vectors (real and imaginary parts)
     */
    public Vector[] getPhaseWeightedPolarizationVectors1(double phase) {
        Vector B1 = new BasicVector(new double[]{0, 0, 0});
        Vector B2 = new BasicVector(new double[]{0, 0, 0});
        double sn = Math.sin(phase);
        double cs = Math.cos(phase);
        //Calculating a sum of two independent polarizations with phase factors
        for (int s = 0; s < 2; s++) {
            sn *= 1 - 2 * s;
            B1 = B1.add(A1[s].multiply(cs).subtract(A2[s].multiply(sn)));
            B2 = B2.add(A1[s].multiply(sn).add(A2[s].multiply(cs)));
        }
        return getOrthogonalPolarizationVectors(B1, B2);
    }

    /**
     * Returning orthogonalized real and imaginary parts of a linear
     * superposition of two independent polarization states
     *
     * @param phase mutual phase of two polarizations
     * @return an array of polarization vectors (real and imaginary parts)
     */
    public Vector[] getPhaseWeightedPolarizationVectors2(double phase) {
        Vector B1 = new BasicVector(new double[]{0, 0, 0});
        Vector B2 = new BasicVector(new double[]{0, 0, 0});
        double[] ff = new double[]{Math.cos(phase) * Math.sqrt(2), Math.sin(phase) * Math.sqrt(2)};
        //Calculating a sum of two independent polarizations with phase factors
        for (int s = 0; s < 2; s++) {
            B1 = B1.add(A1[s].multiply(ff[s]));
            B2 = B2.add(A2[s].multiply(ff[s]));
        }
        return getOrthogonalPolarizationVectors(B1, B2);
    }
    
    /**
     * Transforming coordinates form the electron bunch coordinate system to the laser pulse coordinate system
     * 
     * @param r a vector in the electron coordinate system
     * @return 
     */
    
    public Vector getTransformedCoordinates(Vector r) {
        double sn = getDirection().get(1);
        double cs = getDirection().get(2);
        //Defining new coordinates
        double x1 = r.get(0);
        double y1 = -sn * r.get(2) + cs * r.get(1);
        double z1 = cs * r.get(2) + sn * r.get(1);
        return new BasicVector(new double[]{x1, y1, -z1});
    }
}
