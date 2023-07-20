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
package electronbunch;

import org.la4j.Vector;
import org.la4j.vector.dense.BasicVector;

/**
 *
 * @author Ruslan
 * @version 1.01
 */
public abstract class AbstractElectronBunch implements Cloneable {

    /**
     * The charge of electron
     */
    public static final double E = 1.602e-19;

    /**
     * The mass of electron in MeV
     */
    public static final double mc2 = 0.5109989461;
    /**
     * Mean electron bunch gamma
     */
    protected double gamma = 50 / mc2;
    /**
     * Number of electrons in the bunch
     */
    protected double number = 0.2 / E * 1e-9;
    /**
     * Relative electron bunch energy spread
     */
    protected double delgamma = 1.25e-3;
    /**
     * Electron bunch semi-length, m
     */
    protected double length = 1.5e-3;
    /**
     * Electron transversal bunch emittance in x direction, m*rad
     */
    protected double epsx = 1e-6;
    /**
     * Electron transversal bunch emittance in y direction, m*rad
     */
    protected double epsy = 1e-6;
    /**
     * Electron bunch beta function in x direction, m
     */
    protected double betax = 0.02;
    /**
     * Electron bunch beta function in y direction, m
     */
    protected double betay = 0.02;
    /**
     * Electron bunch shift relative to the laser pulse, m
     */
    protected Vector shift;

    public AbstractElectronBunch() {
        shift = new BasicVector(new double[]{0.0, 0.0, 0.0});
    }

    @Override
    public Object clone() throws CloneNotSupportedException {
        Object tm = super.clone();
        ((AbstractElectronBunch) tm).setShift(new BasicVector(new double[]{this.getShift().get(0),
            this.getShift().get(1), this.getShift().get(2)}));
        return tm;
    }

    /**
     * sets electron's gamma parameter
     *
     * @param gamma electron's gamma
     */
    public void setGamma(double gamma) {
        this.gamma = gamma;
    }

    /**
     * Returning gamma
     *
     * @return gamma
     */
    public double getGamma() {
        return gamma;
    }

    /**
     * @return the number
     */
    public double getNumber() {
        return number;
    }

    /**
     * @param number the number to set
     */
    public void setNumber(double number) {
        this.number = number;
    }

    /**
     * @return the delgamma
     */
    public double getDelGamma() {
        return delgamma;
    }

    /**
     * @param delgamma the delgamma to set
     */
    public void setDelgamma(double delgamma) {
        this.delgamma = delgamma;
    }

    /**
     * @return the length
     */
    public double getLength() {
        return length;
    }

    /**
     * @param length the length to set
     */
    public void setLength(double length) {
        this.length = length;
    }

    /**
     * @return returning the emittance in the x direction
     */
    public double getEpsx() {
        return epsx;
    }

    /**
     * @param epsx setting the emittance in the x direction
     */
    public void setEpsx(double epsx) {
        this.epsx = epsx;
    }

    /**
     * @return returning the emittance in the y direction
     */
    public double getEpsy() {
        return epsy;
    }

    /**
     * @param epsy setting the emittance in the y direction
     */
    public void setEpsy(double epsy) {
        this.epsy = epsy;
    }

    /**
     * @return the betax
     */
    public double getBetax() {
        return betax;
    }

    /**
     * @param betax the betax to set
     */
    public void setBetax(double betax) {
        this.betax = betax;
    }

    /**
     * @return the betay
     */
    public double getBetay() {
        return betay;
    }

    /**
     * @param betay the betay to set
     */
    public void setBetay(double betay) {
        this.betay = betay;
    }

    /**
     * @return the shift
     */
    public Vector getShift() {
        return shift.copy();
    }

    /**
     * @param shift the shift to set
     */
    public void setShift(Vector shift) {
        this.shift = shift;
    }

    /**
     * Setting the width of the electron bunch in x direction
     *
     * @param w width
     */
    public void setxWidth(double w) {
        setBetax(w * w / getEpsx() * getGamma());
    }

    /**
     * Returning the velocity spread of the electron bunch in x direction
     *
     * @return velocity spread in the x direction
     */
    public double getXSpread() {
        return Math.sqrt(getEpsx() / getGamma() / getBetax());
    }

    /**
     * Setting the width of the electron bunch in y direction
     *
     * @param w width
     */
    public void setyWidth(double w) {
        setBetay(w * w / getEpsy() * getGamma());
    }

    /**
     * Returning the velocity spread of the electron bunch in y direction
     *
     * @return velocity spread in the y direction
     */
    public double getYSpread() {
        return Math.sqrt(getEpsy() / getGamma() / getBetay());
    }

    /**
     * Returning the average width of the electron bunch at position z
     *
     * @param z coordinate z
     * @return average width
     */
    public abstract double getWidth(double z);

    /**
     * Returning the average squared width of the electron bunch at position z
     *
     * @param z coordinate z
     * @return average width squared
     */
    public abstract double getWidth2(double z);

    /**
     * Returning the average velocity spread of the electron bunch
     *
     * @return average velocity spread
     */
    public double getSpread() {
        return Math.sqrt(getXSpread() * getYSpread());
    }

    /**
     * The angular distribution of electrons in the bunch
     *
     * @param thetax
     * @param thetay
     * @return
     */
    public abstract double angleDistribution(double thetax, double thetay);

    /**
     * The transversal spatial distribution of electrons in the bunch
     *
     * @param r
     * @return
     */
    public abstract double tSpatialDistribution(Vector r);

    /**
     * The longitudinal spatial distribution of electrons in the bunch
     *
     * @param r
     * @return
     */
    public abstract double lSpatialDistribution(Vector r);

    /**
     * Returning the width of the electron bunch in x direction at position z
     *
     * @param z coordinate z
     * @return width in the x direction
     */
    public abstract double getxWidth(double z);

    /**
     * Returning the width squared of the electron bunch in x direction at
     * position z
     *
     * @param z coordinate z
     * @return width in the x direction squared
     */
    public abstract double getxWidth2(double z);

    /**
     * Returning the width of the electron bunch in y direction at position z
     *
     * @param z coordinate z
     * @return width in the y direction
     */
    public abstract double getyWidth(double z);

    /**
     * Returning the width squared of the electron bunch in y direction at
     * position z
     *
     * @param z coordinate z
     * @return width in the x direction squared
     */
    public abstract double getyWidth2(double z);

    /**
     * Distribution of electrons by factor gamma
     *
     * @param g - gamma factor
     * @return
     */
    public abstract double gammaDistribution(double g);

}
