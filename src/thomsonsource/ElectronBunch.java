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

/**
 * The package for the laser-electron source simulation
 *
 * @author Ruslan Feshchenko
 * @version 2.0
 */
import org.la4j.Vector;
import org.la4j.vector.dense.BasicVector;

/**
 * A class for the electron bunch properties. All properties are in the SI
 * system of units
 *
 * @author Ruslan Feshchenko
 * @version 1.3
 */
public class ElectronBunch implements Cloneable {

    /**
     * Constructor creating vector shift
     */
    public ElectronBunch() {
        shift = new BasicVector(new double[]{0.0, 0.0, 0.0});
    }

    @Override
    public Object clone() throws CloneNotSupportedException {
        Object tm = super.clone();
        ((ElectronBunch) tm).setShift(new BasicVector(new double[]{0.0, 0.0, 0.0}));
        ((ElectronBunch) tm).getShift().set(0, this.getShift().get(0));
        ((ElectronBunch) tm).getShift().set(1, this.getShift().get(1));
        ((ElectronBunch) tm).getShift().set(2, this.getShift().get(2));
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
    public double getDelgamma() {
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
        return shift;
    }

    /**
     * @param shift the shift to set
     */
    public void setShift(Vector shift) {
        this.shift = shift;
    }

    /**
     * Returning the width of the electron bunch in x direction
     *
     * @param z coordinate z
     * @return width in the x direction
     */
    public double getxWidth(double z) {
        return Math.sqrt((getBetax() + z * z / getBetax()) * getEpsx() / getGamma());
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
     * Returning the width squared of the electron bunch in x direction
     *
     * @param z coordinate z
     * @return width in the x direction squared
     */
    public double getxWidth2(double z) {
        return (getBetax() + z * z / getBetax()) * getEpsx() / getGamma();
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
     * Returning the width of the electron bunch in y direction
     *
     * @param z coordinate z
     * @return width in the y direction
     */
    public double getyWidth(double z) {
        return Math.sqrt((getBetay() + z * z / getBetay()) * getEpsy() / getGamma());
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
     * Returning the width squared of the electron bunch in y direction
     *
     * @param z coordinate z
     * @return width in the x direction squared
     */
    public double getyWidth2(double z) {
        return (getBetay() + z * z / getBetay()) * getEpsy() / getGamma();
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
     * Returning the average width of the electron bunch
     *
     * @param z coordinate z
     * @return average width
     */
    public double getWidth(double z) {
        return Math.sqrt(getxWidth(z) * getyWidth(z));
    }

    /**
     * Returning the average squared width of the electron bunch
     *
     * @param z coordinate z
     * @return average width squared
     */
    public double getWidth2(double z) {
        return Math.sqrt(getxWidth2(z) * getyWidth2(z));
    }

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
    public double angleDistribution(double thetax, double thetay) {
        double dpx = getXSpread();
        double dpy = getYSpread();
        return Math.exp(-Math.pow(thetax / dpx, 2) - Math.pow(thetay / dpy, 2))
                / dpx / dpy / Math.PI;
    }

    /**
     * Mean electron bunch gamma
     */
    private double gamma = 100;

    /**
     * Number of electrons in the bunch
     */
    private double number = 1 / E * 1e-9;

    /**
     * Relative electron bunch energy spread
     */
    private double delgamma = 1e-2;

    /**
     * Electron bunch semi-length, m
     */
    private double length = 0.0045;

    /**
     * Electron transversal bunch emittance in x direction, m*rad
     */
    private double epsx = 5e-6;
    
    /**
     * Electron transversal bunch emittance in y direction, m*rad
     */
    private double epsy = 5e-6;

    /**
     * Electron bunch beta function in x direction, m
     */
    private double betax = 0.01;

    /**
     * Electron bunch beta function in y direction, m
     */
    private double betay = 0.01;

    /**
     * Electron bunch shift relative to the laser pulse, m
     */
    private Vector shift;

    /**
     * The charge of electron
     */
    public static final double E = 1.602e-19;

}
