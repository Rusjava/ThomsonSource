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
 * @version 0.7
 */
import org.la4j.vector.Vector;
import org.la4j.vector.dense.BasicVector;

/**
 * A class for the electron bunch properties. All properties are in the SI
 * system of units
 *
 * @author Ruslan Feshchenko
 * @version 1.1
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
        ((ElectronBunch) tm).shift = new BasicVector(new double[]{0.0, 0.0, 0.0});
        ((ElectronBunch) tm).shift.set(0, this.shift.get(0));
        ((ElectronBunch) tm).shift.set(1, this.shift.get(1));
        ((ElectronBunch) tm).shift.set(2, this.shift.get(2));
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
     * Returning the width of the electron bunch in x direction
     *
     * @param z coordinate z
     * @return width in the x direction
     */
    public double getxWidth(double z) {
        return Math.sqrt((betax + z * z / betax) * eps / gamma);
    }

    /**
     * Setting the width of the electron bunch in x direction
     *
     * @param w width
     */
    public void setxWidth(double w) {
        betax = w * w / eps * gamma;
    }

    /**
     * Returning the width squared of the electron bunch in x direction
     *
     * @param z coordinate z
     * @return width in the x direction squared
     */
    public double getxWidth2(double z) {
        return (betax + z * z / betax) * eps / gamma;
    }

    /**
     * Returning the velocity spread of the electron bunch in x direction
     *
     * @return velocity spread in the x direction
     */
    public double getXSpread() {
        return Math.sqrt(eps / gamma / betax);
    }

    /**
     * Returning the width of the electron bunch in y direction
     *
     * @param z coordinate z
     * @return width in the y direction
     */
    public double getyWidth(double z) {
        return Math.sqrt((betay + z * z / betay) * eps / gamma);
    }

    /**
     * Setting the width of the electron bunch in y direction
     *
     * @param w width
     */
    public void setyWidth(double w) {
        betay = w * w / eps * gamma;
    }

    /**
     * Returning the width squared of the electron bunch in y direction
     *
     * @param z coordinate z
     * @return width in the x direction squared
     */
    public double getyWidth2(double z) {
        return (betay + z * z / betay) * eps / gamma;
    }

    /**
     * Returning the velocity spread of the electron bunch in y direction
     *
     * @return velocity spread in the y direction
     */
    public double getYSpread() {
        return Math.sqrt(eps / gamma / betay);
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
    public double number = 1 / 1.6e-10;

    /**
     * Relative electron bunch energy spread
     */
    public double delgamma = 1e-2;

    /**
     * Electron bunch semi-length, m
     */
    public double length = 0.0045;

    /**
     * Electron transversal bunch emittance, m*rad
     */
    public double eps = 5e-6;

    /**
     * Electron bunch beta function in x direction, m
     */
    public double betax = 0.01;

    /**
     * Electron bunch beta function in y direction, m
     */
    public double betay = 0.01;

    /**
     * Electron bunch shift relative to the laser pulse, m
     */
    public Vector shift;
}
