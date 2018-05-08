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
package electronbunch;

/**
 * The package for the laser-electron source simulation
 *
 */
import org.la4j.Vector;

/**
 * A class for properties of the Gaussian electron bunch. All properties are in the SI
 * system of units
 *
 * @author Ruslan Feshchenko
 * @version 3.0
 */
public class GaussianElectronBunch extends AbstractElectronBunch {

    /**
     * Constructor
     */
    public GaussianElectronBunch() {
        super();
    }

    /**
     * The angular distribution of electrons in the Gaussian bunch
     *
     * @param thetax
     * @param thetay
     * @return
     */
    @Override
    public double angleDistribution(double thetax, double thetay) {
        double dpx = getXSpread();
        double dpy = getYSpread();
        return Math.exp(-Math.pow(thetax / dpx, 2) - Math.pow(thetay / dpy, 2))
                / dpx / dpy / Math.PI;
    }

    /**
     * The transversal spatial distribution of electrons in the Gaussian bunch
     *
     * @return
     */
    @Override
    public double tSpatialDistribution(Vector r) {
        double K = Math.pow((r.get(0) - getShift().get(0)), 2) / getxWidth2(r.get(2) - getShift().get(2))
                + Math.pow((r.get(1) - getShift().get(1)), 2) / getyWidth2(r.get(2) - getShift().get(2));
        return Math.exp(-K) / getxWidth(r.get(2) - getShift().get(2))
                / getyWidth(r.get(2) - getShift().get(2)) / Math.PI;
    }

    /**
     * Returning the width of the Gaussian electron bunch in x direction at position z
     *
     * @param z coordinate z
     * @return width in the x direction
     */
    @Override
    public double getxWidth(double z) {
        return Math.sqrt((getBetax() + z * z / getBetax()) * getEpsx() / getGamma());
    }

    /**
     * Returning the width squared of the Gaussian electron bunch in x direction at position z
     *
     * @param z coordinate z
     * @return width in the x direction squared
     */
    @Override
    public double getxWidth2(double z) {
        return (getBetax() + z * z / getBetax()) * getEpsx() / getGamma();
    }

    /**
     * Returning the width of the Gaussian electron bunch in y direction at position z
     *
     * @param z coordinate z
     * @return width in the y direction
     */
    @Override
    public double getyWidth(double z) {
        return Math.sqrt((getBetay() + z * z / getBetay()) * getEpsy() / getGamma());
    }

    /**
     * Returning the width squared of the Gaussian electron bunch in y direction at position z
     *
     * @param z coordinate z
     * @return width in the x direction squared
     */
    @Override
    public double getyWidth2(double z) {
        return (getBetay() + z * z / getBetay()) * getEpsy() / getGamma();
    }

    /**
     * Returning the average width of the Gaussian electron bunch at position z
     *
     * @param z coordinate z
     * @return average width
     */
    @Override
    public double getWidth(double z) {
        return Math.sqrt(getxWidth(z) * getyWidth(z));
    }

    /**
     * Returning the average squared width of the Gaussian electron bunch at position z
     *
     * @param z coordinate z
     * @return average width squared
     */
    @Override
    public double getWidth2(double z) {
        return Math.sqrt(getxWidth2(z) * getyWidth2(z));
    }

    @Override
    public double lSpatialDistribution(Vector r) {
        return Math.exp(-Math.pow(r.get(2) / getLength(), 2)) / getLength() / Math.sqrt(Math.PI);
    }
}
