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
package laserpulse;

import org.la4j.Vector;

/**
 * A class for laser pulse properties. All properties are in the SI system of
 * units
 *
 * @author Ruslan feshchenko
 * @version 3.0
 */
public class LaserPulse extends AbstractLaserPulse {

    /**
     * Constructor
     */
    public LaserPulse() {
        super();
    }

    /**
     * Returns the Gaussian laser bunch width at position z
     *
     * @param z coordinate
     * @return
     */
    @Override
    public double getWidth(double z) {
        return Math.sqrt((getRlength() + z * z / getRlength()) * rk / 2);
    }

    /**
     * Getting the Gaussian laser bunch width squared as function of z
     *
     * @param z coordinate z
     * @return
     */
    @Override
    public double getWidth2(double z) {
        return (getRlength() + z * z / getRlength()) * rk / 2;
    }

    /**
     * The transversal spatial distribution of photons in the Gaussian pulse
     *
     * @param r
     * @return
     */
    @Override
    public double tSpatialDistribution(Vector r) {
        double K = (Math.pow((r.get(0)), 2) + Math.pow((r.get(1)), 2)) / getWidth2(r.get(2));
        return Math.exp(-K) / getWidth2(r.get(2)) / Math.PI;
    }
}
