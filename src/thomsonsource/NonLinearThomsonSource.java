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
package thomsonsource;

import electronbunch.AbstractElectronBunch;
import laserpulse.LaserPulse;
import org.la4j.Vector;

/**
 * The main class containing all physics of LEXG in non-linear case
 * 
 * @version 1.0
 * @author Ruslan Feshchenko
 */
public class NonLinearThomsonSource extends AbstractThomsonSource {

    /**
     * Constructor
     * 
     * @param l
     * @param b
     * @param n - non-linear order number
     */
    /**
     * The non-linear order number
     */
    private int n;
    
    public NonLinearThomsonSource(LaserPulse l, AbstractElectronBunch b, int n) {
        super(l, b);
    }

    @Override
    public double directionFrequencyFluxNoSpread(Vector n, Vector v, double e) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public double[] directionFrequencyPolarizationNoSpread(Vector n, Vector v, double e) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public double directionFlux(Vector n, Vector v) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public double directionEnergy(Vector n, Vector v) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    /**
     * Setting the non-linear order number
     * @return the n
     */
    public int getN() {
        return n;
    }

    /**
     * Getting the non-linear order number
     * @param n the n to set
     */
    public void setN(int n) {
        this.n = n;
    }
    
}
