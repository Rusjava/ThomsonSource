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
import laserpulse.AbstractLaserPulse;
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
    private int ordernumber = 1;

    /**
     * The relative laser pulse intensity
     */
    private double zeta = 2;

    public NonLinearThomsonSource(AbstractLaserPulse l, AbstractElectronBunch b, int ordernumber) {
        super(l, b);
        this.ordernumber = ordernumber;
    }

    @Override
    public double directionFrequencyFluxNoSpread(Vector n, Vector v, double e) {
        double K, th;
        th = (1 - n.innerProduct(v)) * 2;
        K = Math.pow((Math.sqrt(e / lp.getPhotonEnergy() / (1 - e * th / lp.getPhotonEnergy() / 4)) - 2 * eb.getGamma()), 2)
                / 4 / Math.pow(eb.getGamma() * eb.getDelgamma(), 2);
        return getTotalFlux() * e * 3.0 / 64 / Math.PI / Math.sqrt(Math.PI) / eb.getDelgamma() / eb.getGamma() / lp.getPhotonEnergy()
                * Math.sqrt(e / lp.getPhotonEnergy()) * (Math.pow((1 - e * th / lp.getPhotonEnergy() / 2), 2) + 1)
                / Math.sqrt(1 - e * th / lp.getPhotonEnergy() / 4) * Math.exp(-K);
    }

    @Override
    public double[] directionFrequencyPolarizationNoSpread(Vector n, Vector v, double e) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public double directionFlux(Vector n, Vector v) {
        double gamma2, th;
        th = (1 - n.innerProduct(v)) * 2;
        gamma2 = eb.getGamma() * eb.getGamma();
        return getTotalFlux() * 3.0 / 2 / Math.PI * gamma2 * (1 + Math.pow(th * gamma2, 2))
                / Math.pow((1 + gamma2 * th), 4) * getGeometricFactor();
    }

    @Override
    public double directionEnergy(Vector n, Vector v) {
        double mv, M;
        mv = Math.sqrt(1.0 - 1.0 / eb.getGamma() / eb.getGamma());
        M = getZeta() * (1 + n.innerProduct(v)) / 4 / eb.getGamma() / eb.getGamma() / (1 + mv);
        return ordernumber * (1 + mv) * lp.getPhotonEnergy() / (1 - n.innerProduct(v) * mv + M);
    }

    /**
     * Setting the non-linear order number
     *
     * @return the n
     */
    public int getOrdernumber() {
        return ordernumber;
    }

    /**
     * Getting the non-linear order number
     *
     * @param n the n to set
     */
    public void setOrdernumber(int n) {
        this.ordernumber = n;
    }

    /**
     * @return the zeta
     */
    public double getZeta() {
        return zeta;
    }

    /**
     * Setting relative intensity based on the known laser pulser intensity
     */
    public void setZeta() {
        double k = 10e-2 * lp.getPhotonEnergy() / AbstractLaserPulse.HC;
        this.zeta = 16 * Math.PI * 1e3 * lp.getIntensity() / k / k / AS / AS;
    }

}
