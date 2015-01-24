/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package thomsonsource;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.la4j.vector.Vector;
import org.apache.commons.math3.analysis.integration.*;
import org.la4j.vector.Vectors;
import org.la4j.vector.dense.BasicVector;

/**
 * A class for the X-ray source with associated electron bunch and laser pulse
 * @author Ruslan Feshchenko
 * @version 0.3
 */

public class ThompsonSource {
    ThompsonSource (LaserPulse l, ElectronBunch b) {
        this.lp=l;
        this.eb=b;
        calculateTotalFlux ();
        calculateGeometricFactor ();
    }
   
    /**
     * Number of rays exported for Shadow
     */
    public int ray_number=1000;

    /**
     *Number of points in Monte Carlo calculation of the geometric factor
     */
    public int np_gf=5000000;

    /**
     *Number of points in Monte Carlo calculation of the brilliance
     */
    public int ni_bril=30000;

    /**
     * Normalized total flux from the source
     */
    public double totalFlux;

    /**
     * Geometric factor. Assumes values from 0 to 1
     */
    public double gf=1;

    /**
     * Flag - whether or not the electron beam transversal velocity spread 
     * is taken into account
     */
    public boolean espread=false;
    
    private LaserPulse lp;
    private ElectronBunch eb;
    private final static double sigma=6.65e-29; /* Thompson cross-section, m2 */
    
    /**
     * A method calculating normalized total flux
     */
    
    public void calculateTotalFlux () {
       this.totalFlux=sigma*eb.number*lp.getPhotonNumber()*
                lp.fq/Math.PI/Math.sqrt((lp.getWidth2(0.0)+eb.getxWidth2(0.0))*
                        (lp.getWidth2(0.0)+eb.getyWidth2(0.0)));
    }
    
    /**
     * A method calculating the geometric factor
     */
    public void calculateGeometricFactor () {
        Vector iter=new BasicVector (new double []{0.0,0.0,0.0});
        double sum=0, wdx, wdy, len;
        int mult=2;
        wdx=mult*Math.max(eb.getxWidth(0.0)+Math.abs(eb.shift.get(0))/2, lp.getWidth(0.0)+Math.abs(eb.shift.get(0))/2);
        wdy=mult*Math.max(eb.getyWidth(0.0)+Math.abs(eb.shift.get(1))/2, lp.getWidth(0.0)+Math.abs(eb.shift.get(1))/2);
        len=mult*Math.max(eb.length+Math.abs(eb.shift.get(2))/2, lp.length+Math.abs(eb.shift.get(2))/2);
        
        for (int i=0; i<np_gf; i++) {
            iter.set(0, eb.shift.get(0)/2+wdx*(2*Math.random()-1.0));
            iter.set(1, eb.shift.get(1)/2+wdy*(2*Math.random()-1.0));
            iter.set(2, eb.shift.get(2)/2+len*(2*Math.random()-1.0));
            sum+=volumeFlux(iter);
        }
        this.gf=8*wdx*wdy*len*sum/np_gf;
    }
    
    /**
     * A method calculating the flux density in a given direction for a given 
     * X-ray photon energy
     * 
     * @param n
     * @param v
     * @param e
     * @return 
     */
    
    public double directionFrequencyFlux(Vector n, Vector v, double e) {
        if (espread) {
            return directionFrequencyFluxSpread(n, v, e);
        } else {
            return directionFrequencyFluxNoSpread(n, v, e);
        }
    }
    
    /**
     * A method calculating the flux density in a given direction for a given 
     * X-ray photon energy without taking into account electron transversal pulse spread
     * 
     * @param n
     * @param v
     * @param e
     * @return 
     */
    
    public double directionFrequencyFluxNoSpread(Vector n, Vector v, double e) {
        double u, K, th;
        th=(1-n.innerProduct(v))*2;
        K=Math.pow((Math.sqrt(e/lp.getPhotonEnergy()/(1-e*th/lp.getPhotonEnergy()/4))-2*eb.getGamma()),2)
                /4/Math.pow(eb.getGamma()*eb.delgamma,2);
        u=totalFlux*e*3.0/64/Math.PI/Math.sqrt(Math.PI)/eb.delgamma/eb.getGamma()/lp.getPhotonEnergy()*
                Math.sqrt(e/lp.getPhotonEnergy())*(Math.pow((1-e*th/lp.getPhotonEnergy()/2),2)+1)
                /Math.sqrt(1-e*th/lp.getPhotonEnergy()/4)*Math.exp(-K);
        return u;
    }
    
    /**
     * A method calculating the flux density in a given direction for a given 
     * X-ray photon energy taking into account electron transversal pulse spread
     * 
     * @param n
     * @param v
     * @param e
     * @return 
     */
    
    public double directionFrequencyFluxSpread(Vector n, Vector v, double e) {
        double u;
        RombergIntegrator spreadflux=new RombergIntegrator(); 
            UnivariateFunction func=
                        new UnivariateFrequencyFluxSpreadOuter (e, v, n);
            try {
                u=spreadflux.integrate(ni_bril, func, 0.0, 2*Math.PI);
                return u;
            } catch (TooManyEvaluationsException ex) {
                return 0; 
            }
    }    
        class UnivariateFrequencyFluxSpreadOuter implements UnivariateFunction {
            double e;
            Vector n, v0;
            public UnivariateFrequencyFluxSpreadOuter (double e, Vector v0, Vector n) {
                this.e=e;
                this.v0=v0;
                this.n=n;
            }
            @Override
            public double value(double phi) {
                double u;
                RombergIntegrator spreadflux=new RombergIntegrator(); 
                UnivariateFunction func=
                        new UnivariateFrequencyFluxSpreadInner (phi, e, v0, n);
                try {
                    u=spreadflux.integrate(ni_bril, func, 0.0, 3*eb.getSpread());
                    return u/Math.PI/eb.getxSpread()/eb.getySpread();
                } catch (TooManyEvaluationsException ex) {
                    return 0;
                }
            }
        }
        
        class UnivariateFrequencyFluxSpreadInner implements UnivariateFunction {
            double phi, e;
            Vector n, v0;
            public UnivariateFrequencyFluxSpreadInner (double phi, double e, Vector v0, Vector n) {
                this.phi=phi;
                this.e=e;
                this.n=n;
                this.v0=v0;
            }
            @Override
            public double value(double theta) {
                double u, dpx, dpy;
                dpx=eb.getxSpread();
                dpy=eb.getySpread();
                Vector dv;
                Vector v=new BasicVector (new double []{Math.sin(theta)*Math.cos(phi),
                    Math.sin(theta)*Math.sin(phi), Math.cos(theta)});
                dv=v.subtract(v0);
                u=theta*directionFrequencyFluxNoSpread(n, v, e)*
                        Math.exp(-Math.pow(dv.get(0)/dpx, 2)-Math.pow(dv.get(1)/dpy, 2));
                if ((new Double (u)).isNaN()) {
                    return 0;
                } 
                return u;
            }
        }
    
    /**
     * A method calculating the flux density in a given direction for a given 
     * X-ray photon energy for a given volume element
     * 
     * @param r
     * @param n
     * @param v
     * @param e
     * @return 
     */
    
    public double directionFrequencyVolumeFlux(Vector r, Vector n, Vector v, double e) {
        if (espread) {
            return directionFrequencyVolumeFluxSpread(r, n, v, e);
        } else {
            return directionFrequencyVolumeFluxNoSpread(r, n, v, e);
        }
    }
    
    /**
     * A method calculating the flux density in a given direction for a given 
     * X-ray photon energy for a given volume element without taking into account electron transversial pulse spread
     * 
     * @param r
     * @param n
     * @param v
     * @param e
     * @return 
     */
    
    public double directionFrequencyVolumeFluxNoSpread(Vector r, Vector n, Vector v, double e) {
        return directionFrequencyFluxNoSpread(n, v, e)*volumeFlux(r);
    }
    
    /**
     * A method calculating the flux density in a given direction for a given 
     * X-ray photon energy for a given volume element taking into account electron transversial pulse spread
     * 
     * @param r
     * @param n
     * @param v
     * @param e
     * @return 
     */
    
    public double directionFrequencyVolumeFluxSpread(Vector r, Vector n, Vector v, double e) {
        return directionFrequencyFluxSpread(n, v, e)*volumeFlux(r);
    }
    
    /**
     * An auxiliary method calculating volume density of the X-ray source
     * @param r
     * @return 
     */
    
    public double volumeFlux(Vector r) {
        double u, z0, z, z1, x0, x, x1, y0, y, y1, sn, cs, K, len;
        len=Math.sqrt(lp.length*lp.length+
                        eb.length*eb.length);
        sn=lp.direction.get(1);
        cs=lp.direction.get(2);
        x0=eb.shift.get(0);
        y0=eb.shift.get(1);
        z0=eb.shift.get(2);
        x=r.get(0);
        y=r.get(1);
        z=r.get(2);
        x1=x;
        y1=-sn*z+cs*y;
        z1=cs*z+sn*y;
        K=Math.pow((z+z1-z0-lp.delay)/len,2)+
                Math.pow((x-x0),2)/eb.getxWidth2(z-z0)+Math.pow((y-y0),2)/eb.getyWidth2(z-z0)+
                (Math.pow(x1,2)+Math.pow(y1,2))/lp.getWidth2(z1);
        if (K<15) {
            u=2.0/Math.pow(Math.PI, 1.5)*Math.sqrt((lp.getWidth2(0.0)+
                    eb.getxWidth2(0.0))*(lp.getWidth2(0.0)+
                    eb.getyWidth2(0.0)))/len/lp.getWidth2(z1)/eb.getxWidth(z-z0)/eb.getyWidth(z-z0)*Math.exp(-K);
        } else {
            u=0;
        }
        return u;
    }
    
    /**
     * A method giving the flux density in a given direction 
     * 
     * @param n
     * @param v
     * @return 
     */   
    
    public double directionFlux(Vector n, Vector v) {
        double u, gamma2, th;
        th=(1-n.innerProduct(v))*2;
        gamma2=eb.getGamma()*eb.getGamma();
        u=totalFlux*3.0/2/Math.PI*gamma2*(1+Math.pow(th*gamma2,2))/
                Math.pow((1+gamma2*th),4)*gf;
        return u;
    }
    
    /**
     * A method calculating X-ray energy in a given direction 
     * 
     * @param n
     * @param v
     * @return 
     */    
    
    public double directionEnergy(Vector n, Vector v) {
        double mv;
        mv=Math.sqrt(1.0-1.0/eb.getGamma()/eb.getGamma());
        return 2*lp.getPhotonEnergy()/(1-n.innerProduct(v)*mv);
    }
    
    /**
     * A method calculating spectral brilliance in a given direction
     * @param r0
     * @param n
     * @param v
     * @param e
     * @return 
     */
    
    public double directionFrequencyBrilliance(Vector r0, Vector n, Vector v, double e) {
        if (espread) {
            return directionFrequencyBrillianceSpread(r0, n, v, e);
        } else {
            return directionFrequencyBrillianceNoSpread(r0, n, v, e);
        }
    }
    
    /**
     * A method calculating spectral brilliance in a given direction
     * without taking into account electron transversal pulse spread
     * @param r0
     * @param n
     * @param v
     * @param e
     * @return 
     */
    
    public double directionFrequencyBrillianceNoSpread(Vector r0, Vector n, Vector v, double e) {
       double u;
       RombergIntegrator intvolumeflux=new RombergIntegrator(); 
       UnivariateVolumeFlux func=new UnivariateVolumeFlux (r0, n);
       try {
            u=intvolumeflux.integrate(30000, func,
               r0.fold(Vectors.mkEuclideanNormAccumulator())-3*eb.length,
               r0.fold(Vectors.mkEuclideanNormAccumulator())+3*eb.length);
            u=u*directionFrequencyFluxNoSpread(n, v, e);
            return u;
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }  
    }
    
    /**
     * A method calculating spectral brilliance in a given direction
     * taking into account electron transversal pulse spread
     * @param r0
     * @param n
     * @param v
     * @param e
     * @return 
     */
    
    public double directionFrequencyBrillianceSpread(Vector r0, Vector n, Vector v, double e) {
        double u;
        RombergIntegrator intvolumeflux=new RombergIntegrator(); 
        UnivariateVolumeFlux func=new UnivariateVolumeFlux (r0, n);
        try {
            u=intvolumeflux.integrate(ni_bril, func,
               r0.fold(Vectors.mkEuclideanNormAccumulator())-3*eb.length,
               r0.fold(Vectors.mkEuclideanNormAccumulator())+3*eb.length);
            u=u*directionFrequencyFluxSpread(n, v, e);
            return u;
        } catch (TooManyEvaluationsException ex) {
            return 0;
        }  
    }
    
    private class UnivariateVolumeFlux implements UnivariateFunction {
        Vector r0;
        Vector n0;
        public UnivariateVolumeFlux (Vector r0, Vector n0) {
            this.r0=r0;
            this.n0=n0;
        }
        @Override
        public double value(double x) {
            Vector r;
            r=r0.add(n0.multiply(x));
            double y = volumeFlux(r);
            if (n0.get(0)+n0.get(1)+n0.get(2)==0f) {
                throw new LocalException(x);
            }
            return y;
        }
    }
    
    private class LocalException extends RuntimeException {
     // The x value that caused the problem.
        private final double x;

        public LocalException(double x) {
            this.x = x;
        }

        public double getX() {
            return x;
        }
    }
    
    /**
     * Returning a random ray
     * @return an array with ray parameters
     */
    public double [] getRay () {
        double [] ray=new double [7];
        Vector n=new BasicVector(new double []{0.0,0.0,1.0});
        Vector r=new BasicVector(new double []{0.0,0.0,0.0});
        double prob0, prob, TrSpreadRange, ESpreadRange, EMax, EMin, mult=0.5;
        TrSpreadRange=Math.pow(eb.getSpread()*eb.getGamma(), 2);
        prob0=directionFrequencyVolumeFlux(r, n, new BasicVector(new double []{0.0,0.0,1.0}), 
                directionEnergy(n, new BasicVector(new double []{0.0,0.0,1.0})));
        do {
            ray[0]=2*(2*Math.random()-1.0)*Math.max(eb.getxWidth(0.0), lp.getWidth(0.0));
            ray[1]=2*(2*Math.random()-1.0)*Math.max(eb.getyWidth(0.0), lp.getWidth(0.0));
            ray[2]=2*(2*Math.random()-1.0)*Math.max(eb.length, lp.length);
            r.set(0,ray[0]);
            r.set(1,ray[1]);
            r.set(2,ray[2]);
            ray[3]=mult*(2*Math.random()-1.0)/eb.getGamma();
            ray[4]=mult*(2*Math.random()-1.0)/eb.getGamma();
            n.set(0,ray[3]);
            n.set(1,ray[4]);
            n.set(2,1.0);
            n=n.divide(n.fold(Vectors.mkEuclideanNormAccumulator()));
            ray[5]=n.get(2); 
            if (espread) {
                EMax=directionEnergy(n, n);
                ray[6]=(Math.random()*(TrSpreadRange)+1.0-TrSpreadRange)*
                    directionEnergy(n, new BasicVector(new double []{0.0,0.0,1.0}));
                prob=directionFrequencyVolumeFluxSpread(r, n, new BasicVector(new double []{0.0,0.0,1.0}), ray[6])/prob0;
            } else {
                ESpreadRange=2*eb.delgamma/(1+eb.getGamma()*eb.getGamma()*(ray[3]*ray[3]+ray[4]*ray[4]));
                ray[6]=(3*(2*Math.random()-1.0)*ESpreadRange+1.0)*directionEnergy(n, new BasicVector(new double []{0.0,0.0,1.0}));
                prob=directionFrequencyVolumeFluxNoSpread(r, n, new BasicVector(new double []{0.0,0.0,1.0}), ray[6])/prob0;
            }
        } while ( prob < Math.random() || (new Double(prob)).isNaN());    
        return ray;
    }
}
