/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package thomsonsource;

import org.la4j.vector.Vector;
import org.la4j.vector.dense.BasicVector;

/**
 * A class for laser pulse properties. All properties are in the SI system of units
 * @author Ruslan feshchenko
 * @version 1.0
 */

public class LaserPulse implements Cloneable {
    public LaserPulse () {
        this.photonenergy=1.1*1.6e-19;
        this.setPulseEnergy(2.0e-2);
        this.rk=3.201e-26/this.photonenergy;
        this.direction=new BasicVector(new double []{0.0,0.0,1.0});
    }
    
    @Override
    public Object clone() throws CloneNotSupportedException {
        Object tm=super.clone();
        ((LaserPulse)tm).direction=new BasicVector(new double []{0.0,0.0,0.0});
        ((LaserPulse)tm).direction.set(0,this.direction.get(0));
        ((LaserPulse)tm).direction.set(1,this.direction.get(1));
        ((LaserPulse)tm).direction.set(2,this.direction.get(2)); 
        return tm;
    }

    /**
     * Returns the laser bunch width at position z
     * @param z coordinate
     * @return
     */
    public double getWidth(double z) {
        return Math.sqrt((rlength+z*z/rlength)*rk/2);
    }
    
    /**
     * Sets the laser bunch width
     * @param w width
     */
    public void setWidth(double w) {
        rlength=2*w*w/rk;
    }
    
    /**
     * Sets the laser bunch width squared
     * @param z coordinate z
     * @return
     */
    public double getWidth2(double z) {
        return (rlength+z*z/rlength)*rk/2;
    }
    
    /**
     * Sets the energy of laser photons
     * @param e photon energy
     */
    public void setPhotonEnergy (double e) {
        this.photonenergy=e;
        this.rk=3.201e-26/this.photonenergy;
    }
    
    /**
     * Return the energy of laser photons
     * @return
     */
    public double getPhotonEnergy () {
        return this.photonenergy;
    }
    
    /**
     * Sets the number of photons in laser pulse
     * @param n number of photons
     */
    public void setPhotonNumber (double n) {
        this.number=n;
    }
    
    /**
     * Returns the number of photons in laser pulse
     * @return
     */
    public double getPhotonNumber () {
        return this.number;
    }
    
    /**
     * Sets the energy of laser pulse (by setting the photon number)
     * @param e energy of laser pulse
     */
    public final void setPulseEnergy (double e) {
        this.number=e/this.photonenergy;
    }
    
    /**
     * Returns the energy of laser pulse
     * @return
     */
    public double getPulseEnergy () {
        return this.number*this.photonenergy;
    }
    
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
    public double length=4.5e-3;

    /**
     * Mean direction of the laser pulse
     */
    public Vector direction;

    /**
     * Laser pulse Rayleigh length, m 
     */
    public double rlength=5.4e-3; 

    /**
     * Pulse frequency, 1/s 
     */
    public double fq=7.9e7;

    /**
     * Laser pulse delay, m 
     */
    public double delay=0;
    
    private double rk;
}
