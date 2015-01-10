/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package thomsonsource;

import org.la4j.vector.Vector;
import org.la4j.vector.dense.BasicVector;

/**
 * A class for laser pulse properties
 * @author Ruslan feshchenko
 */

public class LaserPulse {
    public LaserPulse () {
        this.photonenergy=1.1*1.6e-19;
        this.setPulseEnergy(2.0e-2);
        this.rk=3.201e-26/this.photonenergy;
        this.direction=new BasicVector(new double []{0.0,0.0,1.0});
    }
    
    public void duplicate (LaserPulse tm) {
        tm.setPhotonEnergy(this.photonenergy);
        tm.setPhotonNumber(this.number);
        tm.length=this.length;
        tm.rlength=this.rlength;
        tm.fq=this.fq;
        tm.delay=this.delay;
        tm.direction.set(0,this.direction.get(0));
        tm.direction.set(1,this.direction.get(1));
        tm.direction.set(2,this.direction.get(2));
    }
    public double getWidth(double z) {
        return Math.sqrt((rlength+z*z/rlength)*rk);
    }
    
    public void setWidth(double w) {
        rlength=w*w/rk;
    }
    
    public double getWidth2(double z) {
        return (rlength+z*z/rlength)*rk;
    }
    
    public void setPhotonEnergy (double e) {
        this.photonenergy=e;
        this.rk=3.201e-26/this.photonenergy;
    }
    
    public double getPhotonEnergy () {
        return this.photonenergy;
    }
    
    public void setPhotonNumber (double n) {
        this.number=n;
    }
    
    public double getPhotonNumber () {
        return this.number;
    }
    
    public void setPulseEnergy (double e) {
        this.number=e/this.photonenergy;
    }
    
    public double getPulseEnergy () {
        return this.number*this.photonenergy;
    }
    
    private double photonenergy; /* Photon energy, J */
    private double number; /* Number of photon in the laser pulse */
    public double length=4.5e-3; /* Laser pulse semi-length, m */  
    public Vector direction; /* Mean direction of the laser pulse */
    public double rlength=2.7e-3; /* Laser pulse Reyley length, m */
    public double fq=7.9e7; /* Pusle frequency, 1/s */
    public double delay=0; /* Laser pulse delay, m */
    private double rk;
}
