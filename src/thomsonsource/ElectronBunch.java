/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package thomsonsource;

/**
 * The package for the laser-electron source simulation
 * @author Ruslan Feshchenko
 * @version 0.3
 */

import org.la4j.vector.Vector;
import org.la4j.vector.dense.BasicVector;

/**
 * A class for the electron bunch properties. All properties are in the SI system of units
 * @author Ruslan Feshchenko
 */
public class ElectronBunch {
    
    /**
     * Constructor creating vector shift
     */
    public ElectronBunch () {
        shift=new BasicVector(new double []{0.0,0.0,0.0});
    }
    
    /**
     * Copies parameters of this object to tm object
     * @param tm object to clone to
     */
    public void duplicate (ElectronBunch tm) {
        tm.gamma=this.gamma;
        tm.number=this.number;
        tm.delgamma=this.delgamma;
        tm.length=this.length;
        tm.eps=this.eps;
        tm.betax=this.betax;
        tm.betay=this.betay;
        tm.shift.set(0,this.shift.get(0));
        tm.shift.set(1,this.shift.get(1));
        tm.shift.set(2,this.shift.get(2));    
    }
    
    /**
     * sets electron's gamma parameter
     * @param gamma electron's gamma
     */
    public void setGamma(double gamma) {
        this.gamma=gamma;
    }
    
    /**
     * Returning gamma
     * @return gamma
     */
    public double getGamma() {
        return gamma;
    }
    
    /**
     * Returning the width of the electron bunch in x direction
     * @param z coordinate z
     * @return width in the x direction
     */
    public double getxWidth(double z) {
        return Math.sqrt((betax+z*z/betax)*eps/gamma);
    }
    
    /**
     * Setting the width of the electron bunch in x direction
     * @param w width
     */
    public void setxWidth(double w) {
        betax=w*w/eps*gamma;
    }
    
    /**
     * Returning the width squared of the electron bunch in x direction
     * @param z coordinate z
     * @return width in the x direction squared
     */
    public double getxWidth2(double z) {
        return (betax+z*z/betax)*eps/gamma;
    }
    
    /**
     * Returning the velocity spread of the electron bunch in x direction
     * @return velocity spread in the x direction
     */
    public double getxSpread() {
        return Math.sqrt(eps/gamma/betax);
    }
    
    /**
     * Returning the width of the electron bunch in y direction
     * @param z coordinate z
     * @return width in the y direction
     */
    public double getyWidth(double z) {
        return Math.sqrt((betay+z*z/betay)*eps/gamma);
    }
    
    /**
     * Setting the width of the electron bunch in y direction
     * @param w width
     */
    public void setyWidth(double w) {
        betay=w*w/eps*gamma;
    }
    
    /**
     * Returning the width squared of the electron bunch in y direction
     * @param z coordinate z
     * @return width in the x direction squared
     */
    public double getyWidth2(double z) {
        return (betay+z*z/betay)*eps/gamma;
    }
    
    /**
     * Returning the velocity spread of the electron bunch in y direction
     * @return velocity spread in the y direction
     */
    public double getySpread() {
        return Math.sqrt(eps/gamma/betay);
    }
    
    /**
     * Returning the average width of the electron bunch
     * @param z coordinate z
     * @return average width
     */
    public double getWidth(double z) {
        return Math.sqrt(getxWidth(z)*getyWidth(z));
    }
    
    /**
     * Returning the average squared width of the electron bunch
     * @param z coordinate z
     * @return average width squared
     */
    public double getWidth2(double z) {
        return Math.sqrt(getxWidth2(z)*getyWidth2(z));
    }
    
    /**
     * Returning the average velocity spread of the electron bunch
     * @return average velocity spread
     */
    public double getSpread() {
        return Math.sqrt(getxSpread()*getySpread());
    }
    
    /**
     * Mean electron bunch gamma
     */
    private double gamma=100;

    /**
     * Number of electrons in the bunch
     */
    public double number=1/1.6e-10;

    /**
     * Relative electron bunch energy spread
     */
    public double delgamma=1e-2;

    /**
     * Electron bunch semi-length, m 
     */
    public double length=0.0045;

    /**
     * Electron transversal bunch emittance, m*rad
     */
    public double eps=5e-6;

    /**
     * Electron bunch beta function in x direction, m
     */
    public double betax=0.01;

    /**
     * Electron bunch beta function in y direction, m
     */
    public double betay=0.01; 

    /**
     * Electron bunch shift relative to the laser pulse, m
     */
    public Vector shift;  
}
