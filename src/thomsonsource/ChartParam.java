/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package thomsonsource;

/**
 *
 * @author Ruslan Feshchenko
 * @version 0.7
 */
public abstract class ChartParam {
    /**
     * Constructor
     */
    public ChartParam () {
            
    }
    
    public ChartParam (double [][] u, double um) {
        this.udata=u;
        this.umax=um;
    }
    
    /**
     * Method generating matrices for Z data and determining the max value
     * 
     * @param xsize
     * @param ysize
     * @param xstep
     * @param ystep
     * @param xoffset
     * @param yoffset
     * @throws java.lang.InterruptedException
     */
    public void setup(int xsize, int ysize, double xstep, double ystep, 
            double xoffset, double yoffset)  throws InterruptedException {
        this.udata=new double[xsize][ysize];
        double x, y;
        for (int j=0; j<xsize; j++) {
            for (int p=0; p<ysize; p++) {
                if (Thread.currentThread().isInterrupted()) {
                    throw new InterruptedException();
                }
                x=xoffset+xstep*(j-xsize/2);
                y=yoffset+ystep*(p-ysize/2);
                this.udata[j][p]=func(x, y);
            }
        } 
        this.umax=this.udata[xsize/2][ysize/2];
        this.xoffset=xoffset;
        this.yoffset=yoffset;
        this.xstep=xstep;
        this.ystep=ystep;
        this.xsize=xsize;
        this.ysize=ysize;
    }
    
    /**
     * Returning the array with data
     * @return 
     */
    public double [][] getudata () {
        return this.udata;
    }
    
    /**
     * Returning the maximal value
     * @return 
     */   
    public double getumax () {
        return this.umax;
    }
    
    /**
     * Returning the x offset
     * @return 
     */
    public double getxoffset () {
        return this.xoffset;
    }
    
     /**
     * Returning the y offset
     * @return 
     */
    public double getyoffset () {
        return this.yoffset;
    }
    
    /**
     * Returning the x step
     * @return 
     */
    public double getxstep () {
        return this.xstep;
    }
    
    /**
     * Returning the y step
     * @return 
     */
    public double getystep () {
        return this.ystep;
    }
    
    /**
     * Returning the x size
     * @return 
     */
    public int getxsize () {
        return this.xsize;
    }
    
    /**
     * Returning the y size
     * @return 
     */
    public int getysize () {
        return this.ysize;
    }
    
    public abstract double func (double x, double y);
    
    private double [][] udata;
    private double umax;
    private double xoffset=0.0;
    private double yoffset=0.0;
    private double xstep;
    private double ystep;
    private int xsize;
    private int ysize;
}
