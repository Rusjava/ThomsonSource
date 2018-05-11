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
package thomsonsourcegui;

/**
 * Class for color chart parameters
 *
 * @author Ruslan Feshchenko
 * @version 1.0
 */
public abstract class ChartParam {

    /**
     * Constructor
     */
    public ChartParam() {

    }

    public ChartParam(double[][] u, double um) {
        this.udata = u;
        this.umax = um;
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
            double xoffset, double yoffset) throws InterruptedException {
        this.xoffset = xoffset;
        this.yoffset = yoffset;
        this.xstep = xstep;
        this.ystep = ystep;
        this.xsize = xsize;
        this.ysize = ysize;
        this.udata = new double[xsize][ysize];
        double x, y;
        for (int j = 0; j < xsize; j++) {
            for (int p = 0; p < ysize; p++) {
                if (Thread.currentThread().isInterrupted()) {
                    throw new InterruptedException();
                }
                x = xoffset + xstep * (j - xsize / 2);
                y = yoffset + ystep * (p - ysize / 2);
                this.udata[j][p] = func(x, y);
                System.out.println(j + " " + p + "\n");
            }
        }
        this.umax = this.udata[xsize / 2][ysize / 2];
    }

    /**
     * Returning the array with data
     *
     * @return
     */
    public double[][] getudata() {
        return this.udata;
    }

    /**
     * Returning the maximal value
     *
     * @return
     */
    public double getumax() {
        return this.umax;
    }

    /**
     * Returning the x offset
     *
     * @return
     */
    public double getxoffset() {
        return this.xoffset;
    }

    /**
     * Returning the y offset
     *
     * @return
     */
    public double getyoffset() {
        return this.yoffset;
    }

    /**
     * Returning the x step
     *
     * @return
     */
    public double getxstep() {
        return this.xstep;
    }

    /**
     * Returning the y step
     *
     * @return
     */
    public double getystep() {
        return this.ystep;
    }

    /**
     * Returning the x size
     *
     * @return
     */
    public int getxsize() {
        return this.xsize;
    }

    /**
     * Returning the y size
     *
     * @return
     */
    public int getysize() {
        return this.ysize;
    }

    public abstract double func(double x, double y);

    private double[][] udata;
    private double umax;
    private double xoffset = 0.0;
    private double yoffset = 0.0;
    private double xstep;
    private double ystep;
    private int xsize;
    private int ysize;
}
