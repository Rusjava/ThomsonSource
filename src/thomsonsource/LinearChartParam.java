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
package thomsonsource;

import java.util.function.Function;
import org.la4j.vector.dense.BasicVector;

/**
 * Class for linear chart parameters
 *
 * @author Ruslan Feshchenko
 * @version 0.1
 */
public class LinearChartParam {

    /**
     * Default constructor
     */
    public LinearChartParam() {

    }

    /**
     * Size of the plot
     */
    protected int size;

    /**
     * Plot step
     */
    protected double step;

    /**
     * Plot offset
     */
    protected double offset;

    /**
     * Plot data
     */
    protected double[] data;

    /**
     * Maximum value
     */
    protected double umax;

    /**
     * Minimum value
     */
    protected double umin;

    /**
     * Function for calculation values
     */
    protected Function<Double, Double> func;

    /**
     * Returning plot size
     *
     * @return
     */
    public int getSize() {
        return size;
    }

    /**
     * returning plot step
     *
     * @return
     */
    public double getStep() {
        return step;
    }

    /**
     * Returning plot offset
     *
     * @return
     */
    public double getOffset() {
        return offset;
    }

    /**
     * returning plot data
     *
     * @return
     */
    public double[] getData() {
        return data;
    }

    /**
     * Returning maximum plot value
     *
     * @return
     */
    public double getUMax() {
        return umax;
    }
    
    /**
     * Returning minimum plot value
     *
     * @return
     */
    public double getUMin() {
        return umin;
    }

    /**
     * Setting up the data based on 2D array
     *
     * @param data 2D data array
     * @param index index of row/column
     * @param row row or column
     * @param size plot size
     * @param step plotting step
     * @param offset plotting offset
     * @throws java.lang.InterruptedException
     */
    public void setup(double[][] data, int index, boolean row,
            int size, double step, double offset) throws InterruptedException {
        if (Thread.currentThread().isInterrupted()) {
            throw new InterruptedException();
        }
        this.data = new double[size];
        for (int i = 0; i < size; i++) {
            this.data[i] = row ? data[i][index] : data[index][i];
        }
        setExtr();
        this.size = size;
        this.step = step;
        this.offset = offset;
    }

    /**
     * Setting up data based on a given function
     *
     * @param f function to get data from
     * @param size
     * @param step
     * @param offset
     * @throws java.lang.InterruptedException
     */
    public void setup(Function<Double, Double> f, int size,
            double step, double offset) throws InterruptedException {
        for (int i = 0; i < size; i++) {
            if (Thread.currentThread().isInterrupted()) {
                throw new InterruptedException();
            }
            this.data[i] = f.apply(step * i + offset);
        }
        setExtr();
        this.size = size;
        this.step = step;
        this.offset = offset;
        this.func = f;
    }

    /**
     * Calculating min and max values of data
     */
    protected void setExtr() {
        this.umax = (new BasicVector(this.data)).max();
        this.umin = (new BasicVector(this.data)).min();
    }
}
