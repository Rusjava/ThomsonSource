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

import java.util.function.BiFunction;
import java.util.stream.Stream;

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
    protected double size;

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
     * Function for calculation values
     */
    protected BiFunction<Double, ThompsonSource, Double> func;

    /**
     * Returning plot size
     *
     * @return
     */
    public double getSize() {
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
     * Setting up the data based on 2D array
     * @param data 2D data array
     * @param index index of row/column
     * @param row row or column
     */
    public void setup(double[][] data, int index, boolean row) {
        int length = row ? data[0].length : data.length;
        for (int i = 0; i < length; i++) {
            this.data[i] = row ? data[index][i] : data[i][index];
        }
        umax=Stream.of(data).max(null).get()[0];
    }
}
