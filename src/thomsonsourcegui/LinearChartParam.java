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

import java.awt.BasicStroke;
import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Formatter;
import java.util.function.Function;
import org.la4j.vector.dense.BasicVector;
import java.util.List;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.DomainOrder;
import org.jfree.data.general.DatasetChangeListener;
import org.jfree.data.general.DatasetGroup;
import org.jfree.data.xy.XYDataset;

/**
 * Class for linear chart parameters
 *
 * @author Ruslan Feshchenko
 * @version 3.02
 */
public class LinearChartParam {

    /**
     * Default constructor
     *
     * @param trfunc - the function producing the plotted values
     */
    public LinearChartParam(List<Function<double[], Double>> trfunc) {
        this.trfunc = trfunc;
    }

    /**
     * Transformation functions
     */
    private final List<Function<double[], Double>> trfunc;

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
    protected double[][] data;

    /**
     * Maximum value
     */
    protected double umax;

    /**
     * Minimum value
     */
    protected double umin;

    /**
     * Functional list for calculation of values
     */
    protected List<Function<Double, Double>> func;

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
    public double[][] getData() {
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
     * Setting up the data based on a given 2D array
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
        this.size = size;
        this.step = step;
        this.offset = offset;
        this.data = new double[1][];
        this.data[0] = new double[size];
        for (int i = 0; i < size; i++) {
            if (Thread.currentThread().isInterrupted()) {
                throw new InterruptedException();
            }
            this.data[0][i] = row ? data[i][index] : data[index][i];
        }
        setExtr();
    }

    /**
     * Setting up the data based on a given list of functions
     *
     * @param f functional array to get data from
     * @param size
     * @param step
     * @param offset
     * @throws java.lang.InterruptedException
     */
    public void setup(List<Function<Double, Double>> f, int size,
            double step, double offset)throws InterruptedException {
        this.size = size;
        this.step = step;
        this.offset = offset;
        this.func = f;
        this.data = new double[f.size()][];
        for (int k = 0; k < f.size(); k++) {
            this.data[k] = new double[size];
        }
        for (int i = 0; i < size; i++) {
            double xp = step * i + offset;
            for (int k = 0; k < f.size(); k++) {
                if (Thread.currentThread().isInterrupted()) {
                    throw new InterruptedException();
                }
                this.data[k][i] = f.get(k).apply(xp);
            }
        }
        setExtr();
    }

    /**
     * Calculating the min and max values of the data
     */
    private void setExtr() {
        int sz = (func != null) ? func.size() : 1;
        double[] umaxt = new double[sz];
        double[] umint = new double[sz];
        double[] tdata = new double[getSize()];
        for (int k = 0; k < sz; k++) {
            for (int i = 0; i < getSize(); i++) {
                tdata[i] = getTransformedData(k, i);
            }
            umaxt[k] = (new BasicVector(tdata)).max();
            umint[k] = (new BasicVector(tdata)).min();
        }
        umax = (new BasicVector(umaxt)).max();
        umin = (new BasicVector(umint)).min();
    }

    /**
     * Getting the data if no function is specified and if specified then using it calculate them
     */
    private double getTransformedData(int k, int i) {
        //If trasformation functions are not defined then just return data
        if (trfunc == null) {
            return getData()[k][i];
        } else {
            //If exist, apply them
            double[] res = new double[data.length];
            for (int s = 0; s < data.length; s++) {
                res[s] = data[s][i];
            }
            return trfunc.get(k).apply(res);
        }
    }

    /**
     * Saves the calculated data into a given text file
     *
     * @param file file to save data to
     * @throws java.io.IOException
     */
    public void save(File file) throws IOException {
        Formatter fm;
        try (PrintWriter pw = new PrintWriter(new FileWriter(file, false))) {
            for (int i = 0; i < getSize(); i++) {
                fm = new Formatter();
                fm.format("%f", i * getStep() + getOffset());
                for (int s = 0; s < getData().length; s++) {
                    fm.format(" %10.8f", getTransformedData(s, i));
                }
                pw.println(fm);
            }
            pw.close();
        } catch (IOException e) {
            throw e;
        }
    }

    /**
     * Creating a line chart based on a dataset
     *
     * @param dataset the value of dataset
     * @param xlabel the value of xlabel
     * @param ylabel the value of ylabel
     * @return returning a JFreeChart
     */
    public static JFreeChart createLineChart(XYDataset dataset, String xlabel, String ylabel) {
        /* X axis */
        NumberAxis xAxis = new NumberAxis(xlabel);
        xAxis.setStandardTickUnits(NumberAxis.createStandardTickUnits());
        xAxis.setLowerMargin(0.0);
        xAxis.setUpperMargin(0.0);
        xAxis.setAutoRangeIncludesZero(false);
        /* Y axis */
        NumberAxis yAxis = new NumberAxis(ylabel);
        yAxis.setStandardTickUnits(NumberAxis.createStandardTickUnits());
        yAxis.setLowerMargin(0.0);
        yAxis.setUpperMargin(0.0);
        yAxis.setAutoRangeIncludesZero(false);
        /* Renderer */
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        for (int i = 0; i < dataset.getSeriesCount(); i++) {
            renderer.setSeriesLinesVisible(i, true);
            renderer.setSeriesShapesVisible(i, false);
            renderer.setSeriesStroke(i, new BasicStroke(2.0F));
        }
        renderer.setSeriesPaint(0, Color.BLUE);
        renderer.setSeriesPaint(1, Color.GREEN);
        renderer.setSeriesPaint(2, Color.MAGENTA);
        renderer.setSeriesPaint(3, Color.BLACK);
        /* Plot creation */
        XYPlot plot = new XYPlot(dataset, xAxis, yAxis, renderer);
        plot.setBackgroundPaint(Color.white);
        plot.setRangeGridlinePaint(Color.black);
        plot.setDomainGridlinePaint(Color.black);
        /* Chart creation */
        JFreeChart chart = new JFreeChart(plot);
        chart.setBackgroundPaint(Color.white);
        return chart;
    }

    /**
     * Creating a dataset based on the data specified in this LinearChartParam
     * type object
     *
     * @param keys the value of keys
     * @return
     */
    public XYDataset createLineDataset(String[] keys) {
        return new XYDataset() {
            @Override
            public int getSeriesCount() {
                return getData().length;
            }

            @Override
            public int getItemCount(int series) {
                return getSize();
            }

            @Override
            public Number getX(int series, int item) {
                return new Double(getXValue(series, item));
            }

            @Override
            public double getXValue(int series, int item) {
                return item * getStep() + getOffset();
            }

            @Override
            public Number getY(int series, int item) {
                return new Double(getYValue(series, item));
            }

            @Override
            public double getYValue(int series, int item) {
                return getTransformedData(series, item);
            }

            @Override
            public void addChangeListener(DatasetChangeListener listener) {
                // ignore - this dataset never changes
            }

            @Override
            public void removeChangeListener(DatasetChangeListener listener) {
                // ignore
            }

            @Override
            public DatasetGroup getGroup() {
                return null;
            }

            @Override
            public void setGroup(DatasetGroup group) {
                // ignore
            }

            @Override
            public Comparable getSeriesKey(int series) {
                return keys[series];
            }

            @Override
            public int indexOf(Comparable seriesKey) {
                return 0;
            }

            @Override
            public DomainOrder getDomainOrder() {
                return DomainOrder.ASCENDING;
            }
        };
    }
}
