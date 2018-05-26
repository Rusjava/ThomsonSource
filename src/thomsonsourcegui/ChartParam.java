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

import java.awt.Color;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.PaintScale;
import org.jfree.chart.renderer.xy.XYBlockRenderer;
import org.jfree.data.DomainOrder;
import org.jfree.data.general.DatasetChangeListener;
import org.jfree.data.general.DatasetGroup;
import org.jfree.data.xy.XYZDataset;
import org.la4j.vector.dense.BasicVector;

/**
 * Class for color chart parameters
 *
 * @author Ruslan Feshchenko
 * @version 1.1
 */
public abstract class ChartParam {

    /**
     * Constructor
     */
    public ChartParam() {
        this.sliderposition = (getxsize() - 1) * 50 / 100;

    }
    /**
     * Constructor which allows one to directly specify the data array
     * @param u data array
     * @param um maximum value of data
     */
    public ChartParam(double[][] u, double um) {
        this.sliderposition = (getxsize() - 1) * 50 / 100;
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
            }
        }
        setExtr();
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
    public final int getxsize() {
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

    /**
     * Calculating the max value of data
     */
    private void setExtr() {
        int sz = udata.length;
        double[] umaxt = new double[sz];
        for (int k = 0; k < sz; k++) {
            umaxt[k] = (new BasicVector(udata[k])).max();
        }
        umax = (new BasicVector(umaxt)).max();
    }

    public abstract double func(double x, double y);

    private double[][] udata; // Data
    private double umax;// Maximum of the data
    private double xoffset = 0.0;
    private double yoffset = 0.0;
    private double xstep;
    private double ystep;
    private int xsize;
    private int ysize;
    private int sliderposition;

    /**
     * Creates a color chart based on a given dataset
     *
     * @param dataset the value of dataset
     * @param xlabel the value of xlabel
     * @param ylabel the value of ylabel
     * @return
     */
    public JFreeChart createChart(XYZDataset dataset, String xlabel, String ylabel) {
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
        XYBlockRenderer renderer = new XYBlockRenderer();
        PaintScale scale = new JetPaintScale(0, this.getumax());
        renderer.setPaintScale(scale);
        renderer.setBlockHeight(this.getystep());
        renderer.setBlockWidth(this.getxstep());
        /* Plot creation */
        XYPlot plot = new XYPlot(dataset, xAxis, yAxis, renderer);
        plot.setBackgroundPaint(Color.white);
        plot.setDomainGridlinesVisible(false);
        plot.setRangeGridlinePaint(Color.black);
        /* Chart creation */
        JFreeChart chart = new JFreeChart(plot);
        chart.removeLegend();
        chart.setBackgroundPaint(Color.white);
        return chart;
    }

    /**
     * Creates a color bar for this chart
     *
     * @param label the value of label
     * @return
     */
    public JFreeChart createColorBar(String label) {
        NumberAxis xAxis = new NumberAxis();
        xAxis.setLowerMargin(0.0);
        xAxis.setUpperMargin(0.0);
        xAxis.setVisible(false);
        NumberAxis yAxis = new NumberAxis(label);
        yAxis.setStandardTickUnits(NumberAxis.createStandardTickUnits());
        yAxis.setLowerMargin(0.0);
        yAxis.setUpperMargin(0.0);
        XYZDataset dataset = new XYZDataset() {
            @Override
            public int getSeriesCount() {
                return 1;
            }

            @Override
            public int getItemCount(int series) {
                return getysize();
            }

            @Override
            public Number getX(int series, int item) {
                return new Double(getXValue(series, item));
            }

            @Override
            public double getXValue(int series, int item) {
                return 0;
            }

            @Override
            public Number getY(int series, int item) {
                return new Double(getYValue(series, item));
            }

            @Override
            public double getYValue(int series, int item) {
                return item * getumax() / (getysize() - 1);
            }

            @Override
            public Number getZ(int series, int item) {
                return new Double(getZValue(series, item));
            }

            @Override
            public double getZValue(int series, int item) {
                return getYValue(series, item);
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
                return "colorbar";
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
        XYBlockRenderer renderer = new XYBlockRenderer();
        PaintScale scale = new JetPaintScale(0, this.getumax());
        renderer.setPaintScale(scale);
        renderer.setBlockHeight((double) this.getumax() / (this.getysize() - 1));
        XYPlot plot = new XYPlot(dataset, xAxis, yAxis, renderer);
        plot.setBackgroundPaint(Color.white);
        plot.setDomainGridlinesVisible(false);
        plot.setRangeGridlinePaint(Color.black);
        JFreeChart chart = new JFreeChart("", plot);
        chart.removeLegend();
        chart.setBackgroundPaint(Color.white);
        return chart;
    }

    /**
     * Returning the current slider position in scale from 0 to 1
     * 
     * @return the position of the slider
     */
    public int getSliderposition() {
        return sliderposition;
    }

    /**
     * Setting the current slider position
     * 
     * @param sliderposition the position of the slider to set in scale from 0 to 100
     */
    public void setSliderposition(int sliderposition) {
        this.sliderposition = (getxsize() - 1) * sliderposition / 100;
    }

    /**
     * Creates a dataset for a color chart
     * 
     * @param linemark flag indicating whether the position of the slider should be shown
     * @return 
     */
    public XYZDataset createDataset(final boolean linemark) {
        return new XYZDataset() {
            @Override
            public int getSeriesCount() {
                return 1;
            }

            @Override
            public int getItemCount(int series) {
                return getxsize() * getysize();
            }

            @Override
            public Number getX(int series, int item) {
                return new Double(getXValue(series, item));
            }

            @Override
            public double getXValue(int series, int item) {
                return (getXindex(series, item) - getxsize() / 2) * getxstep() + getxoffset();
            }

            public int getXindex(int series, int item) {
                return item / getysize();
            }

            @Override
            public Number getY(int series, int item) {
                return new Double(getYValue(series, item));
            }

            @Override
            public double getYValue(int series, int item) {
                return (getYindex(series, item) - getysize() / 2) * getystep() + getyoffset();
            }

            public int getYindex(int series, int item) {
                return item - (item / getysize()) * getysize();
            }

            @Override
            public Number getZ(int series, int item) {
                return new Double(getZValue(series, item));
            }

            @Override
            public double getZValue(int series, int item) {
                int x = getXindex(series, item);
                int y = getYindex(series, item);
                if (!linemark) {
                    return getudata()[x][y];
                } else {
                    if (x == sliderposition) {
                        return getumax() / 2;
                    } else {
                        return getudata()[x][y];
                    }
                }
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
                return "Flux";
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
