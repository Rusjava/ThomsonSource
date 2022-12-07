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

import laserpulse.GaussianLaserPulse;
import electronbunch.GaussianElectronBunch;
import java.text.*;
import java.awt.*;
import java.awt.event.ItemEvent;
import javax.swing.*;
import javax.swing.border.TitledBorder;

import java.io.File;
import java.io.IOException;
import java.util.Locale;
import java.util.concurrent.ExecutionException;
import java.util.IllegalFormatException;
import java.util.concurrent.CancellationException;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.HashMap;
import java.util.Map;
import java.util.List;
import java.util.ArrayList;
import java.util.Date;
import java.util.Enumeration;
import java.util.jar.Manifest;
import java.util.Properties;
import java.util.function.Function;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.renderer.PaintScale;
import org.jfree.chart.renderer.xy.XYBlockRenderer;

import org.la4j.*;
import org.la4j.vector.dense.*;

import static TextUtilities.MyTextUtilities.*;
import electronbunch.AbstractElectronBunch;
import java.net.URL;
import javax.swing.filechooser.FileNameExtensionFilter;
import laserpulse.AbstractLaserPulse;
import shadowfileconverter.ShadowFiles;
import thomsonsource.AbstractThomsonSource;
import thomsonsource.LinearThomsonSource;
import thomsonsource.NonLinearThomsonSource;

/**
 * The GUI for non-linear Thomson source program
 *
 * @author Ruslan Feshchenko
 * @version 3.42
 */
public class ThomsonJFrame extends javax.swing.JFrame {

    /**
     * Creates new form ThomsonJFrame
     */
    public ThomsonJFrame() {
        this.xrayenergyborder = javax.swing.BorderFactory.createTitledBorder(null, "X-ray photon energy",
                javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION);
        this.ebunch = new GaussianElectronBunch();
        this.lpulse = new GaussianLaserPulse();
        this.tsource = new NonLinearThomsonSource(lpulse, ebunch, 1);
        this.tsourcelinear = new LinearThomsonSource(lpulse, ebunch);
        this.xsize = 300;
        this.ysize = 200;
        this.xstep = 20.0 / xsize;
        this.ystep = 20.0 / ysize;
        this.estep = 2000 / xsize;
        this.oldStrings = new HashMap<>();
        this.rayNumberBox = getIntegerFormattedTextField(1000, 1, 1000000);
        this.orderNumberBox = getIntegerFormattedTextField(1, 1, 10);
        this.rayXAngleRangeBox = getDoubleFormattedTextField(0.1, 0.0, 100.0, false);
        this.rayYAngleRangeBox = getDoubleFormattedTextField(0.1, 0.0, 100.0, false);
        this.gfMonteCarloNumberBox = getIntegerFormattedTextField(50000, 1, 100000000);
        this.numericallPrecisionBox = getDoubleFormattedTextField(1e-3, 1e-10, 1e-1, true);
        this.shiftFactorBox = getDoubleFormattedTextField(1.0, 1e-20, 1e20, true);
        this.xSizeBox = getIntegerFormattedTextField(300, 1, 10000);
        this.ySizeBox = getIntegerFormattedTextField(200, 1, 10000);
        this.xRangeBox = getDoubleFormattedTextField(20.0, 0.0, 100.0, false);
        this.yRangeBox = getDoubleFormattedTextField(20.0, 0.0, 100.0, false);
        this.xEnergyRangeBox = getDoubleFormattedTextField(2000.0, 0.0, 20000.0, false);
        this.rayMinEnergyBox = getDoubleFormattedTextField(42.0, 0.0, 10000.0, false);
        this.rayEnergyRangeBox = getDoubleFormattedTextField(6.0, 0.0, 1000.0, false);
        this.threadsNumberBox = getIntegerFormattedTextField(2, 1, 100);
        this.ksi1Box = getDoubleFormattedTextField(1.0, -1.0, 1.0, false);
        this.ksi2Box = getDoubleFormattedTextField(0.0, -1.0, 1.0, false);
        this.ksi3Box = getDoubleFormattedTextField(0.0, -1.0, 1.0, false);
        this.parametersIniFileName = "my.ini";
        this.orderofmagnitude = 10;
        this.normfactor = Math.pow(10, -15 - orderofmagnitude);

        /**
         * Defining polarization transformation functions
         */
        fn = new ArrayList<>();
        fn.add(x -> {
            return (x[0] == 0) ? 0 : x[1] / x[0];
        });
        fn.add(x -> {
            return (x[0] == 0) ? 0 : x[2] / x[0];
        });
        fn.add(x -> {
            return (x[0] == 0) ? 0 : x[3] / x[0];
        });
        fn.add(x -> {
            return (x[0] == 0) ? 0 : Math.sqrt(x[1] * x[1] + x[2] * x[2] + x[3] * x[3]) / x[0];
        });
        /**
         * An auxiliary method giving the flux density in a given direction
         *
         */
        this.fluxdata = new ChartParam() {
            @Override
            public double func(double thetax, double thetay) {
                Vector v, n;
                v = new BasicVector(new double[]{0.0, 0.0, 1.0});
                n = new BasicVector(new double[]{thetax * 1e-3, thetay * 1e-3, 1.0});
                n = n.divide(n.fold(Vectors.mkEuclideanNormAccumulator()));
                return 1e-6 * tsource.directionFlux(n, v) / 1e10;
            }
        };
        fluxdata.setSliderposition(50);
        /**
         * An auxiliary method calculating the flux density in a given direction
         * for a given X-ray photon energy
         *
         */
        this.fluxcrossdata = new ChartParam() {
            @Override
            public double func(double e, double theta) {
                Vector n, v;
                v = new BasicVector(new double[]{0.0, 0.0, 1.0});
                n = new BasicVector(new double[]{hoffset * 1e-3, theta * 1e-3, 1.0});
                n = n.divide(n.fold(Vectors.mkEuclideanNormAccumulator()));
                try {
                    return 1e-9 * tsource.getGeometricFactor() * tsource.directionFrequencyFlux(n, v, null, e * GaussianElectronBunch.E) / 1e10;
                } catch (InterruptedException ex) {
                    return 0;
                }
            }
        };

        /**
         * An auxiliary method calculating X-ray energy in a given direction
         *
         */
        this.xenergydata = new ChartParam() {
            @Override
            public double func(double thetax, double thetay) {
                Vector n, v;
                v = new BasicVector(new double[]{0.0, 0.0, 1.0});
                n = new BasicVector(new double[]{thetax * 1e-3, thetay * 1e-3, 1.0});
                n = n.divide(n.fold(Vectors.mkEuclideanNormAccumulator()));
                return tsource.directionEnergy(n, v) / GaussianElectronBunch.E * 1e-3;
            }
        };
        xenergydata.setSliderposition(50);
        /**
         * Auxiliary object for linear energy chart parameters
         */
        this.xenergycrossdata = new LinearChartParam(null);

        /**
         * Objects for the brilliance calculation
         */
        this.brilForm = new CalcBoxParam(new String[]{"Spectral brilliance"}, null);
        this.brilForm.valueUnitLabels = new String[]{"mrad", "ps", "mm", "mm", "mm mrad", "mm mrad", "mm mrad",
            "mm", "<html>&mu;m</html>", "", "keV", "mrad"};
        this.brilForm.plotLabels = new String[]{"Laser-electron angle, mrad", "Delay, ps", "Z-shift, mm", "beta, mm",
            "eps, mm mrad", "X-eps, mm mrad", "Y-eps, mm mrad", "Reyleigh length, mm", "Waist semi-width, \u03BCm", "\u0394\u03B3/\u03B3",
            "X-ray energy, keV", "Observation angle, mrad"};
        this.brilForm.conversionValues = new double[]{1e-3, 3e-4, 1e-3, 1e-3, 1e-6, 1e-6, 1e-6, 1e-3, 1e-6, 1.0, GaussianElectronBunch.E * 1e3, 1e-3};
        this.brilForm.minValues = new String[]{"0", "0", "0", "10", "0.5", "0.5", "0.5", "0.3", "5", "0.1", "0", "0"};
        this.brilForm.maxValues = new String[]{"50", "100", "10", "50", "5", "5", "5", "3", "50", "1", "100", "5"};
        this.brilForm.savetext = "Choose file to save spectral brilliance data";
        this.brilForm.numberOfItems = 12;

        /**
         * Objects for the non-linear brilliance calculation
         */
        this.brilFormNonLinear = new CalcBoxParam(new String[]{"Spectral brilliance"}, null);
        this.brilFormNonLinear.valueUnitLabels = new String[]{"mrad", "ps", "mm", "mm", "mm mrad", "mm mrad", "mm mrad",
            "mm", "<html>&mu;m</html>", "", "keV", "mrad"};
        this.brilFormNonLinear.plotLabels = new String[]{"Laser-electron angle, mrad", "Delay, ps", "Z-shift, mm", "beta, mm",
            "eps, mm mrad", "X-eps, mm mrad", "Y-eps, mm mrad", "Reyleigh length, mm", "Waist semi-width, \u03BCm", "\u0394\u03B3/\u03B3",
            "X-ray energy, keV", "Observation angle, mrad"};
        this.brilFormNonLinear.conversionValues = new double[]{1e-3, 3e-4, 1e-3, 1e-3, 1e-6, 1e-6, 1e-6, 1e-3, 1e-6, 1.0, GaussianElectronBunch.E * 1e3, 1e-3};
        this.brilFormNonLinear.minValues = new String[]{"0", "0", "0", "10", "0.5", "0.5", "0.5", "0.3", "5", "0.1", "0", "0"};
        this.brilFormNonLinear.maxValues = new String[]{"50", "100", "10", "50", "5", "5", "5", "3", "50", "1", "100", "5"};
        this.brilFormNonLinear.savetext = "Choose file to save spectral brilliance data";
        this.brilFormNonLinear.numberOfItems = 12;

        /**
         * Objects for the exact and approximate GF calculations
         */
        this.gfForm = new CalcBoxParam(new String[]{"Full flux", "Approximate full flux"}, null);
        this.gfForm.valueUnitLabels = new String[]{"mrad", "ps", "mm", "mm", "mm mrad", "mm", "<html>&mu;m</html>"};
        this.gfForm.plotLabels = new String[]{"Angle, mrad", "Delay, ps", "Z-shift, mm", "beta, mm",
            "eps, mm mrad", "Reyleigh length, mm", "Waist semi-width, \u03BCm"};
        this.gfForm.conversionValues = new double[]{1e-3, 3e-4, 1e-3, 1e-3, 1e-6, 1e-3, 1e-6};
        this.gfForm.minValues = new String[]{"0", "0", "0", "10", "0.5", "0.3", "5"};
        this.gfForm.maxValues = new String[]{"50", "100", "10", "50", "5", "3", "50"};
        this.gfForm.savetext = "Choose file to save geometric factor data";
        this.gfForm.numberOfItems = 7;
        this.gfForm.selectedItemIndex = 0;
        this.gfForm.selectedItemIndexClone = 0;

        /**
         * Objects for the linear polarization calculation
         */
        this.polForm = new CalcBoxParam(new String[]{"\u03BE1", "\u03BE2", "\u03BE3", "polarization degree"}, fn);
        this.polForm.valueUnitLabels = new String[]{"mrad", "ps", "mm", "mm", "mm mrad", "mm mrad", "mm mrad",
            "mm", "<html>&mu;m</html>", "", "keV", "mrad"};
        this.polForm.plotLabels = new String[]{"Laser-electron angle, mrad", "Delay, ps", "Z-shift, mm", "beta, mm",
            "eps, mm mrad", "X-eps, mm mrad", "Y-eps, mm mrad", "Reyleigh length, mm", "Waist semi-width, \u03BCm", "\u0394\u03B3/\u03B3",
            "X-ray energy, keV", "Observation angle, mrad"};
        this.polForm.conversionValues = new double[]{1e-3, 3e-4, 1e-3, 1e-3, 1e-6, 1e-6, 1e-6, 1e-3, 1e-6, 1.0, GaussianElectronBunch.E * 1e3, 1e-3};
        this.polForm.minValues = new String[]{"0", "0", "0", "10", "0.5", "0.5", "0.5", "0.3", "5", "0.1", "0", "0"};
        this.polForm.maxValues = new String[]{"50", "100", "10", "50", "10", "10", "10", "10", "100", "0.01", "100", "5"};
        this.polForm.savetext = "Choose file to save polarization data";
        this.polForm.numberOfItems = 12;
        /**
         * Objects for the non-linear polarization calculation
         */
        this.polFormNonLinear = new CalcBoxParam(new String[]{"\u03BE1", "\u03BE2", "\u03BE3", "polarization degree"}, fn);
        this.polFormNonLinear.valueUnitLabels = new String[]{"mrad", "ps", "mm", "mm", "mm mrad", "mm mrad", "mm mrad",
            "mm", "<html>&mu;m</html>", "", "keV", "mrad"};
        this.polFormNonLinear.plotLabels = new String[]{"Laser-electron angle, mrad", "Delay, ps", "Z-shift, mm", "beta, mm",
            "eps, mm mrad", "X-eps, mm mrad", "Y-eps, mm mrad", "Reyleigh length, mm", "Waist semi-width, \u03BCm", "\u0394\u03B3/\u03B3",
            "X-ray energy, keV", "Observation angle, mrad"};
        this.polFormNonLinear.conversionValues = new double[]{1e-3, 3e-4, 1e-3, 1e-3, 1e-6, 1e-6, 1e-6, 1e-3, 1e-6, 1.0, GaussianElectronBunch.E * 1e3, 1e-3};
        this.polFormNonLinear.minValues = new String[]{"0", "0", "0", "10", "0.5", "0.5", "0.5", "0.3", "5", "0.1", "0", "0"};
        this.polFormNonLinear.maxValues = new String[]{"50", "100", "10", "50", "5", "5", "5", "3", "50", "1", "100", "5"};
        this.polFormNonLinear.savetext = "Choose file to save polarization data";
        this.polFormNonLinear.numberOfItems = 12;
        //Setting default number of threads
        this.threadsNumberBox.setValue(new Integer(Runtime.getRuntime().availableProcessors()));

        initComponents();
        // Adding a listerner to the UI manager for skin update
        UIManager.addPropertyChangeListener(e -> {
            SwingUtilities.updateComponentTreeUI(this);
            SwingUtilities.updateComponentTreeUI(gfCalc);
            SwingUtilities.updateComponentTreeUI(brillianceCalc);
            SwingUtilities.updateComponentTreeUI(brillianceCalcNonLinear);
            SwingUtilities.updateComponentTreeUI(polarizationCalc);
            SwingUtilities.updateComponentTreeUI(rayProgressFrame);
        });

        //Get the current path for the java program and trying to open the "my.ini" file  
        java.awt.EventQueue.invokeLater(() -> {
            try {
                pFile = new File(new File(".").getCanonicalPath() + File.separator + parametersIniFileName);
            } catch (IOException ex) {
                //Do nothing;
            }
            if (pFile.exists() && pFile.isFile()) {
                loadParameters(pFile);
            }
            pFile = new File(".");
        });
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        brillianceCalc = new javax.swing.JFrame();
        BrillianceParam = new javax.swing.JPanel();
        BrillianceCalcBox = new javax.swing.JComboBox();
        BrillianceCalcStart = new javax.swing.JButton();
        BrillianceCalcSave = new javax.swing.JButton();
        Brilminvalue = new javax.swing.JTextField();
        Brilminvaluelabel = new javax.swing.JLabel();
        Brilminvalueunitlabel = new javax.swing.JLabel();
        Brilmaxvalueunitlabel = new javax.swing.JLabel();
        Brilmaxvalue = new javax.swing.JTextField();
        Brilmaxvaluelabel = new javax.swing.JLabel();
        jCheckBoxSpread = new javax.swing.JCheckBox();
        BrilProgressBar = new javax.swing.JProgressBar();
        jAngleLabel = new javax.swing.JLabel();
        angleValue = new javax.swing.JTextField();
        angleValueUnitLable = new javax.swing.JLabel();
        jEnergyLabel = new javax.swing.JLabel();
        energyValue = new javax.swing.JTextField();
        energyValueUnitLable = new javax.swing.JLabel();
        BrillianceCalcGraph = new javax.swing.JPanel();
        brillianceCalcNonLinear = new javax.swing.JFrame();
        BrillianceParamNonLinear = new javax.swing.JPanel();
        BrillianceCalcBoxNonLinear = new javax.swing.JComboBox();
        BrillianceCalcStartNonLinear = new javax.swing.JButton();
        BrillianceCalcSaveNonLinear = new javax.swing.JButton();
        BrilminvalueNonLinear = new javax.swing.JTextField();
        BrilminvaluelabelNonLinear = new javax.swing.JLabel();
        BrilminvalueunitlabelNonLinear = new javax.swing.JLabel();
        BrilmaxvalueunitlabelNonLinear = new javax.swing.JLabel();
        BrilmaxvalueNonLinear = new javax.swing.JTextField();
        BrilmaxvaluelabelNonLinear = new javax.swing.JLabel();
        jCheckBoxSpreadNonLinear = new javax.swing.JCheckBox();
        BrilProgressBarNonLinear = new javax.swing.JProgressBar();
        jAngleLabelNonLinear = new javax.swing.JLabel();
        angleValueNonLinear = new javax.swing.JTextField();
        angleValueUnitLableNonLinear = new javax.swing.JLabel();
        jEnergyLabelNonLinear = new javax.swing.JLabel();
        energyValueNonLinear = new javax.swing.JTextField();
        energyValueUnitLableNonLinear = new javax.swing.JLabel();
        OrderNumber = new javax.swing.JTextField();
        OrderNumberLabel = new javax.swing.JLabel();
        BrillianceCalcGraphNonLinear = new javax.swing.JPanel();
        gfCalc = new javax.swing.JFrame();
        GFParam = new javax.swing.JPanel();
        GFCalcBox = new javax.swing.JComboBox();
        GFCalcStart = new javax.swing.JButton();
        GFCalcSave = new javax.swing.JButton();
        GFminvalue = new javax.swing.JTextField();
        GFminvaluelabel = new javax.swing.JLabel();
        GFminvalueunitlabel = new javax.swing.JLabel();
        GFmaxvalueunitlabel = new javax.swing.JLabel();
        GFmaxvalue = new javax.swing.JTextField();
        GFmaxvaluelabel = new javax.swing.JLabel();
        GFProgressBar = new javax.swing.JProgressBar();
        GFValueSelectionBox = new javax.swing.JComboBox();
        GFCalcGraph = new javax.swing.JPanel();
        polarizationCalc = new javax.swing.JFrame();
        polarizationParam = new javax.swing.JPanel();
        polarizationCalcBox = new javax.swing.JComboBox();
        polarizationCalcStart = new javax.swing.JButton();
        polarizationCalcSave = new javax.swing.JButton();
        polminvalue = new javax.swing.JTextField();
        polminvaluelabel = new javax.swing.JLabel();
        polminvalueunitlabel = new javax.swing.JLabel();
        polmaxvalueunitlabel = new javax.swing.JLabel();
        polmaxvalue = new javax.swing.JTextField();
        polmaxvaluelabel = new javax.swing.JLabel();
        jPolCheckBoxSpread = new javax.swing.JCheckBox();
        polProgressBar = new javax.swing.JProgressBar();
        jPolAngleLabel = new javax.swing.JLabel();
        polAngleValue = new javax.swing.JTextField();
        polAngleValueUnitLable = new javax.swing.JLabel();
        jPolEnergyLabel = new javax.swing.JLabel();
        polEnergyValue = new javax.swing.JTextField();
        polEnergyValueUnitLable = new javax.swing.JLabel();
        polarizationCalcGraph = new javax.swing.JPanel();
        polarizationCalcNonLinear = new javax.swing.JFrame();
        polarizationParamcNonLinear = new javax.swing.JPanel();
        polarizationCalcBoxNonLinear = new javax.swing.JComboBox();
        polarizationCalcStartNonLinear = new javax.swing.JButton();
        polarizationCalcSaveNonLinear = new javax.swing.JButton();
        polminvalueNonLinear = new javax.swing.JTextField();
        polminvaluelabelNonLinear = new javax.swing.JLabel();
        polminvalueunitlabelNonLinear = new javax.swing.JLabel();
        polmaxvalueunitlabelNonLinear = new javax.swing.JLabel();
        polmaxvalueNonLinear = new javax.swing.JTextField();
        polmaxvaluelabelNonLinear = new javax.swing.JLabel();
        jPolCheckBoxSpreadNonLinear = new javax.swing.JCheckBox();
        polProgressBarNonLinear = new javax.swing.JProgressBar();
        jPolAngleLabelNonLinear = new javax.swing.JLabel();
        polAngleValueNonLinear = new javax.swing.JTextField();
        polAngleValueUnitLableNonLinear = new javax.swing.JLabel();
        jPolEnergyLabelNonLinear = new javax.swing.JLabel();
        polEnergyValueNonLinear = new javax.swing.JTextField();
        polEnergyValueUnitLableNonLinear = new javax.swing.JLabel();
        OrderNumberNonLinear = new javax.swing.JTextField();
        OrderNumberLabelNonLinear = new javax.swing.JLabel();
        polarizationCalcGraphNonLinear = new javax.swing.JPanel();
        rayProgressFrame = new javax.swing.JFrame();
        jRayProgressBar = new javax.swing.JProgressBar();
        jRayStopButton = new javax.swing.JButton();
        jLabelPartialFlux = new javax.swing.JLabel();
        buttonGroupPolarization = new javax.swing.ButtonGroup();
        buttonGroupSkin = new javax.swing.ButtonGroup();
        jScrollPane1 = new javax.swing.JScrollPane();
        jPanel1 = new javax.swing.JPanel();
        jPanel_el = new javax.swing.JPanel();
        energylabel = new javax.swing.JLabel();
        energyvalue = new javax.swing.JTextField();
        energyunitlabel = new javax.swing.JLabel();
        chargelabel = new javax.swing.JLabel();
        chargevalue = new javax.swing.JTextField();
        chargeunitlabel = new javax.swing.JLabel();
        spreadlabel = new javax.swing.JLabel();
        spreadvalue = new javax.swing.JTextField();
        elengthlabel = new javax.swing.JLabel();
        elengthvalue = new javax.swing.JTextField();
        elengthunitlabel = new javax.swing.JLabel();
        eemitxlabel = new javax.swing.JLabel();
        eemitxvalue = new javax.swing.JTextField();
        eemitxunitlabel = new javax.swing.JLabel();
        ebetaxlabel = new javax.swing.JLabel();
        ebetaxvalue = new javax.swing.JTextField();
        ebetaxunitlabel = new javax.swing.JLabel();
        eemitylabel = new javax.swing.JLabel();
        eemityvalue = new javax.swing.JTextField();
        eemityunitlabel = new javax.swing.JLabel();
        chargeunitlabel1 = new javax.swing.JLabel();
        jPanel_ph = new javax.swing.JPanel();
        phenergylabel = new javax.swing.JLabel();
        phenergyvalue = new javax.swing.JTextField();
        phenergyunitlabel = new javax.swing.JLabel();
        pulseenergylabel = new javax.swing.JLabel();
        pulseenergyvalue = new javax.swing.JTextField();
        pulseenergyunitlabel = new javax.swing.JLabel();
        puslelengthlabel = new javax.swing.JLabel();
        pulselengthvalue = new javax.swing.JTextField();
        pulselengthunitlabel = new javax.swing.JLabel();
        pulserellabel = new javax.swing.JLabel();
        pulsefreqlabel = new javax.swing.JLabel();
        pulsedelaylabel = new javax.swing.JLabel();
        pulserelvalue = new javax.swing.JTextField();
        pulsefreqvalue = new javax.swing.JTextField();
        pulsedelayvalue = new javax.swing.JTextField();
        pulserelunitlable = new javax.swing.JLabel();
        pulsefrequnitlabel = new javax.swing.JLabel();
        pulsedelayunitlabel = new javax.swing.JLabel();
        ebetayunitlabel = new javax.swing.JLabel();
        ebetayvalue = new javax.swing.JTextField();
        ebetaylabel = new javax.swing.JLabel();
        jPanel_exec = new javax.swing.JPanel();
        startbutton = new javax.swing.JButton();
        MainProgressBar = new javax.swing.JProgressBar();
        jTabbedPane1 = new javax.swing.JTabbedPane();
        jPanel_xflux = new javax.swing.JPanel();
        jPanel_xflux_left = new javax.swing.JPanel();
        jPanel_xflux_right = new javax.swing.JPanel();
        jPanel_xenergy = new javax.swing.JPanel();
        jPanel_xenergy_left = new javax.swing.JPanel();
        jPanel_xenergy_right = new javax.swing.JPanel();
        jPanel_slider = new javax.swing.JPanel();
        jSlider_pickup = new javax.swing.JSlider();
        totalFluxLabel = new javax.swing.JLabel();
        totalFluxAngleLabel = new javax.swing.JLabel();
        jPanel_sh = new javax.swing.JPanel();
        eshiftxlabel = new javax.swing.JLabel();
        eshiftylabel = new javax.swing.JLabel();
        eshiftzlabel = new javax.swing.JLabel();
        pulseanglelabel = new javax.swing.JLabel();
        eshiftxvalue = new javax.swing.JTextField();
        eshiftyvalue = new javax.swing.JTextField();
        eshiftzvalue = new javax.swing.JTextField();
        pulseanglevalue = new javax.swing.JTextField();
        eshiftxunitlabel = new javax.swing.JLabel();
        eshiftyunitlabel = new javax.swing.JLabel();
        eshiftzunitlabel = new javax.swing.JLabel();
        pulseangleunitlabel = new javax.swing.JLabel();
        jMenuBarMain = new javax.swing.JMenuBar();
        jMenuFile = new javax.swing.JMenu();
        jMenuItemSaveParam = new javax.swing.JMenuItem();
        jMenuItemLoadParam = new javax.swing.JMenuItem();
        jSeparator1 = new javax.swing.JPopupMenu.Separator();
        jMenuItemExit = new javax.swing.JMenuItem();
        jMenuCalc = new javax.swing.JMenu();
        jMenuItemBrilliance = new javax.swing.JMenuItem();
        jMenuItemBrillianceNonLinear = new javax.swing.JMenuItem();
        jMenuItemPolarization = new javax.swing.JMenuItem();
        jMenuItemPolarizationNonLinear = new javax.swing.JMenuItem();
        jMenuItemGeometricFactor = new javax.swing.JMenuItem();
        jMenuShadow = new javax.swing.JMenu();
        jMenuItemSource = new javax.swing.JMenuItem();
        jMenuItemSourceParam = new javax.swing.JMenuItem();
        jMenuPolarization = new javax.swing.JMenu();
        jRadioButtonMenuItemUnPolarized = new javax.swing.JRadioButtonMenuItem();
        jRadioButtonMenuItemLinearPolarized = new javax.swing.JRadioButtonMenuItem();
        jRadioButtonMenuItemCircularPolarized = new javax.swing.JRadioButtonMenuItem();
        jRadioButtonMenuItemAutoPolarized = new javax.swing.JRadioButtonMenuItem();
        jSeparator3 = new javax.swing.JPopupMenu.Separator();
        jMenuItemConv = new javax.swing.JMenuItem();
        jMenuOptions = new javax.swing.JMenu();
        jMenuItemLaserPolarization = new javax.swing.JMenuItem();
        jSeparator5 = new javax.swing.JPopupMenu.Separator();
        jMenuItemSize = new javax.swing.JMenuItem();
        jMenuItemNumerical = new javax.swing.JMenuItem();
        jSeparator2 = new javax.swing.JPopupMenu.Separator();
        jMenuItemOrderNumber = new javax.swing.JMenuItem();
        jCheckBoxMenuItemSpread = new javax.swing.JCheckBoxMenuItem();
        jSeparator4 = new javax.swing.JPopupMenu.Separator();
        jMenuSkin = new javax.swing.JMenu();
        jRadioButtonMenuDefault = new javax.swing.JRadioButtonMenuItem();
        jRadioButtonMenuSystem = new javax.swing.JRadioButtonMenuItem();
        jRadioButtonMenuNimbus = new javax.swing.JRadioButtonMenuItem();
        jMenuHelp = new javax.swing.JMenu();
        HelpItem = new javax.swing.JMenuItem();
        jMenuItemAbout = new javax.swing.JMenuItem();

        brillianceCalc.setTitle("Linear brilliance box");
        brillianceCalc.setMinimumSize(new java.awt.Dimension(800, 500));

        BrillianceParam.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Plot parameter selection", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION));

        BrillianceCalcBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Laser-electron angle", "Delay", "Z-shift", "Beta function", "Emittance", "X-emittance", "Y-emittance", "Rayleigh length", "Waist semi-width", "Energy spread", "X-ray energy", "Observation angle" }));
        BrillianceCalcBox.setSelectedIndex(10);
        BrillianceCalcBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                BrillianceCalcBoxActionPerformed(evt);
            }
        });

        BrillianceCalcStart.setText("Calculate");
        BrillianceCalcStart.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                BrillianceCalcStartActionPerformed(evt);
            }
        });

        BrillianceCalcSave.setText("Save");
        BrillianceCalcSave.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                BrillianceCalcSaveActionPerformed(evt);
            }
        });

        Brilminvalue.setText("0");
        Brilminvalue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                BrilminvalueFocusLost(evt);
            }
        });
        Brilminvalue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                BrilminvalueActionPerformed(evt);
            }
        });

        Brilminvaluelabel.setText("Min value");

        Brilminvalueunitlabel.setText("keV");

        Brilmaxvalueunitlabel.setText("keV");

        Brilmaxvalue.setText("100");
        Brilmaxvalue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                BrilmaxvalueFocusLost(evt);
            }
        });
        Brilmaxvalue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                BrilmaxvalueActionPerformed(evt);
            }
        });

        Brilmaxvaluelabel.setText("Max value");

        jCheckBoxSpread.setText("Spread");
        jCheckBoxSpread.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jCheckBoxSpreadActionPerformed(evt);
            }
        });

        jAngleLabel.setText("Angle");

        angleValue.setText("0");
        angleValue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                angleValueFocusLost(evt);
            }
        });
        angleValue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                angleValueActionPerformed(evt);
            }
        });

        angleValueUnitLable.setText("mrad");

        jEnergyLabel.setText("Energy");

        energyValue.setText("46");
        energyValue.setEnabled(false);
        energyValue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                energyValueFocusLost(evt);
            }
        });
        energyValue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                energyValueActionPerformed(evt);
            }
        });

        energyValueUnitLable.setText("keV");

        javax.swing.GroupLayout BrillianceParamLayout = new javax.swing.GroupLayout(BrillianceParam);
        BrillianceParam.setLayout(BrillianceParamLayout);
        BrillianceParamLayout.setHorizontalGroup(
            BrillianceParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(BrillianceParamLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(BrillianceParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(BrillianceParamLayout.createSequentialGroup()
                        .addComponent(BrillianceCalcBox, javax.swing.GroupLayout.PREFERRED_SIZE, 372, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(29, 29, 29))
                    .addGroup(BrillianceParamLayout.createSequentialGroup()
                        .addComponent(Brilminvaluelabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(Brilminvalue, javax.swing.GroupLayout.PREFERRED_SIZE, 45, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(Brilminvalueunitlabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(Brilmaxvaluelabel, javax.swing.GroupLayout.PREFERRED_SIZE, 57, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(Brilmaxvalue, javax.swing.GroupLayout.PREFERRED_SIZE, 43, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(Brilmaxvalueunitlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 36, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(jCheckBoxSpread)
                        .addGap(18, 18, 18)))
                .addGroup(BrillianceParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(BrillianceParamLayout.createSequentialGroup()
                        .addComponent(BrillianceCalcStart)
                        .addGap(31, 31, 31)
                        .addComponent(BrillianceCalcSave, javax.swing.GroupLayout.PREFERRED_SIZE, 77, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addComponent(BrilProgressBar, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addGap(18, 18, 18)
                .addGroup(BrillianceParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(BrillianceParamLayout.createSequentialGroup()
                        .addComponent(jAngleLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(angleValue, javax.swing.GroupLayout.PREFERRED_SIZE, 33, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(angleValueUnitLable, javax.swing.GroupLayout.PREFERRED_SIZE, 36, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(BrillianceParamLayout.createSequentialGroup()
                        .addComponent(jEnergyLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(energyValue, javax.swing.GroupLayout.PREFERRED_SIZE, 33, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(energyValueUnitLable, javax.swing.GroupLayout.PREFERRED_SIZE, 36, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap(16, Short.MAX_VALUE))
        );
        BrillianceParamLayout.setVerticalGroup(
            BrillianceParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(BrillianceParamLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(BrillianceParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(BrillianceCalcBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(BrillianceCalcStart)
                    .addComponent(BrillianceCalcSave)
                    .addComponent(jAngleLabel)
                    .addComponent(angleValue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(angleValueUnitLable, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(BrillianceParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(BrillianceParamLayout.createSequentialGroup()
                        .addGroup(BrillianceParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(Brilminvaluelabel)
                            .addComponent(Brilminvalue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(Brilminvalueunitlabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(Brilmaxvaluelabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(Brilmaxvalue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(Brilmaxvalueunitlabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jCheckBoxSpread, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(BrilProgressBar, javax.swing.GroupLayout.DEFAULT_SIZE, 34, Short.MAX_VALUE))
                        .addGap(8, 8, 8))
                    .addGroup(BrillianceParamLayout.createSequentialGroup()
                        .addGroup(BrillianceParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(energyValue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(energyValueUnitLable, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jEnergyLabel))
                        .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))))
        );

        BrillianceCalcGraph.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Spectral brilliance", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION));
        BrillianceCalcGraph.setPreferredSize(new java.awt.Dimension(639, 215));
        BrillianceCalcGraph.setRequestFocusEnabled(false);

        javax.swing.GroupLayout BrillianceCalcGraphLayout = new javax.swing.GroupLayout(BrillianceCalcGraph);
        BrillianceCalcGraph.setLayout(BrillianceCalcGraphLayout);
        BrillianceCalcGraphLayout.setHorizontalGroup(
            BrillianceCalcGraphLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );
        BrillianceCalcGraphLayout.setVerticalGroup(
            BrillianceCalcGraphLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 213, Short.MAX_VALUE)
        );

        javax.swing.GroupLayout brillianceCalcLayout = new javax.swing.GroupLayout(brillianceCalc.getContentPane());
        brillianceCalc.getContentPane().setLayout(brillianceCalcLayout);
        brillianceCalcLayout.setHorizontalGroup(
            brillianceCalcLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(BrillianceParam, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
            .addComponent(BrillianceCalcGraph, javax.swing.GroupLayout.DEFAULT_SIZE, 753, Short.MAX_VALUE)
        );
        brillianceCalcLayout.setVerticalGroup(
            brillianceCalcLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(brillianceCalcLayout.createSequentialGroup()
                .addComponent(BrillianceParam, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(BrillianceCalcGraph, javax.swing.GroupLayout.DEFAULT_SIZE, 236, Short.MAX_VALUE))
        );

        brillianceCalcNonLinear.setTitle("Non-linear brilliance box");
        brillianceCalcNonLinear.setMinimumSize(new java.awt.Dimension(850, 500));

        BrillianceParamNonLinear.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Plot parameter selection", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION));

        BrillianceCalcBoxNonLinear.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Laser-electron angle", "Delay", "Z-shift", "Beta function", "Emittance", "X-emittance", "Y-emittance", "Rayleigh length", "Waist semi-width", "Energy spread", "X-ray energy", "Observation angle" }));
        BrillianceCalcBoxNonLinear.setSelectedIndex(10);
        BrillianceCalcBoxNonLinear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                BrillianceCalcBoxNonLinearActionPerformed(evt);
            }
        });

        BrillianceCalcStartNonLinear.setText("Calculate");
        BrillianceCalcStartNonLinear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                BrillianceCalcStartNonLinearActionPerformed(evt);
            }
        });

        BrillianceCalcSaveNonLinear.setText("Save");
        BrillianceCalcSaveNonLinear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                BrillianceCalcSaveNonLinearActionPerformed(evt);
            }
        });

        BrilminvalueNonLinear.setText("0");
        BrilminvalueNonLinear.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                BrilminvalueNonLinearFocusLost(evt);
            }
        });
        BrilminvalueNonLinear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                BrilminvalueNonLinearActionPerformed(evt);
            }
        });

        BrilminvaluelabelNonLinear.setText("Min value");

        BrilminvalueunitlabelNonLinear.setText("keV");

        BrilmaxvalueunitlabelNonLinear.setText("keV");

        BrilmaxvalueNonLinear.setText("100");
        BrilmaxvalueNonLinear.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                BrilmaxvalueNonLinearFocusLost(evt);
            }
        });
        BrilmaxvalueNonLinear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                BrilmaxvalueNonLinearActionPerformed(evt);
            }
        });

        BrilmaxvaluelabelNonLinear.setText("Max value");

        jCheckBoxSpreadNonLinear.setText("Spread");
        jCheckBoxSpreadNonLinear.setEnabled(false);
        jCheckBoxSpreadNonLinear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jCheckBoxSpreadNonLinearActionPerformed(evt);
            }
        });

        jAngleLabelNonLinear.setText("Angle");

        angleValueNonLinear.setText("0");
        angleValueNonLinear.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                angleValueNonLinearFocusLost(evt);
            }
        });
        angleValueNonLinear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                angleValueNonLinearActionPerformed(evt);
            }
        });

        angleValueUnitLableNonLinear.setText("mrad");

        jEnergyLabelNonLinear.setText("Energy");

        energyValueNonLinear.setText("46");
        energyValueNonLinear.setToolTipText("");
        energyValueNonLinear.setEnabled(false);
        energyValueNonLinear.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                energyValueNonLinearFocusLost(evt);
            }
        });
        energyValueNonLinear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                energyValueNonLinearActionPerformed(evt);
            }
        });

        energyValueUnitLableNonLinear.setText("keV");

        OrderNumber.setText("1");
        OrderNumber.setToolTipText("");
        OrderNumber.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                OrderNumberFocusLost(evt);
            }
        });
        OrderNumber.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                OrderNumberActionPerformed(evt);
            }
        });

        OrderNumberLabel.setText("Order");

        javax.swing.GroupLayout BrillianceParamNonLinearLayout = new javax.swing.GroupLayout(BrillianceParamNonLinear);
        BrillianceParamNonLinear.setLayout(BrillianceParamNonLinearLayout);
        BrillianceParamNonLinearLayout.setHorizontalGroup(
            BrillianceParamNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(BrillianceParamNonLinearLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(BrillianceParamNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(BrillianceParamNonLinearLayout.createSequentialGroup()
                        .addComponent(BrilminvaluelabelNonLinear)
                        .addGap(18, 18, 18)
                        .addComponent(BrilminvalueNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 45, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(BrilminvalueunitlabelNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 40, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(BrilmaxvaluelabelNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 57, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(BrilmaxvalueNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 43, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addComponent(BrillianceCalcBoxNonLinear, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(BrilmaxvalueunitlabelNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 36, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(BrillianceParamNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                    .addGroup(BrillianceParamNonLinearLayout.createSequentialGroup()
                        .addComponent(BrillianceCalcStartNonLinear)
                        .addGap(31, 31, 31)
                        .addComponent(BrillianceCalcSaveNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 77, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addComponent(BrilProgressBarNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 185, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(24, 24, 24)
                .addGroup(BrillianceParamNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(jCheckBoxSpreadNonLinear)
                    .addGroup(BrillianceParamNonLinearLayout.createSequentialGroup()
                        .addComponent(OrderNumberLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(OrderNumber, javax.swing.GroupLayout.PREFERRED_SIZE, 45, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addGap(33, 33, 33)
                .addGroup(BrillianceParamNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jEnergyLabelNonLinear)
                    .addComponent(jAngleLabelNonLinear))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(BrillianceParamNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(BrillianceParamNonLinearLayout.createSequentialGroup()
                        .addComponent(energyValueNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 33, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(energyValueUnitLableNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 36, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(BrillianceParamNonLinearLayout.createSequentialGroup()
                        .addComponent(angleValueNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 33, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(angleValueUnitLableNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 36, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap(24, Short.MAX_VALUE))
        );
        BrillianceParamNonLinearLayout.setVerticalGroup(
            BrillianceParamNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(BrillianceParamNonLinearLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(BrillianceParamNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(BrillianceCalcBoxNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(BrillianceCalcStartNonLinear)
                    .addComponent(BrillianceCalcSaveNonLinear)
                    .addComponent(jAngleLabelNonLinear)
                    .addComponent(angleValueNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(angleValueUnitLableNonLinear, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(OrderNumber, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(OrderNumberLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(BrillianceParamNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(BrillianceParamNonLinearLayout.createSequentialGroup()
                        .addGroup(BrillianceParamNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(BrilProgressBarNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 30, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addGroup(BrillianceParamNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                .addComponent(BrilminvaluelabelNonLinear)
                                .addComponent(BrilminvalueNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addComponent(BrilminvalueunitlabelNonLinear, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addComponent(BrilmaxvaluelabelNonLinear, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addComponent(BrilmaxvalueNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addComponent(BrilmaxvalueunitlabelNonLinear, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                        .addGap(20, 20, 20))
                    .addGroup(BrillianceParamNonLinearLayout.createSequentialGroup()
                        .addGroup(BrillianceParamNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(energyValueNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(energyValueUnitLableNonLinear, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jEnergyLabelNonLinear)
                            .addComponent(jCheckBoxSpreadNonLinear, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addGap(28, 28, 28))))
        );

        BrillianceCalcGraphNonLinear.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Spectral brilliance", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION));
        BrillianceCalcGraphNonLinear.setPreferredSize(new java.awt.Dimension(639, 215));
        BrillianceCalcGraphNonLinear.setRequestFocusEnabled(false);

        javax.swing.GroupLayout BrillianceCalcGraphNonLinearLayout = new javax.swing.GroupLayout(BrillianceCalcGraphNonLinear);
        BrillianceCalcGraphNonLinear.setLayout(BrillianceCalcGraphNonLinearLayout);
        BrillianceCalcGraphNonLinearLayout.setHorizontalGroup(
            BrillianceCalcGraphNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );
        BrillianceCalcGraphNonLinearLayout.setVerticalGroup(
            BrillianceCalcGraphNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 259, Short.MAX_VALUE)
        );

        javax.swing.GroupLayout brillianceCalcNonLinearLayout = new javax.swing.GroupLayout(brillianceCalcNonLinear.getContentPane());
        brillianceCalcNonLinear.getContentPane().setLayout(brillianceCalcNonLinearLayout);
        brillianceCalcNonLinearLayout.setHorizontalGroup(
            brillianceCalcNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(BrillianceParamNonLinear, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
            .addComponent(BrillianceCalcGraphNonLinear, javax.swing.GroupLayout.DEFAULT_SIZE, 810, Short.MAX_VALUE)
        );
        brillianceCalcNonLinearLayout.setVerticalGroup(
            brillianceCalcNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(brillianceCalcNonLinearLayout.createSequentialGroup()
                .addComponent(BrillianceParamNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(BrillianceCalcGraphNonLinear, javax.swing.GroupLayout.DEFAULT_SIZE, 282, Short.MAX_VALUE)
                .addGap(9, 9, 9))
        );

        gfCalc.setTitle("Full flux box");
        gfCalc.setMinimumSize(new java.awt.Dimension(800, 500));

        GFParam.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Plot parameter selection", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION));

        GFCalcBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Laser-electron angle", "Delay", "Z-shift", "Beta function", "Emittance", "Rayleigh length", "Waist semi-width" }));
        GFCalcBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                GFCalcBoxActionPerformed(evt);
            }
        });

        GFCalcStart.setText("Calculate");
        GFCalcStart.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                GFCalcStartActionPerformed(evt);
            }
        });

        GFCalcSave.setText("Save");
        GFCalcSave.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                GFCalcSaveActionPerformed(evt);
            }
        });

        GFminvalue.setText("0");
        GFminvalue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                GFminvalueFocusLost(evt);
            }
        });
        GFminvalue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                GFminvalueActionPerformed(evt);
            }
        });

        GFminvaluelabel.setText("Min value");

        GFminvalueunitlabel.setText("mrad");

        GFmaxvalueunitlabel.setText("mrad");

        GFmaxvalue.setText("50");
        GFmaxvalue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                GFmaxvalueFocusLost(evt);
            }
        });
        GFmaxvalue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                GFmaxvalueActionPerformed(evt);
            }
        });

        GFmaxvaluelabel.setText("Max value");

        GFValueSelectionBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Full flux", "Geopmetric factor" }));
        GFValueSelectionBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                GFValueSelectionBoxActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout GFParamLayout = new javax.swing.GroupLayout(GFParam);
        GFParam.setLayout(GFParamLayout);
        GFParamLayout.setHorizontalGroup(
            GFParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(GFParamLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(GFParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(GFParamLayout.createSequentialGroup()
                        .addComponent(GFminvaluelabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(GFminvalue, javax.swing.GroupLayout.PREFERRED_SIZE, 45, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(GFminvalueunitlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 33, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(13, 13, 13)
                        .addComponent(GFmaxvaluelabel, javax.swing.GroupLayout.PREFERRED_SIZE, 58, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(GFmaxvalue, javax.swing.GroupLayout.PREFERRED_SIZE, 43, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(GFmaxvalueunitlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 33, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(GFValueSelectionBox, 0, 138, Short.MAX_VALUE))
                    .addComponent(GFCalcBox, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(GFParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(GFParamLayout.createSequentialGroup()
                        .addComponent(GFCalcStart)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(GFCalcSave, javax.swing.GroupLayout.PREFERRED_SIZE, 77, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addComponent(GFProgressBar, javax.swing.GroupLayout.PREFERRED_SIZE, 160, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(106, Short.MAX_VALUE))
        );
        GFParamLayout.setVerticalGroup(
            GFParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(GFParamLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(GFParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(GFCalcBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(GFCalcStart)
                    .addComponent(GFCalcSave))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(GFParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(GFProgressBar, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(GFParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(GFminvaluelabel)
                        .addComponent(GFmaxvalueunitlabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(GFmaxvalue)
                        .addComponent(GFmaxvaluelabel)
                        .addComponent(GFminvalue)
                        .addComponent(GFminvalueunitlabel)
                        .addComponent(GFValueSelectionBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap())
        );

        GFCalcGraph.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Full flux", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION));
        GFCalcGraph.setPreferredSize(new java.awt.Dimension(418, 216));
        GFCalcGraph.setRequestFocusEnabled(false);

        javax.swing.GroupLayout GFCalcGraphLayout = new javax.swing.GroupLayout(GFCalcGraph);
        GFCalcGraph.setLayout(GFCalcGraphLayout);
        GFCalcGraphLayout.setHorizontalGroup(
            GFCalcGraphLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );
        GFCalcGraphLayout.setVerticalGroup(
            GFCalcGraphLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 193, Short.MAX_VALUE)
        );

        javax.swing.GroupLayout gfCalcLayout = new javax.swing.GroupLayout(gfCalc.getContentPane());
        gfCalc.getContentPane().setLayout(gfCalcLayout);
        gfCalcLayout.setHorizontalGroup(
            gfCalcLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(GFParam, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
            .addComponent(GFCalcGraph, javax.swing.GroupLayout.DEFAULT_SIZE, 734, Short.MAX_VALUE)
        );
        gfCalcLayout.setVerticalGroup(
            gfCalcLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(gfCalcLayout.createSequentialGroup()
                .addComponent(GFParam, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(GFCalcGraph, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        polarizationCalc.setTitle("Linear polarization box");
        polarizationCalc.setMinimumSize(new java.awt.Dimension(800, 500));

        polarizationParam.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Plot parameter selection", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION));

        polarizationCalcBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Laser-electron angle", "Delay", "Z-shift", "Beta function", "Emittance", "X-emittance", "Y-emittance", "Rayleigh length", "Waist semi-width", "Energy spread", "X-ray energy", "Observation angle" }));
        polarizationCalcBox.setSelectedIndex(10);
        polarizationCalcBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                polarizationCalcBoxActionPerformed(evt);
            }
        });

        polarizationCalcStart.setText("Calculate");
        polarizationCalcStart.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                polarizationCalcStartActionPerformed(evt);
            }
        });

        polarizationCalcSave.setText("Save");
        polarizationCalcSave.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                polarizationCalcSaveActionPerformed(evt);
            }
        });

        polminvalue.setText("0");
        polminvalue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                polminvalueFocusLost(evt);
            }
        });
        polminvalue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                polminvalueActionPerformed(evt);
            }
        });

        polminvaluelabel.setText("Min value");

        polminvalueunitlabel.setText("keV");

        polmaxvalueunitlabel.setText("keV");

        polmaxvalue.setText("100");
        polmaxvalue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                polmaxvalueFocusLost(evt);
            }
        });
        polmaxvalue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                polmaxvalueActionPerformed(evt);
            }
        });

        polmaxvaluelabel.setText("Max value");

        jPolCheckBoxSpread.setText("Spread");
        jPolCheckBoxSpread.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jPolCheckBoxSpreadActionPerformed(evt);
            }
        });

        jPolAngleLabel.setText("Angle");

        polAngleValue.setText("0");
        polAngleValue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                polAngleValueFocusLost(evt);
            }
        });
        polAngleValue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                polAngleValueActionPerformed(evt);
            }
        });

        polAngleValueUnitLable.setText("mrad");

        jPolEnergyLabel.setText("Energy");

        polEnergyValue.setText("46");
        polEnergyValue.setEnabled(false);
        polEnergyValue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                polEnergyValueFocusLost(evt);
            }
        });
        polEnergyValue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                polEnergyValueActionPerformed(evt);
            }
        });

        polEnergyValueUnitLable.setText("keV");

        javax.swing.GroupLayout polarizationParamLayout = new javax.swing.GroupLayout(polarizationParam);
        polarizationParam.setLayout(polarizationParamLayout);
        polarizationParamLayout.setHorizontalGroup(
            polarizationParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(polarizationParamLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(polarizationParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(polarizationParamLayout.createSequentialGroup()
                        .addComponent(polarizationCalcBox, javax.swing.GroupLayout.PREFERRED_SIZE, 372, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(29, 29, 29))
                    .addGroup(polarizationParamLayout.createSequentialGroup()
                        .addComponent(polminvaluelabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(polminvalue, javax.swing.GroupLayout.PREFERRED_SIZE, 45, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(polminvalueunitlabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(polmaxvaluelabel, javax.swing.GroupLayout.PREFERRED_SIZE, 57, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(polmaxvalue, javax.swing.GroupLayout.PREFERRED_SIZE, 43, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(polmaxvalueunitlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 36, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(jPolCheckBoxSpread)
                        .addGap(18, 18, 18)))
                .addGroup(polarizationParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(polarizationParamLayout.createSequentialGroup()
                        .addComponent(polarizationCalcStart)
                        .addGap(31, 31, 31)
                        .addComponent(polarizationCalcSave, javax.swing.GroupLayout.PREFERRED_SIZE, 77, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addComponent(polProgressBar, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addGap(18, 18, 18)
                .addGroup(polarizationParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(polarizationParamLayout.createSequentialGroup()
                        .addComponent(jPolAngleLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(polAngleValue, javax.swing.GroupLayout.PREFERRED_SIZE, 33, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(polAngleValueUnitLable, javax.swing.GroupLayout.PREFERRED_SIZE, 36, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(polarizationParamLayout.createSequentialGroup()
                        .addComponent(jPolEnergyLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(polEnergyValue, javax.swing.GroupLayout.PREFERRED_SIZE, 33, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(polEnergyValueUnitLable, javax.swing.GroupLayout.PREFERRED_SIZE, 36, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap(16, Short.MAX_VALUE))
        );
        polarizationParamLayout.setVerticalGroup(
            polarizationParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(polarizationParamLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(polarizationParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(polarizationCalcBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(polarizationCalcStart)
                    .addComponent(polarizationCalcSave)
                    .addComponent(jPolAngleLabel)
                    .addComponent(polAngleValue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(polAngleValueUnitLable, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(polarizationParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(polarizationParamLayout.createSequentialGroup()
                        .addGroup(polarizationParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(polminvaluelabel)
                            .addComponent(polminvalue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(polminvalueunitlabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(polmaxvaluelabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(polmaxvalue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(polmaxvalueunitlabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jPolCheckBoxSpread, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(polProgressBar, javax.swing.GroupLayout.DEFAULT_SIZE, 34, Short.MAX_VALUE))
                        .addGap(8, 8, 8))
                    .addGroup(polarizationParamLayout.createSequentialGroup()
                        .addGroup(polarizationParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(polEnergyValue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(polEnergyValueUnitLable, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jPolEnergyLabel))
                        .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))))
        );

        polarizationCalcGraph.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Polarization parameters", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION));
        polarizationCalcGraph.setPreferredSize(new java.awt.Dimension(639, 215));
        polarizationCalcGraph.setRequestFocusEnabled(false);

        javax.swing.GroupLayout polarizationCalcGraphLayout = new javax.swing.GroupLayout(polarizationCalcGraph);
        polarizationCalcGraph.setLayout(polarizationCalcGraphLayout);
        polarizationCalcGraphLayout.setHorizontalGroup(
            polarizationCalcGraphLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );
        polarizationCalcGraphLayout.setVerticalGroup(
            polarizationCalcGraphLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 213, Short.MAX_VALUE)
        );

        javax.swing.GroupLayout polarizationCalcLayout = new javax.swing.GroupLayout(polarizationCalc.getContentPane());
        polarizationCalc.getContentPane().setLayout(polarizationCalcLayout);
        polarizationCalcLayout.setHorizontalGroup(
            polarizationCalcLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(polarizationParam, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
            .addComponent(polarizationCalcGraph, javax.swing.GroupLayout.DEFAULT_SIZE, 753, Short.MAX_VALUE)
        );
        polarizationCalcLayout.setVerticalGroup(
            polarizationCalcLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(polarizationCalcLayout.createSequentialGroup()
                .addComponent(polarizationParam, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(polarizationCalcGraph, javax.swing.GroupLayout.DEFAULT_SIZE, 236, Short.MAX_VALUE))
        );

        polarizationCalcNonLinear.setTitle("Non-linear polarization box");
        polarizationCalcNonLinear.setFocusTraversalPolicyProvider(true);
        polarizationCalcNonLinear.setMinimumSize(new java.awt.Dimension(850, 500));

        polarizationParamcNonLinear.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Plot parameter selection", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION));

        polarizationCalcBoxNonLinear.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Laser-electron angle", "Delay", "Z-shift", "Beta function", "Emittance", "X-emittance", "Y-emittance", "Rayleigh length", "Waist semi-width", "Energy spread", "X-ray energy", "Observation angle" }));
        polarizationCalcBoxNonLinear.setSelectedIndex(10);
        polarizationCalcBoxNonLinear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                polarizationCalcBoxNonLinearActionPerformed(evt);
            }
        });

        polarizationCalcStartNonLinear.setText("Calculate");
        polarizationCalcStartNonLinear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                polarizationCalcStartNonLinearActionPerformed(evt);
            }
        });

        polarizationCalcSaveNonLinear.setText("Save");
        polarizationCalcSaveNonLinear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                polarizationCalcSaveNonLinearActionPerformed(evt);
            }
        });

        polminvalueNonLinear.setText("0");
        polminvalueNonLinear.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                polminvalueNonLinearFocusLost(evt);
            }
        });
        polminvalueNonLinear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                polminvalueNonLinearActionPerformed(evt);
            }
        });

        polminvaluelabelNonLinear.setText("Min value");

        polminvalueunitlabelNonLinear.setText("keV");

        polmaxvalueunitlabelNonLinear.setText("keV");

        polmaxvalueNonLinear.setText("100");
        polmaxvalueNonLinear.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                polmaxvalueNonLinearFocusLost(evt);
            }
        });
        polmaxvalueNonLinear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                polmaxvalueNonLinearActionPerformed(evt);
            }
        });

        polmaxvaluelabelNonLinear.setText("Max value");

        jPolCheckBoxSpreadNonLinear.setText("Spread");
        jPolCheckBoxSpreadNonLinear.setEnabled(false);
        jPolCheckBoxSpreadNonLinear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jPolCheckBoxSpreadNonLinearActionPerformed(evt);
            }
        });

        jPolAngleLabelNonLinear.setText("Angle");

        polAngleValueNonLinear.setText("0");
        polAngleValueNonLinear.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                polAngleValueNonLinearFocusLost(evt);
            }
        });
        polAngleValueNonLinear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                polAngleValueNonLinearActionPerformed(evt);
            }
        });

        polAngleValueUnitLableNonLinear.setText("mrad");

        jPolEnergyLabelNonLinear.setText("Energy");

        polEnergyValueNonLinear.setText("46");
        polEnergyValueNonLinear.setToolTipText("");
        polEnergyValueNonLinear.setEnabled(false);
        polEnergyValueNonLinear.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                polEnergyValueNonLinearFocusLost(evt);
            }
        });
        polEnergyValueNonLinear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                polEnergyValueNonLinearActionPerformed(evt);
            }
        });

        polEnergyValueUnitLableNonLinear.setText("keV");

        OrderNumberNonLinear.setText("1");
        OrderNumberNonLinear.setToolTipText("");
        OrderNumberNonLinear.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                OrderNumberNonLinearFocusLost(evt);
            }
        });
        OrderNumberNonLinear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                OrderNumberNonLinearActionPerformed(evt);
            }
        });

        OrderNumberLabelNonLinear.setText("Order");

        javax.swing.GroupLayout polarizationParamcNonLinearLayout = new javax.swing.GroupLayout(polarizationParamcNonLinear);
        polarizationParamcNonLinear.setLayout(polarizationParamcNonLinearLayout);
        polarizationParamcNonLinearLayout.setHorizontalGroup(
            polarizationParamcNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(polarizationParamcNonLinearLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(polarizationParamcNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(polarizationParamcNonLinearLayout.createSequentialGroup()
                        .addComponent(polminvaluelabelNonLinear)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(polminvalueNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 45, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(polminvalueunitlabelNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 54, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(polmaxvaluelabelNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 57, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(polmaxvalueNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 43, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(polmaxvalueunitlabelNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 36, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addComponent(polarizationCalcBoxNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 283, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(polarizationParamcNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(polarizationParamcNonLinearLayout.createSequentialGroup()
                        .addComponent(polarizationCalcStartNonLinear)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(polarizationCalcSaveNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 77, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addComponent(polProgressBarNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 185, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGroup(polarizationParamcNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(polarizationParamcNonLinearLayout.createSequentialGroup()
                        .addGap(20, 20, 20)
                        .addComponent(OrderNumberLabelNonLinear)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(OrderNumberNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 28, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(14, 14, 14)
                        .addComponent(jPolAngleLabelNonLinear)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(polAngleValueNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 33, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(polAngleValueUnitLableNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 36, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(polarizationParamcNonLinearLayout.createSequentialGroup()
                        .addGap(25, 25, 25)
                        .addComponent(jPolCheckBoxSpreadNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 66, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(jPolEnergyLabelNonLinear)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(polEnergyValueNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 33, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(10, 10, 10)
                        .addComponent(polEnergyValueUnitLableNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 36, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap(13, Short.MAX_VALUE))
        );
        polarizationParamcNonLinearLayout.setVerticalGroup(
            polarizationParamcNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(polarizationParamcNonLinearLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(polarizationParamcNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(polarizationCalcBoxNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(polarizationCalcStartNonLinear)
                    .addComponent(polarizationCalcSaveNonLinear)
                    .addComponent(jPolAngleLabelNonLinear)
                    .addComponent(polAngleValueNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(polAngleValueUnitLableNonLinear, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(OrderNumberNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(OrderNumberLabelNonLinear))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(polarizationParamcNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(polarizationParamcNonLinearLayout.createSequentialGroup()
                        .addGroup(polarizationParamcNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(polarizationParamcNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                .addComponent(polminvaluelabelNonLinear)
                                .addComponent(polminvalueNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addComponent(polminvalueunitlabelNonLinear, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addComponent(polmaxvaluelabelNonLinear, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addComponent(polmaxvalueNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addComponent(polmaxvalueunitlabelNonLinear, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                            .addGroup(polarizationParamcNonLinearLayout.createSequentialGroup()
                                .addComponent(polProgressBarNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 30, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addGap(0, 0, Short.MAX_VALUE)))
                        .addGap(8, 8, 8))
                    .addGroup(polarizationParamcNonLinearLayout.createSequentialGroup()
                        .addGroup(polarizationParamcNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(polEnergyValueNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(polEnergyValueUnitLableNonLinear, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jPolEnergyLabelNonLinear)
                            .addComponent(jPolCheckBoxSpreadNonLinear))
                        .addGap(35, 35, 35))))
        );

        polarizationCalcGraphNonLinear.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Polarization parameters", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION));
        polarizationCalcGraphNonLinear.setPreferredSize(new java.awt.Dimension(639, 215));
        polarizationCalcGraphNonLinear.setRequestFocusEnabled(false);

        javax.swing.GroupLayout polarizationCalcGraphNonLinearLayout = new javax.swing.GroupLayout(polarizationCalcGraphNonLinear);
        polarizationCalcGraphNonLinear.setLayout(polarizationCalcGraphNonLinearLayout);
        polarizationCalcGraphNonLinearLayout.setHorizontalGroup(
            polarizationCalcGraphNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );
        polarizationCalcGraphNonLinearLayout.setVerticalGroup(
            polarizationCalcGraphNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 210, Short.MAX_VALUE)
        );

        javax.swing.GroupLayout polarizationCalcNonLinearLayout = new javax.swing.GroupLayout(polarizationCalcNonLinear.getContentPane());
        polarizationCalcNonLinear.getContentPane().setLayout(polarizationCalcNonLinearLayout);
        polarizationCalcNonLinearLayout.setHorizontalGroup(
            polarizationCalcNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(polarizationParamcNonLinear, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
            .addComponent(polarizationCalcGraphNonLinear, javax.swing.GroupLayout.DEFAULT_SIZE, 762, Short.MAX_VALUE)
        );
        polarizationCalcNonLinearLayout.setVerticalGroup(
            polarizationCalcNonLinearLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(polarizationCalcNonLinearLayout.createSequentialGroup()
                .addComponent(polarizationParamcNonLinear, javax.swing.GroupLayout.PREFERRED_SIZE, 113, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(polarizationCalcGraphNonLinear, javax.swing.GroupLayout.DEFAULT_SIZE, 233, Short.MAX_VALUE))
        );

        rayProgressFrame.setTitle("Ray generation progress");
        rayProgressFrame.setAlwaysOnTop(true);
        rayProgressFrame.setMinimumSize(new java.awt.Dimension(400, 100));
        rayProgressFrame.setResizable(false);

        jRayStopButton.setText("Stop");
        jRayStopButton.setToolTipText("");
        jRayStopButton.setEnabled(false);
        jRayStopButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jRayStopButtonActionPerformed(evt);
            }
        });

        jLabelPartialFlux.setText("Flux: ");

        javax.swing.GroupLayout rayProgressFrameLayout = new javax.swing.GroupLayout(rayProgressFrame.getContentPane());
        rayProgressFrame.getContentPane().setLayout(rayProgressFrameLayout);
        rayProgressFrameLayout.setHorizontalGroup(
            rayProgressFrameLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, rayProgressFrameLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(rayProgressFrameLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(jLabelPartialFlux, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(jRayProgressBar, javax.swing.GroupLayout.DEFAULT_SIZE, 273, Short.MAX_VALUE))
                .addGap(18, 18, 18)
                .addComponent(jRayStopButton)
                .addGap(21, 21, 21))
        );
        rayProgressFrameLayout.setVerticalGroup(
            rayProgressFrameLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(rayProgressFrameLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(rayProgressFrameLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(jRayProgressBar, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(jRayStopButton, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 21, Short.MAX_VALUE)
                .addComponent(jLabelPartialFlux)
                .addContainerGap())
        );

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setTitle("TSourceNX");
        setMinimumSize(new java.awt.Dimension(750, 0));
        addMouseMotionListener(new java.awt.event.MouseMotionAdapter() {
            public void mouseMoved(java.awt.event.MouseEvent evt) {
                formMouseMoved(evt);
            }
        });

        jScrollPane1.setHorizontalScrollBarPolicy(javax.swing.ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
        jScrollPane1.setDoubleBuffered(true);
        jScrollPane1.setHorizontalScrollBar(null);
        jScrollPane1.setMinimumSize(new java.awt.Dimension(720, 6));
        jScrollPane1.setPreferredSize(new java.awt.Dimension(780, 633));

        jPanel1.setMinimumSize(new java.awt.Dimension(720, 0));
        jPanel1.setPreferredSize(new java.awt.Dimension(768, 631));

        jPanel_el.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Electron bunch parameters", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION));
        jPanel_el.setToolTipText("Input parameters of the electron bunch");
        jPanel_el.setPreferredSize(new java.awt.Dimension(247, 350));
        jPanel_el.setRequestFocusEnabled(false);
        jPanel_el.setVerifyInputWhenFocusTarget(false);

        energylabel.setText("Electron energy");

        energyvalue.setText("50");
        energyvalue.setToolTipText("");
        energyvalue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                energyvalueFocusLost(evt);
            }
        });
        energyvalue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                energyvalueActionPerformed(evt);
            }
        });

        energyunitlabel.setText("MeV");

        chargelabel.setText("Charge");

        chargevalue.setText("0.2");
        chargevalue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                chargevalueFocusLost(evt);
            }
        });
        chargevalue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                chargevalueActionPerformed(evt);
            }
        });

        chargeunitlabel.setText("nQ");

        spreadlabel.setText("Gamma-spread");

        spreadvalue.setText("0.25");
        spreadvalue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                spreadvalueFocusLost(evt);
            }
        });
        spreadvalue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                spreadvalueActionPerformed(evt);
            }
        });

        elengthlabel.setText("Length");

        elengthvalue.setText("10");
        elengthvalue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                elengthvalueFocusLost(evt);
            }
        });
        elengthvalue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                elengthvalueActionPerformed(evt);
            }
        });

        elengthunitlabel.setText("ps");

        eemitxlabel.setText("X-emittance");

        eemitxvalue.setText("1");
        eemitxvalue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                eemitxvalueFocusLost(evt);
            }
        });
        eemitxvalue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                eemitxvalueActionPerformed(evt);
            }
        });

        eemitxunitlabel.setText("mm*mrad");

        ebetaxlabel.setText("Beta x");

        ebetaxvalue.setText("20");
        ebetaxvalue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                ebetaxvalueFocusLost(evt);
            }
        });
        ebetaxvalue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ebetaxvalueActionPerformed(evt);
            }
        });

        ebetaxunitlabel.setText("mm");

        eemitylabel.setText("Y-emittance");

        eemityvalue.setText("1");
        eemityvalue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                eemityvalueFocusLost(evt);
            }
        });
        eemityvalue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                eemityvalueActionPerformed(evt);
            }
        });

        eemityunitlabel.setText("mm*mrad");

        chargeunitlabel1.setText("%");

        javax.swing.GroupLayout jPanel_elLayout = new javax.swing.GroupLayout(jPanel_el);
        jPanel_el.setLayout(jPanel_elLayout);
        jPanel_elLayout.setHorizontalGroup(
            jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel_elLayout.createSequentialGroup()
                .addGap(6, 6, 6)
                .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel_elLayout.createSequentialGroup()
                        .addComponent(ebetaxlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 62, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(ebetaxvalue, javax.swing.GroupLayout.PREFERRED_SIZE, 41, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel_elLayout.createSequentialGroup()
                        .addGap(0, 0, Short.MAX_VALUE)
                        .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addGroup(jPanel_elLayout.createSequentialGroup()
                                .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(energylabel, javax.swing.GroupLayout.PREFERRED_SIZE, 92, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addComponent(elengthlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 62, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                                        .addComponent(eemitylabel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                        .addComponent(eemitxlabel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                        .addComponent(spreadlabel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, 75, Short.MAX_VALUE)
                                        .addComponent(chargelabel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                                .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addGroup(jPanel_elLayout.createSequentialGroup()
                                        .addGap(20, 20, 20)
                                        .addComponent(energyvalue, javax.swing.GroupLayout.PREFERRED_SIZE, 41, javax.swing.GroupLayout.PREFERRED_SIZE))
                                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel_elLayout.createSequentialGroup()
                                        .addGap(13, 13, 13)
                                        .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                                            .addComponent(spreadvalue, javax.swing.GroupLayout.DEFAULT_SIZE, 41, Short.MAX_VALUE)
                                            .addComponent(elengthvalue)
                                            .addComponent(eemitxvalue)
                                            .addComponent(eemityvalue)))))
                            .addComponent(chargevalue, javax.swing.GroupLayout.PREFERRED_SIZE, 41, javax.swing.GroupLayout.PREFERRED_SIZE))))
                .addGap(18, 18, 18)
                .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel_elLayout.createSequentialGroup()
                        .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(elengthunitlabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addGroup(jPanel_elLayout.createSequentialGroup()
                                .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(eemitxunitlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 59, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addComponent(ebetaxunitlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 48, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addComponent(eemityunitlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 59, javax.swing.GroupLayout.PREFERRED_SIZE))
                                .addGap(0, 0, Short.MAX_VALUE)))
                        .addContainerGap())
                    .addGroup(jPanel_elLayout.createSequentialGroup()
                        .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(chargeunitlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 24, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(energyunitlabel)
                            .addComponent(chargeunitlabel1, javax.swing.GroupLayout.PREFERRED_SIZE, 24, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addGap(0, 0, Short.MAX_VALUE))))
        );
        jPanel_elLayout.setVerticalGroup(
            jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel_elLayout.createSequentialGroup()
                .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(energylabel)
                    .addComponent(energyvalue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(energyunitlabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(chargelabel)
                    .addComponent(chargevalue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(chargeunitlabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(spreadlabel)
                    .addComponent(spreadvalue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(chargeunitlabel1))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(elengthunitlabel)
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(elengthlabel)
                        .addComponent(elengthvalue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(eemitxunitlabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(eemitxvalue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(eemitxlabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(eemitylabel)
                    .addComponent(eemityvalue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(eemityunitlabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addGap(8, 8, 8)
                .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(ebetaxunitlabel)
                    .addComponent(ebetaxvalue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(ebetaxlabel))
                .addGap(12, 12, 12))
        );

        jPanel_ph.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Laser pulse parameters", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION));
        jPanel_ph.setToolTipText("Input parameters of the laser bunch");
        jPanel_ph.setPreferredSize(new java.awt.Dimension(223, 350));

        phenergylabel.setText("Photon energy");

        phenergyvalue.setText("1.204");
        phenergyvalue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                phenergyvalueFocusLost(evt);
            }
        });
        phenergyvalue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                phenergyvalueActionPerformed(evt);
            }
        });

        phenergyunitlabel.setText("eV");

        pulseenergylabel.setText("Pulse energy");

        pulseenergyvalue.setText("100");
        pulseenergyvalue.setToolTipText("");
        pulseenergyvalue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                pulseenergyvalueFocusLost(evt);
            }
        });
        pulseenergyvalue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                pulseenergyvalueActionPerformed(evt);
            }
        });

        pulseenergyunitlabel.setText("mJ");

        puslelengthlabel.setText("Pulse length");

        pulselengthvalue.setText("10");
        pulselengthvalue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                pulselengthvalueFocusLost(evt);
            }
        });
        pulselengthvalue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                pulselengthvalueActionPerformed(evt);
            }
        });

        pulselengthunitlabel.setText("ps");

        pulserellabel.setText("Rayleigh length");

        pulsefreqlabel.setText("Pulse frequency");

        pulsedelaylabel.setText("Delay");

        pulserelvalue.setText("0.175");
        pulserelvalue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                pulserelvalueFocusLost(evt);
            }
        });
        pulserelvalue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                pulserelvalueActionPerformed(evt);
            }
        });

        pulsefreqvalue.setText("1000");
        pulsefreqvalue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                pulsefreqvalueFocusLost(evt);
            }
        });
        pulsefreqvalue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                pulsefreqvalueActionPerformed(evt);
            }
        });

        pulsedelayvalue.setText("0");
        pulsedelayvalue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                pulsedelayvalueFocusLost(evt);
            }
        });
        pulsedelayvalue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                pulsedelayvalueActionPerformed(evt);
            }
        });

        pulserelunitlable.setText("ps");

        pulsefrequnitlabel.setText("mm");

        pulsedelayunitlabel.setText("Hz");

        ebetayunitlabel.setText("mm");

        ebetayvalue.setText("20");
        ebetayvalue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                ebetayvalueFocusLost(evt);
            }
        });
        ebetayvalue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ebetayvalueActionPerformed(evt);
            }
        });

        ebetaylabel.setText("Beta  y");

        javax.swing.GroupLayout jPanel_phLayout = new javax.swing.GroupLayout(jPanel_ph);
        jPanel_ph.setLayout(jPanel_phLayout);
        jPanel_phLayout.setHorizontalGroup(
            jPanel_phLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel_phLayout.createSequentialGroup()
                .addContainerGap(28, Short.MAX_VALUE)
                .addGroup(jPanel_phLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel_phLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                        .addComponent(puslelengthlabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(phenergylabel, javax.swing.GroupLayout.DEFAULT_SIZE, 84, Short.MAX_VALUE)
                        .addComponent(pulserellabel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(pulsefreqlabel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(pulsedelaylabel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                    .addComponent(pulseenergylabel, javax.swing.GroupLayout.PREFERRED_SIZE, 78, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(ebetaylabel, javax.swing.GroupLayout.PREFERRED_SIZE, 62, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel_phLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel_phLayout.createSequentialGroup()
                        .addGroup(jPanel_phLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(pulseenergyvalue, javax.swing.GroupLayout.DEFAULT_SIZE, 41, Short.MAX_VALUE)
                            .addComponent(phenergyvalue, javax.swing.GroupLayout.DEFAULT_SIZE, 41, Short.MAX_VALUE)
                            .addComponent(pulselengthvalue)
                            .addComponent(pulserelvalue)
                            .addComponent(pulsefreqvalue)
                            .addComponent(pulsedelayvalue))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addGroup(jPanel_phLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(phenergyunitlabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(pulseenergyunitlabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(pulselengthunitlabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(pulserelunitlable, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(pulsefrequnitlabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(pulsedelayunitlabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                    .addGroup(jPanel_phLayout.createSequentialGroup()
                        .addComponent(ebetayvalue, javax.swing.GroupLayout.PREFERRED_SIZE, 41, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(ebetayunitlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 34, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap())
        );
        jPanel_phLayout.setVerticalGroup(
            jPanel_phLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel_phLayout.createSequentialGroup()
                .addGroup(jPanel_phLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(phenergylabel)
                    .addComponent(phenergyvalue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(phenergyunitlabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel_phLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(pulseenergyvalue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(pulseenergyunitlabel)
                    .addComponent(pulseenergylabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel_phLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(pulselengthvalue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(puslelengthlabel)
                    .addComponent(pulselengthunitlabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel_phLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(pulserelvalue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(pulsefrequnitlabel)
                    .addComponent(pulserellabel))
                .addGap(8, 8, 8)
                .addGroup(jPanel_phLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(pulsefreqlabel)
                    .addComponent(pulsefreqvalue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(pulsedelayunitlabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel_phLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(pulsedelayvalue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(pulsedelaylabel)
                    .addComponent(pulserelunitlable))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addGroup(jPanel_phLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(ebetayunitlabel)
                    .addComponent(ebetayvalue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(ebetaylabel))
                .addContainerGap())
        );

        jPanel_exec.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Execution", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION));

        startbutton.setAction(startbutton.getAction());
        startbutton.setText("Start");
        startbutton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                startbuttonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout jPanel_execLayout = new javax.swing.GroupLayout(jPanel_exec);
        jPanel_exec.setLayout(jPanel_execLayout);
        jPanel_execLayout.setHorizontalGroup(
            jPanel_execLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel_execLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(startbutton, javax.swing.GroupLayout.PREFERRED_SIZE, 64, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(MainProgressBar, javax.swing.GroupLayout.PREFERRED_SIZE, 91, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        jPanel_execLayout.setVerticalGroup(
            jPanel_execLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel_execLayout.createSequentialGroup()
                .addGroup(jPanel_execLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(startbutton, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(MainProgressBar, javax.swing.GroupLayout.DEFAULT_SIZE, 23, Short.MAX_VALUE))
                .addGap(0, 5, Short.MAX_VALUE))
        );

        jTabbedPane1.setBorder(new javax.swing.border.SoftBevelBorder(javax.swing.border.BevelBorder.RAISED));
        jTabbedPane1.setAutoscrolls(true);
        jTabbedPane1.setDoubleBuffered(true);
        jTabbedPane1.setMinimumSize(new java.awt.Dimension(0, 0));
        jTabbedPane1.setName(""); // NOI18N
        jTabbedPane1.setPreferredSize(new java.awt.Dimension(713, 500));
        jTabbedPane1.setRequestFocusEnabled(false);

        jPanel_xflux.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "X-ray photon flux", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION));
        jPanel_xflux.setAutoscrolls(true);
        jPanel_xflux.setMinimumSize(new java.awt.Dimension(0, 0));

        jPanel_xflux_left.setMinimumSize(new java.awt.Dimension(0, 0));
        jPanel_xflux_left.setName(""); // NOI18N
        jPanel_xflux_left.setPreferredSize(new java.awt.Dimension(362, 318));

        javax.swing.GroupLayout jPanel_xflux_leftLayout = new javax.swing.GroupLayout(jPanel_xflux_left);
        jPanel_xflux_left.setLayout(jPanel_xflux_leftLayout);
        jPanel_xflux_leftLayout.setHorizontalGroup(
            jPanel_xflux_leftLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 343, Short.MAX_VALUE)
        );
        jPanel_xflux_leftLayout.setVerticalGroup(
            jPanel_xflux_leftLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 226, Short.MAX_VALUE)
        );

        jPanel_xflux_right.setMinimumSize(new java.awt.Dimension(0, 0));
        jPanel_xflux_right.setPreferredSize(new java.awt.Dimension(309, 318));
        jPanel_xflux_right.setRequestFocusEnabled(false);
        jPanel_xflux_right.setVerifyInputWhenFocusTarget(false);

        javax.swing.GroupLayout jPanel_xflux_rightLayout = new javax.swing.GroupLayout(jPanel_xflux_right);
        jPanel_xflux_right.setLayout(jPanel_xflux_rightLayout);
        jPanel_xflux_rightLayout.setHorizontalGroup(
            jPanel_xflux_rightLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 299, Short.MAX_VALUE)
        );
        jPanel_xflux_rightLayout.setVerticalGroup(
            jPanel_xflux_rightLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );

        javax.swing.GroupLayout jPanel_xfluxLayout = new javax.swing.GroupLayout(jPanel_xflux);
        jPanel_xflux.setLayout(jPanel_xfluxLayout);
        jPanel_xfluxLayout.setHorizontalGroup(
            jPanel_xfluxLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel_xfluxLayout.createSequentialGroup()
                .addComponent(jPanel_xflux_left, javax.swing.GroupLayout.DEFAULT_SIZE, 343, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jPanel_xflux_right, javax.swing.GroupLayout.DEFAULT_SIZE, 299, Short.MAX_VALUE)
                .addGap(4, 4, 4))
        );
        jPanel_xfluxLayout.setVerticalGroup(
            jPanel_xfluxLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jPanel_xflux_left, javax.swing.GroupLayout.DEFAULT_SIZE, 226, Short.MAX_VALUE)
            .addComponent(jPanel_xflux_right, javax.swing.GroupLayout.DEFAULT_SIZE, 226, Short.MAX_VALUE)
        );

        jTabbedPane1.addTab("Flux", jPanel_xflux);

        jPanel_xenergy.setBorder(xrayenergyborder);
        jPanel_xenergy.setAutoscrolls(true);
        jPanel_xenergy.setMinimumSize(new java.awt.Dimension(0, 0));
        jPanel_xenergy.setPreferredSize(new java.awt.Dimension(704, 352));

        jPanel_xenergy_left.setMinimumSize(new java.awt.Dimension(0, 0));
        jPanel_xenergy_left.setName(""); // NOI18N
        jPanel_xenergy_left.setPreferredSize(new java.awt.Dimension(362, 318));
        jPanel_xenergy_left.setVerifyInputWhenFocusTarget(false);

        javax.swing.GroupLayout jPanel_xenergy_leftLayout = new javax.swing.GroupLayout(jPanel_xenergy_left);
        jPanel_xenergy_left.setLayout(jPanel_xenergy_leftLayout);
        jPanel_xenergy_leftLayout.setHorizontalGroup(
            jPanel_xenergy_leftLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 349, Short.MAX_VALUE)
        );
        jPanel_xenergy_leftLayout.setVerticalGroup(
            jPanel_xenergy_leftLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );

        jPanel_xenergy_right.setMinimumSize(new java.awt.Dimension(0, 0));
        jPanel_xenergy_right.setPreferredSize(new java.awt.Dimension(309, 318));
        jPanel_xenergy_right.setVerifyInputWhenFocusTarget(false);

        javax.swing.GroupLayout jPanel_xenergy_rightLayout = new javax.swing.GroupLayout(jPanel_xenergy_right);
        jPanel_xenergy_right.setLayout(jPanel_xenergy_rightLayout);
        jPanel_xenergy_rightLayout.setHorizontalGroup(
            jPanel_xenergy_rightLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 309, Short.MAX_VALUE)
        );
        jPanel_xenergy_rightLayout.setVerticalGroup(
            jPanel_xenergy_rightLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 249, Short.MAX_VALUE)
        );

        javax.swing.GroupLayout jPanel_xenergyLayout = new javax.swing.GroupLayout(jPanel_xenergy);
        jPanel_xenergy.setLayout(jPanel_xenergyLayout);
        jPanel_xenergyLayout.setHorizontalGroup(
            jPanel_xenergyLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel_xenergyLayout.createSequentialGroup()
                .addComponent(jPanel_xenergy_left, javax.swing.GroupLayout.DEFAULT_SIZE, 349, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jPanel_xenergy_right, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
        );
        jPanel_xenergyLayout.setVerticalGroup(
            jPanel_xenergyLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel_xenergyLayout.createSequentialGroup()
                .addGap(0, 0, 0)
                .addGroup(jPanel_xenergyLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jPanel_xenergy_right, javax.swing.GroupLayout.DEFAULT_SIZE, 249, Short.MAX_VALUE)
                    .addComponent(jPanel_xenergy_left, javax.swing.GroupLayout.DEFAULT_SIZE, 249, Short.MAX_VALUE)))
        );

        jTabbedPane1.addTab("Energy", jPanel_xenergy);

        jPanel_slider.setBorder(javax.swing.BorderFactory.createBevelBorder(javax.swing.border.BevelBorder.RAISED));

        jSlider_pickup.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                jSlider_pickupStateChanged(evt);
            }
        });

        totalFluxLabel.setText("Total flux:");

        totalFluxAngleLabel.setText("Within limits: ");

        javax.swing.GroupLayout jPanel_sliderLayout = new javax.swing.GroupLayout(jPanel_slider);
        jPanel_slider.setLayout(jPanel_sliderLayout);
        jPanel_sliderLayout.setHorizontalGroup(
            jPanel_sliderLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel_sliderLayout.createSequentialGroup()
                .addGap(69, 69, 69)
                .addComponent(jSlider_pickup, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addComponent(totalFluxLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 149, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addComponent(totalFluxAngleLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 183, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        jPanel_sliderLayout.setVerticalGroup(
            jPanel_sliderLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel_sliderLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel_sliderLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(jSlider_pickup, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(jPanel_sliderLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(totalFluxLabel)
                        .addComponent(totalFluxAngleLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        jPanel_sh.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Relative position", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION));
        jPanel_sh.setPreferredSize(new java.awt.Dimension(221, 180));

        eshiftxlabel.setText("X-shift");

        eshiftylabel.setText("Y-shift");

        eshiftzlabel.setText("Z-shift");

        pulseanglelabel.setText("Angle");

        eshiftxvalue.setText("0");
        eshiftxvalue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                eshiftxvalueFocusLost(evt);
            }
        });
        eshiftxvalue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                eshiftxvalueActionPerformed(evt);
            }
        });

        eshiftyvalue.setText("0");
        eshiftyvalue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                eshiftyvalueFocusLost(evt);
            }
        });
        eshiftyvalue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                eshiftyvalueActionPerformed(evt);
            }
        });

        eshiftzvalue.setText("0");
        eshiftzvalue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                eshiftzvalueFocusLost(evt);
            }
        });
        eshiftzvalue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                eshiftzvalueActionPerformed(evt);
            }
        });

        pulseanglevalue.setText("52");
        pulseanglevalue.setToolTipText("");
        pulseanglevalue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                pulseanglevalueFocusLost(evt);
            }
        });
        pulseanglevalue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                pulseanglevalueActionPerformed(evt);
            }
        });

        eshiftxunitlabel.setText("mm");

        eshiftyunitlabel.setText("mm");

        eshiftzunitlabel.setText("mm");

        pulseangleunitlabel.setText("mrad");

        javax.swing.GroupLayout jPanel_shLayout = new javax.swing.GroupLayout(jPanel_sh);
        jPanel_sh.setLayout(jPanel_shLayout);
        jPanel_shLayout.setHorizontalGroup(
            jPanel_shLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel_shLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel_shLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(eshiftxlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 49, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addGroup(jPanel_shLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                        .addComponent(eshiftylabel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, 34, Short.MAX_VALUE)
                        .addComponent(pulseanglelabel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(eshiftzlabel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                .addGap(18, 18, 18)
                .addGroup(jPanel_shLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel_shLayout.createSequentialGroup()
                        .addComponent(eshiftzvalue, javax.swing.GroupLayout.PREFERRED_SIZE, 37, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(eshiftzunitlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 32, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(jPanel_shLayout.createSequentialGroup()
                        .addComponent(eshiftyvalue, javax.swing.GroupLayout.PREFERRED_SIZE, 37, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(eshiftyunitlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 31, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(jPanel_shLayout.createSequentialGroup()
                        .addComponent(pulseanglevalue, javax.swing.GroupLayout.PREFERRED_SIZE, 37, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(pulseangleunitlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 46, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(jPanel_shLayout.createSequentialGroup()
                        .addComponent(eshiftxvalue, javax.swing.GroupLayout.PREFERRED_SIZE, 37, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(eshiftxunitlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 31, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        jPanel_shLayout.setVerticalGroup(
            jPanel_shLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel_shLayout.createSequentialGroup()
                .addGroup(jPanel_shLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(eshiftxlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 20, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(eshiftxvalue)
                    .addComponent(eshiftxunitlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 20, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel_shLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(eshiftyvalue)
                    .addComponent(eshiftyunitlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 20, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(eshiftylabel, javax.swing.GroupLayout.PREFERRED_SIZE, 16, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel_shLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(eshiftzvalue)
                    .addComponent(eshiftzunitlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 20, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(eshiftzlabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel_shLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(pulseanglelabel)
                    .addComponent(pulseanglevalue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(pulseangleunitlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 20, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap())
        );

        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel1Layout.createSequentialGroup()
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(javax.swing.GroupLayout.Alignment.LEADING, jPanel1Layout.createSequentialGroup()
                        .addComponent(jPanel_el, javax.swing.GroupLayout.PREFERRED_SIZE, 258, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jPanel_ph, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(jPanel_sh, javax.swing.GroupLayout.DEFAULT_SIZE, 192, Short.MAX_VALUE)
                            .addComponent(jPanel_exec, javax.swing.GroupLayout.PREFERRED_SIZE, 192, Short.MAX_VALUE))
                        .addGap(0, 0, Short.MAX_VALUE))
                    .addGroup(javax.swing.GroupLayout.Alignment.LEADING, jPanel1Layout.createSequentialGroup()
                        .addContainerGap()
                        .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(jPanel_slider, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jTabbedPane1, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, 675, Short.MAX_VALUE))))
                .addGap(137, 137, 137))
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addComponent(jPanel_sh, javax.swing.GroupLayout.DEFAULT_SIZE, 155, Short.MAX_VALUE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(jPanel_exec, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addComponent(jPanel_ph, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, 213, Short.MAX_VALUE)
                    .addComponent(jPanel_el, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, 213, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jTabbedPane1, javax.swing.GroupLayout.PREFERRED_SIZE, 283, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(jPanel_slider, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addGap(67, 67, 67))
        );

        jTabbedPane1.getAccessibleContext().setAccessibleName("Flux");

        jScrollPane1.setViewportView(jPanel1);

        jMenuFile.setText("File");

        jMenuItemSaveParam.setText("Save parameters");
        jMenuItemSaveParam.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItemSaveParamActionPerformed(evt);
            }
        });
        jMenuFile.add(jMenuItemSaveParam);

        jMenuItemLoadParam.setText("Load parameters");
        jMenuItemLoadParam.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItemLoadParamActionPerformed(evt);
            }
        });
        jMenuFile.add(jMenuItemLoadParam);
        jMenuFile.add(jSeparator1);

        jMenuItemExit.setText("Exit");
        jMenuItemExit.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItemExitActionPerformed(evt);
            }
        });
        jMenuFile.add(jMenuItemExit);

        jMenuBarMain.add(jMenuFile);

        jMenuCalc.setText("Calculations");

        jMenuItemBrilliance.setText("Linear brilliance...");
        jMenuItemBrilliance.setToolTipText("");
        jMenuItemBrilliance.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItemBrillianceActionPerformed(evt);
            }
        });
        jMenuCalc.add(jMenuItemBrilliance);

        jMenuItemBrillianceNonLinear.setText("Non-linear brilliance...");
        jMenuItemBrillianceNonLinear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItemBrillianceNonLinearActionPerformed(evt);
            }
        });
        jMenuCalc.add(jMenuItemBrillianceNonLinear);

        jMenuItemPolarization.setText("Linear polarization...");
        jMenuItemPolarization.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItemPolarizationActionPerformed(evt);
            }
        });
        jMenuCalc.add(jMenuItemPolarization);

        jMenuItemPolarizationNonLinear.setText("Non-linear polarization...");
        jMenuItemPolarizationNonLinear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItemPolarizationNonLinearActionPerformed(evt);
            }
        });
        jMenuCalc.add(jMenuItemPolarizationNonLinear);

        jMenuItemGeometricFactor.setText("Full flux/Geometric factor...");
        jMenuItemGeometricFactor.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItemGeometricFactorActionPerformed(evt);
            }
        });
        jMenuCalc.add(jMenuItemGeometricFactor);

        jMenuBarMain.add(jMenuCalc);

        jMenuShadow.setText("Shadow");

        jMenuItemSource.setText("Generate source...");
        jMenuItemSource.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItemSourceActionPerformed(evt);
            }
        });
        jMenuShadow.add(jMenuItemSource);

        jMenuItemSourceParam.setText("Parameters...");
        jMenuItemSourceParam.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItemSourceParamActionPerformed(evt);
            }
        });
        jMenuShadow.add(jMenuItemSourceParam);

        jMenuPolarization.setText("Polarization...");

        buttonGroupPolarization.add(jRadioButtonMenuItemUnPolarized);
        jRadioButtonMenuItemUnPolarized.setText("Unpolarized");
        jRadioButtonMenuItemUnPolarized.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                jRadioButtonMenuItemUnPolarizedItemStateChanged(evt);
            }
        });
        jMenuPolarization.add(jRadioButtonMenuItemUnPolarized);

        buttonGroupPolarization.add(jRadioButtonMenuItemLinearPolarized);
        jRadioButtonMenuItemLinearPolarized.setSelected(true);
        jRadioButtonMenuItemLinearPolarized.setText("Linear");
        jRadioButtonMenuItemLinearPolarized.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                jRadioButtonMenuItemLinearPolarizedItemStateChanged(evt);
            }
        });
        jMenuPolarization.add(jRadioButtonMenuItemLinearPolarized);

        buttonGroupPolarization.add(jRadioButtonMenuItemCircularPolarized);
        jRadioButtonMenuItemCircularPolarized.setText("Circular");
        jRadioButtonMenuItemCircularPolarized.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                jRadioButtonMenuItemCircularPolarizedItemStateChanged(evt);
            }
        });
        jMenuPolarization.add(jRadioButtonMenuItemCircularPolarized);

        buttonGroupPolarization.add(jRadioButtonMenuItemAutoPolarized);
        jRadioButtonMenuItemAutoPolarized.setText("Automatic");
        jRadioButtonMenuItemAutoPolarized.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                jRadioButtonMenuItemAutoPolarizedItemStateChanged(evt);
            }
        });
        jMenuPolarization.add(jRadioButtonMenuItemAutoPolarized);

        jMenuShadow.add(jMenuPolarization);
        jMenuShadow.add(jSeparator3);

        jMenuItemConv.setText("Converter...");
        jMenuItemConv.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItemConvActionPerformed(evt);
            }
        });
        jMenuShadow.add(jMenuItemConv);

        jMenuBarMain.add(jMenuShadow);

        jMenuOptions.setText("Options");

        jMenuItemLaserPolarization.setText("Laser polarization...");
        jMenuItemLaserPolarization.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItemLaserPolarizationActionPerformed(evt);
            }
        });
        jMenuOptions.add(jMenuItemLaserPolarization);
        jMenuOptions.add(jSeparator5);

        jMenuItemSize.setText("Graphical parameters...");
        jMenuItemSize.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItemSizeActionPerformed(evt);
            }
        });
        jMenuOptions.add(jMenuItemSize);

        jMenuItemNumerical.setText("Numerical parameters...");
        jMenuItemNumerical.setToolTipText("");
        jMenuItemNumerical.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItemNumericalActionPerformed(evt);
            }
        });
        jMenuOptions.add(jMenuItemNumerical);
        jMenuOptions.add(jSeparator2);

        jMenuItemOrderNumber.setText("Order number...");
        jMenuItemOrderNumber.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItemOrderNumberActionPerformed(evt);
            }
        });
        jMenuOptions.add(jMenuItemOrderNumber);

        jCheckBoxMenuItemSpread.setText("Velocity spread");
        jCheckBoxMenuItemSpread.setToolTipText("");
        jCheckBoxMenuItemSpread.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jCheckBoxMenuItemSpreadActionPerformed(evt);
            }
        });
        jMenuOptions.add(jCheckBoxMenuItemSpread);
        jMenuOptions.add(jSeparator4);

        jMenuSkin.setText("Look&Feel...");

        buttonGroupSkin.add(jRadioButtonMenuDefault);
        jRadioButtonMenuDefault.setText("Default");
        jRadioButtonMenuDefault.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                jRadioButtonMenuDefaultItemStateChanged(evt);
            }
        });
        jMenuSkin.add(jRadioButtonMenuDefault);

        buttonGroupSkin.add(jRadioButtonMenuSystem);
        jRadioButtonMenuSystem.setText("System");
        jRadioButtonMenuSystem.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                jRadioButtonMenuSystemItemStateChanged(evt);
            }
        });
        jMenuSkin.add(jRadioButtonMenuSystem);

        buttonGroupSkin.add(jRadioButtonMenuNimbus);
        jRadioButtonMenuNimbus.setSelected(true);
        jRadioButtonMenuNimbus.setText("Nimbus");
        jRadioButtonMenuNimbus.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                jRadioButtonMenuNimbusItemStateChanged(evt);
            }
        });
        jMenuSkin.add(jRadioButtonMenuNimbus);

        jMenuOptions.add(jMenuSkin);

        jMenuBarMain.add(jMenuOptions);

        jMenuHelp.setText("Help");

        HelpItem.setText("Contents...");
        HelpItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                HelpItemActionPerformed(evt);
            }
        });
        jMenuHelp.add(HelpItem);

        jMenuItemAbout.setText("About TSourceNX");
        jMenuItemAbout.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItemAboutActionPerformed(evt);
            }
        });
        jMenuHelp.add(jMenuItemAbout);

        jMenuBarMain.add(jMenuHelp);

        setJMenuBar(jMenuBarMain);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 730, Short.MAX_VALUE)
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 572, Short.MAX_VALUE)
                .addContainerGap(56, Short.MAX_VALUE))
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents
    private AbstractElectronBunch ebunch;
    private AbstractLaserPulse lpulse;
    private AbstractThomsonSource tsource, tsourcelinear, tsourceRayClone = null;

    /**
     * Parameters for the calculation boxes
     */
    class CalcBoxParam {

        public final static double MIN_DIF = 1e-10;
        public String[] valueUnitLabels, plotLabels;
        public String[] minValues, maxValues;
        public int selectedItemIndex = 10, selectedItemIndexClone = 10;
        public int numberOfItems;
        public double minValue = 0, maxValue = 100, minValueClone, maxValueClone;
        public double[] conversionValues;
        public AbstractThomsonSource tsourceclone;
        public boolean espread = false;
        public boolean working = false;
        public SwingWorker<Void, Void> worker;
        public String savetext;
        final private String[] keys;
        public LinearChartParam chartParam;
        public ChartPanel chartPanel = null;
        public JFreeChart chart = null;
        public double angle = 0, angleclone, energy = 46, energyclone;
        public int ordernumber = 1;
        private File file = null;

        /**
         * Constructor
         *
         * @param keys
         */
        public CalcBoxParam(String[] keys, List<Function<double[], Double>> trfunc) {
            super();
            this.keys = keys;
            this.chartParam = new LinearChartParam(trfunc);
        }

        /**
         * Initializing code
         */
        public void initialize(AbstractThomsonSource ts) {
            working = true;
            try {
                tsourceclone = (AbstractThomsonSource) ts.clone();
            } catch (CloneNotSupportedException ex) {
                Logger.getLogger(ThomsonJFrame.class.getName()).log(Level.SEVERE, null, ex);
            }
            tsourceclone.seteSpread(espread);
            minValueClone = minValue;
            maxValueClone = maxValue;
            angleclone = angle;
            energyclone = energy;
            selectedItemIndexClone = selectedItemIndex;
        }

        /**
         * Saving the results into the text file
         */
        public void save() {
            JFileChooser fo = new JFileChooser(file);
            fo.setDialogTitle(savetext);
            int ans = fo.showSaveDialog(null);
            if (ans == JFileChooser.APPROVE_OPTION) {
                file = fo.getSelectedFile();
                if (file.exists()) {
                    int n = JOptionPane.showConfirmDialog(null, "The file already exists. Overwrite?", "Warning",
                            JOptionPane.YES_NO_OPTION, JOptionPane.WARNING_MESSAGE);
                    if (n == JOptionPane.NO_OPTION) {
                        return;
                    }
                }
                try {
                    chartParam.save(file);
                } catch (IOException e) {
                    JOptionPane.showMessageDialog(null, "Error while writing to the file", "Error",
                            JOptionPane.ERROR_MESSAGE);
                }
            }
        }

        /**
         * Updating or creating Chart and ChartPanel
         */
        public void updateGraph(JPanel panel, String label) {
            if (chartPanel == null) {
                /*
                 * Creating chart and plot dataset
                 */
                chart = LinearChartParam.createLineChart(chartParam.createLineDataset(keys), plotLabels[selectedItemIndexClone], label);
                chart.getXYPlot().getRangeAxis().setRange(chartParam.getUMin(),
                        chartParam.getUMax() + MIN_DIF);
                chart.fireChartChanged();
                /**
                 * Creation of the ChartPanel
                 */
                chartPanel = new ChartPanel(chart,
                        (int) (panel.getWidth()), (int) panel.getHeight(), 0, 0,
                        (int) (10 * panel.getWidth()), (int) (10 * panel.getHeight()),
                        false, true, true, true, true, true);
                panel.setLayout(new BorderLayout(10, 10));
                panel.add(chartPanel, BorderLayout.CENTER);
                panel.revalidate();
                panel.repaint();
            } else {
                chart.getXYPlot().getDomainAxis().setRange(minValueClone, maxValueClone);
                chart.getXYPlot().getDomainAxis().setLabel(plotLabels[selectedItemIndexClone]);
                chart.getXYPlot().getRangeAxis().setRange(chartParam.getUMin(),
                        chartParam.getUMax() + MIN_DIF);
                chart.getXYPlot().getRangeAxis().setLabel(label);
                chart.fireChartChanged();
            }
            working = false;
        }

        /**
         * Canceling worker
         */
        public void cancel() {
            working = false;
            worker.cancel(true);
        }

        /**
         * Returning the array with graph keys
         */
        public String[] getKeys() {
            return keys;
        }
    }

    /**
     * Object for the color charts
     */
    class ColorChart {

        private final JFreeChart chart;
        private final JFreeChart colorbarchart;
        private final ChartPanel chartpanel;

        /**
         * Creating a color chart with a colorbar and and attaching it to a
         * JPanel
         *
         * @param data
         * @param xlabel
         * @param ylabel
         * @param colorBarlabel
         * @param jPanel
         * @param fraction
         */
        ColorChart(ChartParam data, String xlabel, String ylabel, String colorBarlabel, JPanel jPanel, double fraction, boolean slider) {
            this.chart = data.createChart(data.createDataset(slider), xlabel, ylabel);
            this.chartpanel = new ChartPanel(chart,
                    (int) (fraction * jPanel.getWidth()), (int) jPanel.getHeight(), 0, 0,
                    (int) (10 * jPanel.getWidth()), (int) (10 * jPanel.getHeight()),
                    false, true, true, true, true, true);
            this.colorbarchart = data.createColorBar(colorBarlabel);
            JPanel fluxcolorbarpanel = new ChartPanel(colorbarchart,
                    (int) ((1 - fraction) * jPanel.getWidth()), (int) jPanel.getHeight(), 0, 0,
                    (int) (10 * jPanel.getWidth()), (int) (10 * jPanel.getHeight()),
                    false, true, true, true, true, true);
            jPanel.setLayout(new BorderLayout(10, 10));
            jPanel.add(chartpanel, BorderLayout.CENTER);
            jPanel.add(fluxcolorbarpanel, BorderLayout.LINE_END);
            jPanel.revalidate();
            jPanel.repaint();
        }

        /**
         * Updating the chart and colorbar
         *
         * @param data
         */
        void fullupdate(ChartParam data) {
            chart.getXYPlot().getDomainAxis().setRangeAboutValue(data.getxoffset(), data.getxsize() * data.getxstep());
            chart.getXYPlot().getRangeAxis().setRangeAboutValue(data.getyoffset(), data.getysize() * data.getystep());
            PaintScale scale = new JetPaintScale(0, data.getumax());
            XYBlockRenderer renderer = ((XYBlockRenderer) chart.getXYPlot().getRenderer());
            renderer.setBlockHeight(data.getystep());
            renderer.setBlockWidth(data.getxstep());
            renderer.setPaintScale(scale);
            chart.fireChartChanged();

            colorbarchart.getXYPlot().getRangeAxis().setRange(0.0, data.getumax());
            renderer = ((XYBlockRenderer) colorbarchart.getXYPlot().getRenderer());
            renderer.setPaintScale(scale);
            renderer.setBlockHeight(data.getumax() / (data.getxsize() - 1));
            colorbarchart.fireChartChanged();
        }

        void update() {
            chart.fireChartChanged();
            colorbarchart.fireChartChanged();
        }

        JFreeChart getchart() {
            return chart;
        }

        ChartPanel getchartpanel() {
            return chartpanel;
        }

        JFreeChart getcolorbarchart() {
            return colorbarchart;
        }
    }

    /**
     * Main program parameters
     */
    /* The size of the graph in x or y direction */
    private int xsize, ysize;

    /* Step in x or y direction */
    private double xstep, ystep, estep, hoffset = 0;

    /* Number of rays exported for Shadow */
    private int numberOfRays = 1000;

    /* The order of magnitude desplayed in linear graphs */
    private int orderofmagnitude;
    private double normfactor;

    private String parametersIniFileName;

    private final ChartParam fluxdata, fluxcrossdata, xenergydata;
    private final LinearChartParam xenergycrossdata;
    private JFreeChart xenergycrosschart = null;
    private final TitledBorder xrayenergyborder;
    private final CalcBoxParam brilForm, brilFormNonLinear, gfForm, polForm, polFormNonLinear;
    private ColorChart fluxChart, fluxCrossChart, xEnergyChart;
    private boolean working = false, rayWorking = false;
    private SwingWorker<Void, Void> mainWorker, rayWorker;
    private Map<JTextField, String> oldStrings;
    JFormattedTextField rayNumberBox, rayXAngleRangeBox, rayYAngleRangeBox, rayMinEnergyBox, rayEnergyRangeBox,
            gfMonteCarloNumberBox, numericallPrecisionBox, shiftFactorBox, xSizeBox, ySizeBox, xRangeBox,
            yRangeBox, xEnergyRangeBox, threadsNumberBox, ksi1Box, ksi2Box, ksi3Box, orderNumberBox;
    List<Function<double[], Double>> fn;
    private File bFile = null, pFile = null;

    private void energyvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_energyvalueActionPerformed
        // TODO add your handling code here:
        ebunch.setGamma(TestValueWithMemory(0, 10000, energyvalue, "50", oldStrings) / 0.512);
    }//GEN-LAST:event_energyvalueActionPerformed

    private void phenergyvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_phenergyvalueActionPerformed
        // TODO add your handling code here:
        lpulse.setPhotonEnergy(TestValueWithMemory(0, 10, phenergyvalue, "1.204", oldStrings) * GaussianElectronBunch.E);
        ((NonLinearThomsonSource) tsource).setsIntensity();
    }//GEN-LAST:event_phenergyvalueActionPerformed

    private void energyvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_energyvalueFocusLost
        // TODO add your handling code here:
        ebunch.setGamma(TestValueWithMemory(0, 10000, energyvalue, "50", oldStrings) / 0.512);
    }//GEN-LAST:event_energyvalueFocusLost

    private void phenergyvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_phenergyvalueFocusLost
        // TODO add your handling code here:
        lpulse.setPhotonEnergy(TestValueWithMemory(0, 10, phenergyvalue, "1.204", oldStrings) * GaussianElectronBunch.E);
    }//GEN-LAST:event_phenergyvalueFocusLost

    private void chargevalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_chargevalueActionPerformed
        // TODO add your handling code here:
        ebunch.setNumber(TestValueWithMemory(0, 10, chargevalue, "0.2", oldStrings) / GaussianElectronBunch.E * 1e-9);
    }//GEN-LAST:event_chargevalueActionPerformed

    private void chargevalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_chargevalueFocusLost
        // TODO add your handling code here:
        ebunch.setNumber(TestValueWithMemory(0, 10, chargevalue, "0.2", oldStrings) / GaussianElectronBunch.E * 1e-9);
    }//GEN-LAST:event_chargevalueFocusLost

    private void spreadvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_spreadvalueActionPerformed
        // TODO add your handling code here:
        ebunch.setDelgamma((double) TestValueWithMemory(0.01, 2, spreadvalue, "0.25", oldStrings) / 200);
    }//GEN-LAST:event_spreadvalueActionPerformed

    private void spreadvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_spreadvalueFocusLost
        // TODO add your handling code here:
        ebunch.setDelgamma((double) TestValueWithMemory(0.01, 2, spreadvalue, "0.25", oldStrings) / 200);
    }//GEN-LAST:event_spreadvalueFocusLost

    private void elengthvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_elengthvalueActionPerformed
        // TODO add your handling code here:
        ebunch.setLength(TestValueWithMemory(0, 1000, elengthvalue, "10", oldStrings) * 3e-4 / 2);
    }//GEN-LAST:event_elengthvalueActionPerformed

    private void elengthvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_elengthvalueFocusLost
        // TODO add your handling code here:
        ebunch.setLength(TestValueWithMemory(0, 1000, elengthvalue, "10", oldStrings) * 3e-4 / 2);
    }//GEN-LAST:event_elengthvalueFocusLost

    private void pulseenergyvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_pulseenergyvalueActionPerformed
        // TODO add your handling code here:
        lpulse.setPulseEnergy(TestValueWithMemory(0, 10000, pulseenergyvalue, "100", oldStrings) * 1e-3);
    }//GEN-LAST:event_pulseenergyvalueActionPerformed

    private void pulseenergyvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_pulseenergyvalueFocusLost
        // TODO add your handling code here:
        lpulse.setPulseEnergy(TestValueWithMemory(0, 10000, pulseenergyvalue, "100", oldStrings) * 1e-3);
    }//GEN-LAST:event_pulseenergyvalueFocusLost

    private void pulselengthvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_pulselengthvalueActionPerformed
        // TODO add your handling code here:
        lpulse.setLength(TestValueWithMemory(0, 1000, pulselengthvalue, "10", oldStrings) * 3e-4 / 2);
    }//GEN-LAST:event_pulselengthvalueActionPerformed

    private void pulselengthvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_pulselengthvalueFocusLost
        // TODO add your handling code here:
        lpulse.setLength(TestValueWithMemory(0, 1000, pulselengthvalue, "10", oldStrings) * 3e-4 / 2);
    }//GEN-LAST:event_pulselengthvalueFocusLost

    private void startbuttonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_startbuttonActionPerformed
        // TODO add your handling code here:  
        if (working) {
            mainWorker.cancel(true);
            return;
        }
        MainProgressBar.setValue(0);
        MainProgressBar.setStringPainted(true);
        jSlider_pickup.setEnabled(false);
        startbutton.setText("Stop ");
        working = true;

        mainWorker = new SwingWorker<Void, Void>() {
            @Override
            protected Void doInBackground() throws Exception {
                try {
                    tsource.calculateLinearTotalFlux();
                    tsource.calculateGeometricFactor();
                    fluxdata.setup(xsize, ysize, xstep, ystep, 0, 0, (x) -> {
                        setStatusBar((int) x / 4);
                    });
                    xenergydata.setup(xsize, ysize, xstep, ystep, 0, 0, (x) -> {
                        setStatusBar((int) (x / 4 + 25));
                    });
                    fluxcrossdata.setup(xsize, ysize, estep, ystep, xenergydata.func(hoffset, 0.0)
                            * 1e3, 0.0, (x) -> {
                                setStatusBar((int) (x / 4 + 50));
                            });
                    xenergycrossdata.setup(xenergydata.getudata(), xenergydata.getSliderposition(),
                            false, ysize, ystep, -ystep * ysize / 2);
                    setStatusBar((int) 100);
                } catch (InterruptedException e) {

                }
                return null;
            }

            @Override
            protected void done() {
                double plotwidth = 0;
                //Checking if there are any errors in the worker thread
                try {
                    get();
                } catch (ExecutionException ex) {
                    Logger.getLogger(ThomsonJFrame.class.getName()).log(Level.SEVERE, null, ex);
                } catch (InterruptedException | CancellationException ex) {

                }
                if (!isCancelled() || fluxChart != null) {
                    //Creating or updating charts
                    if (fluxChart != null) {
                        fluxChart.fullupdate(fluxdata);
                    } else {
                        fluxChart = new ColorChart(fluxdata, "theta_x, mrad", "theta_y, mrad", "mrad\u207B\u00B2\u00B7s\u207B\u00B9\u00B710\u00B9\u2070",
                                jPanel_xflux_left, 0.75, true);
                    }
                    if (fluxCrossChart != null) {
                        fluxCrossChart.fullupdate(fluxcrossdata);
                    } else {
                        fluxCrossChart = new ColorChart(fluxcrossdata, "X-ray energy, eV", "theta_y, mrad",
                                "mrad\u207B\u00B2\u00B7s\u207B\u00B9\u00B70.1%\u00B710\u00B9\u2070", jPanel_xflux_right, 0.73, false);
                    }
                    if (xEnergyChart != null) {
                        xEnergyChart.fullupdate(xenergydata);
                    } else {
                        xEnergyChart = new ColorChart(xenergydata, "theta_x, mrad", "theta_y, mrad", "kev", jPanel_xenergy_left, 0.75, true);
                    }
                    if (xenergycrosschart != null) {
                        xenergycrosschart.getXYPlot().getRangeAxis().setRange(xenergycrossdata.getUMin(), xenergycrossdata.getUMax());
                        xenergycrosschart.getXYPlot().getDomainAxis().setRangeAboutValue(xenergycrossdata.getOffset()
                                + xenergycrossdata.getSize() * xenergycrossdata.getStep() / 2, xenergycrossdata.getSize() * xenergycrossdata.getStep());
                        xenergycrosschart.fireChartChanged();
                    } else {
                        xenergycrosschart = LinearChartParam.createLineChart(xenergycrossdata
                                .createLineDataset(new String[]{"Energy cross section"}), "theta_y, mrad", "Energy, keV");
                        ChartPanel chartpanel = new ChartPanel(xenergycrosschart,
                                (int) (jPanel_xenergy_right.getWidth()), (int) jPanel_xenergy_right.getHeight(), 0, 0,
                                (int) (10 * jPanel_xenergy_right.getWidth()), (int) (10 * jPanel_xenergy_right.getHeight()),
                                false, true, true, true, true, true);
                        jPanel_xenergy_right.setLayout(new BorderLayout(10, 10));
                        jPanel_xenergy_right.add(chartpanel, BorderLayout.CENTER);
                        jPanel_xenergy_right.revalidate();
                        jPanel_xenergy_right.repaint();
                    }
                    fluxChart.getchartpanel().revalidate();
                    fluxChart.getchartpanel().repaint();
                    plotwidth = xEnergyChart.getchartpanel().getChartRenderingInfo().
                            getPlotInfo().getDataArea().getWidth();
                    xrayenergyborder.setTitle("X-ray photon energy" + ". Order: "
                            + (new DecimalFormat("##")).format(((NonLinearThomsonSource) tsource).getOrdernumber()) + ", Max: "
                            + (new DecimalFormat("########.##")).format(xenergydata.getumax()) + " keV");
                    jPanel_xenergy.repaint();
                    totalFluxLabel.setText("Total flux: "
                            + (new DecimalFormat("##.#######")).format(tsource.getLinearTotalFlux() * tsource.getGeometricFactor() * 1e-12)
                            + "\u00B710\u00B9\u00B2\u00B7ph\u00B7s\u207B\u00B9");
                    totalFluxAngleLabel.setText("Within angle: "
                            + (new DecimalFormat("##.#######"))
                                    .format(tsource.calculateAngleLinearTotalFlux(Math.max(xsize * xstep,
                                            ysize * ystep) * 1e-3 / 2) * 1e-12)
                            + "\u00B710\u00B9\u00B2\u00B7ph\u00B7s\u207B\u00B9");
                }
                startbutton.setText("Start");
                jSlider_pickup.setEnabled(true);
                if (plotwidth != 0) {
                    jSlider_pickup.setPreferredSize(new Dimension((int) plotwidth,
                            (int) jSlider_pickup.getSize().getHeight()));
                }
                working = false;
            }

            /**
             * Updating progress bar
             *
             * @param status
             */
            public void setStatusBar(final int status) {
                SwingUtilities.invokeLater(() -> MainProgressBar.setValue(status));
            }
        };
        mainWorker.execute();
    }//GEN-LAST:event_startbuttonActionPerformed

    private void eemitxvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_eemitxvalueFocusLost
        // TODO add your handling code here:
        ebunch.setEpsx(TestValueWithMemory(0.1, 100, eemitxvalue, "1", oldStrings) * 1e-6);
    }//GEN-LAST:event_eemitxvalueFocusLost

    private void eemitxvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_eemitxvalueActionPerformed
        // TODO add your handling code here:
        ebunch.setEpsx(TestValueWithMemory(0.1, 100, eemitxvalue, "1", oldStrings) * 1e-6);
    }//GEN-LAST:event_eemitxvalueActionPerformed

    private void ebetaxvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_ebetaxvalueFocusLost
        // TODO add your handling code here:
        ebunch.setBetax(TestValueWithMemory(0.1, 100, ebetaxvalue, "20", oldStrings) * 1e-3);
    }//GEN-LAST:event_ebetaxvalueFocusLost

    private void ebetaxvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ebetaxvalueActionPerformed
        // TODO add your handling code here:
        ebunch.setBetax(TestValueWithMemory(0.1, 100, ebetaxvalue, "20", oldStrings) * 1e-3);
    }//GEN-LAST:event_ebetaxvalueActionPerformed

    private void jSlider_pickupStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_jSlider_pickupStateChanged
        // Setting the slider position
        JSlider source = (JSlider) evt.getSource();
        fluxdata.setSliderposition(source.getValue());
        xenergydata.setSliderposition(source.getValue());

        //Updating charts when the slider is moved
        if (!source.getValueIsAdjusting()) {
            if (working) {
                return;
            }
            MainProgressBar.setValue(0);
            MainProgressBar.setStringPainted(true);
            startbutton.setText("Stop");
            working = true;

            hoffset = xsize * xstep * (source.getValue() - 50) / 100;
            mainWorker = new SwingWorker<Void, Void>() {
                @Override
                protected Void doInBackground() throws Exception {
                    try {
                        fluxcrossdata.setup(xsize, ysize, estep, ystep, xenergydata.func(hoffset, 0.0) * 1e3,
                                0.0, (x) -> {
                                    setStatusBar((int) x / 2);
                                });
                        xenergycrossdata.setup(xenergydata.getudata(), xenergydata.getSliderposition(),
                                false, ysize, ystep, -ystep * ysize / 2);
                        setStatusBar((int) 100);
                    } catch (InterruptedException e) {

                    }
                    return null;
                }

                @Override
                protected void done() {
                    //Checking if there are any errors in the worker thread
                    try {
                        get();
                    } catch (ExecutionException ex) {
                        Logger.getLogger(ThomsonJFrame.class.getName()).log(Level.SEVERE, null, ex);
                    } catch (InterruptedException | CancellationException ex) {

                    }
                    if (fluxChart != null) {
                        fluxChart.update();
                    }
                    //Updating charts
                    if (fluxCrossChart != null) {
                        fluxCrossChart.fullupdate(fluxcrossdata);
                    }
                    if (xEnergyChart != null) {
                        xEnergyChart.update();
                    }
                    if (xenergycrosschart != null) {
                        xenergycrosschart.getXYPlot().getRangeAxis().setRange(xenergycrossdata.getUMin(), xenergycrossdata.getUMax());
                        xenergycrosschart.fireChartChanged();
                    }
                    startbutton.setText("Start");
                    working = false;
                }

                /**
                 * Updating progress bar
                 *
                 * @param status
                 */
                public void setStatusBar(final int status) {
                    SwingUtilities.invokeLater(() -> MainProgressBar.setValue(status));
                }
            };
            mainWorker.execute();
        }
    }//GEN-LAST:event_jSlider_pickupStateChanged

    private void pulserelvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_pulserelvalueFocusLost
        // TODO add your handling code here:
        lpulse.setRlength(TestValueWithMemory(0.01, 100, pulserelvalue, "0.35", oldStrings) * 1e-3);
    }//GEN-LAST:event_pulserelvalueFocusLost

    private void pulserelvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_pulserelvalueActionPerformed
        // TODO add your handling code here:
        lpulse.setRlength(TestValueWithMemory(0.01, 100, pulserelvalue, "0.35", oldStrings) * 1e-3);
    }//GEN-LAST:event_pulserelvalueActionPerformed

    private void pulsefreqvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_pulsefreqvalueFocusLost
        // TODO add your handling code here:
        lpulse.setFq(TestValueWithMemory(0, 1000, pulsefreqvalue, "1000", oldStrings));
    }//GEN-LAST:event_pulsefreqvalueFocusLost

    private void pulsefreqvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_pulsefreqvalueActionPerformed
        // TODO add your handling code here:
        lpulse.setFq(TestValueWithMemory(0, 1000, pulsefreqvalue, "1000", oldStrings));
    }//GEN-LAST:event_pulsefreqvalueActionPerformed

    private void pulsedelayvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_pulsedelayvalueFocusLost
        // TODO add your handling code here:
        lpulse.setDelay(TestValueWithMemory(0, 1000, pulsedelayvalue, "0", oldStrings) * 3e-4);
    }//GEN-LAST:event_pulsedelayvalueFocusLost

    private void pulsedelayvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_pulsedelayvalueActionPerformed
        // TODO add your handling code here:
        lpulse.setDelay(TestValueWithMemory(0, 1000, pulsedelayvalue, "0", oldStrings) * 3e-4);
    }//GEN-LAST:event_pulsedelayvalueActionPerformed

    private void eshiftxvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_eshiftxvalueFocusLost
        // TODO add your handling code here:
        Double value = TestValueWithMemory(0, 1000, eshiftxvalue, "0", oldStrings) * 1e-3;
        Vector dir = ebunch.getShift();
        dir.set(0, value);
        ebunch.setShift(dir);
    }//GEN-LAST:event_eshiftxvalueFocusLost

    private void eshiftxvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_eshiftxvalueActionPerformed
        // TODO add your handling code here:
        Double value = TestValueWithMemory(0, 1000, eshiftxvalue, "0", oldStrings) * 1e-3;
        Vector dir = ebunch.getShift();
        dir.set(0, value);
        ebunch.setShift(dir);
    }//GEN-LAST:event_eshiftxvalueActionPerformed

    private void eshiftyvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_eshiftyvalueFocusLost
        // TODO add your handling code here:
        Double value = TestValueWithMemory(0, 1000, eshiftyvalue, "0", oldStrings) * 1e-3;
        Vector dir = ebunch.getShift();
        dir.set(1, value);
        ebunch.setShift(dir);
    }//GEN-LAST:event_eshiftyvalueFocusLost

    private void eshiftyvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_eshiftyvalueActionPerformed
        // TODO add your handling code here:
        Double value = TestValueWithMemory(0, 1000, eshiftyvalue, "0", oldStrings) * 1e-3;
        Vector dir = ebunch.getShift();
        dir.set(1, value);
        ebunch.setShift(dir);
    }//GEN-LAST:event_eshiftyvalueActionPerformed

    private void eshiftzvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_eshiftzvalueFocusLost
        // TODO add your handling code here:
        Double value = TestValueWithMemory(0, 1000, eshiftzvalue, "0", oldStrings) * 1e-3;
        Vector dir = ebunch.getShift();
        dir.set(2, value);
        ebunch.setShift(dir);
    }//GEN-LAST:event_eshiftzvalueFocusLost

    private void eshiftzvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_eshiftzvalueActionPerformed
        // TODO add your handling code here:
        Double value = TestValueWithMemory(0, 1000, eshiftzvalue, "0", oldStrings) * 1e-3;
        Vector dir = new BasicVector(new double[3]);
        dir.set(2, value);
        ebunch.setShift(dir);
    }//GEN-LAST:event_eshiftzvalueActionPerformed

    private void pulseanglevalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_pulseanglevalueFocusLost
        // TODO add your handling code here:
        Double value = TestValueWithMemory(0, 100, pulseanglevalue, "52", oldStrings) * 1e-3;
        Vector dir = new BasicVector(new double[3]);
        dir.set(2, Math.cos(value));
        dir.set(1, Math.sin(value));
        lpulse.setDirection(dir);
    }//GEN-LAST:event_pulseanglevalueFocusLost

    private void pulseanglevalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_pulseanglevalueActionPerformed
        // TODO add your handling code here:
        Double value = TestValueWithMemory(0, 100, pulseanglevalue, "52", oldStrings) * 1e-3;
        Vector dir = new BasicVector(new double[3]);
        dir.set(2, Math.cos(value));
        dir.set(1, Math.sin(value));
        lpulse.setDirection(dir);
    }//GEN-LAST:event_pulseanglevalueActionPerformed

    private void jMenuItemBrillianceActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItemBrillianceActionPerformed
        // TODO add your handling code here:  
        brillianceCalc.setVisible(true);
    }//GEN-LAST:event_jMenuItemBrillianceActionPerformed

    private void BrillianceCalcStartActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_BrillianceCalcStartActionPerformed
        // Brilliance box calculations
        // Checking if already running
        if (brilForm.working) {
            brilForm.cancel();
            BrillianceCalcStart.setText("Calculate");
            BrillianceCalcSave.setEnabled(true);
            return;
        }
        BrilProgressBar.setValue(0);
        BrilProgressBar.setStringPainted(true);
        BrillianceCalcStart.setText("Terminate");
        BrillianceCalcSave.setEnabled(false);
        brilForm.initialize(tsourcelinear);

        /**
         * Calculating data array. Using SwingWorker class
         */
        brilForm.worker = new SwingWorker<Void, Void>() {
            @Override
            protected Void doInBackground() throws Exception {
                double step = (brilForm.maxValueClone - brilForm.minValueClone) / (xsize - 1);
                double offset = brilForm.minValueClone;
                List<Function<Double, Double>> func = new ArrayList<>();
                switch (brilForm.selectedItemIndexClone) {
                    case 0:
                        func.add(xp -> {
                            double ang = brilForm.angleclone * 1e-3;
                            double e = brilForm.energyclone * GaussianElectronBunch.E * 1e3;
                            double x = xp * brilForm.conversionValues[brilForm.selectedItemIndexClone];
                            Vector dir = new BasicVector(new double[3]);
                            dir.set(2, Math.cos(x));
                            dir.set(1, Math.sin(x));
                            brilForm.tsourceclone.getLaserPulse().setDirection(dir);
                            brilForm.tsourceclone.calculateLinearTotalFlux();
                            setStatusBar((xp - offset) / step / (xsize - 1));
                            try {
                                return brilForm.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0, 0, 0}),
                                        new BasicVector(new double[]{Math.sin(ang), 0, Math.cos(ang)}), new BasicVector(new double[]{0, 0, 1}),
                                        e) * normfactor;
                            } catch (InterruptedException ex) {
                                Thread.currentThread().interrupt();
                                return 0.0;
                            }
                        });
                        break;
                    case 1:
                        func.add(xp -> {
                            double ang = brilForm.angleclone * 1e-3;
                            double e = brilForm.energyclone * GaussianElectronBunch.E * 1e3;
                            double x = xp * brilForm.conversionValues[brilForm.selectedItemIndexClone];
                            brilForm.tsourceclone.getLaserPulse().setDelay(x);
                            setStatusBar((xp - offset) / step / (xsize - 1));
                            try {
                                return brilForm.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0, 0, 0}),
                                        new BasicVector(new double[]{Math.sin(ang), 0, Math.cos(ang)}), new BasicVector(new double[]{0, 0, 1}),
                                        e) * normfactor;
                            } catch (InterruptedException ex) {
                                Thread.currentThread().interrupt();
                                return 0.0;
                            }
                        });
                        break;
                    case 2:
                        func.add(xp -> {
                            double ang = brilForm.angleclone * 1e-3;
                            double e = brilForm.energyclone * GaussianElectronBunch.E * 1e3;
                            double x = xp * brilForm.conversionValues[brilForm.selectedItemIndexClone];
                            Vector sht = new BasicVector(new double[3]);
                            sht.set(2, x);
                            brilForm.tsourceclone.getElectronBunch().setShift(sht);
                            brilForm.tsourceclone.calculateLinearTotalFlux();
                            setStatusBar((xp - offset) / step / (xsize - 1));
                            try {
                                return brilForm.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0, 0, 0}),
                                        new BasicVector(new double[]{Math.sin(ang), 0, Math.cos(ang)}), new BasicVector(new double[]{0, 0, 1}),
                                        e) * normfactor;
                            } catch (InterruptedException ex) {
                                Thread.currentThread().interrupt();
                                return 0.0;
                            }
                        });
                        break;
                    case 3:
                        func.add(xp -> {
                            double ang = brilForm.angleclone * 1e-3;
                            double e = brilForm.energyclone * GaussianElectronBunch.E * 1e3;
                            double x = xp * brilForm.conversionValues[brilForm.selectedItemIndexClone];
                            brilForm.tsourceclone.getElectronBunch().setBetax(x);
                            brilForm.tsourceclone.getElectronBunch().setBetay(x);
                            brilForm.tsourceclone.calculateLinearTotalFlux();
                            setStatusBar((xp - offset) / step / (xsize - 1));
                            try {
                                return brilForm.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0, 0, 0}),
                                        new BasicVector(new double[]{Math.sin(ang), 0, Math.cos(ang)}), new BasicVector(new double[]{0, 0, 1}),
                                        e) * normfactor;
                            } catch (InterruptedException ex) {
                                Thread.currentThread().interrupt();
                                return 0.0;
                            }
                        });
                        break;
                    case 4:
                        func.add(xp -> {
                            double ang = brilForm.angleclone * 1e-3;
                            double e = brilForm.energyclone * GaussianElectronBunch.E * 1e3;
                            double x = xp * brilForm.conversionValues[brilForm.selectedItemIndexClone];
                            brilForm.tsourceclone.getElectronBunch().setEpsx(x);
                            brilForm.tsourceclone.getElectronBunch().setEpsy(x);
                            brilForm.tsourceclone.calculateLinearTotalFlux();
                            setStatusBar((xp - offset) / step / (xsize - 1));
                            try {
                                return brilForm.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0, 0, 0}),
                                        new BasicVector(new double[]{Math.sin(ang), 0, Math.cos(ang)}), new BasicVector(new double[]{0, 0, 1}),
                                        e) * normfactor;
                            } catch (InterruptedException ex) {
                                Thread.currentThread().interrupt();
                                return 0.0;
                            }
                        });
                        break;
                    case 5:
                        func.add(xp -> {
                            double ang = brilForm.angleclone * 1e-3;
                            double e = brilForm.energyclone * GaussianElectronBunch.E * 1e3;
                            double x = xp * brilForm.conversionValues[brilForm.selectedItemIndexClone];
                            brilForm.tsourceclone.getElectronBunch().setEpsx(x);
                            brilForm.tsourceclone.calculateLinearTotalFlux();
                            setStatusBar((xp - offset) / step / (xsize - 1));
                            try {
                                return brilForm.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0, 0, 0}),
                                        new BasicVector(new double[]{Math.sin(ang), 0, Math.cos(ang)}), new BasicVector(new double[]{0, 0, 1}),
                                        e) * normfactor;
                            } catch (InterruptedException ex) {
                                Thread.currentThread().interrupt();
                                return 0.0;
                            }
                        });
                        break;
                    case 6:
                        func.add(xp -> {
                            double ang = brilForm.angleclone * 1e-3;
                            double e = brilForm.energyclone * GaussianElectronBunch.E * 1e3;
                            double x = xp * brilForm.conversionValues[brilForm.selectedItemIndexClone];
                            brilForm.tsourceclone.getElectronBunch().setEpsy(x);
                            brilForm.tsourceclone.calculateLinearTotalFlux();
                            setStatusBar((xp - offset) / step / (xsize - 1));
                            try {
                                return brilForm.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0, 0, 0}),
                                        new BasicVector(new double[]{Math.sin(ang), 0, Math.cos(ang)}), new BasicVector(new double[]{0, 0, 1}),
                                        e) * normfactor;
                            } catch (InterruptedException ex) {
                                Thread.currentThread().interrupt();
                                return 0.0;
                            }
                        });
                        break;
                    case 7:
                        func.add(xp -> {
                            double ang = brilForm.angleclone * 1e-3;
                            double e = brilForm.energyclone * GaussianElectronBunch.E * 1e3;
                            double x = xp * brilForm.conversionValues[brilForm.selectedItemIndexClone];
                            brilForm.tsourceclone.getLaserPulse().setRlength(x);
                            brilForm.tsourceclone.calculateLinearTotalFlux();
                            setStatusBar((xp - offset) / step / (xsize - 1));
                            try {
                                return brilForm.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0, 0, 0}),
                                        new BasicVector(new double[]{Math.sin(ang), 0, Math.cos(ang)}), new BasicVector(new double[]{0, 0, 1}),
                                        e) * normfactor;
                            } catch (InterruptedException ex) {
                                Thread.currentThread().interrupt();
                                return 0.0;
                            }
                        });
                        break;
                    case 8:
                        func.add(xp -> {
                            double ang = brilForm.angleclone * 1e-3;
                            double e = brilForm.energyclone * GaussianElectronBunch.E * 1e3;
                            double x = xp * brilForm.conversionValues[brilForm.selectedItemIndexClone];
                            brilForm.tsourceclone.getLaserPulse().setWidth(x);
                            brilForm.tsourceclone.getElectronBunch().setxWidth(x);
                            brilForm.tsourceclone.getElectronBunch().setyWidth(x);
                            brilForm.tsourceclone.calculateLinearTotalFlux();
                            setStatusBar((xp - offset) / step / (xsize - 1));
                            try {
                                return brilForm.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0, 0, 0}),
                                        new BasicVector(new double[]{Math.sin(ang), 0, Math.cos(ang)}), new BasicVector(new double[]{0, 0, 1}),
                                        e) * normfactor;
                            } catch (InterruptedException ex) {
                                Thread.currentThread().interrupt();
                                return 0.0;
                            }
                        });
                        break;
                    case 9:
                        func.add(xp -> {
                            double ang = brilForm.angleclone * 1e-3;
                            double e = brilForm.energyclone * GaussianElectronBunch.E * 1e3;
                            double x = xp * brilForm.conversionValues[brilForm.selectedItemIndexClone];
                            brilForm.tsourceclone.getElectronBunch().setDelgamma(x);
                            brilForm.tsourceclone.calculateLinearTotalFlux();
                            setStatusBar((xp - offset) / step / (xsize - 1));
                            try {
                                return brilForm.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0, 0, 0}),
                                        new BasicVector(new double[]{Math.sin(ang), 0, Math.cos(ang)}), new BasicVector(new double[]{0, 0, 1}),
                                        e) * normfactor;
                            } catch (InterruptedException ex) {
                                Thread.currentThread().interrupt();
                                return 0.0;
                            }
                        });
                        break;
                    case 10:
                        func.add(xp -> {
                            double ang = brilForm.angleclone * 1e-3;
                            double e = xp * brilForm.conversionValues[brilForm.selectedItemIndexClone];
                            brilForm.tsourceclone.calculateLinearTotalFlux();
                            setStatusBar((xp - offset) / step / (xsize - 1));
                            try {
                                return brilForm.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0.0, 0.0, 0.0}),
                                        new BasicVector(new double[]{Math.sin(ang), 0, Math.cos(ang)}), new BasicVector(new double[]{0.0, 0.0, 1.0}),
                                        e) * normfactor;
                            } catch (InterruptedException ex) {
                                Thread.currentThread().interrupt();
                                return 0.0;
                            }
                        });
                        break;
                    case 11:
                        func.add(xp -> {
                            double ang = xp * brilForm.conversionValues[brilForm.selectedItemIndexClone];
                            double e = brilForm.energyclone * GaussianElectronBunch.E * 1e3;
                            brilForm.tsourceclone.calculateLinearTotalFlux();
                            setStatusBar((xp - offset) / step / (xsize - 1));
                            try {
                                return brilForm.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0, 0, 0}),
                                        new BasicVector(new double[]{Math.sin(ang), 0, Math.cos(ang)}), new BasicVector(new double[]{0, 0, 1}),
                                        e) * normfactor;
                            } catch (InterruptedException ex) {
                                Thread.currentThread().interrupt();
                                return 0.0;
                            }
                        });
                        break;
                }
                brilForm.chartParam.setup(func, xsize, step, offset);
                return null;
            }

            @Override
            protected void done() {
                try {
                    get();
                } catch (ExecutionException ex) {
                    Logger.getLogger(ThomsonJFrame.class.getName()).log(Level.SEVERE, null, ex);
                } catch (InterruptedException | CancellationException ex) {

                }
                brilForm.updateGraph(BrillianceCalcGraph,
                        "mm\u207B\u00B2\u00B7mrad\u207B\u00B2\u00B7s\u207B\u00B9\u00B70.1%\u00B710\u00B9\u2070");
                BrillianceCalcStart.setText("Calculate");
                BrillianceCalcSave.setEnabled(true);
            }

            /**
             * Updating progress bar
             *
             * @param status
             */
            public void setStatusBar(final double status) {
                SwingUtilities.invokeLater(() -> BrilProgressBar.setValue((int) Math.round(100 * status)));
            }
        };
        brilForm.worker.execute();
    }//GEN-LAST:event_BrillianceCalcStartActionPerformed

    private void BrillianceCalcSaveActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_BrillianceCalcSaveActionPerformed
        // Saving the brilliance plot data
        if (brilForm.chartPanel != null) {
            brilForm.save();
        }
    }//GEN-LAST:event_BrillianceCalcSaveActionPerformed

    private void BrillianceCalcBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_BrillianceCalcBoxActionPerformed
        // Update boxes, labels and other items when parameters selection is performed
        int sInd = BrillianceCalcBox.getSelectedIndex();
        brilForm.selectedItemIndex = sInd;
        Brilminvalueunitlabel.setText(brilForm.valueUnitLabels[sInd]);
        Brilmaxvalueunitlabel.setText(brilForm.valueUnitLabels[sInd]);
        Brilminvalue.setText(brilForm.minValues[sInd]);
        Brilmaxvalue.setText(brilForm.maxValues[sInd]);
        brilForm.minValue = Float.parseFloat(Brilminvalue.getText());
        brilForm.maxValue = Float.parseFloat(Brilmaxvalue.getText());
        switch (sInd) {
            case 10:
                angleValue.setEnabled(true);
                energyValue.setEnabled(false);
                break;
            case 11:
                angleValue.setEnabled(false);
                energyValue.setEnabled(true);
                break;
            default:
                angleValue.setEnabled(true);
                energyValue.setEnabled(true);
        }
    }//GEN-LAST:event_BrillianceCalcBoxActionPerformed

    private void BrilmaxvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_BrilmaxvalueActionPerformed
        // TODO add your handling code here:
        brilForm.maxValue = TestValueWithMemory(0, 10000, Brilmaxvalue, brilForm.maxValues[brilForm.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_BrilmaxvalueActionPerformed

    private void BrilminvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_BrilminvalueActionPerformed
        // TODO add your handling code here:
        brilForm.minValue = TestValueWithMemory(0, 10000, Brilminvalue, brilForm.minValues[brilForm.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_BrilminvalueActionPerformed

    private void BrilminvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_BrilminvalueFocusLost
        // TODO add your handling code here:
        brilForm.minValue = TestValueWithMemory(0, 10000, Brilminvalue, brilForm.minValues[brilForm.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_BrilminvalueFocusLost

    private void BrilmaxvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_BrilmaxvalueFocusLost
        // TODO add your handling code here:
        brilForm.maxValue = TestValueWithMemory(0, 10000, Brilmaxvalue, brilForm.maxValues[brilForm.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_BrilmaxvalueFocusLost

    private void jMenuItemExitActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItemExitActionPerformed
        // TODO add your handling code here:
        System.exit(0);
    }//GEN-LAST:event_jMenuItemExitActionPerformed

    private void jMenuItemGeometricFactorActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItemGeometricFactorActionPerformed
        // TODO add your handling code here:
        gfCalc.setVisible(true);
    }//GEN-LAST:event_jMenuItemGeometricFactorActionPerformed

    private void jMenuItemAboutActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItemAboutActionPerformed
        // Extracting the build date from the MANIFEST.MF file
        Package pk = Package.getPackage("thomsonsource");
        Date dt = new Date();
        DateFormat dtf = new SimpleDateFormat("yyyy-MM-dd");
        try {
            Enumeration<URL> mfs = this.getClass().getClassLoader().getResources("META-INF/MANIFEST.MF");
            while (mfs.hasMoreElements()) {
                Manifest mft = new Manifest(mfs.nextElement().openStream());
                if (mft.getMainAttributes().getValue("Built-Date") != null) {
                    dt = dtf.parse(mft.getMainAttributes().getValue("Built-Date"));
                    break;
                }
            }
        } catch (IOException | ParseException ex) {
            Logger.getLogger(ThomsonJFrame.class.getName()).log(Level.SEVERE, null, ex);
        }
        // Extracting vendor and version from the MANIFEST.MF file and showing up About popup window
        JOptionPane.showMessageDialog(null,
                "<html>" + pk.getImplementationTitle()
                + "<br>Version: " + pk.getImplementationVersion()
                + "<br>Build date: " + DateFormat.getDateInstance(DateFormat.LONG).format(dt)
                + "<br>Author: " + pk.getImplementationVendor()
                + "</html>",
                "About TSourceNX", JOptionPane.INFORMATION_MESSAGE);
    }//GEN-LAST:event_jMenuItemAboutActionPerformed

    private void jMenuItemSaveParamActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItemSaveParamActionPerformed
        // Saving the LEXG parameters into a file
        JFileChooser fo = new JFileChooser(pFile);
        fo.setDialogTitle("Choose a file to save Thompson source parameters");
        fo.setFileFilter(new FileNameExtensionFilter("ini file", "ini"));
        int ans = fo.showSaveDialog(this);
        if (ans == JFileChooser.APPROVE_OPTION) {
            pFile = fo.getSelectedFile();
            if (pFile.exists()) {
                if (!pFile.isFile()) {
                    //If not a file then do nothing
                    return;
                }
                int n = JOptionPane.showConfirmDialog(null, "The file already exists. Overwrite?", "Warning",
                        JOptionPane.YES_NO_OPTION, JOptionPane.WARNING_MESSAGE);
                if (n == JOptionPane.NO_OPTION) {
                    return;
                }
            } else {
                String extension = "";
                int ind = pFile.getName().lastIndexOf('.');
                if (ind > -1) {
                    //If the file has extension then get it
                    extension = pFile.getName().substring(ind + 1);
                }
                //If the extension is not 'ini' then add the 'ini' extension to the file name
                if (!extension.equals("ini")) {
                    try {

                        pFile = (ind != -1) ? new File(pFile.getCanonicalPath().substring(0, ind + 1) + ".ini")
                                : new File(pFile.getCanonicalPath() + ".ini");
                    } catch (IOException ex) {
                        //Do nothing
                    }
                }
            }
            //Saving Thomson source parameters into the file
            try {
                tsource.savePropperties(pFile);
            } catch (IOException e) {
                JOptionPane.showMessageDialog(null, "Error while writing to the file", "Error",
                        JOptionPane.ERROR_MESSAGE);
            }
        }
    }//GEN-LAST:event_jMenuItemSaveParamActionPerformed

    private void jMenuItemLoadParamActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItemLoadParamActionPerformed
        // Loading LEXG parameters from a file
        JFileChooser fo = new JFileChooser(pFile);
        fo.setDialogTitle("Choose a file to load Thompson source parameters from");
        fo.setFileFilter(new FileNameExtensionFilter("ini file", "ini"));
        int ans = fo.showOpenDialog(this);
        if (ans == JFileChooser.APPROVE_OPTION) {
            pFile = fo.getSelectedFile();
            loadParameters(pFile);
        }
    }//GEN-LAST:event_jMenuItemLoadParamActionPerformed

    private void jMenuItemSourceActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItemSourceActionPerformed
        // TODO add your handling code here:
        if (rayWorking) {
            rayProgressFrame.setVisible(true);
            return;
        }
        jRayProgressBar.setStringPainted(true);
        jRayProgressBar.setValue(0);
        jRayStopButton.setEnabled(true);
        try {
            tsourceRayClone = (AbstractThomsonSource) tsourcelinear.clone();
        } catch (CloneNotSupportedException ex) {

        }
        final int rayNumber = tsourceRayClone.getThreadNumber() * (numberOfRays / tsourceRayClone.getThreadNumber());
        tsourceRayClone.calculateLinearTotalFlux();
        rayWorking = true;
        rayWorker = new SwingWorker<Void, Void>() {
            @Override
            protected Void doInBackground() throws Exception {
                // Open a file for rays
                try (ShadowFiles shadowFile = new ShadowFiles(true, true, AbstractThomsonSource.NUMBER_OF_COLUMNS, rayNumber, bFile)) {
                    bFile = shadowFile.getFile();
                    tsourceRayClone.writeRays(shadowFile, numberOfRays, this::setStatusBar);
                }
                return null;
            }

            @Override
            protected void done() {
                jRayStopButton.setEnabled(false);
                jLabelPartialFlux.setText("Flux: " + tsourceRayClone.getPartialFlux() * 1e-8
                        + " 10\u2078 s\u207B\u00B9");
                try {
                    get();
                } catch (InterruptedException | CancellationException e) {
                    /* If the thread is interrupted or cancelled */

                } catch (ExecutionException e) {
                    /* If an exception is thrown during execution */
                    if (e.getCause() instanceof IOException) {
                        JOptionPane.showMessageDialog(null, "Error while writing to the file", "Error",
                                JOptionPane.ERROR_MESSAGE);
                    } else if (e.getCause() instanceof IllegalFormatException) {
                        JOptionPane.showMessageDialog(null, "Format error while writing to the file", "Error",
                                JOptionPane.ERROR_MESSAGE);
                    } else if (e.getCause() instanceof ShadowFiles.FileNotOpenedException) {

                    } else {

                    }
                }
                rayWorking = false;
            }

            /**
             * Updating progress bar
             *
             * @param status
             */
            public void setStatusBar(final int status) {
                SwingUtilities.invokeLater(() -> {
                    if (status != jRayProgressBar.getValue()) {
                        jRayProgressBar.setValue(status);
                    }
                });
            }
        };
        rayProgressFrame.setVisible(true);
        rayWorker.execute();
    }//GEN-LAST:event_jMenuItemSourceActionPerformed

    private void jMenuItemSizeActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItemSizeActionPerformed
        // Obtaining display parameters
        Object[] message = {
            "x-size:", xSizeBox,
            "y-size:", ySizeBox,
            "x-range (mrad):", xRangeBox,
            "y-range (mrad):", yRangeBox,
            "xenergy-range (eV):", xEnergyRangeBox
        };
        int option = JOptionPane.showConfirmDialog(null, message, "Graphical parameters", JOptionPane.OK_CANCEL_OPTION);
        if (option == JOptionPane.OK_OPTION) {
            xsize = (int) xSizeBox.getValue();
            ysize = (int) ySizeBox.getValue();
            xstep = (double) xRangeBox.getValue() / xsize;
            ystep = (double) yRangeBox.getValue() / ysize;
            estep = (double) xEnergyRangeBox.getValue() / xsize;
        }
    }//GEN-LAST:event_jMenuItemSizeActionPerformed

    private void GFCalcBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_GFCalcBoxActionPerformed
        // TODO add your handling code here:
        int sInd = GFCalcBox.getSelectedIndex();
        gfForm.selectedItemIndex = sInd;
        GFminvalueunitlabel.setText(gfForm.valueUnitLabels[sInd]);
        GFmaxvalueunitlabel.setText(gfForm.valueUnitLabels[sInd]);
        GFminvalue.setText(gfForm.minValues[sInd]);
        GFmaxvalue.setText(gfForm.maxValues[sInd]);
        gfForm.minValue = Float.parseFloat(GFminvalue.getText());
        gfForm.maxValue = Float.parseFloat(GFmaxvalue.getText());
    }//GEN-LAST:event_GFCalcBoxActionPerformed

    private void GFCalcStartActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_GFCalcStartActionPerformed
        // TODO add your handling code here:
        // Checking if already running
        if (gfForm.working) {
            gfForm.cancel();
            GFCalcStart.setText("Calculate");
            GFCalcSave.setEnabled(true);
            return;
        }
        GFProgressBar.setValue(0);
        GFProgressBar.setStringPainted(true);
        GFCalcStart.setText("Terminate");
        GFCalcSave.setEnabled(false);
        GFValueSelectionBox.setEnabled(false);
        gfForm.initialize(tsourcelinear);
        /**
         * Calculating data array. Using SwingWorker class
         */
        gfForm.worker = new SwingWorker<Void, Void>() {
            @Override
            protected Void doInBackground() throws Exception {
                final double step = (gfForm.maxValueClone - gfForm.minValueClone) / (xsize - 1);
                final double offset = gfForm.minValueClone;
                List<Function<Double, Double>> func = new ArrayList<>();
                //Precise geometric factor
                func.add(xp -> {
                    double x = xp * gfForm.conversionValues[gfForm.selectedItemIndexClone];
                    switch (gfForm.selectedItemIndexClone) {
                        case 0:
                            Vector dir = new BasicVector(new double[3]);
                            dir.set(2, Math.cos(x));
                            dir.set(1, Math.sin(x));
                            gfForm.tsourceclone.getLaserPulse().setDirection(dir);
                            break;
                        case 1:
                            gfForm.tsourceclone.getLaserPulse().setDelay(x);
                            break;
                        case 2:
                            Vector sht = new BasicVector(new double[3]);
                            sht.set(2, x);
                            gfForm.tsourceclone.getElectronBunch().setShift(sht);
                            break;
                        case 3:
                            gfForm.tsourceclone.getElectronBunch().setBetax(x);
                            gfForm.tsourceclone.getElectronBunch().setBetay(x);
                            break;
                        case 4:
                            gfForm.tsourceclone.getElectronBunch().setEpsx(x);
                            gfForm.tsourceclone.getElectronBunch().setEpsy(x);
                            break;
                        case 5:
                            gfForm.tsourceclone.getLaserPulse().setRlength(x);
                            break;
                        case 6:
                            gfForm.tsourceclone.getLaserPulse().setWidth(x);
                            gfForm.tsourceclone.getElectronBunch().setxWidth(x);
                            gfForm.tsourceclone.getElectronBunch().setyWidth(x);
                            break;
                    }
                    gfForm.tsourceclone.calculateGeometricFactor();
                    gfForm.tsourceclone.calculateLinearTotalFlux();
                    setStatusBar((xp - offset) / step / (xsize - 1));
                    return GFValueSelectionBox.getSelectedIndex() == 1 ? gfForm.tsourceclone.getGeometricFactor()
                            : gfForm.tsourceclone.getGeometricFactor() * gfForm.tsourceclone.getLinearTotalFlux() * 1e-12;
                });
                //Approximate geopmetric factor
                func.add(xp -> {
                    double x = xp * gfForm.conversionValues[gfForm.selectedItemIndexClone];
                    switch (gfForm.selectedItemIndexClone) {
                        case 0:
                            Vector dir = new BasicVector(new double[3]);
                            dir.set(2, Math.cos(x));
                            dir.set(1, Math.sin(x));
                            gfForm.tsourceclone.getLaserPulse().setDirection(dir);
                            break;
                        case 1:
                            gfForm.tsourceclone.getLaserPulse().setDelay(x);
                            break;
                        case 2:
                            Vector sht = new BasicVector(new double[3]);
                            sht.set(2, x);
                            gfForm.tsourceclone.getElectronBunch().setShift(sht);
                            break;
                        case 3:
                            gfForm.tsourceclone.getElectronBunch().setBetax(x);
                            gfForm.tsourceclone.getElectronBunch().setBetay(x);
                            break;
                        case 4:
                            gfForm.tsourceclone.getElectronBunch().setEpsx(x);
                            gfForm.tsourceclone.getElectronBunch().setEpsy(x);
                            break;
                        case 5:
                            gfForm.tsourceclone.getLaserPulse().setRlength(x);
                            break;
                        case 6:
                            gfForm.tsourceclone.getLaserPulse().setWidth(x);
                            gfForm.tsourceclone.getElectronBunch().setxWidth(x);
                            gfForm.tsourceclone.getElectronBunch().setyWidth(x);
                            break;
                    }
                    gfForm.tsourceclone.calculateLinearTotalFlux();
                    return GFValueSelectionBox.getSelectedIndex() == 1 ? gfForm.tsourceclone.getApproxGeometricFactor()
                            : gfForm.tsourceclone.getApproxGeometricFactor() * gfForm.tsourceclone.getLinearTotalFlux() * 1e-15;
                });
                gfForm.chartParam.setup(func, xsize, step, offset);
                return null;
            }

            @Override
            protected void done() {
                try {
                    get();
                } catch (ExecutionException ex) {
                    Logger.getLogger(ThomsonJFrame.class.getName()).log(Level.SEVERE, null, ex);
                } catch (InterruptedException | CancellationException ex) {

                }

                if (GFValueSelectionBox.getSelectedIndex() == 1) {
                    gfForm.updateGraph(GFCalcGraph, " ");
                    gfForm.getKeys()[0] = "Geometric factor";
                    gfForm.getKeys()[1] = "Approximate geometric factor";
                } else {
                    gfForm.updateGraph(GFCalcGraph, "ph/s\u00B710\u00B9\u00B2");
                    gfForm.getKeys()[0] = "Full flux";
                    gfForm.getKeys()[1] = "Approximate full flux";
                }
                GFCalcStart.setText("Calculate");
                GFCalcSave.setEnabled(true);
                GFValueSelectionBox.setEnabled(true);
            }

            /**
             * Updating progress bar
             *
             * @param status
             */
            public void setStatusBar(final double status) {
                SwingUtilities.invokeLater(() -> GFProgressBar.setValue((int) Math.round(100 * status)));
            }
        };
        gfForm.worker.execute();
    }//GEN-LAST:event_GFCalcStartActionPerformed

    private void GFCalcSaveActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_GFCalcSaveActionPerformed
        // TODO add your handling code here:
        if (gfForm.chartPanel != null) {
            gfForm.save();
        }
    }//GEN-LAST:event_GFCalcSaveActionPerformed

    private void GFminvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_GFminvalueFocusLost
        // TODO add your handling code here:
        gfForm.minValue = TestValueWithMemory(0, 10000, GFminvalue, gfForm.minValues[gfForm.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_GFminvalueFocusLost

    private void GFminvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_GFminvalueActionPerformed
        // TODO add your handling code here:
        gfForm.minValue = TestValueWithMemory(0, 10000, GFminvalue, gfForm.minValues[gfForm.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_GFminvalueActionPerformed

    private void GFmaxvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_GFmaxvalueFocusLost
        // TODO add your handling code here:
        gfForm.maxValue = TestValueWithMemory(0, 10000, GFmaxvalue, gfForm.maxValues[gfForm.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_GFmaxvalueFocusLost

    private void GFmaxvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_GFmaxvalueActionPerformed
        // TODO add your handling code here:
        gfForm.maxValue = TestValueWithMemory(0, 10000, GFmaxvalue, gfForm.maxValues[gfForm.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_GFmaxvalueActionPerformed

    private void jCheckBoxSpreadActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jCheckBoxSpreadActionPerformed
        // Cheking spread check box
        brilForm.espread = jCheckBoxSpread.isSelected();
    }//GEN-LAST:event_jCheckBoxSpreadActionPerformed

    private void HelpItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_HelpItemActionPerformed
        // Displaying help
        JTextPane textArea = new JTextPane();
        //Reading the HTML help file
        try {
            textArea.setPage(ThomsonJFrame.class.
                    getResource("/thomsonsourcehelp/thomsonsourcehelp.html"));
        } catch (IOException e) {
            JOptionPane.showMessageDialog(null, "The help file does not exist!", "Error",
                    JOptionPane.ERROR_MESSAGE);
            return;
        }
        textArea.setPreferredSize(new Dimension(600, 400));
        textArea.setEditable(false);

        JScrollPane scrollPane = new JScrollPane();
        scrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
        scrollPane.getViewport().add(textArea, BorderLayout.CENTER);
        Object[] message = {
            "Program description", scrollPane
        };
        JOptionPane.showMessageDialog(null, message, "Help", JOptionPane.INFORMATION_MESSAGE);
    }//GEN-LAST:event_HelpItemActionPerformed

    private void ebetayvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_ebetayvalueFocusLost
        // TODO add your handling code here:
        ebunch.setBetay(TestValueWithMemory(0.1, 100, ebetayvalue, "20", oldStrings) * 1e-3);
    }//GEN-LAST:event_ebetayvalueFocusLost

    private void ebetayvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ebetayvalueActionPerformed
        // TODO add your handling code here:
        ebunch.setBetay(TestValueWithMemory(0.1, 100, ebetayvalue, "20", oldStrings) * 1e-3);
    }//GEN-LAST:event_ebetayvalueActionPerformed

    private void jMenuItemNumericalActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItemNumericalActionPerformed
        // Dispalying a window to enter numerical parameters
        Object[] message = {
            "<html>Number of points in the Monte Carlo<br/> calculation of the geometric factor:</html>", gfMonteCarloNumberBox,
            "<html>Relative precision of <br/> the numerical integration in<br/> calculations of the brilliance and polarization:</html>", numericallPrecisionBox,
            "<html>Multiplication factor for numerical shift<br/> in integrals:</html>", shiftFactorBox,
            "<html>Number of used threads:</html>", threadsNumberBox
        };
        int option = JOptionPane.showConfirmDialog(null, message, "Shadow parameters", JOptionPane.OK_CANCEL_OPTION);
        if (option == JOptionPane.OK_OPTION) {
            //Non-linear source
            tsource.setPrecision((double) numericallPrecisionBox.getValue());
            tsource.setShiftfactor((double) shiftFactorBox.getValue());
            tsource.setNpGeometricFactor((int) gfMonteCarloNumberBox.getValue());
            tsource.setThreadNumber((int) threadsNumberBox.getValue());
            //Linear source
            tsourcelinear.setPrecision((double) numericallPrecisionBox.getValue());
            tsourcelinear.setShiftfactor((double) shiftFactorBox.getValue());
            tsourcelinear.setNpGeometricFactor((int) gfMonteCarloNumberBox.getValue());
            tsourcelinear.setThreadNumber((int) threadsNumberBox.getValue());
            //ChatParam
            fluxdata.setThreadNumber((int) threadsNumberBox.getValue());
            fluxcrossdata.setThreadNumber((int) threadsNumberBox.getValue());
            xenergydata.setThreadNumber((int) threadsNumberBox.getValue());
        }
    }//GEN-LAST:event_jMenuItemNumericalActionPerformed

    private void jMenuItemSourceParamActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItemSourceParamActionPerformed
        // TODO add your handling code here:
        Object[] message = {
            "Numner of rays:", rayNumberBox,
            "X-range, mrad", rayXAngleRangeBox,
            "Y-range, mrad", rayYAngleRangeBox,
            "Minimal energy, keV", rayMinEnergyBox,
            "Energy range, keV", rayEnergyRangeBox
        };
        int option = JOptionPane.showConfirmDialog(null, message, "Shadow parameters", JOptionPane.OK_CANCEL_OPTION);
        if (option == JOptionPane.OK_OPTION) {
            double eMin = (double) rayMinEnergyBox.getValue() * GaussianElectronBunch.E * 1e3;
            numberOfRays = (int) rayNumberBox.getValue();
            //Setting linear and non-linear sources
            tsource.setRayRanges((double) rayXAngleRangeBox.getValue() * 1e-3,
                    (double) rayYAngleRangeBox.getValue() * 1e-3, eMin,
                    eMin + (double) rayEnergyRangeBox.getValue() * GaussianElectronBunch.E * 1e3);
            tsourcelinear.setRayRanges((double) rayXAngleRangeBox.getValue() * 1e-3,
                    (double) rayYAngleRangeBox.getValue() * 1e-3, eMin,
                    eMin + (double) rayEnergyRangeBox.getValue() * GaussianElectronBunch.E * 1e3);
        }
    }//GEN-LAST:event_jMenuItemSourceParamActionPerformed

    private void jCheckBoxMenuItemSpreadActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jCheckBoxMenuItemSpreadActionPerformed
        // TODO add your handling code here:
        tsource.seteSpread(jCheckBoxMenuItemSpread.isSelected());
        tsourcelinear.seteSpread(jCheckBoxMenuItemSpread.isSelected());
    }//GEN-LAST:event_jCheckBoxMenuItemSpreadActionPerformed

    private void jRayStopButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jRayStopButtonActionPerformed
        // TODO add your handling code here:
        rayWorker.cancel(true);
    }//GEN-LAST:event_jRayStopButtonActionPerformed

    private void formMouseMoved(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_formMouseMoved
        // TODO add your handling code here:
        if (fluxChart != null) {
            double plotwidth = fluxChart.getchartpanel().getChartRenderingInfo().
                    getPlotInfo().getDataArea().getWidth();
            if (plotwidth != 0) {
                jSlider_pickup.setPreferredSize(new Dimension((int) plotwidth,
                        (int) jSlider_pickup.getSize().getHeight()));
            }
        }
    }//GEN-LAST:event_formMouseMoved

    private void angleValueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_angleValueFocusLost
        // TODO add your handling code here:
        brilForm.angle = TestValueWithMemory(0, 100, angleValue, "0", oldStrings);
    }//GEN-LAST:event_angleValueFocusLost

    private void angleValueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_angleValueActionPerformed
        // TODO add your handling code here:
        brilForm.angle = TestValueWithMemory(0, 100, angleValue, "0", oldStrings);
    }//GEN-LAST:event_angleValueActionPerformed

    private void energyValueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_energyValueFocusLost
        // TODO add your handling code here:
        brilForm.energy = TestValueWithMemory(20, 10000, energyValue, "46", oldStrings);
    }//GEN-LAST:event_energyValueFocusLost

    private void energyValueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_energyValueActionPerformed
        // TODO add your handling code here:
        brilForm.energy = TestValueWithMemory(20, 10000, energyValue, "46", oldStrings);
    }//GEN-LAST:event_energyValueActionPerformed

    private void jMenuItemConvActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItemConvActionPerformed
        try {
            // TODO add your handling code here:
            Runtime.getRuntime().exec("javaw.exe -jar lib\\ShadowFileConverter.jar");
        } catch (IOException ex) {
            Logger.getLogger(ThomsonJFrame.class.getName()).log(Level.SEVERE, null, ex);
        }
    }//GEN-LAST:event_jMenuItemConvActionPerformed

    private void jRadioButtonMenuItemUnPolarizedItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_jRadioButtonMenuItemUnPolarizedItemStateChanged
        // Selecting no polarization
        pRadioButtons();
    }//GEN-LAST:event_jRadioButtonMenuItemUnPolarizedItemStateChanged

    private void jRadioButtonMenuItemLinearPolarizedItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_jRadioButtonMenuItemLinearPolarizedItemStateChanged
        // Selecting s-polarization
        pRadioButtons();
    }//GEN-LAST:event_jRadioButtonMenuItemLinearPolarizedItemStateChanged

    private void jRadioButtonMenuItemCircularPolarizedItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_jRadioButtonMenuItemCircularPolarizedItemStateChanged
        // Selecting p-polarization
        pRadioButtons();
    }//GEN-LAST:event_jRadioButtonMenuItemCircularPolarizedItemStateChanged

    private void jRadioButtonMenuDefaultItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_jRadioButtonMenuDefaultItemStateChanged
        // Selecting default look&feel
        if (evt.getStateChange() == ItemEvent.SELECTED) {
            try {
                UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
            } catch (ClassNotFoundException | InstantiationException | IllegalAccessException | UnsupportedLookAndFeelException ex) {
                Logger.getLogger(ThomsonJFrame.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }//GEN-LAST:event_jRadioButtonMenuDefaultItemStateChanged

    private void jRadioButtonMenuSystemItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_jRadioButtonMenuSystemItemStateChanged
        // Selecting system look&feel
        if (evt.getStateChange() == ItemEvent.SELECTED) {
            try {
                UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
            } catch (ClassNotFoundException | InstantiationException | IllegalAccessException | UnsupportedLookAndFeelException ex) {
                Logger.getLogger(ThomsonJFrame.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }//GEN-LAST:event_jRadioButtonMenuSystemItemStateChanged

    private void jRadioButtonMenuNimbusItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_jRadioButtonMenuNimbusItemStateChanged
        // Selecting nimbus look&feel
        if (evt.getStateChange() == ItemEvent.SELECTED) {
            try {
                UIManager.setLookAndFeel("javax.swing.plaf.nimbus.NimbusLookAndFeel");
            } catch (ClassNotFoundException | InstantiationException | IllegalAccessException | UnsupportedLookAndFeelException ex) {
                Logger.getLogger(ThomsonJFrame.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }//GEN-LAST:event_jRadioButtonMenuNimbusItemStateChanged

    private void jMenuItemLaserPolarizationActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItemLaserPolarizationActionPerformed
        // Setting up laser polarization state
        String warning = "";
        double p2;
        JPanel panel = new JPanel();
        panel.add(new JLabel(""));
        ((JLabel) panel.getComponents()[0]).setForeground(Color.red);
        Object[] message = {
            "ksi1:", ksi1Box,
            "ksi2:", ksi2Box,
            "ksi3:", ksi3Box,
            panel
        };
        int option;
        do {
            ((JLabel) panel.getComponents()[0]).setText(warning);
            option = JOptionPane.showConfirmDialog(null, message, "Laser light polarization", JOptionPane.OK_CANCEL_OPTION);
            p2 = Math.pow((double) ksi1Box.getValue(), 2) + Math.pow((double) ksi2Box.getValue(), 2)
                    + Math.pow((double) ksi3Box.getValue(), 2);
            warning = p2 > 1 ? "The sum of squares of ksi1, ksi2 and ksi3 must be not exceed unity!" : "";
        } while (p2 > 1 && option == JOptionPane.OK_OPTION);
        if (option == JOptionPane.OK_OPTION) {
            lpulse.setPolarization((Double) ksi1Box.getValue(),
                    (Double) ksi2Box.getValue(), (Double) ksi3Box.getValue());
        }
    }//GEN-LAST:event_jMenuItemLaserPolarizationActionPerformed

    private void jRadioButtonMenuItemAutoPolarizedItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_jRadioButtonMenuItemAutoPolarizedItemStateChanged
        // Selecting automatic polarization
        pRadioButtons();
    }//GEN-LAST:event_jRadioButtonMenuItemAutoPolarizedItemStateChanged

    private void polarizationCalcBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_polarizationCalcBoxActionPerformed
        // Update boxes, labels and other items when parameters selection is performed
        int sInd = polarizationCalcBox.getSelectedIndex();
        polForm.selectedItemIndex = sInd;
        polminvalueunitlabel.setText(polForm.valueUnitLabels[sInd]);
        polmaxvalueunitlabel.setText(polForm.valueUnitLabels[sInd]);
        polminvalue.setText(polForm.minValues[sInd]);
        polmaxvalue.setText(polForm.maxValues[sInd]);
        polForm.minValue = Float.parseFloat(polminvalue.getText());
        polForm.maxValue = Float.parseFloat(polmaxvalue.getText());
        switch (sInd) {
            case 10:
                polAngleValue.setEnabled(true);
                polEnergyValue.setEnabled(false);
                break;
            case 11:
                polAngleValue.setEnabled(false);
                polEnergyValue.setEnabled(true);
                break;
            default:
                polAngleValue.setEnabled(true);
                polEnergyValue.setEnabled(true);
        }
    }//GEN-LAST:event_polarizationCalcBoxActionPerformed

    private void polarizationCalcStartActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_polarizationCalcStartActionPerformed
        //Polarization box calculations
        // Checking if already running
        if (polForm.working) {
            polForm.cancel();
            polarizationCalcStart.setText("Calculate");
            polarizationCalcSave.setEnabled(true);
            return;
        }
        polProgressBar.setValue(0);
        polProgressBar.setStringPainted(true);
        polarizationCalcStart.setText("Terminate");
        polarizationCalcSave.setEnabled(false);
        polForm.initialize(tsourcelinear);
        /**
         * Calculating data array. Using SwingWorker class
         */
        polForm.worker = new SwingWorker<Void, Void>() {
            @Override
            protected Void doInBackground() throws Exception {
                double step = (polForm.maxValueClone - polForm.minValueClone) / (xsize - 1);
                double offset = polForm.minValueClone;
                //A list of functions calculating intensity and polarization
                List<Function<Double, Double>> func = new ArrayList<>();
                //Itterating over polarization parameters
                for (int i = 0; i < AbstractThomsonSource.NUMBER_OF_POL_PARAM; i++) {
                    final int[] ia = new int[]{i};
                    func.add(xp -> {
                        double ang, e, x;
                        x = xp * polForm.conversionValues[polForm.selectedItemIndexClone];
                        ang = polForm.angleclone * 1e-3;
                        e = polForm.energyclone * GaussianElectronBunch.E * 1e3;
                        switch (polForm.selectedItemIndexClone) {
                            case 0:
                                Vector dir = new BasicVector(new double[3]);
                                dir.set(2, Math.cos(x));
                                dir.set(1, Math.sin(x));
                                polForm.tsourceclone.getLaserPulse().setDirection(dir);
                                break;
                            case 1:
                                polForm.tsourceclone.getLaserPulse().setDelay(x);
                                break;
                            case 2:
                                Vector sht = new BasicVector(new double[3]);
                                sht.set(2, x);
                                polForm.tsourceclone.getElectronBunch().setShift(sht);
                                break;
                            case 3:
                                polForm.tsourceclone.getElectronBunch().setBetax(x);
                                polForm.tsourceclone.getElectronBunch().setBetay(x);
                                break;
                            case 4:
                                polForm.tsourceclone.getElectronBunch().setEpsx(x);
                                polForm.tsourceclone.getElectronBunch().setEpsy(x);
                                break;
                            case 5:
                                polForm.tsourceclone.getElectronBunch().setEpsx(x);
                                break;
                            case 6:
                                polForm.tsourceclone.getElectronBunch().setEpsy(x);
                                break;
                            case 7:
                                polForm.tsourceclone.getLaserPulse().setRlength(x);
                                break;
                            case 8:
                                polForm.tsourceclone.getLaserPulse().setWidth(x);
                                polForm.tsourceclone.getElectronBunch().setxWidth(x);
                                polForm.tsourceclone.getElectronBunch().setyWidth(x);
                                break;
                            case 9:
                                polForm.tsourceclone.getElectronBunch().setDelgamma(x);
                                break;
                            case 10:
                                ang = polForm.angleclone * 1e-3;
                                e = xp * polForm.conversionValues[polForm.selectedItemIndexClone];
                                break;
                            case 11:
                                ang = xp * polForm.conversionValues[polForm.selectedItemIndexClone];
                                break;
                        }
                        setStatusBar((xp - offset) / step / (xsize - 1));
                        polForm.tsourceclone.calculateLinearTotalFlux();
                        //Calculating and returning an intensity multiplied Stocks parameter
                        try {
                            return polForm.tsourceclone.directionFrequencyPolarization(new BasicVector(new double[]{Math.sin(ang), 0, Math.cos(ang)}),
                                    new BasicVector(new double[]{0, 0, 1}), null, e, ia[0]);
                        } catch (InterruptedException ex) {
                            Thread.currentThread().interrupt();
                            return 0.0;
                        }
                    });
                }
                polForm.chartParam.setup(func, xsize, step, offset);
                return null;
            }

            @Override
            protected void done() {
                try {
                    get();
                } catch (ExecutionException ex) {
                    Logger.getLogger(ThomsonJFrame.class.getName()).log(Level.SEVERE, null, ex);
                } catch (InterruptedException | CancellationException ex) {

                }
                polForm.updateGraph(polarizationCalcGraph, "Polarization parameters");
                polarizationCalcStart.setText("Calculate");
                polarizationCalcSave.setEnabled(true);
            }

            /**
             * Updating progress bar
             *
             * @param status
             */
            public void setStatusBar(final double status) {
                SwingUtilities.invokeLater(() -> polProgressBar.setValue((int) Math.round(100 * status)));
            }
        };
        polForm.worker.execute();
    }//GEN-LAST:event_polarizationCalcStartActionPerformed

    private void polarizationCalcSaveActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_polarizationCalcSaveActionPerformed
        // Saving the brilliance plot data
        if (polForm.chartPanel != null) {
            polForm.save();
        }
    }//GEN-LAST:event_polarizationCalcSaveActionPerformed

    private void polminvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_polminvalueFocusLost
        // TODO add your handling code here:
        polForm.minValue = TestValueWithMemory(0, 10000, polminvalue, polForm.minValues[polForm.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_polminvalueFocusLost

    private void polminvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_polminvalueActionPerformed
        // TODO add your handling code here:
        polForm.minValue = TestValueWithMemory(0, 10000, polminvalue, polForm.minValues[polForm.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_polminvalueActionPerformed

    private void polmaxvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_polmaxvalueFocusLost
        // TODO add your handling code here:
        polForm.maxValue = TestValueWithMemory(0, 10000, polmaxvalue, polForm.maxValues[polForm.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_polmaxvalueFocusLost

    private void polmaxvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_polmaxvalueActionPerformed
        // TODO add your handling code here:
        polForm.maxValue = TestValueWithMemory(0, 10000, polmaxvalue, polForm.maxValues[polForm.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_polmaxvalueActionPerformed

    private void jPolCheckBoxSpreadActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jPolCheckBoxSpreadActionPerformed
        // Cheking spread check box
        polForm.espread = jPolCheckBoxSpread.isSelected();
    }//GEN-LAST:event_jPolCheckBoxSpreadActionPerformed

    private void polAngleValueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_polAngleValueFocusLost
        // Setting scattering angle in polarization form:
        polForm.angle = TestValueWithMemory(0, 100, polAngleValue, "0", oldStrings);
    }//GEN-LAST:event_polAngleValueFocusLost

    private void polAngleValueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_polAngleValueActionPerformed
        // Setting scattering angle in polarization form:
        polForm.angle = TestValueWithMemory(0, 100, polAngleValue, "0", oldStrings);
    }//GEN-LAST:event_polAngleValueActionPerformed

    private void polEnergyValueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_polEnergyValueFocusLost
        //Setting energy in polarization form:
        polForm.energy = TestValueWithMemory(20, 10000, polEnergyValue, "46", oldStrings);
    }//GEN-LAST:event_polEnergyValueFocusLost

    private void polEnergyValueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_polEnergyValueActionPerformed
        //Setting energy in polarization form:
        polForm.energy = TestValueWithMemory(20, 10000, polEnergyValue, "46", oldStrings);
    }//GEN-LAST:event_polEnergyValueActionPerformed

    private void jMenuItemPolarizationActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItemPolarizationActionPerformed
        // TODO add your handling code here:
        polarizationCalc.setVisible(true);
    }//GEN-LAST:event_jMenuItemPolarizationActionPerformed

    private void eemityvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_eemityvalueFocusLost
        // TODO add your handling code here:
        ebunch.setEpsy(TestValueWithMemory(0.1, 100, eemityvalue, "1", oldStrings) * 1e-6);
    }//GEN-LAST:event_eemityvalueFocusLost

    private void eemityvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_eemityvalueActionPerformed
        // TODO add your handling code here:
        ebunch.setEpsy(TestValueWithMemory(0.1, 100, eemityvalue, "1", oldStrings) * 1e-6);
    }//GEN-LAST:event_eemityvalueActionPerformed

    private void GFValueSelectionBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_GFValueSelectionBoxActionPerformed
        // Selecting either the full flux of geometric factor:
        if (GFValueSelectionBox.getSelectedIndex() == 0) {
            gfCalc.setTitle("Full flux box");
            ((TitledBorder) GFCalcGraph.getBorder()).setTitle("Full flux");
        } else {
            gfCalc.setTitle("Geometric factor box");
            ((TitledBorder) GFCalcGraph.getBorder()).setTitle("Geometric factor");
        }
        GFCalcGraph.repaint();
    }//GEN-LAST:event_GFValueSelectionBoxActionPerformed

    private void BrillianceCalcBoxNonLinearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_BrillianceCalcBoxNonLinearActionPerformed
        // Update boxes, labels and other items when parameters selection is performed
        int sInd = BrillianceCalcBoxNonLinear.getSelectedIndex();
        brilFormNonLinear.selectedItemIndex = sInd;
        BrilminvalueunitlabelNonLinear.setText(brilFormNonLinear.valueUnitLabels[sInd]);
        BrilmaxvalueunitlabelNonLinear.setText(brilFormNonLinear.valueUnitLabels[sInd]);
        BrilminvalueNonLinear.setText(brilFormNonLinear.minValues[sInd]);
        BrilmaxvalueNonLinear.setText(brilFormNonLinear.maxValues[sInd]);
        brilFormNonLinear.minValue = Float.parseFloat(BrilminvalueNonLinear.getText());
        brilFormNonLinear.maxValue = Float.parseFloat(BrilmaxvalueNonLinear.getText());
        switch (sInd) {
            case 10:
                angleValueNonLinear.setEnabled(true);
                energyValueNonLinear.setEnabled(false);
                break;
            case 11:
                angleValueNonLinear.setEnabled(false);
                energyValueNonLinear.setEnabled(true);
                break;
            default:
                angleValueNonLinear.setEnabled(true);
                energyValueNonLinear.setEnabled(true);
        }
    }//GEN-LAST:event_BrillianceCalcBoxNonLinearActionPerformed

    private void BrillianceCalcStartNonLinearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_BrillianceCalcStartNonLinearActionPerformed
        // Non-linear brilliance box calculations
        // Checking if already running
        if (brilFormNonLinear.working) {
            brilFormNonLinear.cancel();
            BrillianceCalcStartNonLinear.setText("Calculate");
            BrillianceCalcSaveNonLinear.setEnabled(true);
            return;
        }
        BrilProgressBarNonLinear.setValue(0);
        BrilProgressBarNonLinear.setStringPainted(true);
        BrillianceCalcStartNonLinear.setText("Terminate");
        BrillianceCalcSaveNonLinear.setEnabled(false);
        brilFormNonLinear.initialize(tsource);
        ((NonLinearThomsonSource) brilFormNonLinear.tsourceclone).setOrdernumber(brilFormNonLinear.ordernumber);
        /**
         * Calculating data array. Using SwingWorker class
         */
        brilFormNonLinear.worker = new SwingWorker<Void, Void>() {
            @Override
            protected Void doInBackground() throws Exception {
                double step = (brilFormNonLinear.maxValueClone - brilFormNonLinear.minValueClone) / (xsize - 1);
                double offset = brilFormNonLinear.minValueClone;
                List<Function<Double, Double>> func = new ArrayList<>();
                switch (brilFormNonLinear.selectedItemIndexClone) {
                    case 0:
                        func.add(xp -> {
                            double ang = brilFormNonLinear.angleclone * 1e-3;
                            double e = brilFormNonLinear.energyclone * GaussianElectronBunch.E * 1e3;
                            double x = xp * brilFormNonLinear.conversionValues[brilFormNonLinear.selectedItemIndexClone];
                            Vector dir = new BasicVector(new double[3]);
                            dir.set(2, Math.cos(x));
                            dir.set(1, Math.sin(x));
                            brilFormNonLinear.tsourceclone.getLaserPulse().setDirection(dir);
                            brilFormNonLinear.tsourceclone.calculateLinearTotalFlux();
                            setStatusBar((xp - offset) / step / (xsize - 1));
                            try {
                                return brilFormNonLinear.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0, 0, 0}),
                                        new BasicVector(new double[]{Math.sin(ang), 0, Math.cos(ang)}), new BasicVector(new double[]{0, 0, 1}),
                                        e) * normfactor;
                            } catch (InterruptedException ex) {
                                Thread.currentThread().interrupt();
                                return 0.0;
                            }
                        });
                        break;
                    case 1:
                        func.add(xp -> {
                            double ang = brilFormNonLinear.angleclone * 1e-3;
                            double e = brilFormNonLinear.energyclone * GaussianElectronBunch.E * 1e3;
                            double x = xp * brilFormNonLinear.conversionValues[brilFormNonLinear.selectedItemIndexClone];
                            brilFormNonLinear.tsourceclone.getLaserPulse().setDelay(x);
                            setStatusBar((xp - offset) / step / (xsize - 1));
                            try {
                                return brilFormNonLinear.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0, 0, 0}),
                                        new BasicVector(new double[]{Math.sin(ang), 0, Math.cos(ang)}), new BasicVector(new double[]{0, 0, 1}),
                                        e) * normfactor;
                            } catch (InterruptedException ex) {
                                Thread.currentThread().interrupt();
                                return 0.0;
                            }
                        });
                        break;
                    case 2:
                        func.add(xp -> {
                            double ang = brilFormNonLinear.angleclone * 1e-3;
                            double e = brilFormNonLinear.energyclone * GaussianElectronBunch.E * 1e3;
                            double x = xp * brilFormNonLinear.conversionValues[brilFormNonLinear.selectedItemIndexClone];
                            Vector sht = new BasicVector(new double[3]);
                            sht.set(2, x);
                            brilFormNonLinear.tsourceclone.getElectronBunch().setShift(sht);
                            brilFormNonLinear.tsourceclone.calculateLinearTotalFlux();
                            setStatusBar((xp - offset) / step / (xsize - 1));
                            try {
                                return brilFormNonLinear.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0, 0, 0}),
                                        new BasicVector(new double[]{Math.sin(ang), 0, Math.cos(ang)}), new BasicVector(new double[]{0, 0, 1}),
                                        e) * normfactor;
                            } catch (InterruptedException ex) {
                                Thread.currentThread().interrupt();
                                return 0.0;
                            }
                        });
                        break;
                    case 3:
                        func.add(xp -> {
                            double ang = brilFormNonLinear.angleclone * 1e-3;
                            double e = brilFormNonLinear.energyclone * GaussianElectronBunch.E * 1e3;
                            double x = xp * brilFormNonLinear.conversionValues[brilFormNonLinear.selectedItemIndexClone];
                            brilFormNonLinear.tsourceclone.getElectronBunch().setBetax(x);
                            brilFormNonLinear.tsourceclone.getElectronBunch().setBetay(x);
                            brilFormNonLinear.tsourceclone.calculateLinearTotalFlux();
                            setStatusBar((xp - offset) / step / (xsize - 1));
                            try {
                                return brilFormNonLinear.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0, 0, 0}),
                                        new BasicVector(new double[]{Math.sin(ang), 0, Math.cos(ang)}), new BasicVector(new double[]{0, 0, 1}),
                                        e) * normfactor;
                            } catch (InterruptedException ex) {
                                Thread.currentThread().interrupt();
                                return 0.0;
                            }
                        });
                        break;
                    case 4:
                        func.add(xp -> {
                            double ang = brilFormNonLinear.angleclone * 1e-3;
                            double e = brilFormNonLinear.energyclone * GaussianElectronBunch.E * 1e3;
                            double x = xp * brilFormNonLinear.conversionValues[brilFormNonLinear.selectedItemIndexClone];
                            brilFormNonLinear.tsourceclone.getElectronBunch().setEpsx(x);
                            brilFormNonLinear.tsourceclone.getElectronBunch().setEpsy(x);
                            brilFormNonLinear.tsourceclone.calculateLinearTotalFlux();
                            setStatusBar((xp - offset) / step / (xsize - 1));
                            try {
                                return brilFormNonLinear.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0, 0, 0}),
                                        new BasicVector(new double[]{Math.sin(ang), 0, Math.cos(ang)}), new BasicVector(new double[]{0, 0, 1}),
                                        e) * normfactor;
                            } catch (InterruptedException ex) {
                                Thread.currentThread().interrupt();
                                return 0.0;
                            }
                        });
                        break;
                    case 5:
                        func.add(xp -> {
                            double ang = brilFormNonLinear.angleclone * 1e-3;
                            double e = brilFormNonLinear.energyclone * GaussianElectronBunch.E * 1e3;
                            double x = xp * brilFormNonLinear.conversionValues[brilFormNonLinear.selectedItemIndexClone];
                            brilFormNonLinear.tsourceclone.getElectronBunch().setEpsx(x);
                            brilFormNonLinear.tsourceclone.calculateLinearTotalFlux();
                            setStatusBar((xp - offset) / step / (xsize - 1));
                            try {
                                return brilFormNonLinear.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0, 0, 0}),
                                        new BasicVector(new double[]{Math.sin(ang), 0, Math.cos(ang)}), new BasicVector(new double[]{0, 0, 1}),
                                        e) * normfactor;
                            } catch (InterruptedException ex) {
                                Thread.currentThread().interrupt();
                                return 0.0;
                            }
                        });
                        break;
                    case 6:
                        func.add(xp -> {
                            double ang = brilFormNonLinear.angleclone * 1e-3;
                            double e = brilFormNonLinear.energyclone * GaussianElectronBunch.E * 1e3;
                            double x = xp * brilFormNonLinear.conversionValues[brilFormNonLinear.selectedItemIndexClone];
                            brilFormNonLinear.tsourceclone.getElectronBunch().setEpsy(x);
                            brilFormNonLinear.tsourceclone.calculateLinearTotalFlux();
                            setStatusBar((xp - offset) / step / (xsize - 1));
                            try {
                                return brilFormNonLinear.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0, 0, 0}),
                                        new BasicVector(new double[]{Math.sin(ang), 0, Math.cos(ang)}), new BasicVector(new double[]{0, 0, 1}),
                                        e) * normfactor;
                            } catch (InterruptedException ex) {
                                Thread.currentThread().interrupt();
                                return 0.0;
                            }
                        });
                        break;
                    case 7:
                        func.add(xp -> {
                            double ang = brilFormNonLinear.angleclone * 1e-3;
                            double e = brilFormNonLinear.energyclone * GaussianElectronBunch.E * 1e3;
                            double x = xp * brilFormNonLinear.conversionValues[brilFormNonLinear.selectedItemIndexClone];
                            brilFormNonLinear.tsourceclone.getLaserPulse().setRlength(x);
                            brilFormNonLinear.tsourceclone.calculateLinearTotalFlux();
                            setStatusBar((xp - offset) / step / (xsize - 1));
                            try {
                                return brilFormNonLinear.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0, 0, 0}),
                                        new BasicVector(new double[]{Math.sin(ang), 0, Math.cos(ang)}), new BasicVector(new double[]{0, 0, 1}),
                                        e) * normfactor;
                            } catch (InterruptedException ex) {
                                Thread.currentThread().interrupt();
                                return 0.0;
                            }
                        });
                        break;
                    case 8:
                        func.add(xp -> {
                            double ang = brilFormNonLinear.angleclone * 1e-3;
                            double e = brilFormNonLinear.energyclone * GaussianElectronBunch.E * 1e3;
                            double x = xp * brilFormNonLinear.conversionValues[brilFormNonLinear.selectedItemIndexClone];
                            brilFormNonLinear.tsourceclone.getLaserPulse().setWidth(x);
                            brilFormNonLinear.tsourceclone.getElectronBunch().setxWidth(x);
                            brilFormNonLinear.tsourceclone.getElectronBunch().setyWidth(x);
                            brilFormNonLinear.tsourceclone.calculateLinearTotalFlux();
                            setStatusBar((xp - offset) / step / (xsize - 1));
                            try {
                                return brilFormNonLinear.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0, 0, 0}),
                                        new BasicVector(new double[]{Math.sin(ang), 0, Math.cos(ang)}), new BasicVector(new double[]{0, 0, 1}),
                                        e) * normfactor;
                            } catch (InterruptedException ex) {
                                Thread.currentThread().interrupt();
                                return 0.0;
                            }
                        });
                        break;
                    case 9:
                        func.add(xp -> {
                            double ang = brilFormNonLinear.angleclone * 1e-3;
                            double e = brilFormNonLinear.energyclone * GaussianElectronBunch.E * 1e3;
                            double x = xp * brilFormNonLinear.conversionValues[brilFormNonLinear.selectedItemIndexClone];
                            brilFormNonLinear.tsourceclone.getElectronBunch().setDelgamma(x);
                            brilFormNonLinear.tsourceclone.calculateLinearTotalFlux();
                            setStatusBar((xp - offset) / step / (xsize - 1));
                            try {
                                return brilFormNonLinear.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0, 0, 0}),
                                        new BasicVector(new double[]{Math.sin(ang), 0, Math.cos(ang)}), new BasicVector(new double[]{0, 0, 1}),
                                        e) * normfactor;
                            } catch (InterruptedException ex) {
                                Thread.currentThread().interrupt();
                                return 0.0;
                            }
                        });
                        break;
                    case 10:
                        func.add(xp -> {
                            double ang = brilFormNonLinear.angleclone * 1e-3;
                            double e = xp * brilFormNonLinear.conversionValues[brilFormNonLinear.selectedItemIndexClone];
                            brilFormNonLinear.tsourceclone.calculateLinearTotalFlux();
                            setStatusBar((xp - offset) / step / (xsize - 1));
                            try {
                                return brilFormNonLinear.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0.0, 0.0, 0.0}),
                                        new BasicVector(new double[]{Math.sin(ang), 0, Math.cos(ang)}), new BasicVector(new double[]{0.0, 0.0, 1.0}),
                                        e) * normfactor;
                            } catch (InterruptedException ex) {
                                Thread.currentThread().interrupt();
                                return 0.0;
                            }
                        });
                        break;
                    case 11:
                        func.add(xp -> {
                            double ang = xp * brilFormNonLinear.conversionValues[brilFormNonLinear.selectedItemIndexClone];
                            double e = brilFormNonLinear.energyclone * GaussianElectronBunch.E * 1e3;
                            brilFormNonLinear.tsourceclone.calculateLinearTotalFlux();
                            setStatusBar((xp - offset) / step / (xsize - 1));
                            try {
                                return brilFormNonLinear.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0, 0, 0}),
                                        new BasicVector(new double[]{Math.sin(ang), 0, Math.cos(ang)}), new BasicVector(new double[]{0, 0, 1}),
                                        e) * normfactor;
                            } catch (InterruptedException ex) {
                                Thread.currentThread().interrupt();
                                return 0.0;
                            }
                        });
                        break;
                }
                brilFormNonLinear.chartParam.setup(func, xsize, step, offset);
                return null;
            }

            @Override
            protected void done() {
                try {
                    get();
                } catch (ExecutionException ex) {
                    Logger.getLogger(ThomsonJFrame.class.getName()).log(Level.SEVERE, null, ex);
                } catch (InterruptedException | CancellationException ex) {

                }
                brilFormNonLinear.updateGraph(BrillianceCalcGraphNonLinear,
                        "mm\u207B\u00B2\u00B7mrad\u207B\u00B2\u00B7s\u207B\u00B9\u00B70.1%\u00B710\u00B9\u2070");
                BrillianceCalcStartNonLinear.setText("Calculate");
                BrillianceCalcSaveNonLinear.setEnabled(true);
            }

            /**
             * Updating progress bar
             *
             * @param status
             */
            public void setStatusBar(final double status) {
                SwingUtilities.invokeLater(() -> BrilProgressBarNonLinear.setValue((int) Math.round(100 * status)));
            }
        };
        brilFormNonLinear.worker.execute();
    }//GEN-LAST:event_BrillianceCalcStartNonLinearActionPerformed

    private void BrillianceCalcSaveNonLinearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_BrillianceCalcSaveNonLinearActionPerformed
        // Saving the brilliance plot data
        if (brilFormNonLinear.chartPanel != null) {
            brilFormNonLinear.save();
        }
    }//GEN-LAST:event_BrillianceCalcSaveNonLinearActionPerformed

    private void BrilminvalueNonLinearFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_BrilminvalueNonLinearFocusLost
        // TODO add your handling code here:
        brilFormNonLinear.minValue = TestValueWithMemory(0, 10000, BrilminvalueNonLinear, brilFormNonLinear.minValues[brilFormNonLinear.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_BrilminvalueNonLinearFocusLost

    private void BrilminvalueNonLinearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_BrilminvalueNonLinearActionPerformed
        // TODO add your handling code here:
        brilFormNonLinear.minValue = TestValueWithMemory(0, 10000, BrilminvalueNonLinear, brilFormNonLinear.minValues[brilFormNonLinear.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_BrilminvalueNonLinearActionPerformed

    private void BrilmaxvalueNonLinearFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_BrilmaxvalueNonLinearFocusLost
        // TODO add your handling code here:
        brilFormNonLinear.maxValue = TestValueWithMemory(0, 10000, BrilmaxvalueNonLinear, brilFormNonLinear.maxValues[brilFormNonLinear.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_BrilmaxvalueNonLinearFocusLost

    private void BrilmaxvalueNonLinearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_BrilmaxvalueNonLinearActionPerformed
        // TODO add your handling code here:
        brilFormNonLinear.maxValue = TestValueWithMemory(0, 10000, BrilmaxvalueNonLinear, brilFormNonLinear.maxValues[brilFormNonLinear.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_BrilmaxvalueNonLinearActionPerformed

    private void jCheckBoxSpreadNonLinearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jCheckBoxSpreadNonLinearActionPerformed
        // TODO add your handling code here:
        brilFormNonLinear.espread = jCheckBoxSpreadNonLinear.isSelected();
    }//GEN-LAST:event_jCheckBoxSpreadNonLinearActionPerformed

    private void angleValueNonLinearFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_angleValueNonLinearFocusLost
        // TODO add your handling code here:
        brilFormNonLinear.angle = TestValueWithMemory(0, 100, angleValueNonLinear, "0", oldStrings);
    }//GEN-LAST:event_angleValueNonLinearFocusLost

    private void angleValueNonLinearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_angleValueNonLinearActionPerformed
        // TODO add your handling code here:
        brilFormNonLinear.angle = TestValueWithMemory(0, 100, angleValueNonLinear, "0", oldStrings);
    }//GEN-LAST:event_angleValueNonLinearActionPerformed

    private void energyValueNonLinearFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_energyValueNonLinearFocusLost
        // TODO add your handling code here:
        brilFormNonLinear.energy = TestValueWithMemory(20, 10000, energyValueNonLinear, "46", oldStrings);
    }//GEN-LAST:event_energyValueNonLinearFocusLost

    private void energyValueNonLinearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_energyValueNonLinearActionPerformed
        // TODO add your handling code here:
        brilFormNonLinear.energy = TestValueWithMemory(20, 10000, energyValueNonLinear, "46", oldStrings);
    }//GEN-LAST:event_energyValueNonLinearActionPerformed

    private void OrderNumberActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_OrderNumberActionPerformed
        // TODO add your handling code here:
        brilFormNonLinear.ordernumber = (int) Math.floor(TestValueWithMemory(1, 10, OrderNumber, "1", oldStrings));
    }//GEN-LAST:event_OrderNumberActionPerformed

    private void OrderNumberFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_OrderNumberFocusLost
        // TODO add your handling code here:
        brilFormNonLinear.ordernumber = (int) Math.floor(TestValueWithMemory(1, 10, OrderNumber, "1", oldStrings));
    }//GEN-LAST:event_OrderNumberFocusLost

    private void jMenuItemBrillianceNonLinearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItemBrillianceNonLinearActionPerformed
        // TODO add your handling code here:
        brillianceCalcNonLinear.setVisible(true);
    }//GEN-LAST:event_jMenuItemBrillianceNonLinearActionPerformed

    private void jMenuItemOrderNumberActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItemOrderNumberActionPerformed
        // Order input
        Object[] message = {
            "Order number:", orderNumberBox
        };
        int option = JOptionPane.showConfirmDialog(null, message, "Non-linear parameters", JOptionPane.OK_CANCEL_OPTION);
        if (option == JOptionPane.OK_OPTION) {
            ((NonLinearThomsonSource) tsource).setOrdernumber((int) orderNumberBox.getValue());
        }
    }//GEN-LAST:event_jMenuItemOrderNumberActionPerformed

    private void polarizationCalcBoxNonLinearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_polarizationCalcBoxNonLinearActionPerformed
        // Update boxes, labels and other items when parameters selection is performed
        int sInd = polarizationCalcBoxNonLinear.getSelectedIndex();
        polFormNonLinear.selectedItemIndex = sInd;
        polminvalueunitlabelNonLinear.setText(polFormNonLinear.valueUnitLabels[sInd]);
        polmaxvalueunitlabelNonLinear.setText(polFormNonLinear.valueUnitLabels[sInd]);
        polminvalueNonLinear.setText(polFormNonLinear.minValues[sInd]);
        polmaxvalueNonLinear.setText(polFormNonLinear.maxValues[sInd]);
        polFormNonLinear.minValue = Float.parseFloat(polminvalueNonLinear.getText());
        polFormNonLinear.maxValue = Float.parseFloat(polmaxvalueNonLinear.getText());
        switch (sInd) {
            case 10:
                polAngleValueNonLinear.setEnabled(true);
                polEnergyValueNonLinear.setEnabled(false);
                break;
            case 11:
                polAngleValueNonLinear.setEnabled(false);
                polEnergyValueNonLinear.setEnabled(true);
                break;
            default:
                polAngleValueNonLinear.setEnabled(true);
                polEnergyValueNonLinear.setEnabled(true);
        }
    }//GEN-LAST:event_polarizationCalcBoxNonLinearActionPerformed

    private void polarizationCalcStartNonLinearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_polarizationCalcStartNonLinearActionPerformed
        //Polarization box calculations
        // Checking if already running
        if (polFormNonLinear.working) {
            polFormNonLinear.cancel();
            polarizationCalcStartNonLinear.setText("Calculate");
            polarizationCalcSaveNonLinear.setEnabled(true);
            return;
        }
        polProgressBarNonLinear.setValue(0);
        polProgressBarNonLinear.setStringPainted(true);
        polarizationCalcStartNonLinear.setText("Terminate");
        polarizationCalcSaveNonLinear.setEnabled(false);
        polFormNonLinear.initialize(tsource);

        /**
         * Calculating data array. Using SwingWorker class
         */
        polFormNonLinear.worker = new SwingWorker<Void, Void>() {
            @Override
            protected Void doInBackground() throws Exception {
                double step = (polFormNonLinear.maxValueClone - polFormNonLinear.minValueClone) / (xsize - 1);
                double offset = polFormNonLinear.minValueClone;
                //A list of functions calculating intensity and polarization
                List<Function<Double, Double>> func = new ArrayList<>();
                //Itterating over polarization parameters
                for (int i = 0; i < AbstractThomsonSource.NUMBER_OF_POL_PARAM; i++) {
                    final int[] ia = new int[]{i};
                    func.add(xp -> {
                        double ang, e, x;
                        x = xp * polFormNonLinear.conversionValues[polFormNonLinear.selectedItemIndexClone];
                        ang = polFormNonLinear.angleclone * 1e-3;
                        e = polFormNonLinear.energyclone * GaussianElectronBunch.E * 1e3;
                        switch (polFormNonLinear.selectedItemIndexClone) {
                            case 0:
                                Vector dir = new BasicVector(new double[3]);
                                dir.set(2, Math.cos(x));
                                dir.set(1, Math.sin(x));
                                polFormNonLinear.tsourceclone.getLaserPulse().setDirection(dir);
                                break;
                            case 1:
                                polFormNonLinear.tsourceclone.getLaserPulse().setDelay(x);
                                break;
                            case 2:
                                Vector sht = new BasicVector(new double[3]);
                                sht.set(2, x);
                                polFormNonLinear.tsourceclone.getElectronBunch().setShift(sht);
                                break;
                            case 3:
                                polFormNonLinear.tsourceclone.getElectronBunch().setBetax(x);
                                polFormNonLinear.tsourceclone.getElectronBunch().setBetay(x);
                                break;
                            case 4:
                                polFormNonLinear.tsourceclone.getElectronBunch().setEpsx(x);
                                polFormNonLinear.tsourceclone.getElectronBunch().setEpsy(x);
                                break;
                            case 5:
                                polFormNonLinear.tsourceclone.getElectronBunch().setEpsx(x);
                                break;
                            case 6:
                                polFormNonLinear.tsourceclone.getElectronBunch().setEpsy(x);
                                break;
                            case 7:
                                polFormNonLinear.tsourceclone.getLaserPulse().setRlength(x);
                                break;
                            case 8:
                                polFormNonLinear.tsourceclone.getLaserPulse().setWidth(x);
                                polFormNonLinear.tsourceclone.getElectronBunch().setxWidth(x);
                                polFormNonLinear.tsourceclone.getElectronBunch().setyWidth(x);
                                break;
                            case 9:
                                polFormNonLinear.tsourceclone.getElectronBunch().setDelgamma(x);
                                break;
                            case 10:
                                ang = polFormNonLinear.angleclone * 1e-3;
                                e = xp * polFormNonLinear.conversionValues[polFormNonLinear.selectedItemIndexClone];
                                break;
                            case 11:
                                ang = xp * polFormNonLinear.conversionValues[polFormNonLinear.selectedItemIndexClone];
                                break;
                        }
                        setStatusBar((xp - offset) / step / (xsize - 1));
                        polFormNonLinear.tsourceclone.calculateLinearTotalFlux();
                        //Calculating and returning the intensity multiplied Stocks parameter
                        try {
                            return polFormNonLinear.tsourceclone.directionFrequencyPolarization(new BasicVector(new double[]{Math.sin(ang),
                                0, Math.cos(ang)}), new BasicVector(new double[]{0, 0, 1}), null, e, ia[0]);
                        } catch (InterruptedException ex) {
                            Thread.currentThread().interrupt();
                            return 0.0;
                        }
                    });
                }
                polFormNonLinear.chartParam.setup(func, xsize, step, offset);
                return null;
            }

            @Override
            protected void done() {
                try {
                    get();
                } catch (ExecutionException ex) {
                    Logger.getLogger(ThomsonJFrame.class.getName()).log(Level.SEVERE, null, ex);
                } catch (InterruptedException | CancellationException ex) {

                }
                polFormNonLinear.updateGraph(polarizationCalcGraphNonLinear, "Polarization parameters");
                polarizationCalcStartNonLinear.setText("Calculate");
                polarizationCalcSaveNonLinear.setEnabled(true);
            }

            /**
             * Updating progress bar
             *
             * @param status
             */
            public void setStatusBar(final double status) {
                SwingUtilities.invokeLater(() -> polProgressBarNonLinear.setValue((int) Math.round(100 * status)));
            }
        };
        polFormNonLinear.worker.execute();
    }//GEN-LAST:event_polarizationCalcStartNonLinearActionPerformed

    private void polarizationCalcSaveNonLinearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_polarizationCalcSaveNonLinearActionPerformed
        // Saving the brilliance plot data
        if (polFormNonLinear.chartPanel != null) {
            polFormNonLinear.save();
        }
    }//GEN-LAST:event_polarizationCalcSaveNonLinearActionPerformed

    private void polminvalueNonLinearFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_polminvalueNonLinearFocusLost
        // TODO add your handling code here:
        polFormNonLinear.minValue = TestValueWithMemory(0, 10000, polminvalueNonLinear,
                polFormNonLinear.minValues[polFormNonLinear.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_polminvalueNonLinearFocusLost

    private void polminvalueNonLinearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_polminvalueNonLinearActionPerformed
        // TODO add your handling code here:
        polFormNonLinear.minValue = TestValueWithMemory(0, 10000, polminvalueNonLinear,
                polFormNonLinear.minValues[polFormNonLinear.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_polminvalueNonLinearActionPerformed

    private void polmaxvalueNonLinearFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_polmaxvalueNonLinearFocusLost
        // TODO add your handling code here:
        polFormNonLinear.maxValue = TestValueWithMemory(0, 10000, polmaxvalueNonLinear,
                polFormNonLinear.maxValues[polFormNonLinear.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_polmaxvalueNonLinearFocusLost

    private void polmaxvalueNonLinearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_polmaxvalueNonLinearActionPerformed
        // TODO add your handling code here:
        polFormNonLinear.maxValue = TestValueWithMemory(0, 10000, polmaxvalueNonLinear,
                polFormNonLinear.maxValues[polFormNonLinear.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_polmaxvalueNonLinearActionPerformed

    private void jPolCheckBoxSpreadNonLinearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jPolCheckBoxSpreadNonLinearActionPerformed
        // Cheking spread check box
        polFormNonLinear.espread = jPolCheckBoxSpreadNonLinear.isSelected();
    }//GEN-LAST:event_jPolCheckBoxSpreadNonLinearActionPerformed

    private void polAngleValueNonLinearFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_polAngleValueNonLinearFocusLost
        // Setting scattering angle in polarization form:
        polFormNonLinear.angle = TestValueWithMemory(0, 100, polAngleValueNonLinear, "0", oldStrings);
    }//GEN-LAST:event_polAngleValueNonLinearFocusLost

    private void polAngleValueNonLinearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_polAngleValueNonLinearActionPerformed
        // Setting scattering angle in polarization form:
        polFormNonLinear.angle = TestValueWithMemory(0, 100, polAngleValueNonLinear, "0", oldStrings);
    }//GEN-LAST:event_polAngleValueNonLinearActionPerformed

    private void polEnergyValueNonLinearFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_polEnergyValueNonLinearFocusLost
        //Setting energy in polarization form:
        polFormNonLinear.energy = TestValueWithMemory(20, 10000, polEnergyValueNonLinear, "46", oldStrings);
    }//GEN-LAST:event_polEnergyValueNonLinearFocusLost

    private void polEnergyValueNonLinearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_polEnergyValueNonLinearActionPerformed
        //Setting energy in polarization form:
        polFormNonLinear.energy = TestValueWithMemory(20, 10000, polEnergyValueNonLinear, "46", oldStrings);
    }//GEN-LAST:event_polEnergyValueNonLinearActionPerformed

    private void jMenuItemPolarizationNonLinearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItemPolarizationNonLinearActionPerformed
        // TODO add your handling code here:
        polarizationCalcNonLinear.setVisible(true);
    }//GEN-LAST:event_jMenuItemPolarizationNonLinearActionPerformed

    private void OrderNumberNonLinearFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_OrderNumberNonLinearFocusLost
        // TODO add your handling code here:
        polFormNonLinear.ordernumber = (int) Math.floor(TestValueWithMemory(1, 10, OrderNumberNonLinear, "1", oldStrings));
    }//GEN-LAST:event_OrderNumberNonLinearFocusLost

    private void OrderNumberNonLinearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_OrderNumberNonLinearActionPerformed
        // TODO add your handling code here:
        polFormNonLinear.ordernumber = (int) Math.floor(TestValueWithMemory(1, 10, OrderNumberNonLinear, "1", oldStrings));
    }//GEN-LAST:event_OrderNumberNonLinearActionPerformed

    /*
     * Setting up polarization of X-ray radiation for linear and non-linear sources
     */
    private void pRadioButtons() {
        if (jRadioButtonMenuItemUnPolarized.isSelected()) {
            tsource.setPolarization(new double[]{0, 0, 0});
            tsourcelinear.setPolarization(new double[]{0, 0, 0});
        } else if (jRadioButtonMenuItemLinearPolarized.isSelected()) {
            tsource.setPolarization(new double[]{1, 0, 0});
            tsourcelinear.setPolarization(new double[]{1, 0, 0});
        } else if (jRadioButtonMenuItemCircularPolarized.isSelected()) {
            tsource.setPolarization(new double[]{0, -1, 0});
            tsourcelinear.setPolarization(new double[]{0, -1, 0});
        } else {
            tsource.setPolarization(null);
            tsourcelinear.setPolarization(null);
        }
    }

    /**
     * + * Loading the Thomson parameters form a file + * + * @param fl +
     */
    private void loadParameters(File fl) {
        Properties prop;
        try {
            prop = tsource.loadProperties(fl);
        } catch (IOException e) {
            JOptionPane.showMessageDialog(null, "Error while reading from the file!", "Error",
                    JOptionPane.ERROR_MESSAGE);
            return;
        } catch (NumberFormatException e) {
            JOptionPane.showMessageDialog(null, "Error in the parameter file!", "Error",
                    JOptionPane.ERROR_MESSAGE);
            return;
        }
        energyvalue.setText(prop.getProperty(tsource.getParamNames()[0], "0"));
        chargevalue.setText(prop.getProperty(tsource.getParamNames()[1], "0"));
        spreadvalue.setText(prop.getProperty(tsource.getParamNames()[2], "0"));
        elengthvalue.setText(prop.getProperty(tsource.getParamNames()[3], "0"));
        eemitxvalue.setText(prop.getProperty(tsource.getParamNames()[4], "0"));
        eemityvalue.setText(prop.getProperty(tsource.getParamNames()[5], "0"));
        ebetaxvalue.setText(prop.getProperty(tsource.getParamNames()[6], "0"));
        ebetayvalue.setText(prop.getProperty(tsource.getParamNames()[7], "0"));
        phenergyvalue.setText(prop.getProperty(tsource.getParamNames()[8], "0"));
        pulseenergyvalue.setText(prop.getProperty(tsource.getParamNames()[9], "0"));
        pulselengthvalue.setText(prop.getProperty(tsource.getParamNames()[10], "0"));
        pulserelvalue.setText(prop.getProperty(tsource.getParamNames()[11], "0"));
        pulsefreqvalue.setText(prop.getProperty(tsource.getParamNames()[12], "0"));
        pulsedelayvalue.setText(prop.getProperty(tsource.getParamNames()[13], "0"));
        eshiftxvalue.setText(prop.getProperty(tsource.getParamNames()[14], "0"));
        eshiftyvalue.setText(prop.getProperty(tsource.getParamNames()[15], "0"));
        eshiftzvalue.setText(prop.getProperty(tsource.getParamNames()[16], "0"));
        pulseanglevalue.setText(Double.toString(Math.acos(lpulse.getDirection().get(2)) * 1e3));
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
            //javax.swing.UIManager.setLookAndFeel( UIManager.getSystemLookAndFeelClassName());
        } catch (ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger(ThomsonJFrame.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(ThomsonJFrame.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(ThomsonJFrame.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(ThomsonJFrame.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>
        Locale.setDefault(new Locale("en", "US"));
        /* Create and display the form */
        java.awt.EventQueue.invokeLater(() -> {
            new ThomsonJFrame().setVisible(true);
        });
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JProgressBar BrilProgressBar;
    private javax.swing.JProgressBar BrilProgressBarNonLinear;
    private javax.swing.JComboBox BrillianceCalcBox;
    private javax.swing.JComboBox BrillianceCalcBoxNonLinear;
    private javax.swing.JPanel BrillianceCalcGraph;
    private javax.swing.JPanel BrillianceCalcGraphNonLinear;
    private javax.swing.JButton BrillianceCalcSave;
    private javax.swing.JButton BrillianceCalcSaveNonLinear;
    private javax.swing.JButton BrillianceCalcStart;
    private javax.swing.JButton BrillianceCalcStartNonLinear;
    private javax.swing.JPanel BrillianceParam;
    private javax.swing.JPanel BrillianceParamNonLinear;
    private javax.swing.JTextField Brilmaxvalue;
    private javax.swing.JTextField BrilmaxvalueNonLinear;
    private javax.swing.JLabel Brilmaxvaluelabel;
    private javax.swing.JLabel BrilmaxvaluelabelNonLinear;
    private javax.swing.JLabel Brilmaxvalueunitlabel;
    private javax.swing.JLabel BrilmaxvalueunitlabelNonLinear;
    private javax.swing.JTextField Brilminvalue;
    private javax.swing.JTextField BrilminvalueNonLinear;
    private javax.swing.JLabel Brilminvaluelabel;
    private javax.swing.JLabel BrilminvaluelabelNonLinear;
    private javax.swing.JLabel Brilminvalueunitlabel;
    private javax.swing.JLabel BrilminvalueunitlabelNonLinear;
    private javax.swing.JComboBox GFCalcBox;
    private javax.swing.JPanel GFCalcGraph;
    private javax.swing.JButton GFCalcSave;
    private javax.swing.JButton GFCalcStart;
    private javax.swing.JPanel GFParam;
    private javax.swing.JProgressBar GFProgressBar;
    private javax.swing.JComboBox GFValueSelectionBox;
    private javax.swing.JTextField GFmaxvalue;
    private javax.swing.JLabel GFmaxvaluelabel;
    private javax.swing.JLabel GFmaxvalueunitlabel;
    private javax.swing.JTextField GFminvalue;
    private javax.swing.JLabel GFminvaluelabel;
    private javax.swing.JLabel GFminvalueunitlabel;
    private javax.swing.JMenuItem HelpItem;
    private javax.swing.JProgressBar MainProgressBar;
    private javax.swing.JTextField OrderNumber;
    private javax.swing.JLabel OrderNumberLabel;
    private javax.swing.JLabel OrderNumberLabelNonLinear;
    private javax.swing.JTextField OrderNumberNonLinear;
    private javax.swing.JTextField angleValue;
    private javax.swing.JTextField angleValueNonLinear;
    private javax.swing.JLabel angleValueUnitLable;
    private javax.swing.JLabel angleValueUnitLableNonLinear;
    private javax.swing.JFrame brillianceCalc;
    private javax.swing.JFrame brillianceCalcNonLinear;
    private javax.swing.ButtonGroup buttonGroupPolarization;
    private javax.swing.ButtonGroup buttonGroupSkin;
    private javax.swing.JLabel chargelabel;
    private javax.swing.JLabel chargeunitlabel;
    private javax.swing.JLabel chargeunitlabel1;
    private javax.swing.JTextField chargevalue;
    private javax.swing.JLabel ebetaxlabel;
    private javax.swing.JLabel ebetaxunitlabel;
    private javax.swing.JTextField ebetaxvalue;
    private javax.swing.JLabel ebetaylabel;
    private javax.swing.JLabel ebetayunitlabel;
    private javax.swing.JTextField ebetayvalue;
    private javax.swing.JLabel eemitxlabel;
    private javax.swing.JLabel eemitxunitlabel;
    private javax.swing.JTextField eemitxvalue;
    private javax.swing.JLabel eemitylabel;
    private javax.swing.JLabel eemityunitlabel;
    private javax.swing.JTextField eemityvalue;
    private javax.swing.JLabel elengthlabel;
    private javax.swing.JLabel elengthunitlabel;
    private javax.swing.JTextField elengthvalue;
    private javax.swing.JTextField energyValue;
    private javax.swing.JTextField energyValueNonLinear;
    private javax.swing.JLabel energyValueUnitLable;
    private javax.swing.JLabel energyValueUnitLableNonLinear;
    private javax.swing.JLabel energylabel;
    private javax.swing.JLabel energyunitlabel;
    private javax.swing.JTextField energyvalue;
    private javax.swing.JLabel eshiftxlabel;
    private javax.swing.JLabel eshiftxunitlabel;
    private javax.swing.JTextField eshiftxvalue;
    private javax.swing.JLabel eshiftylabel;
    private javax.swing.JLabel eshiftyunitlabel;
    private javax.swing.JTextField eshiftyvalue;
    private javax.swing.JLabel eshiftzlabel;
    private javax.swing.JLabel eshiftzunitlabel;
    private javax.swing.JTextField eshiftzvalue;
    private javax.swing.JFrame gfCalc;
    private javax.swing.JLabel jAngleLabel;
    private javax.swing.JLabel jAngleLabelNonLinear;
    private javax.swing.JCheckBoxMenuItem jCheckBoxMenuItemSpread;
    private javax.swing.JCheckBox jCheckBoxSpread;
    private javax.swing.JCheckBox jCheckBoxSpreadNonLinear;
    private javax.swing.JLabel jEnergyLabel;
    private javax.swing.JLabel jEnergyLabelNonLinear;
    private javax.swing.JLabel jLabelPartialFlux;
    private javax.swing.JMenuBar jMenuBarMain;
    private javax.swing.JMenu jMenuCalc;
    private javax.swing.JMenu jMenuFile;
    private javax.swing.JMenu jMenuHelp;
    private javax.swing.JMenuItem jMenuItemAbout;
    private javax.swing.JMenuItem jMenuItemBrilliance;
    private javax.swing.JMenuItem jMenuItemBrillianceNonLinear;
    private javax.swing.JMenuItem jMenuItemConv;
    private javax.swing.JMenuItem jMenuItemExit;
    private javax.swing.JMenuItem jMenuItemGeometricFactor;
    private javax.swing.JMenuItem jMenuItemLaserPolarization;
    private javax.swing.JMenuItem jMenuItemLoadParam;
    private javax.swing.JMenuItem jMenuItemNumerical;
    private javax.swing.JMenuItem jMenuItemOrderNumber;
    private javax.swing.JMenuItem jMenuItemPolarization;
    private javax.swing.JMenuItem jMenuItemPolarizationNonLinear;
    private javax.swing.JMenuItem jMenuItemSaveParam;
    private javax.swing.JMenuItem jMenuItemSize;
    private javax.swing.JMenuItem jMenuItemSource;
    private javax.swing.JMenuItem jMenuItemSourceParam;
    private javax.swing.JMenu jMenuOptions;
    private javax.swing.JMenu jMenuPolarization;
    private javax.swing.JMenu jMenuShadow;
    private javax.swing.JMenu jMenuSkin;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel jPanel_el;
    private javax.swing.JPanel jPanel_exec;
    private javax.swing.JPanel jPanel_ph;
    private javax.swing.JPanel jPanel_sh;
    private javax.swing.JPanel jPanel_slider;
    private javax.swing.JPanel jPanel_xenergy;
    private javax.swing.JPanel jPanel_xenergy_left;
    private javax.swing.JPanel jPanel_xenergy_right;
    private javax.swing.JPanel jPanel_xflux;
    private javax.swing.JPanel jPanel_xflux_left;
    private javax.swing.JPanel jPanel_xflux_right;
    private javax.swing.JLabel jPolAngleLabel;
    private javax.swing.JLabel jPolAngleLabelNonLinear;
    private javax.swing.JCheckBox jPolCheckBoxSpread;
    private javax.swing.JCheckBox jPolCheckBoxSpreadNonLinear;
    private javax.swing.JLabel jPolEnergyLabel;
    private javax.swing.JLabel jPolEnergyLabelNonLinear;
    private javax.swing.JRadioButtonMenuItem jRadioButtonMenuDefault;
    private javax.swing.JRadioButtonMenuItem jRadioButtonMenuItemAutoPolarized;
    private javax.swing.JRadioButtonMenuItem jRadioButtonMenuItemCircularPolarized;
    private javax.swing.JRadioButtonMenuItem jRadioButtonMenuItemLinearPolarized;
    private javax.swing.JRadioButtonMenuItem jRadioButtonMenuItemUnPolarized;
    private javax.swing.JRadioButtonMenuItem jRadioButtonMenuNimbus;
    private javax.swing.JRadioButtonMenuItem jRadioButtonMenuSystem;
    private javax.swing.JProgressBar jRayProgressBar;
    private javax.swing.JButton jRayStopButton;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JPopupMenu.Separator jSeparator1;
    private javax.swing.JPopupMenu.Separator jSeparator2;
    private javax.swing.JPopupMenu.Separator jSeparator3;
    private javax.swing.JPopupMenu.Separator jSeparator4;
    private javax.swing.JPopupMenu.Separator jSeparator5;
    private javax.swing.JSlider jSlider_pickup;
    private javax.swing.JTabbedPane jTabbedPane1;
    private javax.swing.JLabel phenergylabel;
    private javax.swing.JLabel phenergyunitlabel;
    private javax.swing.JTextField phenergyvalue;
    private javax.swing.JTextField polAngleValue;
    private javax.swing.JTextField polAngleValueNonLinear;
    private javax.swing.JLabel polAngleValueUnitLable;
    private javax.swing.JLabel polAngleValueUnitLableNonLinear;
    private javax.swing.JTextField polEnergyValue;
    private javax.swing.JTextField polEnergyValueNonLinear;
    private javax.swing.JLabel polEnergyValueUnitLable;
    private javax.swing.JLabel polEnergyValueUnitLableNonLinear;
    private javax.swing.JProgressBar polProgressBar;
    private javax.swing.JProgressBar polProgressBarNonLinear;
    private javax.swing.JFrame polarizationCalc;
    private javax.swing.JComboBox polarizationCalcBox;
    private javax.swing.JComboBox polarizationCalcBoxNonLinear;
    private javax.swing.JPanel polarizationCalcGraph;
    private javax.swing.JPanel polarizationCalcGraphNonLinear;
    private javax.swing.JFrame polarizationCalcNonLinear;
    private javax.swing.JButton polarizationCalcSave;
    private javax.swing.JButton polarizationCalcSaveNonLinear;
    private javax.swing.JButton polarizationCalcStart;
    private javax.swing.JButton polarizationCalcStartNonLinear;
    private javax.swing.JPanel polarizationParam;
    private javax.swing.JPanel polarizationParamcNonLinear;
    private javax.swing.JTextField polmaxvalue;
    private javax.swing.JTextField polmaxvalueNonLinear;
    private javax.swing.JLabel polmaxvaluelabel;
    private javax.swing.JLabel polmaxvaluelabelNonLinear;
    private javax.swing.JLabel polmaxvalueunitlabel;
    private javax.swing.JLabel polmaxvalueunitlabelNonLinear;
    private javax.swing.JTextField polminvalue;
    private javax.swing.JTextField polminvalueNonLinear;
    private javax.swing.JLabel polminvaluelabel;
    private javax.swing.JLabel polminvaluelabelNonLinear;
    private javax.swing.JLabel polminvalueunitlabel;
    private javax.swing.JLabel polminvalueunitlabelNonLinear;
    private javax.swing.JLabel pulseanglelabel;
    private javax.swing.JLabel pulseangleunitlabel;
    private javax.swing.JTextField pulseanglevalue;
    private javax.swing.JLabel pulsedelaylabel;
    private javax.swing.JLabel pulsedelayunitlabel;
    private javax.swing.JTextField pulsedelayvalue;
    private javax.swing.JLabel pulseenergylabel;
    private javax.swing.JLabel pulseenergyunitlabel;
    private javax.swing.JTextField pulseenergyvalue;
    private javax.swing.JLabel pulsefreqlabel;
    private javax.swing.JLabel pulsefrequnitlabel;
    private javax.swing.JTextField pulsefreqvalue;
    private javax.swing.JLabel pulselengthunitlabel;
    private javax.swing.JTextField pulselengthvalue;
    private javax.swing.JLabel pulserellabel;
    private javax.swing.JLabel pulserelunitlable;
    private javax.swing.JTextField pulserelvalue;
    private javax.swing.JLabel puslelengthlabel;
    private javax.swing.JFrame rayProgressFrame;
    private javax.swing.JLabel spreadlabel;
    private javax.swing.JTextField spreadvalue;
    private javax.swing.JButton startbutton;
    private javax.swing.JLabel totalFluxAngleLabel;
    private javax.swing.JLabel totalFluxLabel;
    // End of variables declaration//GEN-END:variables
}
