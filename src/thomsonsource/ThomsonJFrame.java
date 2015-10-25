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

import java.text.*;
import java.awt.*;
import java.awt.event.ItemEvent;
import javax.swing.*;
import javax.swing.border.TitledBorder;

import java.util.Formatter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Locale;
import java.util.concurrent.ExecutionException;
import java.util.IllegalFormatException;
import java.util.concurrent.CancellationException;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.HashMap;
import java.util.Map;
import java.util.Date;
import java.util.Enumeration;
import java.util.jar.Manifest;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.DomainOrder;
import org.jfree.data.general.DatasetChangeListener;
import org.jfree.data.general.DatasetGroup;
import org.jfree.data.xy.XYZDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.chart.renderer.PaintScale;
import org.jfree.chart.renderer.xy.XYBlockRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;

import org.la4j.vector.*;
import org.la4j.vector.dense.*;

import static TextUtilities.MyTextUtilities.*;
import java.net.URL;
import shadowfileconverter.ShadowFiles;

/**
 *
 * @author Ruslan Feshchenko
 * @version 1.7.2
 */
public class ThomsonJFrame extends javax.swing.JFrame {

    /**
     * Creates new form ThomsonJFrame
     */
    public ThomsonJFrame() {
        this.xrayenergyborder = javax.swing.BorderFactory.createTitledBorder(null, "X-ray photon energy",
                javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION);
        this.ebunch = new ElectronBunch();
        this.lpulse = new LaserPulse();
        this.tsource = new ThompsonSource(lpulse, ebunch);
        this.xsize = 300;
        this.ysize = 200;
        this.xstep = 20.0 / xsize;
        this.ystep = 20.0 / ysize;
        this.estep = 2000 / xsize;
        this.oldStrings = new HashMap<>();
        rayNumberBox = getIntegerFormattedTextField(1000, 1, 1000000);
        rayXAngleRangeBox = getDoubleFormattedTextField(0.3, 0.0, 100.0, false);
        rayYAngleRangeBox = getDoubleFormattedTextField(0.3, 0.0, 100.0, false);
        gfMonteCarloNumberBox = getIntegerFormattedTextField(5000000, 1, 100000000);
        brilPrecisionBox = getDoubleFormattedTextField(1e-4, 1e-10, 1e-1, true);
        xSizeBox = getIntegerFormattedTextField(300, 1, 10000);
        ySizeBox = getIntegerFormattedTextField(200, 1, 10000);
        xRangeBox = getDoubleFormattedTextField(20.0, 0.0, 100.0, false);
        yRangeBox = getDoubleFormattedTextField(20.0, 0.0, 100.0, false);
        xEnergyRangeBox = getDoubleFormattedTextField(2000.0, 0.0, 20000.0, false);
        rayMinEnergyBox = getDoubleFormattedTextField(36.0, 0.0, 100.0, false);
        rayEnergyRangeBox = getDoubleFormattedTextField(10.0, 0.0, 100.0, false);
        threadsNumberBox = getIntegerFormattedTextField(2, 1, 100);
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
                return 1e-9 * tsource.getGeometricFactor() * tsource.directionFrequencyFlux(n, v, e * ElectronBunch.E) / 1e10;
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
                return tsource.directionEnergy(n, v) / ElectronBunch.E * 1e-3;
            }
        };
        /**
         * Auxiliary object for linear energy chart parameters
         */
        this.xenergycrossdata = new LinearChartParam();

        /**
         * Objects for the brilliance calculation
         */
        this.brilForm = new CalcBoxParam("Brilliance");
        this.brilForm.valueUnitLabels = new String[]{"mrad", "ps", "mm", "mm", "mm mrad",
            "mm", "<html>&mu;m</html>", "", "keV", "mrad"};
        this.brilForm.plotLabels = new String[]{"Laser-electron angle, mrad", "Delay, ps", "Z-shift, mm", "beta, mm",
            "eps, mm mrad", "Reyleigh length, mm", "Waist semi-width, \u03BCm", "\u0394\u03B3/\u03B3",
            "X-ray energy, keV", "Observation angle, mrad"};
        this.brilForm.conversionValues = new double[]{1e-3, 3e-4, 1e-3, 1e-3, 1e-6, 1e-3, 1e-6, 1.0, ElectronBunch.E * 1e3, 1e-3};
        this.brilForm.minValues = new String[]{"0", "0", "0", "10", "3", "5.4", "20", "0.001", "30", "0"};
        this.brilForm.maxValues = new String[]{"35", "100", "10", "50", "10", "10", "100", "0.01", "46", "5"};
        this.brilForm.savetext = "Choose file to save spectral brilliance data";
        this.brilForm.numberOfItems = 10;

        /**
         * Objects for the GF calculation
         */
        this.gfForm = new CalcBoxParam("Geometric factor");
        this.gfForm.valueUnitLabels = new String[]{"mrad", "ps", "mm", "mm", "mm mrad", "mm", "<html>&mu;m</html>"};
        this.gfForm.plotLabels = new String[]{"Angle, mrad", "Delay, ps", "Z-shift, mm", "beta, mm",
            "eps, mm mrad", "Reyleigh length, mm", "Waist semi-width, \u03BCm"};
        this.gfForm.conversionValues = new double[]{1e-3, 3e-4, 1e-3, 1e-3, 1e-6, 1e-3, 1e-6};
        this.gfForm.minValues = new String[]{"0", "0", "0", "10", "3", "2.7", "20"};
        this.gfForm.maxValues = new String[]{"35", "100", "10", "50", "10", "10", "100"};
        this.gfForm.savetext = "Choose file to save geometric factor data";
        this.gfForm.numberOfItems = 7;

        this.paramNames = new String[]{"Electron energy, MeV", "Electron bunch charge nQ",
            "Electron bunch relative energy spread", "Electron bunch length, ps",
            "Emittance, mm*mrad", "Beta-x function, mm", "Beta-y function, mm", "Photon energy, eV",
            "Pulse energy, mJ", "Laser pulse length, ps", "Rayleigh length, mm",
            "Pulse frequency, MHz", "Delay, ps", "X-shift, mm",
            "Y-shift, mm", "Z-shift, mm", "Laser-electron angle, mrad"};
        this.threadsNumberBox.setValue(new Integer(Runtime.getRuntime().availableProcessors()));

        initComponents();
        // Adding skin menua items to their button group
        buttonGroupSkin.add(jRadioButtonMenuDefault);
        buttonGroupSkin.add(jRadioButtonMenuSystem);
        buttonGroupSkin.add(jRadioButtonMenuNimbus);
        // Adding a listerner of the UI manager
        UIManager.addPropertyChangeListener(e -> {
            SwingUtilities.updateComponentTreeUI(this);
            SwingUtilities.updateComponentTreeUI(gfCalc);
            SwingUtilities.updateComponentTreeUI(brillianceCalc);
            SwingUtilities.updateComponentTreeUI(rayProgressFrame);
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
        GFCalcGraph = new javax.swing.JPanel();
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
        eemitlabel = new javax.swing.JLabel();
        eemitvalue = new javax.swing.JTextField();
        eemitunitlabel = new javax.swing.JLabel();
        ebetaxlabel = new javax.swing.JLabel();
        ebetaxvalue = new javax.swing.JTextField();
        ebetaxunitlabel = new javax.swing.JLabel();
        ebetaylabel = new javax.swing.JLabel();
        ebetayvalue = new javax.swing.JTextField();
        ebetayunitlabel = new javax.swing.JLabel();
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
        jMenuItemGeometricFactor = new javax.swing.JMenuItem();
        jMenuShadow = new javax.swing.JMenu();
        jMenuItemSource = new javax.swing.JMenuItem();
        jMenuItemSourceParam = new javax.swing.JMenuItem();
        jMenuPolarization = new javax.swing.JMenu();
        jRadioButtonMenuItemUnPolarized = new javax.swing.JRadioButtonMenuItem();
        jRadioButtonMenuItemSPolarized = new javax.swing.JRadioButtonMenuItem();
        jRadioButtonMenuItemPPolarized = new javax.swing.JRadioButtonMenuItem();
        jSeparator3 = new javax.swing.JPopupMenu.Separator();
        jMenuItemConv = new javax.swing.JMenuItem();
        jMenuOptions = new javax.swing.JMenu();
        jMenuItemSize = new javax.swing.JMenuItem();
        jMenuItemNumerical = new javax.swing.JMenuItem();
        jSeparator2 = new javax.swing.JPopupMenu.Separator();
        jCheckBoxMenuItemSpread = new javax.swing.JCheckBoxMenuItem();
        jSeparator4 = new javax.swing.JPopupMenu.Separator();
        jMenuSkin = new javax.swing.JMenu();
        jRadioButtonMenuDefault = new javax.swing.JRadioButtonMenuItem();
        jRadioButtonMenuSystem = new javax.swing.JRadioButtonMenuItem();
        jRadioButtonMenuNimbus = new javax.swing.JRadioButtonMenuItem();
        jMenuHelp = new javax.swing.JMenu();
        HelpItem = new javax.swing.JMenuItem();
        jMenuItemAbout = new javax.swing.JMenuItem();

        brillianceCalc.setTitle("Brilliance box");
        brillianceCalc.setMinimumSize(new java.awt.Dimension(760, 313));

        BrillianceParam.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Plot parameter selection", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION));

        BrillianceCalcBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Laser-electron angle", "Delay", "Z-shift", "Beta function", "Emittance", "Rayleigh length", "Waist semi-width", "Energy spread", "X-ray energy", "Observation angle" }));
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

        Brilminvalueunitlabel.setText("mrad");

        Brilmaxvalueunitlabel.setText("mrad");

        Brilmaxvalue.setText("35");
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

        energyValue.setText("44");
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

        energyValueUnitLable.setText("kev");

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

        gfCalc.setTitle("Geometric factor box");
        gfCalc.setMinimumSize(new java.awt.Dimension(586, 313));

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

        GFmaxvalue.setText("35");
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

        javax.swing.GroupLayout GFParamLayout = new javax.swing.GroupLayout(GFParam);
        GFParam.setLayout(GFParamLayout);
        GFParamLayout.setHorizontalGroup(
            GFParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(GFParamLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(GFParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(GFParamLayout.createSequentialGroup()
                        .addComponent(GFCalcBox, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addGap(42, 42, 42))
                    .addGroup(GFParamLayout.createSequentialGroup()
                        .addComponent(GFminvaluelabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(GFminvalue, javax.swing.GroupLayout.PREFERRED_SIZE, 45, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(GFminvalueunitlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 54, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(GFmaxvaluelabel, javax.swing.GroupLayout.PREFERRED_SIZE, 58, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(GFmaxvalue, javax.swing.GroupLayout.PREFERRED_SIZE, 43, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(12, 12, 12)
                        .addComponent(GFmaxvalueunitlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 56, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                .addGroup(GFParamLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(GFParamLayout.createSequentialGroup()
                        .addComponent(GFCalcStart)
                        .addGap(18, 18, 18)
                        .addComponent(GFCalcSave, javax.swing.GroupLayout.PREFERRED_SIZE, 77, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addComponent(GFProgressBar, javax.swing.GroupLayout.PREFERRED_SIZE, 172, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(51, Short.MAX_VALUE))
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
                        .addComponent(GFminvalueunitlabel)))
                .addContainerGap())
        );

        GFCalcGraph.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Geometric factor", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION));
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
            .addComponent(GFCalcGraph, javax.swing.GroupLayout.DEFAULT_SIZE, 586, Short.MAX_VALUE)
        );
        gfCalcLayout.setVerticalGroup(
            gfCalcLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(gfCalcLayout.createSequentialGroup()
                .addComponent(GFParam, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(GFCalcGraph, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
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
        setTitle("TSource");
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

        energyvalue.setText("51.2");
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

        chargevalue.setText("1");
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

        spreadlabel.setText("Spread");

        spreadvalue.setText("0.01");
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

        elengthvalue.setText("30");
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

        eemitlabel.setText("Emittance");

        eemitvalue.setText("5");
        eemitvalue.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                eemitvalueFocusLost(evt);
            }
        });
        eemitvalue.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                eemitvalueActionPerformed(evt);
            }
        });

        eemitunitlabel.setText("mm*mrad");

        ebetaxlabel.setText("Beta x");

        ebetaxvalue.setText("10");
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

        ebetaylabel.setText("Beta  y");

        ebetayvalue.setText("10");
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

        ebetayunitlabel.setText("mm");

        javax.swing.GroupLayout jPanel_elLayout = new javax.swing.GroupLayout(jPanel_el);
        jPanel_el.setLayout(jPanel_elLayout);
        jPanel_elLayout.setHorizontalGroup(
            jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel_elLayout.createSequentialGroup()
                .addGap(6, 6, 6)
                .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(spreadlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 71, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(energylabel, javax.swing.GroupLayout.PREFERRED_SIZE, 92, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(chargelabel, javax.swing.GroupLayout.PREFERRED_SIZE, 75, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(elengthlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 62, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(eemitlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 62, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(ebetaxlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 62, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(ebetaylabel, javax.swing.GroupLayout.PREFERRED_SIZE, 62, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel_elLayout.createSequentialGroup()
                        .addComponent(ebetaxvalue, javax.swing.GroupLayout.PREFERRED_SIZE, 41, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(ebetaxunitlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 48, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addContainerGap(30, Short.MAX_VALUE))
                    .addGroup(jPanel_elLayout.createSequentialGroup()
                        .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                                .addComponent(eemitvalue)
                                .addComponent(chargevalue)
                                .addComponent(energyvalue)
                                .addComponent(spreadvalue, javax.swing.GroupLayout.DEFAULT_SIZE, 41, Short.MAX_VALUE)
                                .addComponent(elengthvalue))
                            .addComponent(ebetayvalue, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, 41, Short.MAX_VALUE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(jPanel_elLayout.createSequentialGroup()
                                .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(elengthunitlabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                    .addGroup(jPanel_elLayout.createSequentialGroup()
                                        .addComponent(eemitunitlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 59, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addGap(0, 13, Short.MAX_VALUE)))
                                .addContainerGap())
                            .addGroup(jPanel_elLayout.createSequentialGroup()
                                .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(energyunitlabel)
                                    .addComponent(chargeunitlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 24, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addComponent(ebetayunitlabel, javax.swing.GroupLayout.PREFERRED_SIZE, 48, javax.swing.GroupLayout.PREFERRED_SIZE))
                                .addGap(0, 0, Short.MAX_VALUE))))))
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
                    .addComponent(spreadvalue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(elengthunitlabel)
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(elengthlabel)
                        .addComponent(elengthvalue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(eemitunitlabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(eemitvalue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(eemitlabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(ebetaxlabel)
                    .addComponent(ebetaxvalue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(ebetaxunitlabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel_elLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(ebetayvalue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(ebetaylabel)
                    .addComponent(ebetayunitlabel)))
        );

        jPanel_ph.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Laser pulse parameters", javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.DEFAULT_POSITION));
        jPanel_ph.setToolTipText("Input parameters of the laser bunch");
        jPanel_ph.setPreferredSize(new java.awt.Dimension(223, 350));

        phenergylabel.setText("Photon energy");

        phenergyvalue.setText("1.1");
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

        pulseenergyvalue.setText("20");
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

        pulselengthvalue.setText("30");
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

        pulserelvalue.setText("2.7");
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

        pulsefreqvalue.setText("79");
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

        pulsedelayunitlabel.setText("MHz");

        javax.swing.GroupLayout jPanel_phLayout = new javax.swing.GroupLayout(jPanel_ph);
        jPanel_ph.setLayout(jPanel_phLayout);
        jPanel_phLayout.setHorizontalGroup(
            jPanel_phLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel_phLayout.createSequentialGroup()
                .addContainerGap(28, Short.MAX_VALUE)
                .addGroup(jPanel_phLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel_phLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                        .addComponent(puslelengthlabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(phenergylabel, javax.swing.GroupLayout.DEFAULT_SIZE, 84, Short.MAX_VALUE)
                        .addComponent(pulserellabel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(pulsefreqlabel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(pulsedelaylabel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                    .addComponent(pulseenergylabel, javax.swing.GroupLayout.PREFERRED_SIZE, 78, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
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
                    .addComponent(pulsedelayunitlabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addGap(24, 24, 24))
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
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
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
            .addGap(0, 338, Short.MAX_VALUE)
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
            .addGap(0, 293, Short.MAX_VALUE)
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
                .addComponent(jPanel_xflux_left, javax.swing.GroupLayout.DEFAULT_SIZE, 338, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jPanel_xflux_right, javax.swing.GroupLayout.DEFAULT_SIZE, 293, Short.MAX_VALUE)
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
            .addGap(0, 338, Short.MAX_VALUE)
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
                .addComponent(jPanel_xenergy_left, javax.swing.GroupLayout.DEFAULT_SIZE, 338, Short.MAX_VALUE)
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

        javax.swing.GroupLayout jPanel_sliderLayout = new javax.swing.GroupLayout(jPanel_slider);
        jPanel_slider.setLayout(jPanel_sliderLayout);
        jPanel_sliderLayout.setHorizontalGroup(
            jPanel_sliderLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel_sliderLayout.createSequentialGroup()
                .addGap(69, 69, 69)
                .addComponent(jSlider_pickup, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(94, 94, 94)
                .addComponent(totalFluxLabel, javax.swing.GroupLayout.PREFERRED_SIZE, 240, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        jPanel_sliderLayout.setVerticalGroup(
            jPanel_sliderLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel_sliderLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel_sliderLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(jSlider_pickup, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(totalFluxLabel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
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

        pulseanglevalue.setText("0");
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
                        .addComponent(jPanel_el, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
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
                            .addComponent(jTabbedPane1, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))))
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

        jMenuItemBrilliance.setText("Brilliance...");
        jMenuItemBrilliance.setToolTipText("");
        jMenuItemBrilliance.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItemBrillianceActionPerformed(evt);
            }
        });
        jMenuCalc.add(jMenuItemBrilliance);

        jMenuItemGeometricFactor.setText("Geometric factor...");
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
        jRadioButtonMenuItemUnPolarized.setSelected(true);
        jRadioButtonMenuItemUnPolarized.setText("Unpolarized");
        jRadioButtonMenuItemUnPolarized.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                jRadioButtonMenuItemUnPolarizedItemStateChanged(evt);
            }
        });
        jMenuPolarization.add(jRadioButtonMenuItemUnPolarized);

        buttonGroupPolarization.add(jRadioButtonMenuItemSPolarized);
        jRadioButtonMenuItemSPolarized.setText("S-polarized");
        jRadioButtonMenuItemSPolarized.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                jRadioButtonMenuItemSPolarizedItemStateChanged(evt);
            }
        });
        jMenuPolarization.add(jRadioButtonMenuItemSPolarized);

        buttonGroupPolarization.add(jRadioButtonMenuItemPPolarized);
        jRadioButtonMenuItemPPolarized.setText("P-polarized");
        jRadioButtonMenuItemPPolarized.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                jRadioButtonMenuItemPPolarizedItemStateChanged(evt);
            }
        });
        jMenuPolarization.add(jRadioButtonMenuItemPPolarized);

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

        jMenuItemSize.setText("Graph parameters...");
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

        jRadioButtonMenuDefault.setText("Default");
        jRadioButtonMenuDefault.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                jRadioButtonMenuDefaultItemStateChanged(evt);
            }
        });
        jMenuSkin.add(jRadioButtonMenuDefault);

        jRadioButtonMenuSystem.setText("System");
        jRadioButtonMenuSystem.addItemListener(new java.awt.event.ItemListener() {
            public void itemStateChanged(java.awt.event.ItemEvent evt) {
                jRadioButtonMenuSystemItemStateChanged(evt);
            }
        });
        jMenuSkin.add(jRadioButtonMenuSystem);

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

        jMenuItemAbout.setText("About TSource");
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
    private ElectronBunch ebunch;
    private LaserPulse lpulse;
    private ThompsonSource tsource, tsourceRayClone = null;

    /**
     * Parameters for the calculation boxes
     */
    class CalcBoxParam {

        public String[] valueUnitLabels, plotLabels;
        public String[] minValues, maxValues;
        public int selectedItemIndex, selectedItemIndexClone;
        public int numberOfItems;
        public double minValue, maxValue = 35, minValueClone, maxValueClone;
        public double[] conversionValues;
        public ElectronBunch ebunchclone;
        public LaserPulse lpulseclone;
        public ThompsonSource tsourceclone;
        public boolean espread = false;
        public boolean working = false;
        public SwingWorker<Void, Void> worker;
        public String savetext;
        final private String key;
        public LinearChartParam chartParam;
        ChartPanel chartPanel = null;
        JFreeChart chart = null;
        double angle = 0, angleclone, energy = 44, energyclone;
        private File file = null;

        public CalcBoxParam(String key) {
            super();
            this.key = key;
            this.chartParam = new LinearChartParam();
        }

        public void initialize() {
            working = true;
            try {
                tsourceclone = (ThompsonSource) tsource.clone();
            } catch (CloneNotSupportedException ex) {
                Logger.getLogger(ThomsonJFrame.class.getName()).log(Level.SEVERE, null, ex);
            }
            ebunchclone = tsourceclone.getElectronBunch();
            lpulseclone = tsourceclone.getLaserPulse();
            tsourceclone.seteSpread(espread);
            minValueClone = minValue;
            maxValueClone = maxValue;
            angleclone = angle;
            energyclone = energy;
            selectedItemIndexClone = selectedItemIndex;
        }
        /*
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
                Formatter fm;
                try (PrintWriter pw = new PrintWriter(new FileWriter(file, false))) {
                    for (int i = 0; i < chartParam.getSize(); i++) {
                        fm = new Formatter();
                        fm.format("%f %f", new Double(i * chartParam.getStep()
                                + chartParam.getOffset()), new Double(chartParam.getData()[i]));
                        pw.println(fm);
                    }
                    pw.close();
                } catch (IOException e) {
                    JOptionPane.showMessageDialog(null, "Error while writing to the file", "Error",
                            JOptionPane.ERROR_MESSAGE);
                }
            }
        }
        /*
         * Updatin or creating chart and chartpanel
         */

        public void updateGraph(JPanel panel, String label) {
            if (chartPanel == null) {
                /**
                 * Creating chart and plot dataset
                 */
                chart = createLineChart(createLineDataset(chartParam), plotLabels[selectedItemIndexClone], label);
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
                chart.getXYPlot().getRangeAxis().setRange(chartParam.getUMin(), chartParam.getUMax());
                chart.fireChartChanged();
            }
            working = false;
        }
        /*
         * Canceling worker
         */

        public void cancel() {
            working = false;
            worker.cancel(true);
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
            this.chart = createChart(createDataset(data, slider), data, xlabel, ylabel);
            this.chartpanel = new ChartPanel(chart,
                    (int) (fraction * jPanel.getWidth()), (int) jPanel.getHeight(), 0, 0,
                    (int) (10 * jPanel.getWidth()), (int) (10 * jPanel.getHeight()),
                    false, true, true, true, true, true);
            this.colorbarchart = createColorBar(data, colorBarlabel);
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
    private int xsize, ysize, sliderposition = 50; /* The size of the graph in x or y direction */

    private double xstep, ystep, estep, hoffset = 0; /* Step in x or y direction */

    private int numberOfRays = 1000; /* Number of rays exported for Shadow */

    private final ChartParam fluxdata, fluxcrossdata, xenergydata;
    private final LinearChartParam xenergycrossdata;
    private JFreeChart xenergycrosschart = null;

    private final TitledBorder xrayenergyborder;
    private final String[] paramNames;

    private final CalcBoxParam brilForm, gfForm;
    private ColorChart fluxChart, fluxCrossChart, xEnergyChart;
    private boolean working = false, rayWorking = false;
    private SwingWorker<Void, Void> mainWorker, rayWorker;
    private Map<JTextField, String> oldStrings;
    JFormattedTextField rayNumberBox, rayXAngleRangeBox, rayYAngleRangeBox, rayMinEnergyBox, rayEnergyRangeBox,
            gfMonteCarloNumberBox, brilPrecisionBox, xSizeBox, ySizeBox, xRangeBox,
            yRangeBox, xEnergyRangeBox, threadsNumberBox;

    private File bFile = null, pFile = null;

    private void energyvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_energyvalueActionPerformed
        // TODO add your handling code here:
        ebunch.setGamma(TestValueWithMemory(0, 100, energyvalue, "51.2", oldStrings) / 0.512);
    }//GEN-LAST:event_energyvalueActionPerformed

    private void phenergyvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_phenergyvalueActionPerformed
        // TODO add your handling code here:
        lpulse.setPhotonEnergy(TestValueWithMemory(0, 10, phenergyvalue, "1.1", oldStrings) * ElectronBunch.E);
    }//GEN-LAST:event_phenergyvalueActionPerformed

    private void energyvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_energyvalueFocusLost
        // TODO add your handling code here:
        ebunch.setGamma(TestValueWithMemory(0, 100, energyvalue, "51.2", oldStrings) / 0.512);
    }//GEN-LAST:event_energyvalueFocusLost

    private void phenergyvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_phenergyvalueFocusLost
        // TODO add your handling code here:
        lpulse.setPhotonEnergy(TestValueWithMemory(0, 10, phenergyvalue, "1.1", oldStrings) * ElectronBunch.E);
    }//GEN-LAST:event_phenergyvalueFocusLost

    private void chargevalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_chargevalueActionPerformed
        // TODO add your handling code here:
        ebunch.setNumber(TestValueWithMemory(0, 10, chargevalue, "1", oldStrings) / ElectronBunch.E * 1e-9);
    }//GEN-LAST:event_chargevalueActionPerformed

    private void chargevalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_chargevalueFocusLost
        // TODO add your handling code here:
        ebunch.setNumber(TestValueWithMemory(0, 10, chargevalue, "1", oldStrings) / ElectronBunch.E * 1e-9);
    }//GEN-LAST:event_chargevalueFocusLost

    private void spreadvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_spreadvalueActionPerformed
        // TODO add your handling code here:
        ebunch.setDelgamma((double) TestValueWithMemory(0.0001, 0.1, spreadvalue, "0.01", oldStrings));
    }//GEN-LAST:event_spreadvalueActionPerformed

    private void spreadvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_spreadvalueFocusLost
        // TODO add your handling code here:
        ebunch.setDelgamma((double) TestValueWithMemory(0.0001, 0.1, spreadvalue, "0.01", oldStrings));
    }//GEN-LAST:event_spreadvalueFocusLost

    private void elengthvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_elengthvalueActionPerformed
        // TODO add your handling code here:
        ebunch.setLength(TestValueWithMemory(0, 1000, elengthvalue, "30", oldStrings) * 3e-4 / 2);
    }//GEN-LAST:event_elengthvalueActionPerformed

    private void elengthvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_elengthvalueFocusLost
        // TODO add your handling code here:
        ebunch.setLength(TestValueWithMemory(0, 1000, elengthvalue, "30", oldStrings) * 3e-4 / 2);
    }//GEN-LAST:event_elengthvalueFocusLost

    private void pulseenergyvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_pulseenergyvalueActionPerformed
        // TODO add your handling code here:
        lpulse.setPulseEnergy(TestValueWithMemory(0, 1000, pulseenergyvalue, "20", oldStrings) * 1e-3);
    }//GEN-LAST:event_pulseenergyvalueActionPerformed

    private void pulseenergyvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_pulseenergyvalueFocusLost
        // TODO add your handling code here:
        lpulse.setPulseEnergy(TestValueWithMemory(0, 1000, pulseenergyvalue, "20", oldStrings) * 1e-3);
    }//GEN-LAST:event_pulseenergyvalueFocusLost

    private void pulselengthvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_pulselengthvalueActionPerformed
        // TODO add your handling code here:
        lpulse.setLength(TestValueWithMemory(0, 1000, pulselengthvalue, "30", oldStrings) * 3e-4 / 2);
    }//GEN-LAST:event_pulselengthvalueActionPerformed

    private void pulselengthvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_pulselengthvalueFocusLost
        // TODO add your handling code here:
        lpulse.setLength(TestValueWithMemory(0, 1000, pulselengthvalue, "30", oldStrings) * 3e-4 / 2);
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
                    tsource.calculateTotalFlux();
                    tsource.calculateGeometricFactor();
                    fluxdata.setup(xsize, ysize, xstep, ystep, 0, 0);
                    setStatusBar((int) 100 / 4);
                    xenergydata.setup(xsize, ysize, xstep, ystep, 0, 0);
                    setStatusBar((int) 100 * 2 / 4);
                    fluxcrossdata.setup(xsize, ysize, estep, ystep, xenergydata.func(hoffset, 0.0) * 1e3, 0.0);
                    setStatusBar((int) 100 * 3 / 4);
                    xenergycrossdata.setup(xenergydata.getudata(),
                            (int) (xenergydata.getxsize() - 1) * sliderposition / 100,
                            false, ysize, ystep, -ystep * ysize / 2);
                    setStatusBar((int) 100);
                } catch (InterruptedException e) {

                }
                return null;
            }

            @Override
            protected void done() {
                double plotwidth = 0;
                if (!isCancelled() || fluxChart != null) {

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
                                "mrad\u207B\u00B2\u00B7s\u207B\u00B9\u00B70.1%\u00B710\u00B9\u2070", jPanel_xflux_right, 0.75, false);
                    }

                    if (xEnergyChart != null) {
                        xEnergyChart.fullupdate(xenergydata);
                    } else {
                        xEnergyChart = new ColorChart(xenergydata, "theta_x, mrad", "theta_y, mrad", "kev", jPanel_xenergy_left, 0.75, true);
                    }

                    if (xenergycrosschart != null) {
                        xenergycrosschart.getXYPlot().getRangeAxis().
                                setRange(xenergycrossdata.getData()[0], xenergycrossdata.getUMax());
                        xenergycrosschart.getXYPlot().getDomainAxis().
                                setRangeAboutValue(xenergycrossdata.getOffset()
                                        + xenergycrossdata.getSize() * xenergycrossdata.getStep() / 2, xenergycrossdata.getSize() * xenergycrossdata.getStep());
                        xenergycrosschart.fireChartChanged();
                    } else {
                        xenergycrosschart = createLineChart(createLineDataset(xenergycrossdata), "theta_y, mrad", "Energy, keV");
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
                    plotwidth = fluxChart.getchartpanel().getChartRenderingInfo().
                            getPlotInfo().getDataArea().getWidth();
                    xrayenergyborder.setTitle("X-ray photon energy" + ". Max: " + (new DecimalFormat("########.##")).format(xenergydata.getumax()) + " keV");
                    totalFluxLabel.setText("Total flux: "
                            + (new DecimalFormat("########.##")).format(tsource.getTotalFlux() * tsource.getGeometricFactor() * 1e-15)
                            + "\u00B710\u00B9\u2075\u00B7ph\u00B7s\u207B\u00B9");
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

    private void eemitvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_eemitvalueFocusLost
        // TODO add your handling code here:
        ebunch.setEps(TestValueWithMemory(0.1, 100, eemitvalue, "5", oldStrings) * 1e-6);
    }//GEN-LAST:event_eemitvalueFocusLost

    private void eemitvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_eemitvalueActionPerformed
        // TODO add your handling code here:
        ebunch.setEps(TestValueWithMemory(0.1, 100, eemitvalue, "5", oldStrings) * 1e-6);
    }//GEN-LAST:event_eemitvalueActionPerformed

    private void ebetaxvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_ebetaxvalueFocusLost
        // TODO add your handling code here:
        ebunch.setBetax(TestValueWithMemory(0.1, 100, ebetaxvalue, "10", oldStrings) * 1e-3);
    }//GEN-LAST:event_ebetaxvalueFocusLost

    private void ebetaxvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ebetaxvalueActionPerformed
        // TODO add your handling code here:
        ebunch.setBetax(TestValueWithMemory(0.1, 100, ebetaxvalue, "10", oldStrings) * 1e-3);
    }//GEN-LAST:event_ebetaxvalueActionPerformed

    private void jSlider_pickupStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_jSlider_pickupStateChanged
        // TODO add your handling code here:
        JSlider source = (JSlider) evt.getSource();
        if (!source.getValueIsAdjusting()) {
            if (working) {
                return;
            }
            MainProgressBar.setValue(0);
            MainProgressBar.setStringPainted(true);
            startbutton.setText("Stop");
            working = true;

            sliderposition = (int) source.getValue();
            hoffset = xsize * xstep * (sliderposition - 50) / 100;
            mainWorker = new SwingWorker<Void, Void>() {
                @Override
                protected Void doInBackground() throws Exception {
                    try {
                        fluxcrossdata.setup(xsize, ysize, estep, ystep, xenergydata.func(hoffset, 0.0) * 1e3, 0.0);
                        setStatusBar((int) 100);
                        xenergycrossdata.setup(xenergydata.getudata(),
                                (int) (xenergydata.getxsize() - 1) * sliderposition / 100,
                                false, ysize, ystep, -ystep * ysize / 2);
                    } catch (InterruptedException e) {

                    }
                    return null;
                }

                @Override
                protected void done() {
                    if (fluxChart != null) {
                        fluxChart.update();
                    }

                    if (fluxCrossChart != null) {
                        fluxCrossChart.fullupdate(fluxcrossdata);
                    }

                    if (xEnergyChart != null) {
                        xEnergyChart.update();
                    }

                    if (xenergycrosschart != null) {
                        xenergycrosschart.getXYPlot().getRangeAxis().
                                setRange(xenergycrossdata.getData()[0], xenergycrossdata.getUMax());
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
                    SwingUtilities.invokeLater(() -> {
                        MainProgressBar.setValue(status);
                    });
                }
            };
            mainWorker.execute();
        }
    }//GEN-LAST:event_jSlider_pickupStateChanged

    private void pulserelvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_pulserelvalueFocusLost
        // TODO add your handling code here:
        lpulse.setRlength(TestValueWithMemory(0.01, 100, pulserelvalue, "2.7", oldStrings) * 1e-3);
    }//GEN-LAST:event_pulserelvalueFocusLost

    private void pulserelvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_pulserelvalueActionPerformed
        // TODO add your handling code here:
        lpulse.setRlength(TestValueWithMemory(0.01, 100, pulserelvalue, "2.7", oldStrings) * 1e-3);
    }//GEN-LAST:event_pulserelvalueActionPerformed

    private void pulsefreqvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_pulsefreqvalueFocusLost
        // TODO add your handling code here:
        lpulse.setFq(TestValueWithMemory(0, 1000, pulsefreqvalue, "79", oldStrings) * 1e6);
    }//GEN-LAST:event_pulsefreqvalueFocusLost

    private void pulsefreqvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_pulsefreqvalueActionPerformed
        // TODO add your handling code here:
        lpulse.setFq(TestValueWithMemory(0, 1000, pulsefreqvalue, "79", oldStrings) * 1e6);
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
        ebunch.getShift().set(0, TestValueWithMemory(0, 10, eshiftxvalue, "0", oldStrings) * 1e-3);
    }//GEN-LAST:event_eshiftxvalueFocusLost

    private void eshiftxvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_eshiftxvalueActionPerformed
        // TODO add your handling code here:
        ebunch.getShift().set(0, TestValueWithMemory(0, 10, eshiftxvalue, "0", oldStrings) * 1e-3);
    }//GEN-LAST:event_eshiftxvalueActionPerformed

    private void eshiftyvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_eshiftyvalueFocusLost
        // TODO add your handling code here:
        ebunch.getShift().set(1, TestValueWithMemory(0, 10, eshiftyvalue, "0", oldStrings) * 1e-3);
    }//GEN-LAST:event_eshiftyvalueFocusLost

    private void eshiftyvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_eshiftyvalueActionPerformed
        // TODO add your handling code here:
        ebunch.getShift().set(1, TestValueWithMemory(0, 10, eshiftyvalue, "0", oldStrings) * 1e-3);
    }//GEN-LAST:event_eshiftyvalueActionPerformed

    private void eshiftzvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_eshiftzvalueFocusLost
        // TODO add your handling code here:
        ebunch.getShift().set(2, TestValueWithMemory(0, 1000, eshiftzvalue, "0", oldStrings) * 1e-3);
    }//GEN-LAST:event_eshiftzvalueFocusLost

    private void eshiftzvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_eshiftzvalueActionPerformed
        // TODO add your handling code here:
        ebunch.getShift().set(2, TestValueWithMemory(0, 1000, eshiftzvalue, "0", oldStrings) * 1e-3);
    }//GEN-LAST:event_eshiftzvalueActionPerformed

    private void pulseanglevalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_pulseanglevalueFocusLost
        // TODO add your handling code here:
        Double value = TestValueWithMemory(0, 300, pulseanglevalue, "0", oldStrings) * 1e-3;
        lpulse.getDirection().set(2, Math.cos(value));
        lpulse.getDirection().set(1, Math.sin(value));
    }//GEN-LAST:event_pulseanglevalueFocusLost

    private void pulseanglevalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_pulseanglevalueActionPerformed
        // TODO add your handling code here:
        Double value = TestValueWithMemory(0, 300, pulseanglevalue, "0", oldStrings) * 1e-3;
        lpulse.getDirection().set(2, Math.cos(value));
        lpulse.getDirection().set(1, Math.sin(value));
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
        brilForm.initialize();

        /**
         * Calculating data array. Using SwingWorker class
         */
        brilForm.worker = new SwingWorker<Void, Void>() {
            @Override
            protected Void doInBackground() throws Exception {
                double step = (brilForm.maxValueClone - brilForm.minValueClone) / (xsize - 1);
                double offset = brilForm.minValueClone;
                switch (brilForm.selectedItemIndexClone) {
                    case 0:
                        brilForm.chartParam.setup(xp -> {
                            double ang = brilForm.angleclone * 1e-3;
                            double e = xenergydata.func(ang * 1e3, 0.0) * ElectronBunch.E * 1e3;
                            double x = xp * brilForm.conversionValues[brilForm.selectedItemIndexClone];
                            brilForm.lpulseclone.getDirection().set(2, Math.cos(x));
                            brilForm.lpulseclone.getDirection().set(1, Math.sin(x));
                            brilForm.tsourceclone.calculateTotalFlux();
                            setStatusBar((int) (100 * (xp - offset) / step / (xsize - 1)));
                            return brilForm.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0.0, 0.0, 0.0}),
                                    new BasicVector(new double[]{Math.sin(ang), 0.0, Math.cos(ang)}), new BasicVector(new double[]{0.0, 0.0, 1.0}),
                                    e) * 1e-15 * 1e-13;
                        }, xsize, step, offset);
                        break;
                    case 1:
                        brilForm.chartParam.setup(xp -> {
                            double ang = brilForm.angleclone * 1e-3;
                            double e = xenergydata.func(ang * 1e3, 0.0) * ElectronBunch.E * 1e3;
                            double x = xp * brilForm.conversionValues[brilForm.selectedItemIndexClone];
                            brilForm.lpulseclone.setDelay(x);
                            setStatusBar((int) (100 * (xp - offset) / step / (xsize - 1)));
                            return brilForm.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0.0, 0.0, 0.0}),
                                    new BasicVector(new double[]{Math.sin(ang), 0.0, Math.cos(ang)}), new BasicVector(new double[]{0.0, 0.0, 1.0}),
                                    e) * 1e-15 * 1e-13;
                        }, xsize, step, offset);
                        break;
                    case 2:
                        brilForm.chartParam.setup(xp -> {
                            double ang = brilForm.angleclone * 1e-3;
                            double e = xenergydata.func(ang * 1e3, 0.0) * ElectronBunch.E * 1e3;
                            double x = xp * brilForm.conversionValues[brilForm.selectedItemIndexClone];
                            brilForm.ebunchclone.getShift().set(2, x);
                            brilForm.tsourceclone.calculateTotalFlux();
                            setStatusBar((int) (100 * (xp - offset) / step / (xsize - 1)));
                            return brilForm.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0.0, 0.0, 0.0}),
                                    new BasicVector(new double[]{Math.sin(ang), 0.0, Math.cos(ang)}), new BasicVector(new double[]{0.0, 0.0, 1.0}),
                                    e) * 1e-15 * 1e-13;
                        }, xsize, step, offset);
                        break;
                    case 3:
                        brilForm.chartParam.setup(xp -> {
                            double ang = brilForm.angleclone * 1e-3;
                            double e = xenergydata.func(ang * 1e3, 0.0) * ElectronBunch.E * 1e3;
                            double x = xp * brilForm.conversionValues[brilForm.selectedItemIndexClone];
                            brilForm.ebunchclone.setBetax(x);
                            brilForm.ebunchclone.setBetay(x);
                            brilForm.tsourceclone.calculateTotalFlux();
                            setStatusBar((int) (100 * (xp - offset) / step / (xsize - 1)));
                            return brilForm.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0.0, 0.0, 0.0}),
                                    new BasicVector(new double[]{Math.sin(ang), 0.0, Math.cos(ang)}), new BasicVector(new double[]{0.0, 0.0, 1.0}),
                                    e) * 1e-15 * 1e-13;
                        }, xsize, step, offset);
                        break;
                    case 4:
                        brilForm.chartParam.setup(xp -> {
                            double ang = brilForm.angleclone * 1e-3;
                            double e = xenergydata.func(ang * 1e3, 0.0) * ElectronBunch.E * 1e3;
                            double x = xp * brilForm.conversionValues[brilForm.selectedItemIndexClone];
                            brilForm.ebunchclone.setEps(x);
                            brilForm.tsourceclone.calculateTotalFlux();
                            setStatusBar((int) (100 * (xp - offset) / step / (xsize - 1)));
                            return brilForm.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0.0, 0.0, 0.0}),
                                    new BasicVector(new double[]{Math.sin(ang), 0.0, Math.cos(ang)}), new BasicVector(new double[]{0.0, 0.0, 1.0}),
                                    e) * 1e-15 * 1e-13;
                        }, xsize, step, offset);
                        break;
                    case 5:
                        brilForm.chartParam.setup(xp -> {
                            double ang = brilForm.angleclone * 1e-3;
                            double e = xenergydata.func(ang * 1e3, 0.0) * ElectronBunch.E * 1e3;
                            double x = xp * brilForm.conversionValues[brilForm.selectedItemIndexClone];
                            brilForm.lpulseclone.setRlength(x);
                            brilForm.tsourceclone.calculateTotalFlux();
                            setStatusBar((int) (100 * (xp - offset) / step / (xsize - 1)));
                            return brilForm.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0.0, 0.0, 0.0}),
                                    new BasicVector(new double[]{Math.sin(ang), 0.0, Math.cos(ang)}), new BasicVector(new double[]{0.0, 0.0, 1.0}),
                                    e) * 1e-15 * 1e-13;
                        }, xsize, step, offset);
                        break;
                    case 6:
                        brilForm.chartParam.setup(xp -> {
                            double ang = brilForm.angleclone * 1e-3;
                            double e = xenergydata.func(ang * 1e3, 0.0) * ElectronBunch.E * 1e3;
                            double x = xp * brilForm.conversionValues[brilForm.selectedItemIndexClone];
                            brilForm.lpulseclone.setWidth(x);
                            brilForm.ebunchclone.setxWidth(x);
                            brilForm.ebunchclone.setyWidth(x);
                            brilForm.tsourceclone.calculateTotalFlux();
                            setStatusBar((int) (100 * (xp - offset) / step / (xsize - 1)));
                            return brilForm.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0.0, 0.0, 0.0}),
                                    new BasicVector(new double[]{Math.sin(ang), 0.0, Math.cos(ang)}), new BasicVector(new double[]{0.0, 0.0, 1.0}),
                                    e) * 1e-15 * 1e-13;
                        }, xsize, step, offset);
                        break;
                    case 7:
                        brilForm.chartParam.setup(xp -> {
                            double ang = brilForm.angleclone * 1e-3;
                            double e = xenergydata.func(ang * 1e3, 0.0) * ElectronBunch.E * 1e3;
                            double x = xp * brilForm.conversionValues[brilForm.selectedItemIndexClone];
                            brilForm.ebunchclone.setDelgamma(x);
                            brilForm.tsourceclone.calculateTotalFlux();
                            setStatusBar((int) (100 * (xp - offset) / step / (xsize - 1)));
                            return brilForm.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0.0, 0.0, 0.0}),
                                    new BasicVector(new double[]{Math.sin(ang), 0.0, Math.cos(ang)}), new BasicVector(new double[]{0.0, 0.0, 1.0}),
                                    e) * 1e-15 * 1e-13;
                        }, xsize, step, offset);
                        break;
                    case 8:
                        brilForm.chartParam.setup(xp -> {
                            double ang = brilForm.angleclone * 1e-3;
                            double e = xp * brilForm.conversionValues[brilForm.selectedItemIndexClone];
                            brilForm.tsourceclone.calculateTotalFlux();
                            setStatusBar((int) (100 * (xp - offset) / step / (xsize - 1)));
                            return brilForm.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0.0, 0.0, 0.0}),
                                    new BasicVector(new double[]{Math.sin(ang), 0.0, Math.cos(ang)}), new BasicVector(new double[]{0.0, 0.0, 1.0}),
                                    e) * 1e-15 * 1e-13;
                        }, xsize, step, offset);
                        break;
                    case 9:
                        brilForm.chartParam.setup(xp -> {
                            double ang = xp * brilForm.conversionValues[brilForm.selectedItemIndexClone];
                            double e = brilForm.energyclone * ElectronBunch.E * 1e3;
                            brilForm.tsourceclone.calculateTotalFlux();
                            setStatusBar((int) (100 * (xp - offset) / step / (xsize - 1)));
                            return brilForm.tsourceclone.directionFrequencyBrilliance(new BasicVector(new double[]{0.0, 0.0, 0.0}),
                                    new BasicVector(new double[]{Math.sin(ang), 0.0, Math.cos(ang)}), new BasicVector(new double[]{0.0, 0.0, 1.0}),
                                    e) * 1e-15 * 1e-13;
                        }, xsize, step, offset);
                        break;
                }
                return null;
            }

            @Override
            protected void done() {
                brilForm.updateGraph(BrillianceCalcGraph,
                        "mm\u207B\u00B2\u00B7mrad\u207B\u00B2\u00B7s\u207B\u00B9\u00B70.1%\u00B710\u00B9\u00B3");
                BrillianceCalcStart.setText("Calculate");
                BrillianceCalcSave.setEnabled(true);
            }

            /**
             * Updating progress bar
             *
             * @param status
             */
            public void setStatusBar(final int status) {
                SwingUtilities.invokeLater(() -> BrilProgressBar.setValue(status));
            }
        };
        brilForm.worker.execute();
    }//GEN-LAST:event_BrillianceCalcStartActionPerformed

    private void BrillianceCalcSaveActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_BrillianceCalcSaveActionPerformed
        // TODO add your handling code here:
        if (brilForm.chartPanel != null) {
            brilForm.save();
        }
    }//GEN-LAST:event_BrillianceCalcSaveActionPerformed

    private void BrillianceCalcBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_BrillianceCalcBoxActionPerformed
        // TODO add your handling code here:
        int sInd = BrillianceCalcBox.getSelectedIndex();
        brilForm.selectedItemIndex = sInd;
        Brilminvalueunitlabel.setText(brilForm.valueUnitLabels[sInd]);
        Brilmaxvalueunitlabel.setText(brilForm.valueUnitLabels[sInd]);
        Brilminvalue.setText(brilForm.minValues[sInd]);
        Brilmaxvalue.setText(brilForm.maxValues[sInd]);
        brilForm.minValue = Float.parseFloat(Brilminvalue.getText());
        brilForm.maxValue = Float.parseFloat(Brilmaxvalue.getText());
        if (sInd == 9) {
            angleValue.setEnabled(false);
            energyValue.setEnabled(true);
        } else {
            angleValue.setEnabled(true);
            energyValue.setEnabled(false);
        }
    }//GEN-LAST:event_BrillianceCalcBoxActionPerformed

    private void BrilmaxvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_BrilmaxvalueActionPerformed
        // TODO add your handling code here:
        brilForm.maxValue = TestValueWithMemory(0, 1000, Brilmaxvalue, brilForm.maxValues[brilForm.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_BrilmaxvalueActionPerformed

    private void BrilminvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_BrilminvalueActionPerformed
        // TODO add your handling code here:
        brilForm.minValue = TestValueWithMemory(0, 1000, Brilminvalue, brilForm.minValues[brilForm.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_BrilminvalueActionPerformed

    private void BrilminvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_BrilminvalueFocusLost
        // TODO add your handling code here:
        brilForm.minValue = TestValueWithMemory(0, 1000, Brilminvalue, brilForm.minValues[brilForm.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_BrilminvalueFocusLost

    private void BrilmaxvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_BrilmaxvalueFocusLost
        // TODO add your handling code here:
        brilForm.maxValue = TestValueWithMemory(0, 1000, Brilmaxvalue, brilForm.maxValues[brilForm.selectedItemIndex], oldStrings);
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
                "About", JOptionPane.INFORMATION_MESSAGE);
    }//GEN-LAST:event_jMenuItemAboutActionPerformed

    private void jMenuItemSaveParamActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItemSaveParamActionPerformed
        // TODO add your handling code here:
        JFileChooser fo = new JFileChooser(pFile);
        fo.setDialogTitle("Choose file to save Thompson source parameters");
        int ans = fo.showSaveDialog(this);
        if (ans == JFileChooser.APPROVE_OPTION) {
            pFile = fo.getSelectedFile();
            if (pFile.exists()) {
                int n = JOptionPane.showConfirmDialog(null, "The file already exists. Overwrite?", "Warning",
                        JOptionPane.YES_NO_OPTION, JOptionPane.WARNING_MESSAGE);
                if (n == JOptionPane.NO_OPTION) {
                    return;
                }
            }
            Formatter fm;
            try (PrintWriter pw = new PrintWriter(new FileWriter(pFile, false))) {
                fm = new Formatter();
                fm.format("%s %.2f", paramNames[0] + ": ", ebunch.getGamma() * 0.512);
                pw.println(fm);
                fm = new Formatter();
                fm.format("%s %.2f", paramNames[1] + ": ", ebunch.getNumber() * ElectronBunch.E * 1e9);
                pw.println(fm);
                fm = new Formatter();
                fm.format("%s %.2f", paramNames[2] + ": ", ebunch.getDelgamma());
                pw.println(fm);
                fm = new Formatter();
                fm.format("%s %.2f", paramNames[3] + ": ", ebunch.getLength() * 2 / 3e-4);
                pw.println(fm);
                fm = new Formatter();
                fm.format("%s %.2f", paramNames[4] + ": ", ebunch.getEps() * 1e6);
                pw.println(fm);
                fm = new Formatter();
                fm.format("%s %.2f", paramNames[5] + ": ", ebunch.getBetax() * 1e3);
                pw.println(fm);
                fm = new Formatter();
                fm.format("%s %.2f", paramNames[6] + ": ", ebunch.getBetay() * 1e3);
                pw.println(fm);
                fm = new Formatter();
                fm.format("%s %.2f", paramNames[7] + ": ", lpulse.getPhotonEnergy() / ElectronBunch.E);
                pw.println(fm);
                fm = new Formatter();
                fm.format("%s %.2f", paramNames[8] + ": ", lpulse.getPulseEnergy() * 1e3);
                pw.println(fm);
                fm = new Formatter();
                fm.format("%s %.2f", paramNames[9] + ": ", lpulse.getLength() * 2 / 3e-4);
                pw.println(fm);
                fm = new Formatter();
                fm.format("%s %.2f", paramNames[10] + ": ", lpulse.getRlength() * 1e3);
                pw.println(fm);
                fm = new Formatter();
                fm.format("%s %.2f", paramNames[11] + ": ", lpulse.getFq() * 1e-6);
                pw.println(fm);
                fm = new Formatter();
                fm.format("%s %.2f", paramNames[12] + ": ", lpulse.getDelay() / 3e-4);
                pw.println(fm);
                fm = new Formatter();
                fm.format("%s %.2f", paramNames[13] + ": ", ebunch.getShift().get(0) * 1e3);
                pw.println(fm);
                fm = new Formatter();
                fm.format("%s %.2f", paramNames[14] + ": ", ebunch.getShift().get(1) * 1e3);
                pw.println(fm);
                fm = new Formatter();
                fm.format("%s %.2f", paramNames[15] + ": ", ebunch.getShift().get(2) * 1e3);
                pw.println(fm);
                fm = new Formatter();
                fm.format("%s %.2f", paramNames[16] + ": ", Math.acos(lpulse.getDirection().get(2)) * 1e3);
                pw.println(fm);
                pw.close();
            } catch (IOException e) {
                JOptionPane.showMessageDialog(null, "Error while writing to the file", "Error",
                        JOptionPane.ERROR_MESSAGE);
            }
        }
    }//GEN-LAST:event_jMenuItemSaveParamActionPerformed

    private void jMenuItemLoadParamActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItemLoadParamActionPerformed
        // TODO add your handling code here:
        JFileChooser fo = new JFileChooser();
        fo.setDialogTitle("Choose file to load Thompson source parameters from");
        int ans = fo.showOpenDialog(this);
        if (ans == JFileChooser.APPROVE_OPTION) {
            pFile = fo.getSelectedFile();
            ArrayList<String> inputList = new ArrayList<>();
            try (BufferedReader pr = new BufferedReader(new FileReader(pFile))) {
                String ts;
                do {
                    ts = pr.readLine();
                    if (ts != null) {
                        inputList.add(ts);
                    }
                } while (ts != null);
                pr.close();
            } catch (IOException e) {
                JOptionPane.showMessageDialog(null, "Error while reading from the file", "Error",
                        JOptionPane.ERROR_MESSAGE);
            }
            Iterator<String> itt = inputList.iterator();
            while (itt.hasNext()) {
                String ts = itt.next();
                try {
                    for (int i = 0; i < paramNames.length; i++) {
                        if (ts.contains(paramNames[i])) {
                            String tss = ts.substring(ts.indexOf(':') + 1);
                            switch (i) {
                                case 0:
                                    ebunch.setGamma(Float.parseFloat(tss) / 0.512);
                                    energyvalue.setText(tss);
                                    break;
                                case 1:
                                    ebunch.setNumber(Float.parseFloat(tss) / ElectronBunch.E * 1e-9);
                                    chargevalue.setText(tss);
                                    break;
                                case 2:
                                    ebunch.setDelgamma(Float.parseFloat(tss));
                                    spreadvalue.setText(tss);
                                    break;
                                case 3:
                                    ebunch.setLength(Float.parseFloat(tss) * 3e-4 / 2);
                                    elengthvalue.setText(tss);
                                    break;
                                case 4:
                                    ebunch.setEps(Float.parseFloat(tss) * 1e-6);
                                    eemitvalue.setText(tss);
                                    break;
                                case 5:
                                    ebunch.setBetax(Float.parseFloat(tss) * 1e-3);
                                    ebetaxvalue.setText(tss);
                                    break;
                                case 6:
                                    ebunch.setBetay(Float.parseFloat(tss) * 1e-3);
                                    ebetayvalue.setText(tss);
                                    break;
                                case 7:
                                    lpulse.setPhotonEnergy(Float.parseFloat(tss) * ElectronBunch.E);
                                    phenergyvalue.setText(tss);
                                    break;
                                case 8:
                                    lpulse.setPulseEnergy(Float.parseFloat(tss) * 1e-3);
                                    pulseenergyvalue.setText(tss);
                                    break;
                                case 9:
                                    lpulse.setLength(Float.parseFloat(tss) * 3e-4 / 2);
                                    pulselengthvalue.setText(tss);
                                    break;
                                case 10:
                                    lpulse.setRlength(Float.parseFloat(tss) * 1e-3);
                                    pulserelvalue.setText(tss);
                                    break;
                                case 11:
                                    lpulse.setFq(Float.parseFloat(tss) * 1e6);
                                    pulsefreqvalue.setText(tss);
                                    break;
                                case 12:
                                    lpulse.setDelay(Float.parseFloat(tss) * 3e-4);
                                    pulsedelayvalue.setText(tss);
                                    break;
                                case 13:
                                    ebunch.getShift().set(0, Float.parseFloat(tss) * 1e-3);
                                    eshiftxvalue.setText(tss);
                                    break;
                                case 14:
                                    ebunch.getShift().set(1, Float.parseFloat(tss) * 1e-3);
                                    eshiftyvalue.setText(tss);
                                    break;
                                case 15:
                                    ebunch.getShift().set(2, Float.parseFloat(tss) * 1e-3);
                                    eshiftzvalue.setText(tss);
                                    break;
                                case 16:
                                    lpulse.getDirection().set(2, Math.cos(Float.parseFloat(tss) * 1e-3));
                                    lpulse.getDirection().set(1, Math.sin(Float.parseFloat(tss) * 1e-3));
                                    pulseanglevalue.setText(tss);
                                    break;
                            }
                        }
                    }
                } catch (NumberFormatException e) {
                    JOptionPane.showMessageDialog(null, "Error in the file", "Error",
                            JOptionPane.ERROR_MESSAGE);
                }
            }
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
        final int number = numberOfRays;
        try {
            tsourceRayClone = (ThompsonSource) tsource.clone();
        } catch (CloneNotSupportedException ex) {

        }
        tsourceRayClone.calculateTotalFlux();
        rayWorking = true;
        rayWorker = new SwingWorker<Void, Void>() {
            @Override
            protected Void doInBackground() throws Exception {
                try (ShadowFiles shadowFile = new ShadowFiles(true, true, ThompsonSource.NUMBER_OF_COLUMNS, number, bFile)) {
                    bFile = shadowFile.getFile();
                    SwingUtilities.invokeLater(() -> rayProgressFrame.setVisible(true));
                    for (int i = 0; i < number; i++) {
                        if (isCancelled()) {
                            break;
                        }
                        //Getting a ray
                        double[] ray = tsourceRayClone.getRay();
                        //Units conversions
                        ray[0] *= 1e2;
                        ray[1] *= 1e2;
                        ray[2] *= 1e2;
                        ray[10] *= 1e-2 / LaserPulse.HC;
                        ray[11] = i;
                        shadowFile.write(ray);
                        setStatusBar((int) 100 * (i + 1) / number);
                    }
                } catch (Exception ex) {
                    throw ex;
                }
                return null;
            }

            @Override
            protected void done() {
                jRayStopButton.setEnabled(false);
                jLabelPartialFlux.setText("Flux: " + tsourceRayClone.getPartialFlux()
                        / tsourceRayClone.getCounter() * 1e-12 + " 10\u00B9\u00B2 s\u207B\u00B9");
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
                SwingUtilities.invokeLater(() -> jRayProgressBar.setValue(status));
            }
        };
        rayWorker.execute();
    }//GEN-LAST:event_jMenuItemSourceActionPerformed

    private void jMenuItemSizeActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItemSizeActionPerformed
        // TODO add your handling code here:
        Object[] message = {
            "x-size:", xSizeBox,
            "y-size:", ySizeBox,
            "x-range (mrad):", xRangeBox,
            "y-range (mrad):", yRangeBox,
            "xenergy-range (eV):", xEnergyRangeBox
        };
        int option = JOptionPane.showConfirmDialog(null, message, "Graph parameters", JOptionPane.OK_CANCEL_OPTION);
        if (option == JOptionPane.OK_OPTION) {
            xsize = (int) xSizeBox.getValue();
            ysize = (int) ySizeBox.getValue();
            xstep = (int) xRangeBox.getValue() / xsize;
            ystep = (int) yRangeBox.getValue() / ysize;
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
        gfForm.initialize();
        /**
         * Calculating data array. Using SwingWorker class
         */
        gfForm.worker = new SwingWorker<Void, Void>() {
            @Override
            protected Void doInBackground() throws Exception {
                double step = (gfForm.maxValueClone - gfForm.minValueClone) / (xsize - 1);
                double offset = gfForm.minValueClone;
                Long nt = System.nanoTime();
                switch (gfForm.selectedItemIndexClone) {
                    case 0:
                        gfForm.chartParam.setup(xp -> {
                            double x = xp * gfForm.conversionValues[gfForm.selectedItemIndexClone];
                            gfForm.lpulseclone.getDirection().set(2, Math.cos(x));
                            gfForm.lpulseclone.getDirection().set(1, Math.sin(x));
                            gfForm.tsourceclone.calculateGeometricFactor();
                            setStatusBar((int) (100 * (xp - offset) / step / (xsize - 1)));
                            return gfForm.tsourceclone.getGeometricFactor();
                        }, xsize, step, offset);
                        break;
                    case 1:
                        gfForm.chartParam.setup(xp -> {
                            double x = xp * gfForm.conversionValues[gfForm.selectedItemIndexClone];
                            gfForm.lpulseclone.setDelay(x);
                            gfForm.tsourceclone.calculateGeometricFactor();
                            setStatusBar((int) (100 * (xp - offset) / step / (xsize - 1)));
                            return gfForm.tsourceclone.getGeometricFactor();
                        }, xsize, step, offset);
                        break;
                    case 2:
                        gfForm.chartParam.setup(xp -> {
                            double x = xp * gfForm.conversionValues[gfForm.selectedItemIndexClone];
                            gfForm.ebunchclone.getShift().set(2, x);
                            gfForm.tsourceclone.calculateGeometricFactor();
                            setStatusBar((int) (100 * (xp - offset) / step / (xsize - 1)));
                            return gfForm.tsourceclone.getGeometricFactor();
                        }, xsize, step, offset);
                        break;
                    case 3:
                        gfForm.chartParam.setup(xp -> {
                            double x = xp * gfForm.conversionValues[gfForm.selectedItemIndexClone];
                            gfForm.ebunchclone.setBetax(x);
                            gfForm.ebunchclone.setBetay(x);
                            gfForm.tsourceclone.calculateGeometricFactor();
                            setStatusBar((int) (100 * (xp - offset) / step / (xsize - 1)));
                            return gfForm.tsourceclone.getGeometricFactor();
                        }, xsize, step, offset);
                        break;
                    case 4:
                        gfForm.chartParam.setup(xp -> {
                            double x = xp * gfForm.conversionValues[gfForm.selectedItemIndexClone];
                            gfForm.ebunchclone.setEps(x);
                            gfForm.tsourceclone.calculateGeometricFactor();
                            setStatusBar((int) (100 * (xp - offset) / step / (xsize - 1)));
                            return gfForm.tsourceclone.getGeometricFactor();
                        }, xsize, step, offset);
                        break;
                    case 5:
                        gfForm.chartParam.setup(xp -> {
                            double x = xp * gfForm.conversionValues[gfForm.selectedItemIndexClone];
                            gfForm.lpulseclone.setRlength(x);
                            gfForm.tsourceclone.calculateGeometricFactor();
                            setStatusBar((int) (100 * (xp - offset) / step / (xsize - 1)));
                            return gfForm.tsourceclone.getGeometricFactor();
                        }, xsize, step, offset);
                        break;
                    case 6:
                        gfForm.chartParam.setup(xp -> {
                            double x = xp * gfForm.conversionValues[gfForm.selectedItemIndexClone];
                            gfForm.lpulseclone.setWidth(x);
                            gfForm.ebunchclone.setxWidth(x);
                            gfForm.ebunchclone.setyWidth(x);
                            gfForm.tsourceclone.calculateGeometricFactor();
                            setStatusBar((int) (100 * (xp - offset) / step / (xsize - 1)));
                            return gfForm.tsourceclone.getGeometricFactor();
                        }, xsize, step, offset);
                        break;
                }
                System.out.println(System.nanoTime() - nt);
                return null;
            }

            @Override
            protected void done() {
                gfForm.updateGraph(GFCalcGraph, "");
                GFCalcStart.setText("Calculate");
                GFCalcSave.setEnabled(true);
            }

            /**
             * Updating progress bar
             *
             * @param status
             */
            public void setStatusBar(final int status) {
                SwingUtilities.invokeLater(() -> GFProgressBar.setValue(status));
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
        gfForm.minValue = TestValueWithMemory(0, 1000, GFminvalue, gfForm.minValues[gfForm.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_GFminvalueFocusLost

    private void GFminvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_GFminvalueActionPerformed
        // TODO add your handling code here:
        gfForm.minValue = TestValueWithMemory(0, 1000, GFminvalue, gfForm.minValues[gfForm.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_GFminvalueActionPerformed

    private void GFmaxvalueFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_GFmaxvalueFocusLost
        // TODO add your handling code here:
        gfForm.maxValue = TestValueWithMemory(0, 1000, GFmaxvalue, gfForm.maxValues[gfForm.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_GFmaxvalueFocusLost

    private void GFmaxvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_GFmaxvalueActionPerformed
        // TODO add your handling code here:
        gfForm.maxValue = TestValueWithMemory(0, 1000, GFmaxvalue, gfForm.maxValues[gfForm.selectedItemIndex], oldStrings);
    }//GEN-LAST:event_GFmaxvalueActionPerformed

    private void jCheckBoxSpreadActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jCheckBoxSpreadActionPerformed
        // TODO add your handling code here:
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
        ebunch.setBetay(TestValueWithMemory(1, 100, ebetayvalue, "10", oldStrings) * 1e-3);
    }//GEN-LAST:event_ebetayvalueFocusLost

    private void ebetayvalueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ebetayvalueActionPerformed
        // TODO add your handling code here:
        ebunch.setBetay(TestValueWithMemory(1, 100, ebetayvalue, "10", oldStrings) * 1e-3);
    }//GEN-LAST:event_ebetayvalueActionPerformed

    private void jMenuItemNumericalActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItemNumericalActionPerformed
        // Dispalying a window to enter numerical parameters
        Object[] message = {
            "<html>Number of points in Monte Carlo<br/> calculation of the geometric factor:</html>", gfMonteCarloNumberBox,
            "<html>Relative precision of <br/> the numerical integration in<br/> calculations of the brilliance:</html>", brilPrecisionBox,
            "<html>Number of threads</html>", threadsNumberBox
        };
        int option = JOptionPane.showConfirmDialog(null, message, "Shadow parameters", JOptionPane.OK_CANCEL_OPTION);
        if (option == JOptionPane.OK_OPTION) {
            tsource.setNpGeometricFactor((int) gfMonteCarloNumberBox.getValue());
            tsource.setPrecision((double) brilPrecisionBox.getValue());
            tsource.setThreadNumber((int) threadsNumberBox.getValue());
        }
    }//GEN-LAST:event_jMenuItemNumericalActionPerformed

    private void jMenuItemSourceParamActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItemSourceParamActionPerformed
        // TODO add your handling code here:
        Object[] message = {
            "Numner of rays:", rayNumberBox,
            "X-range, mrad", rayXAngleRangeBox,
            "Y-range, mrad", rayYAngleRangeBox,
            "Minimal energy, kev", rayMinEnergyBox,
            "Energy range, kev", rayEnergyRangeBox
        };
        int option = JOptionPane.showConfirmDialog(null, message, "Shadow parameters", JOptionPane.OK_CANCEL_OPTION);
        if (option == JOptionPane.OK_OPTION) {
            double eMin = (double) rayMinEnergyBox.getValue() * ElectronBunch.E * 1e3;
            numberOfRays = (int) rayNumberBox.getValue();
            tsource.setRayRanges((double) rayXAngleRangeBox.getValue() * 1e-3,
                    (double) rayYAngleRangeBox.getValue() * 1e-3, eMin,
                    eMin + (double) rayEnergyRangeBox.getValue() * ElectronBunch.E * 1e3);
        }
    }//GEN-LAST:event_jMenuItemSourceParamActionPerformed

    private void jCheckBoxMenuItemSpreadActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jCheckBoxMenuItemSpreadActionPerformed
        // TODO add your handling code here:
        tsource.seteSpread(jCheckBoxMenuItemSpread.isSelected());
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
        brilForm.energy = TestValueWithMemory(20, 100, energyValue, "44", oldStrings);
    }//GEN-LAST:event_energyValueFocusLost

    private void energyValueActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_energyValueActionPerformed
        // TODO add your handling code here:
        brilForm.energy = TestValueWithMemory(20, 100, energyValue, "44", oldStrings);
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

    private void jRadioButtonMenuItemSPolarizedItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_jRadioButtonMenuItemSPolarizedItemStateChanged
        // Selecting s-polarization
        pRadioButtons();
    }//GEN-LAST:event_jRadioButtonMenuItemSPolarizedItemStateChanged

    private void jRadioButtonMenuItemPPolarizedItemStateChanged(java.awt.event.ItemEvent evt) {//GEN-FIRST:event_jRadioButtonMenuItemPPolarizedItemStateChanged
        // Selecting p-polarization
        pRadioButtons();
    }//GEN-LAST:event_jRadioButtonMenuItemPPolarizedItemStateChanged

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
    /*
     * Setting up polarization of X-ray radiation
     */

    private void pRadioButtons() {
        if (jRadioButtonMenuItemUnPolarized.isSelected()) {
            tsource.setPolarization(0, 0, 0);
        } else if (jRadioButtonMenuItemSPolarized.isSelected()) {
            tsource.setPolarization(-1, 0, 0);
        } else {
            tsource.setPolarization(1, 0, 0);
        }
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

    private JFreeChart createChart(XYZDataset dataset, ChartParam data, String xlabel, String ylabel) {
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
        PaintScale scale = new JetPaintScale(0, data.getumax());
        renderer.setPaintScale(scale);
        renderer.setBlockHeight(data.getystep());
        renderer.setBlockWidth(data.getxstep());
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

    private JFreeChart createColorBar(final ChartParam data, String label) {
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
                return data.getysize();
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
                return item * data.getumax() / (data.getysize() - 1);
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
        PaintScale scale = new JetPaintScale(0, data.getumax());
        renderer.setPaintScale(scale);
        renderer.setBlockHeight((double) data.getumax() / (data.getysize() - 1));
        XYPlot plot = new XYPlot(dataset, xAxis, yAxis, renderer);
        plot.setBackgroundPaint(Color.white);
        plot.setDomainGridlinesVisible(false);
        plot.setRangeGridlinePaint(Color.black);
        JFreeChart chart = new JFreeChart("", plot);
        chart.removeLegend();
        chart.setBackgroundPaint(Color.white);
        return chart;
    }

    private XYZDataset createDataset(final ChartParam data, final boolean linemark) {
        return new XYZDataset() {
            @Override
            public int getSeriesCount() {
                return 1;
            }

            @Override
            public int getItemCount(int series) {
                return data.getxsize() * data.getysize();
            }

            @Override
            public Number getX(int series, int item) {
                return new Double(getXValue(series, item));
            }

            @Override
            public double getXValue(int series, int item) {
                return (getXindex(series, item) - data.getxsize() / 2) * data.getxstep() + data.getxoffset();
            }

            public int getXindex(int series, int item) {
                return item / data.getysize();
            }

            @Override
            public Number getY(int series, int item) {
                return new Double(getYValue(series, item));
            }

            @Override
            public double getYValue(int series, int item) {
                return (getYindex(series, item) - data.getysize() / 2) * data.getystep() + data.getyoffset();
            }

            public int getYindex(int series, int item) {
                return item - (item / data.getysize()) * data.getysize();
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
                    return data.getudata()[x][y];
                } else {
                    if (x == (int) (data.getxsize() - 1) * sliderposition / 100) {
                        return data.getumax() / 2;
                    } else {
                        return data.getudata()[x][y];
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

    private JFreeChart createLineChart(XYDataset dataset, String xlabel, String ylabel) {
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
        renderer.setSeriesLinesVisible(0, true);
        renderer.setSeriesShapesVisible(0, false);
        renderer.setSeriesStroke(0, new BasicStroke(2.0f));
        renderer.setSeriesPaint(0, Color.GREEN);
        /* Plot creation */
        XYPlot plot = new XYPlot(dataset, xAxis, yAxis, renderer);
        plot.setBackgroundPaint(Color.white);
        plot.setRangeGridlinePaint(Color.black);
        plot.setDomainGridlinePaint(Color.black);
        /* Chart creation */
        JFreeChart chart = new JFreeChart(plot);
        chart.removeLegend();
        chart.setBackgroundPaint(Color.white);
        return chart;
    }

    private XYDataset createLineDataset(final LinearChartParam data) {
        return new XYDataset() {
            @Override
            public int getSeriesCount() {
                return 1;
            }

            @Override
            public int getItemCount(int series) {
                return data.getSize();
            }

            @Override
            public Number getX(int series, int item) {
                return new Double(getXValue(series, item));
            }

            @Override
            public double getXValue(int series, int item) {
                return item * data.getStep() + data.getOffset();
            }

            @Override
            public Number getY(int series, int item) {
                return new Double(getYValue(series, item));
            }

            @Override
            public double getYValue(int series, int item) {
                return data.getData()[item];
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
                return "EnergyCrossSection";
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

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JProgressBar BrilProgressBar;
    private javax.swing.JComboBox BrillianceCalcBox;
    private javax.swing.JPanel BrillianceCalcGraph;
    private javax.swing.JButton BrillianceCalcSave;
    private javax.swing.JButton BrillianceCalcStart;
    private javax.swing.JPanel BrillianceParam;
    private javax.swing.JTextField Brilmaxvalue;
    private javax.swing.JLabel Brilmaxvaluelabel;
    private javax.swing.JLabel Brilmaxvalueunitlabel;
    private javax.swing.JTextField Brilminvalue;
    private javax.swing.JLabel Brilminvaluelabel;
    private javax.swing.JLabel Brilminvalueunitlabel;
    private javax.swing.JComboBox GFCalcBox;
    private javax.swing.JPanel GFCalcGraph;
    private javax.swing.JButton GFCalcSave;
    private javax.swing.JButton GFCalcStart;
    private javax.swing.JPanel GFParam;
    private javax.swing.JProgressBar GFProgressBar;
    private javax.swing.JTextField GFmaxvalue;
    private javax.swing.JLabel GFmaxvaluelabel;
    private javax.swing.JLabel GFmaxvalueunitlabel;
    private javax.swing.JTextField GFminvalue;
    private javax.swing.JLabel GFminvaluelabel;
    private javax.swing.JLabel GFminvalueunitlabel;
    private javax.swing.JMenuItem HelpItem;
    private javax.swing.JProgressBar MainProgressBar;
    private javax.swing.JTextField angleValue;
    private javax.swing.JLabel angleValueUnitLable;
    private javax.swing.JFrame brillianceCalc;
    private javax.swing.ButtonGroup buttonGroupPolarization;
    private javax.swing.ButtonGroup buttonGroupSkin;
    private javax.swing.JLabel chargelabel;
    private javax.swing.JLabel chargeunitlabel;
    private javax.swing.JTextField chargevalue;
    private javax.swing.JLabel ebetaxlabel;
    private javax.swing.JLabel ebetaxunitlabel;
    private javax.swing.JTextField ebetaxvalue;
    private javax.swing.JLabel ebetaylabel;
    private javax.swing.JLabel ebetayunitlabel;
    private javax.swing.JTextField ebetayvalue;
    private javax.swing.JLabel eemitlabel;
    private javax.swing.JLabel eemitunitlabel;
    private javax.swing.JTextField eemitvalue;
    private javax.swing.JLabel elengthlabel;
    private javax.swing.JLabel elengthunitlabel;
    private javax.swing.JTextField elengthvalue;
    private javax.swing.JTextField energyValue;
    private javax.swing.JLabel energyValueUnitLable;
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
    private javax.swing.JCheckBoxMenuItem jCheckBoxMenuItemSpread;
    private javax.swing.JCheckBox jCheckBoxSpread;
    private javax.swing.JLabel jEnergyLabel;
    private javax.swing.JLabel jLabelPartialFlux;
    private javax.swing.JMenuBar jMenuBarMain;
    private javax.swing.JMenu jMenuCalc;
    private javax.swing.JMenu jMenuFile;
    private javax.swing.JMenu jMenuHelp;
    private javax.swing.JMenuItem jMenuItemAbout;
    private javax.swing.JMenuItem jMenuItemBrilliance;
    private javax.swing.JMenuItem jMenuItemConv;
    private javax.swing.JMenuItem jMenuItemExit;
    private javax.swing.JMenuItem jMenuItemGeometricFactor;
    private javax.swing.JMenuItem jMenuItemLoadParam;
    private javax.swing.JMenuItem jMenuItemNumerical;
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
    private javax.swing.JRadioButtonMenuItem jRadioButtonMenuDefault;
    private javax.swing.JRadioButtonMenuItem jRadioButtonMenuItemPPolarized;
    private javax.swing.JRadioButtonMenuItem jRadioButtonMenuItemSPolarized;
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
    private javax.swing.JSlider jSlider_pickup;
    private javax.swing.JTabbedPane jTabbedPane1;
    private javax.swing.JLabel phenergylabel;
    private javax.swing.JLabel phenergyunitlabel;
    private javax.swing.JTextField phenergyvalue;
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
    private javax.swing.JLabel totalFluxLabel;
    // End of variables declaration//GEN-END:variables
}
