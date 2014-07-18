/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ui;

import IAClasses.IsoGaussian;
import ParticleTracking.Particle;
import ParticleTracking.ParticleArray;
import ParticleTracking.Timelapse_Analysis;
import ParticleTracking.UserVariables;
import UIClasses.UIMethods;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.ImageCanvas;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;
import ij.process.TypeConverter;
import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import javax.swing.DefaultBoundedRangeModel;

/**
 *
 * @author David Barry <david.barry at cancer.org.uk>
 */
public class UserInterface extends javax.swing.JDialog {

    private final Timelapse_Analysis analyser;
    private final ImagePlus imp;
    private final String title;
    private boolean wasOKed = false;

    /**
     * Creates new form UserInterface
     */
    public UserInterface(java.awt.Frame parent, boolean modal, String title, Timelapse_Analysis analyser) {
        super(parent, modal);
        this.title = title;
        this.analyser = analyser;
        imp = new ImagePlus("", analyser.getStack().getProcessor(1));
        initComponents();
        UIMethods.centreDialog(this);
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {
        java.awt.GridBagConstraints gridBagConstraints;

        jPanel1 = new javax.swing.JPanel();
        jLabel1 = new javax.swing.JLabel();
        c1Label = new javax.swing.JLabel();
        c2Label = new javax.swing.JLabel();
        c1ComboBox = new javax.swing.JComboBox();
        c2ComboBox = new javax.swing.JComboBox();
        spatResLabel = new javax.swing.JLabel();
        timeResLabel = new javax.swing.JLabel();
        minTrajLengthLabel = new javax.swing.JLabel();
        maxTrajStepLabel = new javax.swing.JLabel();
        chan1MaxThreshLabel = new javax.swing.JLabel();
        chan2MaxThreshLabel = new javax.swing.JLabel();
        spatResTextField = new javax.swing.JTextField();
        timeResTextField = new javax.swing.JTextField();
        minTrajLengthTextField = new javax.swing.JTextField();
        maxTrajStepTextField = new javax.swing.JTextField();
        chan1MaxThreshTextField = new javax.swing.JTextField();
        chan2MaxThreshTextField = new javax.swing.JTextField();
        colocaliseToggleButton = new javax.swing.JToggleButton();
        preprocessToggleButton = new javax.swing.JToggleButton();
        curveFitTolLabel = new javax.swing.JLabel();
        curveFitTolTextField = new javax.swing.JTextField();
        nMaxTextField = new javax.swing.JTextField();
        nMaxLabel = new javax.swing.JLabel();
        jPanel2 = new javax.swing.JPanel();
        canvas1 = new ImageCanvas(imp);
        previewScrollBar = new javax.swing.JScrollBar();
        previewTextField = new javax.swing.JTextField();
        previewToggleButton = new javax.swing.JToggleButton();
        jPanel3 = new javax.swing.JPanel();
        okButton = new javax.swing.JButton();
        cancelButton = new javax.swing.JButton();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle(title);
        getContentPane().setLayout(new java.awt.GridBagLayout());

        jPanel1.setBorder(javax.swing.BorderFactory.createEtchedBorder());
        jPanel1.setLayout(new java.awt.GridBagLayout());

        jLabel1.setText("Channel 2 will be co-localised with channel1:");
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 0;
        gridBagConstraints.gridwidth = 2;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_START;
        gridBagConstraints.weightx = 1.0;
        gridBagConstraints.weighty = 0.09;
        gridBagConstraints.insets = new java.awt.Insets(10, 10, 0, 10);
        jPanel1.add(jLabel1, gridBagConstraints);

        c1Label.setText("Channel 1:");
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 1;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_START;
        gridBagConstraints.weightx = 0.5;
        gridBagConstraints.weighty = 0.09;
        gridBagConstraints.insets = new java.awt.Insets(0, 10, 0, 0);
        jPanel1.add(c1Label, gridBagConstraints);

        c2Label.setText("Channel 2:");
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 2;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_START;
        gridBagConstraints.weightx = 0.5;
        gridBagConstraints.weighty = 0.09;
        gridBagConstraints.insets = new java.awt.Insets(0, 10, 0, 0);
        jPanel1.add(c2Label, gridBagConstraints);

        c1ComboBox.setModel(new javax.swing.DefaultComboBoxModel(UserVariables.channels));
        c1ComboBox.setSelectedIndex(UserVariables.getC1Index());
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 1;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_END;
        gridBagConstraints.weightx = 0.5;
        gridBagConstraints.weighty = 0.09;
        gridBagConstraints.insets = new java.awt.Insets(0, 0, 0, 10);
        jPanel1.add(c1ComboBox, gridBagConstraints);

        c2ComboBox.setModel(new javax.swing.DefaultComboBoxModel(UserVariables.channels));
        c2ComboBox.setSelectedIndex(UserVariables.getC2Index());
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 2;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_END;
        gridBagConstraints.weightx = 0.5;
        gridBagConstraints.weighty = 0.09;
        gridBagConstraints.insets = new java.awt.Insets(0, 0, 0, 10);
        jPanel1.add(c2ComboBox, gridBagConstraints);

        spatResLabel.setText("Spatial resolution ("+IJ.micronSymbol+"m/pixel):");
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 3;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_START;
        gridBagConstraints.weightx = 0.5;
        gridBagConstraints.weighty = 0.09;
        gridBagConstraints.insets = new java.awt.Insets(0, 10, 0, 0);
        jPanel1.add(spatResLabel, gridBagConstraints);

        timeResLabel.setText("Frames per second:");
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 4;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_START;
        gridBagConstraints.weightx = 0.5;
        gridBagConstraints.weighty = 0.09;
        gridBagConstraints.insets = new java.awt.Insets(0, 10, 0, 0);
        jPanel1.add(timeResLabel, gridBagConstraints);

        minTrajLengthLabel.setText("Minimum trajectory length (s):");
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 5;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_START;
        gridBagConstraints.weightx = 0.5;
        gridBagConstraints.weighty = 0.09;
        gridBagConstraints.insets = new java.awt.Insets(0, 10, 0, 0);
        jPanel1.add(minTrajLengthLabel, gridBagConstraints);

        maxTrajStepLabel.setText("Maximum linking distance:");
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 6;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_START;
        gridBagConstraints.weightx = 0.5;
        gridBagConstraints.weighty = 0.09;
        gridBagConstraints.insets = new java.awt.Insets(0, 10, 0, 0);
        jPanel1.add(maxTrajStepLabel, gridBagConstraints);

        chan1MaxThreshLabel.setText("Minimum peak size (C1):");
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 7;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_START;
        gridBagConstraints.weightx = 0.5;
        gridBagConstraints.weighty = 0.09;
        gridBagConstraints.insets = new java.awt.Insets(0, 10, 0, 0);
        jPanel1.add(chan1MaxThreshLabel, gridBagConstraints);

        chan2MaxThreshLabel.setText("Minimum peak size (C2):");
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 8;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_START;
        gridBagConstraints.weightx = 0.5;
        gridBagConstraints.weighty = 0.09;
        gridBagConstraints.insets = new java.awt.Insets(0, 10, 0, 0);
        jPanel1.add(chan2MaxThreshLabel, gridBagConstraints);

        spatResTextField.setText(String.valueOf(UserVariables.getSpatialRes()));
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 3;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_END;
        gridBagConstraints.weightx = 0.5;
        gridBagConstraints.weighty = 0.09;
        gridBagConstraints.insets = new java.awt.Insets(0, 0, 0, 10);
        jPanel1.add(spatResTextField, gridBagConstraints);

        timeResTextField.setText(String.valueOf(UserVariables.getTimeRes()));
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 4;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_END;
        gridBagConstraints.weightx = 0.5;
        gridBagConstraints.weighty = 0.09;
        gridBagConstraints.insets = new java.awt.Insets(0, 0, 0, 10);
        jPanel1.add(timeResTextField, gridBagConstraints);

        minTrajLengthTextField.setText(String.valueOf(UserVariables.getMinTrajLength()));
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 5;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_END;
        gridBagConstraints.weightx = 0.5;
        gridBagConstraints.weighty = 0.09;
        gridBagConstraints.insets = new java.awt.Insets(0, 0, 0, 10);
        jPanel1.add(minTrajLengthTextField, gridBagConstraints);

        maxTrajStepTextField.setText(String.valueOf(UserVariables.getTrajMaxStep()));
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 6;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_END;
        gridBagConstraints.weightx = 0.5;
        gridBagConstraints.weighty = 0.09;
        gridBagConstraints.insets = new java.awt.Insets(0, 0, 0, 10);
        jPanel1.add(maxTrajStepTextField, gridBagConstraints);

        chan1MaxThreshTextField.setText(String.valueOf(UserVariables.getChan1MaxThresh()));
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 7;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_END;
        gridBagConstraints.weightx = 0.5;
        gridBagConstraints.weighty = 0.09;
        gridBagConstraints.insets = new java.awt.Insets(0, 0, 0, 10);
        jPanel1.add(chan1MaxThreshTextField, gridBagConstraints);

        chan2MaxThreshTextField.setText(String.valueOf(UserVariables.getChan2MaxThresh()));
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 8;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_END;
        gridBagConstraints.weightx = 0.5;
        gridBagConstraints.weighty = 0.09;
        gridBagConstraints.insets = new java.awt.Insets(0, 0, 0, 10);
        jPanel1.add(chan2MaxThreshTextField, gridBagConstraints);

        colocaliseToggleButton.setText("Co-Localised Only");
        colocaliseToggleButton.setSelected(UserVariables.isColocal());
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 11;
        gridBagConstraints.gridwidth = 2;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.weightx = 1.0;
        gridBagConstraints.weighty = 0.09;
        gridBagConstraints.insets = new java.awt.Insets(0, 10, 0, 10);
        jPanel1.add(colocaliseToggleButton, gridBagConstraints);

        preprocessToggleButton.setText("Pre-Process Images");
        preprocessToggleButton.setSelected(UserVariables.isPreProcess());
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 12;
        gridBagConstraints.gridwidth = 2;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.weightx = 1.0;
        gridBagConstraints.weighty = 0.09;
        gridBagConstraints.insets = new java.awt.Insets(0, 10, 10, 10);
        jPanel1.add(preprocessToggleButton, gridBagConstraints);

        curveFitTolLabel.setText("Curve fit tolerance:");
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 9;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_START;
        gridBagConstraints.weightx = 0.5;
        gridBagConstraints.weighty = 0.09;
        gridBagConstraints.insets = new java.awt.Insets(0, 10, 0, 0);
        jPanel1.add(curveFitTolLabel, gridBagConstraints);

        curveFitTolTextField.setText(String.valueOf(UserVariables.getCurveFitTol()));
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 9;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_END;
        gridBagConstraints.weightx = 0.5;
        gridBagConstraints.weighty = 0.09;
        gridBagConstraints.insets = new java.awt.Insets(0, 0, 0, 10);
        jPanel1.add(curveFitTolTextField, gridBagConstraints);

        nMaxTextField.setText(String.valueOf(UserVariables.getnMax()));
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 10;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_END;
        gridBagConstraints.weightx = 0.5;
        gridBagConstraints.weighty = 0.09;
        gridBagConstraints.insets = new java.awt.Insets(0, 0, 0, 10);
        jPanel1.add(nMaxTextField, gridBagConstraints);

        nMaxLabel.setText("Nmax:");
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 10;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_START;
        gridBagConstraints.weightx = 0.5;
        gridBagConstraints.weighty = 0.09;
        gridBagConstraints.insets = new java.awt.Insets(0, 10, 0, 0);
        jPanel1.add(nMaxLabel, gridBagConstraints);

        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 0;
        gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.NORTHWEST;
        gridBagConstraints.weightx = 0.2;
        gridBagConstraints.weighty = 0.8;
        getContentPane().add(jPanel1, gridBagConstraints);

        jPanel2.setBorder(javax.swing.BorderFactory.createEtchedBorder());
        jPanel2.setLayout(new java.awt.GridBagLayout());

        canvas1.setPreferredSize(new java.awt.Dimension(
            imp.getProcessor().getWidth(), imp.getProcessor().getHeight()));
    gridBagConstraints = new java.awt.GridBagConstraints();
    gridBagConstraints.gridx = 0;
    gridBagConstraints.gridy = 0;
    gridBagConstraints.gridwidth = 2;
    gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
    gridBagConstraints.weightx = 1.0;
    gridBagConstraints.weighty = 0.8;
    gridBagConstraints.insets = new java.awt.Insets(10, 10, 0, 10);
    jPanel2.add(canvas1, gridBagConstraints);

    previewScrollBar.setOrientation(javax.swing.JScrollBar.HORIZONTAL);
    previewScrollBar.setModel(new DefaultBoundedRangeModel(1, 0, 1, analyser.getStack().getSize()));
    previewScrollBar.setEnabled(previewToggleButton.isSelected());
    previewScrollBar.addAdjustmentListener(new java.awt.event.AdjustmentListener() {
        public void adjustmentValueChanged(java.awt.event.AdjustmentEvent evt) {
            previewScrollBarAdjustmentValueChanged(evt);
        }
    });
    gridBagConstraints = new java.awt.GridBagConstraints();
    gridBagConstraints.gridx = 0;
    gridBagConstraints.gridy = 2;
    gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
    gridBagConstraints.anchor = java.awt.GridBagConstraints.WEST;
    gridBagConstraints.weightx = 0.8;
    gridBagConstraints.weighty = 0.1;
    gridBagConstraints.insets = new java.awt.Insets(0, 10, 10, 0);
    jPanel2.add(previewScrollBar, gridBagConstraints);

    previewTextField.setText(String.valueOf(previewScrollBar.getValue()));
    previewTextField.setEditable(false);
    previewTextField.setEnabled(previewToggleButton.isSelected());
    gridBagConstraints = new java.awt.GridBagConstraints();
    gridBagConstraints.gridx = 1;
    gridBagConstraints.gridy = 2;
    gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
    gridBagConstraints.anchor = java.awt.GridBagConstraints.EAST;
    gridBagConstraints.weightx = 0.2;
    gridBagConstraints.weighty = 0.1;
    gridBagConstraints.insets = new java.awt.Insets(0, 0, 10, 10);
    jPanel2.add(previewTextField, gridBagConstraints);

    previewToggleButton.setText("Preview");
    previewToggleButton.addActionListener(new java.awt.event.ActionListener() {
        public void actionPerformed(java.awt.event.ActionEvent evt) {
            previewToggleButtonActionPerformed(evt);
        }
    });
    gridBagConstraints = new java.awt.GridBagConstraints();
    gridBagConstraints.gridx = 0;
    gridBagConstraints.gridy = 1;
    gridBagConstraints.gridwidth = 2;
    gridBagConstraints.weighty = 0.1;
    gridBagConstraints.insets = new java.awt.Insets(10, 10, 10, 10);
    jPanel2.add(previewToggleButton, gridBagConstraints);

    gridBagConstraints = new java.awt.GridBagConstraints();
    gridBagConstraints.gridx = 1;
    gridBagConstraints.gridy = 0;
    gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
    gridBagConstraints.anchor = java.awt.GridBagConstraints.NORTHWEST;
    gridBagConstraints.weightx = 0.8;
    gridBagConstraints.weighty = 0.8;
    getContentPane().add(jPanel2, gridBagConstraints);

    jPanel3.setBorder(javax.swing.BorderFactory.createEtchedBorder());
    jPanel3.setLayout(new java.awt.GridBagLayout());

    okButton.setText("Run");
    okButton.addActionListener(new java.awt.event.ActionListener() {
        public void actionPerformed(java.awt.event.ActionEvent evt) {
            okButtonActionPerformed(evt);
        }
    });
    gridBagConstraints = new java.awt.GridBagConstraints();
    gridBagConstraints.gridx = 0;
    gridBagConstraints.gridy = 0;
    gridBagConstraints.anchor = java.awt.GridBagConstraints.NORTHWEST;
    jPanel3.add(okButton, gridBagConstraints);

    cancelButton.setText("Cancel");
    cancelButton.addActionListener(new java.awt.event.ActionListener() {
        public void actionPerformed(java.awt.event.ActionEvent evt) {
            cancelButtonActionPerformed(evt);
        }
    });
    gridBagConstraints = new java.awt.GridBagConstraints();
    gridBagConstraints.gridx = 1;
    gridBagConstraints.gridy = 0;
    gridBagConstraints.anchor = java.awt.GridBagConstraints.NORTHWEST;
    jPanel3.add(cancelButton, gridBagConstraints);

    gridBagConstraints = new java.awt.GridBagConstraints();
    gridBagConstraints.gridx = 0;
    gridBagConstraints.gridy = 1;
    gridBagConstraints.gridwidth = 2;
    gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
    gridBagConstraints.weighty = 0.2;
    getContentPane().add(jPanel3, gridBagConstraints);

    pack();
    }// </editor-fold>//GEN-END:initComponents

    private void cancelButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_cancelButtonActionPerformed
        this.dispose();
        wasOKed = false;
    }//GEN-LAST:event_cancelButtonActionPerformed

    private void okButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_okButtonActionPerformed
        if (!setVariables()) {
            return;
        }
        this.dispose();
        wasOKed = true;
    }//GEN-LAST:event_okButtonActionPerformed

    private void previewScrollBarAdjustmentValueChanged(java.awt.event.AdjustmentEvent evt) {//GEN-FIRST:event_previewScrollBarAdjustmentValueChanged
        previewTextField.setText(String.valueOf(previewScrollBar.getValue()));
        if (previewScrollBar.getValueIsAdjusting() || !setVariables()) {
            return;
        }
        viewDetections();
    }//GEN-LAST:event_previewScrollBarAdjustmentValueChanged

    private void previewToggleButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_previewToggleButtonActionPerformed
        previewScrollBar.setEnabled(previewToggleButton.isSelected());
        previewTextField.setEnabled(previewToggleButton.isSelected());
        previewScrollBarAdjustmentValueChanged(null);
    }//GEN-LAST:event_previewToggleButtonActionPerformed

    boolean setVariables() {
        try {
            UserVariables.setChan1MaxThresh(Double.parseDouble(chan1MaxThreshTextField.getText()));
            UserVariables.setChan2MaxThresh(Double.parseDouble(chan2MaxThreshTextField.getText()));
            UserVariables.setMinTrajLength(Double.parseDouble(minTrajLengthTextField.getText()));
            UserVariables.setSpatialRes(Double.parseDouble(spatResTextField.getText()));
            UserVariables.setTimeRes(Double.parseDouble(timeResTextField.getText()));
            UserVariables.setTrajMaxStep(Double.parseDouble(maxTrajStepTextField.getText()));
            UserVariables.setCurveFitTol(Double.parseDouble(curveFitTolTextField.getText()));
            UserVariables.setColocal(colocaliseToggleButton.isSelected());
            UserVariables.setPreProcess(preprocessToggleButton.isSelected());
            UserVariables.setC1Index(c1ComboBox.getSelectedIndex());
            UserVariables.setC2Index(c2ComboBox.getSelectedIndex());
            UserVariables.setnMax(Integer.parseInt(nMaxTextField.getText()));
        } catch (NumberFormatException e) {
            IJ.error("Number formatting error " + e.toString());
            return false;
        }
        return true;
    }

    public void viewDetections() {
        analyser.calcParticleRadius(UserVariables.getSpatialRes());
        ParticleArray detections = analyser.findParticles(1.0, false, previewScrollBar.getValue()-1, previewScrollBar.getValue()-1, 0.0);
        ImageProcessor output;
        int slice = previewScrollBar.getValue();
        ImageStack stack = analyser.getStack();
        boolean monoChrome = analyser.isMonoChrome();
        if (analyser.isMonoChrome()) {
            output = (new TypeConverter((stack.getProcessor(slice)).duplicate(), true)).convertToByte();
        } else {
            output = (new TypeConverter((stack.getProcessor(slice)).duplicate(), true)).convertToRGB();
        }
        double mag = 1.0 / UIMethods.getMagnification(output, canvas1);
        double sr = 1.0 / Double.parseDouble(spatResTextField.getText());
//        int radius = (int)Math.round(sr);
        int radius = analyser.getXyPartRad();
        IsoGaussian c1, c2;
        ArrayList<Particle> particles = detections.getLevel(0);
        Color c1Color = !monoChrome ? analyser.getDrawColor(c1ComboBox.getSelectedIndex()) : Color.white;
        Color c2Color = !monoChrome ? analyser.getDrawColor(c2ComboBox.getSelectedIndex()) : Color.white;
        output.setLineWidth(1);
        for (Particle particle : particles) {
            c1 = particle.getC1Gaussian();
            c2 = particle.getC2Gaussian();
            drawDetections(output, (int) (Math.round(sr * c1.getX())), (int) (Math.round(sr * c1.getY())),
                    radius, c1.getFit() > UserVariables.getCurveFitTol(), c1Color);
            if (c2 != null) {
                drawDetections(output, (int) (Math.round(sr * c2.getX())),
                        (int) (Math.round(sr * c2.getY())), radius,
                        false, c2Color);
            }
        }
        imp.setProcessor("", output);
        ((ImageCanvas) canvas1).setMagnification(mag);

        canvas1.repaint();
    }

    public void drawDetections(ImageProcessor image, int x, int y, int radius,
            boolean drawOval, Color colour) {
        image.setColor(colour);
        if (drawOval) {
            image.drawOval((x - radius), (y - radius), 2 * radius, 2 * radius);
        } else {
            image.drawRect((x - radius), (y - radius), 2 * radius, 2 * radius);
        }
    }

    public byte[] getPixels(int channel, int frame) {
        ColorProcessor processor = (ColorProcessor) analyser.getStack().getProcessor(frame + 1);
        int width = processor.getWidth(), height = processor.getHeight();
        int size = width * height;
        byte redPix[] = new byte[size], greenPix[] = new byte[size],
                bluePix[] = new byte[size], emptyPix[] = new byte[size];
        Arrays.fill(emptyPix, (byte) 0);
        processor.getRGB(redPix, greenPix, bluePix);
        switch (channel) {
            case UserVariables.RED:
                return redPix;
            case UserVariables.GREEN:
                return greenPix;
            case UserVariables.BLUE:
                return bluePix;
            default:
                return emptyPix;
        }
    }

    public boolean isWasOKed() {
        return wasOKed;
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JComboBox c1ComboBox;
    private javax.swing.JLabel c1Label;
    private javax.swing.JComboBox c2ComboBox;
    private javax.swing.JLabel c2Label;
    private javax.swing.JButton cancelButton;
    private java.awt.Canvas canvas1;
    private javax.swing.JLabel chan1MaxThreshLabel;
    private javax.swing.JTextField chan1MaxThreshTextField;
    private javax.swing.JLabel chan2MaxThreshLabel;
    private javax.swing.JTextField chan2MaxThreshTextField;
    private javax.swing.JToggleButton colocaliseToggleButton;
    private javax.swing.JLabel curveFitTolLabel;
    private javax.swing.JTextField curveFitTolTextField;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JPanel jPanel3;
    private javax.swing.JLabel maxTrajStepLabel;
    private javax.swing.JTextField maxTrajStepTextField;
    private javax.swing.JLabel minTrajLengthLabel;
    private javax.swing.JTextField minTrajLengthTextField;
    private javax.swing.JLabel nMaxLabel;
    private javax.swing.JTextField nMaxTextField;
    private javax.swing.JButton okButton;
    private javax.swing.JToggleButton preprocessToggleButton;
    private javax.swing.JScrollBar previewScrollBar;
    private javax.swing.JTextField previewTextField;
    private javax.swing.JToggleButton previewToggleButton;
    private javax.swing.JLabel spatResLabel;
    private javax.swing.JTextField spatResTextField;
    private javax.swing.JLabel timeResLabel;
    private javax.swing.JTextField timeResTextField;
    // End of variables declaration//GEN-END:variables
}
