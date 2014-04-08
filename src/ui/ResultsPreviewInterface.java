/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ui;

import ParticleTracking.ParticleTrajectory;
import ParticleTracking.Timelapse_Analysis;
import ParticleTracking.UserVariables;
import UIClasses.UIMethods;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.ImageCanvas;
import java.util.ArrayList;
import javax.swing.DefaultBoundedRangeModel;

/**
 *
 * @author David Barry <david.barry at cancer.org.uk>
 */
public class ResultsPreviewInterface extends javax.swing.JDialog {

    String title;
    Timelapse_Analysis analyser;
    ImagePlus imp;
    ImageStack stack;
    ArrayList<ParticleTrajectory> trajectories;
    boolean wasOKed;

    /**
     * Creates new form ResultsPreviewer
     */
    public ResultsPreviewInterface(java.awt.Frame parent, boolean modal, String title, Timelapse_Analysis analyser) {
        super(parent, modal);
        this.title = title;
        this.analyser = analyser;
        stack = analyser.getStack();
        imp = new ImagePlus("", stack.getProcessor(1));
        trajectories = analyser.getTrajectories();
        initComponents();
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

        canvas = new ImageCanvas(imp);
        trajScrollBar = new javax.swing.JScrollBar();
        trajTextField = new javax.swing.JTextField();
        removeButton = new javax.swing.JButton();
        jPanel1 = new javax.swing.JPanel();
        okButton = new javax.swing.JButton();
        cancelButton = new javax.swing.JButton();
        imageScrollBar = new javax.swing.JScrollBar();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle(title);
        getContentPane().setLayout(new java.awt.GridBagLayout());

        canvas.setPreferredSize(new java.awt.Dimension(256, 256));
        canvas.addComponentListener(new java.awt.event.ComponentAdapter() {
            public void componentResized(java.awt.event.ComponentEvent evt) {
                canvasComponentResized(evt);
            }
        });
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 0;
        gridBagConstraints.gridwidth = 2;
        gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints.weightx = 0.6;
        gridBagConstraints.weighty = 0.6;
        gridBagConstraints.insets = new java.awt.Insets(10, 10, 0, 10);
        getContentPane().add(canvas, gridBagConstraints);

        trajScrollBar.setOrientation(javax.swing.JScrollBar.HORIZONTAL);
        trajScrollBar.setModel(new DefaultBoundedRangeModel(1, 0, 0, trajectories.size()-1));
        trajScrollBar.addAdjustmentListener(new java.awt.event.AdjustmentListener() {
            public void adjustmentValueChanged(java.awt.event.AdjustmentEvent evt) {
                trajScrollBarAdjustmentValueChanged(evt);
            }
        });
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 3;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.WEST;
        gridBagConstraints.weightx = 0.8;
        gridBagConstraints.weighty = 0.1;
        gridBagConstraints.insets = new java.awt.Insets(10, 10, 10, 10);
        getContentPane().add(trajScrollBar, gridBagConstraints);

        trajTextField.setEditable(false);
        trajTextField.setText(String.valueOf(trajScrollBar.getValue()));
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 3;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.EAST;
        gridBagConstraints.weightx = 0.2;
        gridBagConstraints.weighty = 0.1;
        gridBagConstraints.insets = new java.awt.Insets(10, 10, 10, 10);
        getContentPane().add(trajTextField, gridBagConstraints);

        removeButton.setText("Remove");
        removeButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                removeButtonActionPerformed(evt);
            }
        });
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 2;
        gridBagConstraints.gridwidth = 2;
        gridBagConstraints.weighty = 0.1;
        gridBagConstraints.insets = new java.awt.Insets(10, 10, 10, 10);
        getContentPane().add(removeButton, gridBagConstraints);

        jPanel1.setLayout(new java.awt.GridBagLayout());

        okButton.setText("OK");
        okButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                okButtonActionPerformed(evt);
            }
        });
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 0;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.NORTHWEST;
        jPanel1.add(okButton, gridBagConstraints);

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
        jPanel1.add(cancelButton, gridBagConstraints);

        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 4;
        gridBagConstraints.gridwidth = 2;
        getContentPane().add(jPanel1, gridBagConstraints);

        imageScrollBar.setOrientation(javax.swing.JScrollBar.HORIZONTAL);
        imageScrollBar.setModel(new DefaultBoundedRangeModel(1, 0, 1, stack.getSize()));
        imageScrollBar.addAdjustmentListener(new java.awt.event.AdjustmentListener() {
            public void adjustmentValueChanged(java.awt.event.AdjustmentEvent evt) {
                imageScrollBarAdjustmentValueChanged(evt);
            }
        });
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 1;
        gridBagConstraints.gridwidth = 2;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.weighty = 0.1;
        gridBagConstraints.insets = new java.awt.Insets(0, 10, 10, 10);
        getContentPane().add(imageScrollBar, gridBagConstraints);

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void imageScrollBarAdjustmentValueChanged(java.awt.event.AdjustmentEvent evt) {//GEN-FIRST:event_imageScrollBarAdjustmentValueChanged
        imp.setProcessor(stack.getProcessor(imageScrollBar.getValue()));
        canvas.repaint();
    }//GEN-LAST:event_imageScrollBarAdjustmentValueChanged

    private void trajScrollBarAdjustmentValueChanged(java.awt.event.AdjustmentEvent evt) {//GEN-FIRST:event_trajScrollBarAdjustmentValueChanged
        stack = analyser.mapTrajectories(analyser.getStack(), trajectories, 1.0, UserVariables.getSpatialRes(), UserVariables.getMinTrajLength(), UserVariables.getTimeRes(), true, (int) Math.round(1.0 / UserVariables.getSpatialRes()), trajScrollBar.getValue(), trajScrollBar.getValue());
        imageScrollBarAdjustmentValueChanged(null);
    }//GEN-LAST:event_trajScrollBarAdjustmentValueChanged

    private void removeButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_removeButtonActionPerformed
//        analyser.removeTrajectory(trajScrollBar.getValue());
    }//GEN-LAST:event_removeButtonActionPerformed

    private void okButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_okButtonActionPerformed
        this.dispose();
        wasOKed = true;
    }//GEN-LAST:event_okButtonActionPerformed

    private void cancelButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_cancelButtonActionPerformed
        this.dispose();
        wasOKed = false;
    }//GEN-LAST:event_cancelButtonActionPerformed

    private void canvasComponentResized(java.awt.event.ComponentEvent evt) {//GEN-FIRST:event_canvasComponentResized
        ((ImageCanvas) canvas).setMagnification(1.0 / UIMethods.getMagnification(imp.getProcessor(), canvas));
    }//GEN-LAST:event_canvasComponentResized

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton cancelButton;
    private java.awt.Canvas canvas;
    private javax.swing.JScrollBar imageScrollBar;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JButton okButton;
    private javax.swing.JButton removeButton;
    private javax.swing.JScrollBar trajScrollBar;
    private javax.swing.JTextField trajTextField;
    // End of variables declaration//GEN-END:variables
}
