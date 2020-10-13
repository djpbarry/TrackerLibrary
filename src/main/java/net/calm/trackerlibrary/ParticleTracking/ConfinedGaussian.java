/*
 * Copyright (C) 2017 Dave Barry <david.barry at crick.ac.uk>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
package net.calm.trackerlibrary.ParticleTracking;

import ij.gui.Roi;
import java.awt.geom.Rectangle2D;
import java.util.Random;

/**
 *
 * @author Dave Barry <david.barry at crick.ac.uk>
 */
public class ConfinedGaussian extends MotileGaussian {

    private Roi roi;
    private double boundsRadius;
    private double roiPostIncX, roiPosIncY;

    public ConfinedGaussian(double x0, double y0, double a, double xsig, double ysig,
            double fit, double sens, boolean persistent, boolean changeState, double boundsRadius) {
        super(x0, y0, a, xsig, ysig, fit, sens, persistent, changeState, 0.001, 0.0);
        this.boundsRadius = boundsRadius;
        this.roi = new Roi(x0 - boundsRadius, y0 - boundsRadius, 2.0 * boundsRadius + 1.0, 2.0 * boundsRadius + 1.0);
        Random r = new Random();
        roiPostIncX = 0.133 / 10.0;
        roiPosIncY = 0.133 / 10.0;
    }

    public void updatePosition(double weighting) {
        this.x += rad * Math.cos(theta);
        this.y += rad * Math.sin(theta);
        Rectangle2D.Double bounds = roi.getFloatBounds();
        if (this.x < bounds.x) {
            this.x = bounds.x;
        }
        if (this.x > bounds.x + bounds.width) {
            this.x = bounds.x + bounds.width;
        }
        if (this.y < bounds.y) {
            this.y = bounds.y;
        }
        if (this.y > bounds.y + bounds.height) {
            this.y = bounds.y + bounds.height;
        }
        roi = new Roi(bounds.x + roiPostIncX, bounds.y + roiPosIncY, bounds.width, bounds.height);
    }
}
