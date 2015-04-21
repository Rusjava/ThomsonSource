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

import java.awt.Color;
import java.awt.Paint;
import java.io.Serializable;
import org.jfree.util.PublicCloneable;
import org.jfree.chart.renderer.PaintScale;

/**
 * A jet paint scale 
 * @author Ruslan Feshchenko
 * @version 1.0
 */

  public class JetPaintScale 
          implements PaintScale, PublicCloneable, Serializable {
 
       /** The lower bound. */
      private double lowerBound;
      
       /** The upper bound. */
      private double upperBound;
      
       public JetPaintScale() {
           this(0.0, 1.0);
      }
       
       public JetPaintScale(double lowerBound, double upperBound) {
          if (lowerBound >= upperBound) {
              throw new IllegalArgumentException(
                       "Requires lowerBound < upperBound.");
           }
           this.lowerBound = lowerBound;
           this.upperBound = upperBound;
       }
       
      @Override
       public double getLowerBound() {
           return this.lowerBound;
       }
   
      @Override
      public double getUpperBound() {
          return this.upperBound;
      }
  
      @Override
      public Paint getPaint(double value) {
          double v = Math.min(Math.max(value, this.lowerBound), this.upperBound);
          float h = (float)(1.0f-(v - this.lowerBound) / (this.upperBound 
                  - this.lowerBound));
          return Color.getHSBColor(h, 1f, 0.8f);
      }
      
      @Override
      public boolean equals(Object obj) {
          if (obj == this) {
              return true;
          }
          
          JetPaintScale that = (JetPaintScale) obj;
          
          if (!(obj instanceof JetPaintScale) | (this.lowerBound != that.lowerBound) | (this.upperBound != that.upperBound)) {
              return false;
          }
          
          return true;    
     }
 
      @Override
      public Object clone() throws CloneNotSupportedException {
          return super.clone();
      }
      
  }
