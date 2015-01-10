/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
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
          double v = Math.max(value, this.lowerBound);
          v = Math.min(v, this.upperBound);
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
