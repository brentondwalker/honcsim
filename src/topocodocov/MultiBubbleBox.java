/*
 * MultiBubbleBox.java
 *
 * Copyright (C) 2015 Brenton Walker
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package topocodocov;
/*
 * Experiment with arrangements of lots of bubbleboxes
 * 
 */

import java.util.Vector;

import m4rjni.Mzd;

public class MultiBubbleBox extends DrainageExperiment {

	/*
	 * class data
	 */
	static final int vsDimension = 2;
	
	Mzd v1 = null;
	Mzd v2 = null;
	Mzd v3 = null;
	
	// edge length of outer square
	double so = 0.0;
	
	// the horizontal/vertical offset to the inner square
	double aa = 0.0;
	
	// fundamental cell size of bubble box, or something relayed.
	// actually radius+spc is the fundamental cell
	double spc = 0.0;
	
	// edge length of inner square
	double si = 0.0;
	
	/**
	 * Constructor
	 * 
	 */
	public MultiBubbleBox(int m, int n) {
		// call CoverageExperiment constructor
		super(vsDimension, 1.0);
		
		// allocate the grid of DPoints
		points = new Vector<DPoint>(9);
		
		v1 = new Mzd(this.basis[0]);
		v2 = new Mzd(this.basis[1]);
		v3 = Mzd.add(v1, v2);
		
		so = radius - 0.000001;
		aa = radius*(1.0-1.0/Math.sqrt(2.0)) - 0.00001;
		si = so - 2.0*aa;
		//spc = 1.0/Math.sqrt(2.0)-0.00001;
		spc = Math.sqrt(3.0)/2.0-0.00001;
		
		for (int i=0; i<m; i++) {
			for (int j=0; j<n; j++) {
				placeBubbleBox(j*(so+spc)-i*so/2, i*(so+spc)+j*so/2, v1, v2);
				if ((i<(m-1)) && (j<(n-1))) {
					placeInnerBubbleBox(j*(so+spc)-i*so/2+so, i*(so+spc)+j*so/2+so);
				}
				//placeBubbleBox(so+spc, so/2, v1, v2);
				//placeBubbleBox(-so/2, so+spc, v1, v2);
				//placeBubbleBox(so+spc-so/2, so/2+so+spc, v1, v2);
				//placeInnerBubbleBox(so,so);
			}
		}
		
		// compute the neighbor set
		// in our world the underlying set of points doesn't change,
		// so we compute the neighbor set in the constructor
		this.computeNeighborSet();
	}
	
	
	/**
	 * Place a lone inner bubble box, rotated by 30 degrees counterclockwise.
	 * The x,y passed in should be to the lower left corner of the outer square of the bubble box.
	 * This method will compute its own offset into the outer box.
	 * 
	 * @param x
	 * @param y
	 */
	public void placeInnerBubbleBox(double x, double y) {
		// the offset to the lower corner of a tilted inner square
		double ax = aa*Math.sqrt(2.0)*Math.cos(Math.toRadians(75));
		double ay = aa*Math.sqrt(2.0)*Math.sin(Math.toRadians(75));
		
		points.add(new DPoint(x+ax, y+ay, vsDimension));
		points.add(new DPoint(x+ax+si*Math.cos(Math.toRadians(30)), y+ay+si*Math.sin(Math.toRadians(30)), vsDimension));
		points.add(new DPoint(x+ax+si*(Math.cos(Math.toRadians(30))+Math.cos(Math.toRadians(120))), y+ay+si*(Math.sin(Math.toRadians(30))+Math.sin(Math.toRadians(120))), vsDimension));
		points.add(new DPoint(x+ax+si*Math.cos(Math.toRadians(120)), y+ay+si*Math.sin(Math.toRadians(120)), vsDimension));
	}
	
	
	/**
	 * Place a bubble box at the coordinates indicated, and put vectors v1 and v2 
	 * in the enventory of opposite corners of the outer box.
	 * 
	 * @param x
	 * @param y
	 * @param v1
	 * @param v2
	 */
	public void placeBubbleBox(double x, double y, Mzd v1, Mzd v2) {
		// add the points on the inner square
		points.add(new DPoint(x+aa, y+aa, vsDimension));
		points.add(new DPoint(x+aa, y+so-aa, vsDimension));
		points.add(new DPoint(x+so-aa, y+so-aa, vsDimension));
		points.add(new DPoint(x+so-aa, y+aa, vsDimension));
		
		// add the two outer points with vector e1
		points.add(new DPoint(x, y, vsDimension));
		points.lastElement().addInventoryVector(v1);
		points.add(new DPoint(x+so, y+so, vsDimension));
		points.lastElement().addInventoryVector(v1);
		
		// add the other two outer points
		points.add(new DPoint(x, y+so, vsDimension));
		points.lastElement().addInventoryVector(v2);
		points.add(new DPoint(x+so, y, vsDimension));
		points.lastElement().addInventoryVector(v2);
	}
	
	
	/* **************************************
	 *
	 * Main Routine
	 *
	 * **************************************
	 */
	public static void main(String[] args) {
		
		if (args.length != 4) {
			System.out.println("usage: MultiBubbleBox <grid width> <grid height> <vsDim> <vectors-per-node>\n");
			System.exit(0);
		}
		int gridWidth = Integer.parseInt(args[0]);
		int gridHeight = Integer.parseInt(args[1]);
		
		MultiBubbleBox tbb = new MultiBubbleBox(gridWidth,gridHeight);
	    
		tbb.buildRipsComplex();
		tbb.buildCoverageComplex(tbb.v1);
	    tbb.computeHomology();
	    tbb.drawComplex();
	    
		tbb.buildCoverageComplex(tbb.v2);
	    tbb.computeHomology();
	    tbb.drawComplex();
	    
		tbb.buildCoverageComplex(tbb.v3);
	    tbb.computeHomology();
	    tbb.drawComplex();
	}
	
}
