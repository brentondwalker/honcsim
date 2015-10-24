/*
 * TrippleBubbleBox.java
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
 * The example of three connected bubble boxes where RC(v) has
 * a hole for every possible v in V
 */

import java.util.Vector;

import m4rjni.Mzd;

public class TrippleBubbleBox extends CoverageExperiment{
	/*
	 * class data
	 */
	static final int vsDimension = 2;
	
	Mzd v1 = null;
	Mzd v2 = null;
	Mzd v3 = null;
	
	/**
	 * Constructor
	 * 
	 */
	public TrippleBubbleBox() {
		// call CoverageExperiment constructor
		super(vsDimension, 1.0);
		
		// allocate the grid of DPoints
		points = new Vector<DPoint>(9);
		
		v1 = new Mzd(this.basis[0]);
		v2 = new Mzd(this.basis[1]);
		v3 = Mzd.add(v1, v2);
		
		placeBubbleBox(0, 0, v1, v2);
		placeBubbleBox(0.5, 1.75, v3, v1);
		placeBubbleBox(1.75, -0.5, v2, v3);
		
		// compute the neighbor set
		// in our world the underlying set of points doesn't change,
		// so we compute the neighbor set in the constructor
		this.computeNeighborSet();
	}

	private void placeBubbleBox(double x, double y, Mzd v1, Mzd v2) {
		// add the points on the inner square
		double s = 1.0/Math.sqrt(2.0)/2.0;
		double s2 = s-0.11;
		points.add(new DPoint(x+s2, y+s2, vsDimension));
		points.add(new DPoint(x+s2, y-s2, vsDimension));
		points.add(new DPoint(x-s2, y-s2, vsDimension));
		points.add(new DPoint(x-s2, y+s2, vsDimension));
		
		// add the two outer points with vector e1
		//double s2p = s + 0.12;
		double s2p = 0.5 - 0.0000001;
		points.add(new DPoint(x+s2p, y+s2p, vsDimension));
		points.lastElement().addInventoryVector(v1);
		points.add(new DPoint(x-s2p, y-s2p, vsDimension));
		points.lastElement().addInventoryVector(v1);
		
		// add the other two outer points
		points.add(new DPoint(x-s2p, y+s2p, vsDimension));
		points.lastElement().addInventoryVector(v2);
		points.add(new DPoint(x+s2p, y-s2p, vsDimension));
		points.lastElement().addInventoryVector(v2);
	}
	
	/* **************************************
	 *
	 * Main Routine
	 *
	 * **************************************
	 */
	public static void main(String[] args) {
		
		TrippleBubbleBox tbb = new TrippleBubbleBox();
	    
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
