/*
 * BubbleBox.java
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
 * The bubble box example
 */

import java.util.*;

import m4rjni.Mzd;
//import java.lang.Math.*;

public class BubbleBox extends CoverageExperiment {

	/*
	 * class data
	 */
	static final int vsDimension = 2;
	
	/**
	 * Constructor
	 * 
	 * @param width
	 * @param height
	 * @param spacing
	 */
	public BubbleBox() {
		// call CoverageExperiment constructor
		super(vsDimension, 1.0);
		
		// allocate the grid of DPoints
		points = new Vector<DPoint>(8);
		
		// add the points on the inner square
		double s2 = 1.0/Math.sqrt(2.0)/2.0;
		points.add(new DPoint(s2, s2, vsDimension));
		points.add(new DPoint(s2, -s2, vsDimension));
		points.add(new DPoint(-s2, -s2, vsDimension));
		points.add(new DPoint(-s2, s2, vsDimension));
		
		// add the two outer points with vector e1
		double s2p = s2 + 0.1;
		Mzd[] basis = Mzd.standardBasis(vsDimension);
		DPoint p = null;
		p = new DPoint(s2p, s2p, vsDimension);
		p.addInventoryVector(basis[0]);
		points.add(p);
		p = new DPoint(-s2p, -s2p, vsDimension);
		p.addInventoryVector(basis[0]);
		points.add(p);
		
		// add the other two outer points
		p = new DPoint(-s2p, s2p, vsDimension);
		p.addInventoryVector(basis[1]);
		points.add(p);
		p = new DPoint(s2p, -s2p, vsDimension);
		p.addInventoryVector(basis[1]);
		points.add(p);
		
		// compute the neighbor set
		this.computeNeighborSet();
		
		// build the Rips complex
		// In our case the rips complex is a javaplex object, and javaplex uses streams
		this.buildRipsComplex();
		System.out.println("Built Rips complex with "+ripsComplexStream.getSize()+" faces");
		
		// construct the standard basis for the full vector space V
		Mzd V = new Mzd(vsDimension, vsDimension);
		for (int i=0; i<vsDimension; i++) {
			V.writeBit(i, i, 1);
		}
		
		// build the filtered RC-->R complex
		//this.buildCoverageRipsComplex(V);
		//System.out.println("Built filtered coverage-->Rips complex with "+coverageRipsComplexStream.getSize()+" faces");
		
		// build the coverage complex
		this.buildCoverageComplex(V);
		System.out.println("Built coverage complex with "+coverageComplexStream.getSize()+" faces");
		
	}

	/* **************************************
	 *
	 * Main Routine
	 *
	 * **************************************
	 */
	public static void main(String[] args) {
		
		BubbleBox bb = new BubbleBox();
	    
	    //Vector<ExplicitSimplexStream> cycleStreams = g.computeHomology();
	    bb.computeHomology();
	    
	    //Vector<ExplicitSimplexStream> coverageRipsCycleStreams = g.computePersistentHomology();
	    
	    // view the complexes
	    //for (ExplicitSimplexStream str : cycleStreams) {
	    //	g.drawComplex(str);
	    //}
	    bb.drawComplex();
	}
	
}
