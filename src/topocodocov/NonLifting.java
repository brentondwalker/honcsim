/*
 * NonLifting.java
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
 * counterexample to the general existence of liftings in the coverage complex
 */

import java.util.*;

import edu.stanford.math.plex4.streams.impl.ExplicitSimplexStream;

public class NonLifting extends CoverageExperiment {

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
	public NonLifting() {
		// call CoverageExperiment constructor
		super(vsDimension, 1.0);

		// allocate the storage object for DPoints
		points = new Vector<DPoint>(8);
		
		// inner tetrahedron
		double t = 0.3;
		for (int i=0; i<4; i++) {
			DPoint p = new DPoint(t*Math.cos(i*Math.PI/2.0), t*Math.sin(i*Math.PI/2.0), vsDimension);
			points.add(p);
		}
		
		// outer helper points
		double tz = 0.8;
		for (int i=0; i<4; i++) {
			DPoint p = new DPoint(tz*Math.cos(i*Math.PI/2.0), tz*Math.sin(i*Math.PI/2.0), vsDimension);
			p.addInventoryVector(this.basis[i<2 ? 0 : 1]);
			points.add(p);
		}		
		
		// compute the neighbor set
		this.computeNeighborSet();
		
		// build the Rips complex
		// In our case the rips complex is a javaplex object, and javaplex uses streams
		this.buildRipsComplex();
		System.out.println("Built Rips complex with "+ripsComplexStream.getSize()+" faces");
		
		
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
		
	    NonLifting nl = new NonLifting();
	    
	    //Vector<ExplicitSimplexStream> cycleStreams = g.computeHomology();
	    nl.computeHomology();
	    
	    nl.buildCoverageRipsComplex(nl.V);
	    //Vector<ExplicitSimplexStream> coverageRipsCycleStreams = nl.computePersistentHomology();
	    nl.computePersistentHomology();
	    
	    // view the complexes
	    //for (ExplicitSimplexStream str : cycleStreams) {
	    //	g.drawComplex(str);
	    //}
	    nl.drawComplex();
	}
	
	
	
}
