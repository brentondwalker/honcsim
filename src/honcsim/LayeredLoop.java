/*
 * LayeredLoop.java
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

package honcsim;
/*
 * example of a layered loop. A case where we have a contractible boundary
 * But H1(RC(V)) != 0
 */

import java.util.*;

import edu.stanford.math.plex4.streams.impl.ExplicitSimplexStream;

public class LayeredLoop extends CoverageExperiment {
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
	public LayeredLoop() {
		// call CoverageExperiment constructor
		super(vsDimension, 1.0);
		
		// allocate the grid of DPoints
		points = new Vector<DPoint>(15);
		DPoint p = null;
		
		// the loop edge [X Y]
		points.add(new DPoint(0.5, 0.0, vsDimension));
		points.lastElement().elevatePoint(true);
		points.add(new DPoint(-0.5, 0.0, vsDimension));
		points.lastElement().elevatePoint(true);
		
		// the cross edge that /is/ part of the cover [A B]
		points.add(new DPoint(0.0, 0.5, vsDimension));
		points.add(new DPoint(0.0, -0.5, vsDimension));
		
		// XY helper vertices
		points.add(new DPoint(0.0, Math.sqrt(3)/2, vsDimension));
		points.lastElement().addInventoryVector(this.basis[0]);
		points.add(new DPoint(0.0, -Math.sqrt(3)/2, vsDimension));
		points.lastElement().addInventoryVector(this.basis[1]);
		
		// AB helper vertices
		points.add(new DPoint(Math.sqrt(3)/2, 0.0, vsDimension));
		points.lastElement().addInventoryVectors(this.basis);
		points.add(new DPoint(-Math.sqrt(3)/2, 0.0, vsDimension));
		points.lastElement().addInventoryVectors(this.basis);
		
		// outer rectangle
		points.add(new DPoint(Math.sqrt(3)/2, Math.sqrt(3)/2, vsDimension));
		points.lastElement().addInventoryVectors(this.basis);
		points.add(new DPoint(-Math.sqrt(3)/2, Math.sqrt(3)/2, vsDimension));
		points.lastElement().addInventoryVectors(this.basis);
		points.add(new DPoint(-Math.sqrt(3)/2, -Math.sqrt(3)/2, vsDimension));
		points.lastElement().addInventoryVectors(this.basis);
		points.add(new DPoint(Math.sqrt(3)/2, -Math.sqrt(3)/2, vsDimension));
		points.lastElement().addInventoryVectors(this.basis);
		
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
		
		LayeredLoop ll = new LayeredLoop();
	    
	    //Vector<ExplicitSimplexStream> cycleStreams = g.computeHomology();
	    ll.computeHomology();
	    
	    ll.buildCoverageRipsComplex(ll.V);
	    //Vector<ExplicitSimplexStream> coverageRipsCycleStreams = ll.computePersistentHomology();
	    ll.computePersistentHomology();
	    
	    // view the complexes
	    //for (ExplicitSimplexStream str : cycleStreams) {
	    //	g.drawComplex(str);
	    //}
	    ll.drawComplex();
	}
	
	
}
