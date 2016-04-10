/*
 * LayeredSimplex.java
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
 * The example with a 2-simplex overlapping a fully-covering 2-chain
 * where only the vertices connect to the cover
 * 
 */

import java.util.*;

import edu.stanford.math.plex4.streams.impl.ExplicitSimplexStream;

public class LayeredSimplex extends CoverageExperiment {

	/*
	 * class data
	 */
	static final int vsDimension = 3;
	
	/**
	 * Constructor
	 * 
	 * @param width
	 * @param height
	 * @param spacing
	 */
	public LayeredSimplex() {
		// call CoverageExperiment constructor
		super(vsDimension, 1.0);
		
		// allocate the grid of DPoints
		points = new Vector<DPoint>(9);
		
		// the elevated triangle
		double r1 = Math.sqrt(3)/3-0.001;
		for (int i=0; i<3; i++) {
			points.add(new DPoint(r1*Math.cos(2*Math.PI*i/3), r1*Math.sin(2*Math.PI*i/3), vsDimension));
			points.lastElement().addInventoryVector(this.basis[i]);
			points.lastElement().elevatePoint(true);
		}

		// the overlapping triangle that is part of the covering 2-chain
		for (int i=0; i<3; i++) {
			points.add(new DPoint(r1*Math.cos(2*Math.PI*(i+0.5)/3), r1*Math.sin(2*Math.PI*(i+0.5)/3), vsDimension));
			points.lastElement().addInventoryVector(this.basis[i]);
		}
		
		// outer triangles
		double r2 = 2*Math.sqrt(3)/3-0.001;
		// the overlapping triangle that is part of the covering 2-chain
		for (int i=0; i<3; i++) {
			points.add(new DPoint(r2*Math.cos(2*Math.PI*i/3), r2*Math.sin(2*Math.PI*i/3), vsDimension));
			points.lastElement().addInventoryVector(this.basis[(i+1)%vsDimension]);
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
		
		LayeredSimplex ls = new LayeredSimplex();
	    
	    //Vector<ExplicitSimplexStream> cycleStreams = g.computeHomology();
	    ls.computeHomology();
	    
	    ls.buildCoverageRipsComplex(ls.V);
	    //Vector<ExplicitSimplexStream> coverageRipsCycleStreams = ls.computePersistentHomology();
	    ls.computePersistentHomology();
	    
	    // view the complexes
	    //for (ExplicitSimplexStream str : cycleStreams) {
	    //	g.drawComplex(str);
	    //}
	    ls.drawComplex();
	}
	
	
}
