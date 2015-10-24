/*
 * Wormhole.java
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
 * build and study the wormhole complex
 */

import java.util.*;

import edu.stanford.math.plex4.streams.impl.ExplicitSimplexStream;


public class Wormhole extends CoverageExperiment {

	/*
	 * class data
	 */
	static final int vsDimension = 1;
	
	/**
	 * Constructor
	 * 
	 * @param width
	 * @param height
	 * @param spacing
	 */
	public Wormhole() {
		// call CoverageExperiment constructor
		super(vsDimension, 1.0);

		// allocate the grid of DPoints
		points = new Vector<DPoint>(15);
		
		double R=1.0;
		
		// empty "red" triangle
		//double t = 0.5;
		double t = R/2.0;
		/*
		for (int i=0; i<3; i++) {
			DPoint p = new DPoint(t*Math.cos(2.0*i*Math.PI/3), t*Math.sin(2.0*i*Math.PI/3), vsDimension);
			points.add(p);
		}
		*/
		
		
		// the inner cover triangles
		//double rt = 1.4*t;
		double rt = Math.sqrt(R*R-3.0*t*t/4.0) - t/2.0 + 0.0001;
		for (int i=0; i<3; i++) {
			DPoint p = new DPoint(rt*Math.cos(2.0*i*Math.PI/3), rt*Math.sin(2.0*i*Math.PI/3), vsDimension);
			p.addInventoryVector(this.basis[0]);
			points.add(p);
		}
		
		// distance of cover triangle to helper vertices
		double re = t; // *0.9;
		// angle to perturb the points by, relative to the central angle of the "red" points
		double th = 0.15;
		for (int i=0; i<3; i++) {
			DPoint p = new DPoint(re*Math.cos(2.0*i*Math.PI/3 + th), re*Math.sin(2.0*i*Math.PI/3 + th), vsDimension);
			points.add(p);
			p = new DPoint(re*Math.cos(2.0*i*Math.PI/3 - th), re*Math.sin(2.0*i*Math.PI/3 - th), vsDimension);
			points.add(p);
		}

		// red support points
		//double rr = 1.1;
		//double rr = 5.0*R/4.0;
		//double rr = Math.sqrt(R*R-3.0*t*t/4.0) + t/2.0;
		double rr = Math.sqrt(R*R-3.0*rt*rt/4.0) + rt/2.0 - 0.0001;
		for (int i=0; i<3; i++) {
			DPoint p = new DPoint(rr*Math.cos(Math.PI/3.0 + 2.0*i*Math.PI/3), rr*Math.sin(Math.PI/3.0 + 2.0*i*Math.PI/3), vsDimension);
			p.addInventoryVector(this.basis[0]);
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
		
		Wormhole wrm = new Wormhole();
	    
	    //Vector<ExplicitSimplexStream> cycleStreams = g.computeHomology();
	    wrm.computeHomology();
	    
	    wrm.buildCoverageRipsComplex(wrm.V);
	    //Vector<ExplicitSimplexStream> coverageRipsCycleStreams = wrm.computePersistentHomology();
	    wrm.computePersistentHomology();

	    // view the complexes
	    //for (ExplicitSimplexStream str : cycleStreams) {
	    //	g.drawComplex(str);
	    //}
	    wrm.drawComplex();
	}
	
	
	
}
