/*
 * SimplexBubbleBox.java
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
 * Apply drainage algs to bubble box
 * 
 * This one is simpler.  No outer ring.
 * How can I be so stupid??
 */

import java.util.*;
import java.lang.Thread;

import m4rjni.Mzd;

public class SimpleBubbleBox extends DrainageExperiment {
	
	Mzd V = null;
	
	/**
	 * Constructor
	 * 
	 * @param width
	 * @param height
	 * @param spacing
	 */
	public SimpleBubbleBox(int dim) {
		// call CoverageExperiment constructor
		super(dim, 1.0);

		// allocate the grid of DPoints
		points = new Vector<DPoint>(8);

		// add the points on the inner square
		double s2 = 1.0/Math.sqrt(2.0)/2.0;
		points.add(new DPoint(s2, s2, DrainageExperiment.vsDimension));
		points.add(new DPoint(s2, -s2, DrainageExperiment.vsDimension));
		points.add(new DPoint(-s2, -s2, DrainageExperiment.vsDimension));
		points.add(new DPoint(-s2, s2, DrainageExperiment.vsDimension));

		// compute the neighbor set
		this.computeNeighborSet();
		
		// construct the standard basis for the full vector space V
		V = Mzd.identityMatrix(vsDimension);

		// give the inner points random bases
		points.get(0).addInventoryVector(this.basis[0]);		//points.get(0).addInventoryVector(this.basis[1]);
		points.get(2).addInventoryVector(this.basis[1]);		//points.get(2).addInventoryVector(this.basis[0]);
		points.get(1).addInventoryVector(this.basis[0]);		//points.get(1).addInventoryVector(this.basis[1]);
		points.get(3).addInventoryVector(this.basis[1]);		//points.get(3).addInventoryVector(this.basis[0]);
		
	}
	
	
	/* **************************************
	 *
	 * Main Routine
	 *
	 * **************************************
	 */
	public static void main(String[] args) {

		if (args.length != 1) {
			System.out.println("usage: FivePointRingDrainage <vsDimension>");
			System.exit(0);
		}
		
		int dim = Integer.parseInt(args[0]);
		
		SimpleBubbleBox bb = new SimpleBubbleBox(dim);
		bb.buildRipsComplex();
		System.out.println("Size of DComplex: "+bb.maximalSimplices.simplices.size());
		
		//Vector<ExplicitSimplexStream> cycleStreams = g.computeHomology();
		//bb.computeHomology();

		//Vector<ExplicitSimplexStream> coverageRipsCycleStreams = g.computePersistentHomology();

		// view the complexes
		//for (ExplicitSimplexStream str : cycleStreams) {
		//	g.drawComplex(str);
		//}
		//bb.drawComplex();
		
		
		if (false) {
			int numDrained = 1;
			int passnum = 0;
			Vector<Integer> inventorySizes = new Vector<Integer>();
			while (numDrained > 0) {
				System.out.println("\n*************** Pass # "+passnum+" **************");
				inventorySizes.add(bb.totalInventorySize());
				numDrained = 0;
				for (DPoint p : bb.points) {
					numDrained += bb.drainVertex(p, Math.max(2,BubbleBoxDrainage.vsDimension/10), true);
				}
				System.out.println("\n**** pass #"+passnum+": drained "+numDrained+" vectors total\n\n");
				passnum++;
			}
			System.out.println("Total inventory change:");
			for (Integer x : inventorySizes) {
				System.out.println(x);
			}
		}		

		bb.buildCoverageComplex(bb.V);
		bb.drawComplex();
		
		
		for (int i=0; i<bb.basis.length; i++) {
			bb.buildCoverageComplex(bb.basis[i]);
			bb.drawComplex();
		}
		
		
		// figure out the largest vector space that does cover the whole complex
		Mzd maxVS = bb.maximalSimplices.vsCover();
		System.out.println("maximum vs covering full complex:");
		maxVS.print();

		// look at the vector space(s) covering the two sheets of the central bubble
		DComplex sheet1 = new DComplex(dim);
		DPoint[] verts012 = { bb.points.get(0) , bb.points.get(1) , bb.points.get(2) };
		DSimplex ds012 = new DSimplex(Arrays.asList(verts012));
		sheet1.addSimplex(ds012);
		DPoint[] verts023 = { bb.points.get(0) , bb.points.get(2) , bb.points.get(3) };
		DSimplex ds023 = new DSimplex(Arrays.asList(verts023));
		sheet1.addSimplex(ds023);
		Mzd sheet1vs = sheet1.vsCover();
		System.out.println("maximum vs covering for sheet1:");
		sheet1vs.print();

		DComplex sheet2 = new DComplex(dim);
		DPoint[] verts013 = { bb.points.get(0) , bb.points.get(1) , bb.points.get(3) };
		DSimplex ds013 = new DSimplex(Arrays.asList(verts013));
		sheet2.addSimplex(ds013);
		DPoint[] verts123 = { bb.points.get(1) , bb.points.get(2) , bb.points.get(3) };
		DSimplex ds123 = new DSimplex(Arrays.asList(verts123));
		sheet2.addSimplex(ds123);
		Mzd sheet2vs = sheet2.vsCover();
		System.out.println("maximum vs covering for sheet2:");
		sheet2vs.print();
		
	}

}
