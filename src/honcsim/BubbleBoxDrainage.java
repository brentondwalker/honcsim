/*
 * BubbleBoxDrainage.java
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
 * Apply drainage algs to bubble box
 */

import java.util.*;
import java.lang.Thread;

import m4rjni.Mzd;

public class BubbleBoxDrainage extends DrainageExperiment {
	
	Mzd V = null;
	
	/**
	 * Constructor
	 * 
	 * @param width
	 * @param height
	 * @param spacing
	 */
	public BubbleBoxDrainage(int dim) {
		// call CoverageExperiment constructor
		super(dim, 1.0);
		
		// allocate the grid of DPoints
		points = new Vector<DPoint>(8);
		
		// add the points on the outer square
		double s2 = 1.0/Math.sqrt(2.0)/2.0;
		double s2p = s2 + 0.1;
		points.add(new DPoint(s2p, s2p, vsDimension));
		points.add(new DPoint(-s2p, -s2p, vsDimension));
		points.add(new DPoint(-s2p, s2p, vsDimension));
		points.add(new DPoint(s2p, -s2p, vsDimension));
		
		// add the points on the inner square
		points.add(new DPoint(s2, s2, vsDimension));
		points.add(new DPoint(s2, -s2, vsDimension));
		points.add(new DPoint(-s2, -s2, vsDimension));
		points.add(new DPoint(-s2, s2, vsDimension));
		
		// compute the neighbor set
		this.computeNeighborSet();
		
		// construct the standard basis for the full vector space V
		V = Mzd.identityMatrix(vsDimension);
		
		// give everything a full basis to start
		/*
		for (DPoint p : points) {
			p.addInventoryVectors(this.basis);
		}
		/* */

		// give the outer points structured bases
		/*
		points.get(0).addInventoryVector(this.basis[1]);
		points.get(2).addInventoryVector(this.basis[0]);
		points.get(1).addInventoryVector(this.basis[1]);
		points.get(3).addInventoryVector(this.basis[0]);
		 /* */
		
		// give the outer points random half-bases
		/*
		for (int i=0; i<4; i++) {
			for (int j=0; j<(1+BubbleBoxDrainage.vsDimension/2); j++) {
				points.get(i).addRandomInventoryVector();
			}
		}
		/* */

		
		// give the inner points structured bases
		//*
		points.get(4).addInventoryVector(this.basis[0]);		//points.get(0).addInventoryVector(this.basis[1]);
		points.get(5).addInventoryVector(this.basis[1]);		//points.get(2).addInventoryVector(this.basis[0]);
		points.get(6).addInventoryVector(this.basis[0]);		//points.get(1).addInventoryVector(this.basis[1]);
		points.get(7).addInventoryVector(this.basis[1]);		//points.get(3).addInventoryVector(this.basis[0]);
		 /* */

		// give the inner points random half-bases
		/*
		for (int i=4; i<8; i++) {
			for (int j=0; j<(BubbleBoxDrainage.vsDimension/2); j++) {
				points.get(i).addRandomInventoryVector();
			}
		}
		/* */
		
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
		
		BubbleBoxDrainage bb = new BubbleBoxDrainage(dim);
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

		/*
		bb.drawComplex();
		for (int i=0; i<bb.basis.length; i++) {
			bb.buildCoverageComplex(bb.basis[i]);
			bb.drawComplex();
		}
		*/
		
		System.out.println("\n\nAll points\n----------------------------------------");
		for (DPoint p : bb.points) {
			System.out.println(p+":");
			p.reducedBasis.print();
		}
		
		System.out.println("\n\nAll 2-maximal simplices\n----------------------------------------");
		// figure out which simplex isn't covered here...
		for (DSimplex ds : bb.maximalSimplices.simplices) {
			System.out.println(ds+"\t rank="+ds.rank);
		}
		
		// figure out the largest vector space that does cover the whole complex
		Mzd maxVS = bb.maximalSimplices.vsCover();
		System.out.println("\nmaximum vs covering full complex:");
		maxVS.print();

		// look at the vector space(s) covering the two sheets of the central bubble
		DComplex sheet1 = new DComplex(dim);
		DPoint[] verts456 = { bb.points.get(4) , bb.points.get(5) , bb.points.get(6) };
		DSimplex ds456 = new DSimplex(Arrays.asList(verts456));
		sheet1.addSimplex(ds456);
		DPoint[] verts467 = { bb.points.get(4) , bb.points.get(6) , bb.points.get(7) };
		DSimplex ds467 = new DSimplex(Arrays.asList(verts467));
		sheet1.addSimplex(ds467);
		Mzd sheet1vs = sheet1.vsCover();
		System.out.println("\nmaximum vs covering for sheet1:");
		sheet1vs.print();

		DComplex sheet2 = new DComplex(dim);
		DPoint[] verts457 = { bb.points.get(4) , bb.points.get(5) , bb.points.get(7) };
		DSimplex ds457 = new DSimplex(Arrays.asList(verts457));
		sheet2.addSimplex(ds457);
		DPoint[] verts567 = { bb.points.get(5) , bb.points.get(6) , bb.points.get(7) };
		DSimplex ds567 = new DSimplex(Arrays.asList(verts567));
		sheet2.addSimplex(ds567);
		Mzd sheet2vs = sheet2.vsCover();
		System.out.println("\nmaximum vs covering for sheet2:");
		sheet2vs.print();
		
	}

}
