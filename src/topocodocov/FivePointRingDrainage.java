/*
 * FivePointRingDrainage.java
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
 * Simple experiment with the basis drainage code
 */

import java.util.*;

public class FivePointRingDrainage extends DrainageExperiment {
	
	public FivePointRingDrainage(int dim) {
		super(dim, 1.0);
		
		points = new Vector<DPoint>(5);
		double RADIUS = 0.8;
		for (int i=0; i<5; i++) {
			points.add(new DPoint(RADIUS*Math.cos(i*2.0*Math.PI/5.0), RADIUS*Math.sin(i*2.0*Math.PI/5.0), FivePointRingDrainage.vsDimension));
		}
		
		this.computeNeighborSet();
		
		// give everything a full basis to start
		/*
		for (DPoint p : points) {
			p.addInventoryVectors(this.basis);
		}
		*/
		
		// give everything a random basis (probably) to start
		for (DPoint p : points) {
			for (int i=0; i<p.vsDimension/2; i++) {
				p.addRandomInventoryVector();
			}
		}
		
		// give everything less than a full basis
		/*
		for (int i=0; i<points.size(); i++) {
			points.get(i).addInventoryVector(this.basis[i%2]);
		}
		points.get(0).addInventoryVector(this.basis[1]);
		*/
		
		// give everything a random set of vectors, dpeending on the size of its neighborhood
		/*
		for (int i=0; i<points.size(); i++) {
			for (int j=0; j<(vsDimension/points.get(i).nbrs.size()); j++) {
				points.get(i).addRandomInventoryVector();
			}
		}
		*/
		
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
		
		FivePointRingDrainage fprd = new FivePointRingDrainage(dim);
		
		fprd.buildRipsComplex();
		fprd.buildCoverageComplex(fprd.V);
		
		/*
		for (DPoint p : fprd.points) {
			fprd.drainVertex(p,0);
		}
		
		System.out.println("\n\n\n  2nd pass\n\n\n");
		
		for (DPoint p : fprd.points) {
			fprd.drainVertex(p,0);
		}
		*/
		
		int numDrained = 1;
		int passnum = 0;
		Vector<Integer> inventorySizes = new Vector<Integer>();
		while (numDrained > 0) {
			System.out.println("\n*************** Pass # "+passnum+" **************");
			inventorySizes.add(fprd.totalInventorySize());
			numDrained = 0;
			for (DPoint p : fprd.points) {
				numDrained += fprd.drainVertex(p,Math.max(1, fprd.vsDimension/10), true);
			}
			System.out.println("\n**** pass #"+passnum+": drained "+numDrained+" vectors total\n\n");
			passnum++;
		}
		System.out.println("Total inventory change:");
		for (Integer x : inventorySizes) {
			System.out.println(x);
		}
	}
	
}
