/*
 * FencedCoverageSquare.java
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
 * Kind of like CoverageGrid, but this time have the points inside the box
 * be randomly distributed.
 */


import java.util.Vector;

public class FencedCoverageSquare extends CoverageExperiment {

	/*
	 * class data
	 */
	int gridHeight = 0;
	int gridWidth = 0;
	double gridSpacing = 0.0;
	int numPoints = 0;
	int inventorySize = 0;

	// add a parameter to keep the 'uniformly random' points from being too close together
	double minSpacing = this.radius/3.0;
	
	
	/**
	 * Constructor
	 * 
	 */
	public FencedCoverageSquare(int height, int width, double spacing, int vsDimension, int npoints, int inventorySize) {
		// call CoverageExperiment constructor
		super(vsDimension, 1.0);

		this.vsDimension = vsDimension;
		this.gridHeight = height;
		this.gridWidth = width;
		this.gridSpacing = spacing;
		this.numPoints = npoints;
		this.inventorySize = inventorySize;
		
		// allocate the grid of DPoints
		points = new Vector<DPoint>(2*gridHeight + 2*gridWidth - 2 + numPoints);
		
		placeFullBasisFenceGrid(gridHeight, gridWidth, gridSpacing);
		
		placeRandomInteriorNodes(null, (gridHeight-1)*gridSpacing, (gridWidth-1)*gridSpacing, numPoints, inventorySize);
		
		// compute the neighbor set
		// in our world the underlying set of points doesn't change,
		// so we compute the neighbor set in the constructor
		this.computeNeighborSet();
		
		
		
	}
	
	
	private void placeRandomInteriorNodes(double[] origin, double height, double width, int n, int numVectors) {
		double r2 = this.minSpacing*this.minSpacing;
		
		// if no origin is supplied assume 0,0
		if (origin==null) {
			origin = new double[]{ 0.0 , 0.0 };
		}
		for (int i=0; i<n; i++) {
			// try not to put the random points too close together
			// give it 10 tries and then just place the point wherever
			// note that it's no longer an iid uniform distribution now
			// also this is REALLY inefficient, but this part not really worth optimizing
			double xx = 0.0;  double yy = 0.0;
			for (int t=0; t<100; t++) {
				xx = origin[0]+Math.random()*width;
				yy = origin[1]+Math.random()*height;
				boolean tooClose = false;
				for (DPoint p : points) {
					if ((xx-p.x)*(xx-p.x)+(yy-p.y)*(yy-p.y) < r2) {
						tooClose = true;
						break;
					}
				}
				if (! tooClose) {
					break;
				}
				if (t==99) {
					System.out.println("used up all the trial placements!!");
				}
			}
			points.add(new DPoint(xx , yy, vsDimension));
			for (int j=0; j<numVectors; j++) {
				points.lastElement().addRandomInventoryVector();
			}
		}
	}
	
	
	private void placeFullBasisFenceGrid(int height, int width, double spacing) {
		// add top and bottom rows
		for (int i=0; i<width; i++) {
			points.add(new DPoint(0, i*spacing, vsDimension));
			points.lastElement().addInventoryVectors(this.basis);
			points.add(new DPoint((height-1)*spacing, i*spacing, vsDimension));
			points.lastElement().addInventoryVectors(this.basis);
		}
		
		// add sides
		for (int i=1; i<(height-1); i++) {
			points.add(new DPoint(i*spacing, 0, vsDimension));
			points.lastElement().addInventoryVectors(this.basis);
			points.add(new DPoint(i*spacing, (width-1)*spacing, vsDimension));
			points.lastElement().addInventoryVectors(this.basis);
		}
	}
	
	
	
	/* **************************************
	 *
	 * Main Routine
	 *
	 * **************************************
	 */
	public static void main(String[] args) {
		if (args.length != 6) {
			System.out.println("usage: FencedCoverageSquare <grid width> <grid height> <grid spacing> <vsDim> <numPoints> <vectors-per-node>\n");
			System.exit(0);
		}
		int gridWidth = Integer.parseInt(args[0]);
		int gridHeight = Integer.parseInt(args[1]);
		double gridSpacing = Double.parseDouble(args[2]);
		int vsDim = Integer.parseInt(args[3]);
		int numPoints = Integer.parseInt(args[4]);
		int inventorySize = Integer.parseInt(args[5]);
		
		FencedCoverageSquare g = new FencedCoverageSquare(gridWidth, gridHeight, gridSpacing, vsDim, numPoints,  inventorySize);

		g.buildRipsComplex();
		g.buildCoverageComplex(g.V);
		g.buildCoverageRipsComplex(g.V);		
	    //g.computeHomology();
		g.computePersistentHomology();
	    g.drawComplex();
	}
	
}
