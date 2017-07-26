/*
 * CoverageGrid.java
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
 * Entry point for a class modeling a grid layout of DPoints
 * Parameters TBD
 * 
 * 
 * Note: the coverage complex /might/ be natural as a ConditionalFlagComplexStream
 * 
 * ./simrun CoverageGrid 10 10 0.6 0.1 10 3
 * ./simrun CoverageGrid 10 10 0.6 0.1 10 4
 * 
 * 
 */

import java.util.*;
import m4rjni.Mzd;
import edu.stanford.math.plex4.api.Plex4;
import edu.stanford.math.plex4.homology.barcodes.*;
import edu.stanford.math.plex4.homology.chain_basis.*;
import edu.stanford.math.plex4.homology.interfaces.*;
import edu.stanford.math.plex4.autogen.homology.*;
import edu.stanford.math.plex4.streams.impl.*;
import edu.stanford.math.primitivelib.autogen.formal_sum.*;
import edu.stanford.math.primitivelib.algebraic.impl.*;
import edu.stanford.math.plex_viewer.*;
import edu.stanford.math.plex4.homology.chain_basis.Simplex;
import edu.stanford.math.plex4.homology.chain_basis.SimplexComparator;


public class CoverageGrid extends CoverageExperiment {
	
 
	/*
	 * class data
	 */
	int gridHeight = 0;
	int gridWidth = 0;
	double gridSpacing = 0.0;

	// extra data structure for keeping track of DPoints
	DPoint [][] pGrid = null;
	
	// randomly leave some points out of the Rips complex
	double missingPointProbability = 0.05;
	
	// number of vectors to put at each interior node
	int inventorySize = 1;
	
	/**
	 * Constructor
	 * 
	 * @param width
	 * @param height
	 * @param spacing
	 */
	public CoverageGrid(int width, int height, double spacing, double prob, int D, int inventorySize) {
		// call CoverageExperiment constructor
		super(D, 1.0);

		this.gridWidth = width;
		this.gridHeight = height;
		this.gridSpacing = spacing;
		this.missingPointProbability = prob;
		this.inventorySize = inventorySize;
		
		// allocate the grid of DPoints
		points = new Vector<DPoint>(width*height);
		pGrid = new DPoint[width][];
		for (int i=0; i<width; i++) {
			pGrid[i] = new DPoint[height];
		}

		for (int j=0; j<height; j++) {
		    for (int i=0; i<width; i++) {
				// exclude interior points with some probability
				// just so we can get a more interesting Rips complex.
				if (i==0 || j==0 || i==(width-1) || j==(height-1) || (Math.random() >= missingPointProbability)) {
				    pGrid[i][j] = new DPoint(j*spacing, i*spacing, vsDimension);
				    points.add(pGrid[i][j]);
				} else {
				    System.out.println("Excluding point "+i+","+j);
				    pGrid[i][j] = null;
				}
		    }
		}
		
		// fill in a full basis all around the perimeter
		this.perimeterFullBasis();
		
		// compute the neighbor set
		this.computeNeighborSet();
		
		// build the Rips complex
		// In our case the rips complex is a javaplex object, and javaplex uses streams
		this.buildRipsComplex();
		System.out.println("Built Rips complex with "+ripsComplexStream.getSize()+" faces");
		
		// add a bunch of random vectors to the points' inventories
		for (int i=0; i<points.size(); i++) {
			for (int j=0; j<inventorySize; j++) {
				points.get(i).addRandomInventoryVector();
			
			}
		}

		// build the filtered RC-->R complex
		this.buildCoverageRipsComplex(V);
		System.out.println("Built filtered coverage-->Rips complex with "+coverageRipsComplexStream.getSize()+" faces");

		// build the coverage complex
		this.buildCoverageComplex(V);
		System.out.println("Built coverage complex with "+coverageComplexStream.getSize()+" faces");
		
	}

	
	/**
	 * put a full basis into the inventory of each point around the boundary
	 */
	public void perimeterFullBasis() {
		for (int i=0; i<gridWidth; i++) {
			pGrid[i][0].addInventoryVectors(this.basis);
			//System.out.println("full basis on point "+pGrid[i][0].index);
			pGrid[i][gridHeight-1].addInventoryVectors(this.basis);
			//System.out.println("full basis on point "+pGrid[i][gridHeight-1].index);
		}
		for (int i=1; i<(gridHeight-1); i++) {
			pGrid[0][i].addInventoryVectors(this.basis);
			System.out.println("full basis on point "+pGrid[0][i].index);
			pGrid[gridWidth-1][i].addInventoryVectors(this.basis);
			System.out.println("full basis on point "+pGrid[gridWidth-1][i].index);
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
			System.out.println("usage: CoverageGrid <grid width> <grid height> <grid spacing> <exclusion prob> <vsDim> <vectors-per-node>\n");
			System.exit(0);
		}
		int gridWidth = Integer.parseInt(args[0]);
		int gridHeight = Integer.parseInt(args[1]);
		double gridSpacing = Double.parseDouble(args[2]);
		double exclusionProb = Double.parseDouble(args[3]);		
		int vsDim = Integer.parseInt(args[4]);
		int inventorySize = Integer.parseInt(args[5]);
		
		myCoverageGrid g = new myCoverageGrid(gridWidth, gridHeight, gridSpacing, exclusionProb, vsDim, inventorySize);
	    g.computePersistentHomology();
	    g.reComputeCoverageComplex();
	}
}
