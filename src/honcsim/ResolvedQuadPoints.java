/*
 * ResolvedQuadPoints.java
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
 * An example resolving all the (interior) vertices and no higher dim simplices
 * 
 * Could be pretty easily modified to resolve the boundary simplices too.
 * 
 */

import java.util.*;

import m4rjni.Mzd;
import edu.stanford.math.plex4.streams.impl.ExplicitSimplexStream;

public class ResolvedQuadPoints extends CoverageExperiment {
    
    static int vsDimension = 18;
    
    int gridHeight = 0;
    int gridWidth = 0;
    double gridSpacing = 0.0;
    
    // extra data structure for keeping track of DPoints
    DPoint [][] pGrid = null;
    
    /**
     * Constructor
     * 
     * @param width
     * @param height
     * @param spacing
     */
    public ResolvedQuadPoints(int width, int height, double spacing) {
	super(vsDimension, 1.0);
	
	this.gridWidth = width;
	this.gridHeight = height;
	this.gridSpacing = spacing;

	// allocate the grid of DPoints
	points = new Vector<DPoint>(width*height);
	pGrid = new DPoint[width][];
	for (int i=0; i<width; i++) {
	    pGrid[i] = new DPoint[height];
	}
	
	for (int j=0; j<height; j++) {
	    for (int i=0; i<width; i++) {
		pGrid[i][j] = new DPoint(j*spacing, i*spacing, vsDimension);
		points.add(pGrid[i][j]);
	    }
	}
	
	// put the appropriate inventories
	for (int i=0; i<width; i++) {
	    for (int j=0; j<height; j++) {
		pGrid[i][j].addInventoryVector(this.basis[3*(i%3) + (j%3)]);
		System.out.println("Putting vector "+(3*(i%3) + (j%3))+" at point "+i+","+j);
	    }
	}
	
	for (int j=0; j<height; j++) {
	    for (int i=0; i<width; i++) {
		pGrid[i][j].addInventoryVector(this.basis[9 + 3*(j%3) + (i%3)]);
		System.out.println("Putting vector "+(9 + 3*(j%3) + (i%3))+" at point "+i+","+j);
	    }
	}
	
	// compute the neighbor set
	this.computeNeighborSet();
    }
    
    
    /* **************************************
     *
     * Main Routine
     *
     * **************************************
     */
    public static void main(String[] args) {

	ResolvedQuadPoints rqp = new ResolvedQuadPoints(10, 10, 0.6);

	
	rqp.buildRipsComplex();
	rqp.buildCoverageComplex(rqp.V);
	rqp.computeHomology();
	rqp.drawComplex();
	
	Mzd randomVector = new Mzd(1,ResolvedQuadPoints.vsDimension);
	randomVector.randomize();
	System.out.println("Random vector: ");
	randomVector.print();
	rqp.buildCoverageComplex(randomVector);
	rqp.computeHomology();
	rqp.drawComplex();
	
	rqp.buildCoverageComplex(rqp.basis[0]);
	rqp.computeHomology();
	rqp.drawComplex();
	
    }


}
