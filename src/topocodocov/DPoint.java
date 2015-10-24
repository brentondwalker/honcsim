/*
 * Dpoint.java
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
 * class that will represent a point holding data
 */

import java.util.*;

import m4rjni.Mzd;


public class DPoint {

	/*
	 * static variable to keep track of global point index
	 */
	private static int indexCounter = 0;
	
	/*
	 * class variables
	 */
	
	// global unique index for this DPoint
	int index = indexCounter++;

	// location
	double x, y;
	
	// dimension of the vector space we're dealing with
	// need this to construct the m4ri objects below
	int vsDimension = 0;
	
	// vector to keep track of my Rips neighbors
	HashSet<DPoint> nbrs = null;

	// vector to track my Rips neighbors who have index larger than mine
	// useful for complex construction
	HashSet<DPoint> nbrsUp = null;

	// array of sets containing the 2-maximal (or whatever) cofaces of this point
	HashSet<DSimplex> ripsCofaces = null;
	
	// m4ri objects.
	// M is the inventory of vectors
	Vector<Mzd> M = null;
	// reducedBasis is a square reduced matrix representing the vector
	//space spanned by the vectors in M
	Mzd reducedBasis = null;
	// the current rank of the reduced basis
	int rank = 0;

	// whether or not to draw this point slightly above the rest of the complex
	protected boolean elevated = false;
	
	/**
	 * Constructor
	 * 
	 * @param x
	 * @param y
	 * @param vsDimension
	 */
	public DPoint(double x, double y, int vsDimension) {
		this.x = x;
		this.y = y;
		this.vsDimension = vsDimension;
		nbrs = new HashSet<DPoint>();
		nbrsUp = new HashSet<DPoint>();
		ripsCofaces = new HashSet<DSimplex>();
		M = new Vector<Mzd>();
		reducedBasis = new Mzd(vsDimension,vsDimension);
	}

	
	/**
	 * Add a vector to the inventory.
	 * Input vector must have the right dimensions (1xD)
	 * 
	 * @param v
	 * @return 1 if rank was increased, 0 if rank stayed the same
	 */
	public int addInventoryVector(Mzd v) {
		if (v.getNcols()!=vsDimension || v.getNrows()!=1) {
			System.out.println("ERROR: addInventoryVector() - vectors must have the right dimension");
			return 0;
		}
		//System.out.println("\n--------------\naddInventoryVector():");
		//v.print();
		//System.out.println("reducedBasis:");
		//reducedBasis.print();
		
		M.add(new Mzd(v));
		
		// update the reduced basis matrix
		int oldrank = rank;
		if (rank < vsDimension) {
			Mzd.copyRow(reducedBasis, rank, v, 0);
			rank = reducedBasis.echelonize(false);
		}
		//System.out.println("new rank="+rank);
		//System.out.println("reducedBasis:");
		//reducedBasis.print();

		if (oldrank != rank) {
			return 1;
		}
		return 0;
	}
	

	/**
	 * Add an array of vectors to the inventory.
	 * Separate function for this just for convenience.
	 * 
	 * @param vv
	 * @return
	 */
	public int addInventoryVectors(Mzd[] vv) {
		if (vv.length==0) { return 0; }
		if (vv[0].getNcols()!=vsDimension) {
			System.out.println("ERROR: addInventoryVectors() - vectors must have the right dimension");
			return 0;
		}
		
		int dRank = 0;
		for (int i=0; i<vv.length; i++) {
			dRank += addInventoryVector(new Mzd(vv[i]));
		}
		
		return dRank;
	}
	
	
	/**
	 * Create a random vector in the vector space and try adding it to the inventory.
	 * Doesn't include the zero vector.  If it pulls the zero vector it will keep
	 * trying until it gets one non-zero.
	 * 
	 * @return 1 if rank was increased, 0 if rank stayed the same
	 */
	public int addRandomInventoryVector() {
	    Mzd v = new Mzd(1,vsDimension);
	    v.randomize();
	    while (v.isZero()) {
	    	v.randomize();
	    }
	    return addInventoryVector(v);
	}
	
	
	/**
	 * Clear out the vectors in the inventory
	 */
	public void clearInventory() {
		for (Mzd v : M) {
			v.destroy();
		}
		M.clear();
		reducedBasis.destroy();
		reducedBasis = new Mzd(vsDimension,vsDimension);
		this.rank = 0;
	}

	
	/**
	 * Set the point's inventory from a basis matrix, where the rows contain the vectors
	 * 
	 * @param BB
	 */
	public void setInventory(Mzd BB) {
		//System.out.println("setInventory()");
		clearInventory();
		if (BB.getNcols() != this.vsDimension) {
			System.out.println("ERROR: setInventory() passed in a matrix of vectors of the wrong dimension");
			return;
		}
		for (int i=0; i<BB.getNrows(); i++) {
			Mzd v = new Mzd(1,this.vsDimension);
			Mzd.copyRow(v, 0, BB, i);
			if (v.isZero()) {
				v.destroy();
			} else {
				M.add(v);
			}
		}
		for (int i=0; i<BB.getNrows(); i++) {
			if (this.rank == vsDimension) { break; }
			Mzd.copyRow(reducedBasis, rank, BB, i);
			this.rank = reducedBasis.echelonize(false);
		}
		//System.out.println("\t...done.");
	}
	
	
	/**
	 * If set to true the point will be drawn slightly above the rest of the complex.
	 * It doesn't affect the location of the point for the purposes of neighbor computation.
	 * 
	 * @param x
	 */
	public void elevatePoint(boolean x) {
		this.elevated = x;
	}
	
	
	/**
	 * Check if this point should be drawn slightly above the rest of the complex.
	 * 
	 * @return true if the point is elevated.
	 */
	public boolean isElevated() {
		return this.elevated;
	}
	
	
	/**
	 * toString method.
	 */
	public String toString() {
		String s = "DP."+index;
		return s;
	}
	
	
	
	public String print() {
		String s = "DP."+index+":("+x+","+y+")";
		return s;
	}
	
	
}
