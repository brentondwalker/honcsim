/*
 * DSimplex.java
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
 * Class representing a simplex of DPoint.
 * Should provide some utilities for keeping track of the common neighbors and vector spaces spanned.
 * The plan right now is to create these just for Rips-maximal 0 or 1-simplices, or 2-simplices
 * 
 * The point of keeping exclusive neighbors is to support vector drainage:
 * For example suppose we have a point X with a coface sigma.
 * Let nb(sigma) be the neighbors common to all points of sigma.
 *   (this is essentially the union of all maximal simplices containing sigma)
 * Let Ui = span(nb(sigma)\X)
 * Then essentially Ui is the vector space that covers sigma without any help from X
 * We'll use Ui to decide how to "drain" vectors from the inventory of X
 * 
 * TODO: DSimplex and DPoint do not coordinate well.  For example adding vertices to a DSimplex
 *       does not check if the DPoint are neighbors of each other.
 */

import java.util.*;

import m4rjni.Mzd;


public class DSimplex {
	
    // the dimension of the vector space
    // this may be set to 0 if no vertices are added to the simplex
    int vsDimension = 0;
    
	// the vertices of the simplex
	//public Vector<DPoint> vertices = null;
	public HashSet<DPoint> vertices;
	
	// the neighbors common to all vertices, including the vertices themselves
	public HashSet<DPoint> neighbors = null;
	
	// the neighbors common to all vertices, excluding the vertices themselves
	public HashSet<DPoint> exclusiveNeighbors = null;
	
	// m4ri objects.
	//
	// reducedBasis is a square reduced matrix representing the vector
	// space spanned by the vectors in M
	Mzd reducedBasis = null;
	int rank = 0;

	// reducedBasis is a square reduced matrix representing the vector
	// space spanned by the vectors in M
	Mzd exclusiveReducedBasis = null;
	int exclusiveRank = 0;

	
	/**
	 * Blank constructor
	 */
	public DSimplex() {
		vertices = new HashSet<DPoint>();
		neighbors = new HashSet<DPoint>();
		exclusiveNeighbors = new HashSet<DPoint>();
	}


	/**
	 * constructor with a List of DPoint provided
	 * 
	 * TODO: ensure that all the vertices added to a simplex consider themselves
	 *       neighbors of each other
	 * 
	 * TODO: fill in the ripsCofaces of any vertex added to a simplex
	 * 
	 * @param verts
	 */
	public DSimplex(List<DPoint> verts) {
		if (verts==null) {
			System.out.println("ERROR: DSimplex constructor called with a null array.");
			System.exit(0);
		}
		
		// record the dimension of the vector space, and make sure it's the
		// same for all vertices
		vsDimension = verts.get(0).vsDimension;
		for (DPoint v : verts) {
		    if (v.vsDimension != vsDimension) {
		        throw(new IllegalArgumentException("dimensions of vector spaces at the vertices are not the same"));
		    }
		}
		
		vertices = new HashSet<DPoint>(verts);
		
		if (verts.size()==0) { return; }
		
		// compute the common neighbors
		for (DPoint p : vertices) {
			if (exclusiveNeighbors == null) {
				exclusiveNeighbors = new HashSet<DPoint>(p.nbrs);
			} else {
				exclusiveNeighbors.retainAll(p.nbrs);
				if (exclusiveNeighbors.size()==0) { break; }
			}
		}
		neighbors = new HashSet<DPoint>(exclusiveNeighbors);
		for (DPoint p : vertices) {
			neighbors.add(p);
		}
		
		computeReducedBasis();
		computeExclusiveReducedBasis();
		
		//System.out.println(this);
	}
	

	/**
	 * we want a DSimplex to be uniquely identified by its vertices
	 */
	public int hashCode() {
		return this.vertices.hashCode();
	}
	
	
	/**
	 * Compute a reduced basis for the space covering this DSimplex, excluding the
	 * contribution from the set of vertices passed in.
	 * 
	 * If the set is null or empty it does no work and just returns a copy of
	 * the current reduced basis.
	 * 
	 * It's ok if the set contains vertices not actually in this simplex.  This
	 * method will just ignore them.
	 * 
	 * @param excludedPoints
	 * @return
	 */
	public Mzd exclusiveBasis(HashSet<DPoint> excludedPoints) {
		if (excludedPoints==null || excludedPoints.size()==0) {
			System.out.println("WARNING: exclusiveBasis() received null or zero.  Returning reducedBasis.");
			return new Mzd(this.reducedBasis);
		}
		Mzd EB = new Mzd(vsDimension,vsDimension);
		int ebRank = 0;
		for (DPoint p : neighbors) {
			if (! excludedPoints.contains(p)) {
				for (Mzd v : p.M) {
					// only try adding if we aren't already at full rank
					if (ebRank < DrainageExperiment.vsDimension) {
						Mzd.copyRow(EB, ebRank, v, 0);
						ebRank = EB.echelonize(false);
					}
				}
			}
		}
		return EB;
	}
	
	
	/**
	 * Same as method above, but allow the caller to pass just a single DPoint.
	 * 
	 * @param p
	 * @return
	 */
	public Mzd exclusiveBasis(DPoint p) {
		HashSet<DPoint> hsp = new HashSet<DPoint>();
		hsp.add(p);
		return exclusiveBasis(hsp);
	}	
	
	
	/**
	 * (re)build the reduced basis
	 * this could be optimized a little by adding the vectors more than one at a time,
	 * but we also want to keep this matrix square.... ah well
	 */
	private void computeReducedBasis() {
        if (reducedBasis != null) {
            reducedBasis.destroy();
        }
	    if (vertices.size() == 0) {
	        vsDimension = 0;
	        rank = 0;
	        reducedBasis = null;
	        return;
	    }
		reducedBasis = new Mzd(vsDimension,vsDimension);
		for (DPoint p : vertices) {
			for (Mzd v : p.M) {
				// only try adding if we aren't already at full rank
				if (rank < vsDimension) {
					Mzd.copyRow(reducedBasis, rank, v, 0);
					rank = reducedBasis.echelonize(false);
				}
			}
		}
	}


	/**
	 * (re)build the exclusive reduced basis
	 */
	private void computeExclusiveReducedBasis() {
        if (exclusiveReducedBasis != null) {
            exclusiveReducedBasis.destroy();
        }
	    if (vertices.size() == 0) {
            vsDimension = 0;
	        rank = 0;
	        reducedBasis = null;
	        return;
	    }
		exclusiveReducedBasis = new Mzd(vsDimension,vsDimension);
		for (DPoint p : exclusiveNeighbors) {
			for (Mzd v : p.M) {
				if (exclusiveRank < vsDimension) {
					Mzd.copyRow(exclusiveReducedBasis, exclusiveRank, v, 0);
					exclusiveRank = exclusiveReducedBasis.echelonize(false);
				}
			}
		}
	}
	

	/**
	 * Add a list of vertices to the DSimplex
	 * Updates the common neighbors HashSets
	 * 
	 * @param verts
	 */
	public void addVertices(List<DPoint> verts) {
		if (verts==null) { return; }
		if (vertices.size() == 0) {
		    vsDimension = verts.get(0).vsDimension;
		}
		for (DPoint p : verts) {
		    if (p.vsDimension != vsDimension) {
                throw(new IllegalArgumentException("dimensions of vector spaces at the vertices are not the same"));
            }
		}
		
		// in case we are starting from an empty simplex
		if (vertices.size() == 0) {
		    exclusiveNeighbors.addAll(verts.get(0).nbrs);
		    neighbors.addAll(exclusiveNeighbors);
		}
		
		// compute the intersection of all the neighbor sets
		// effectively finds the maximal coface of the simplex
		for (DPoint p : verts) {
			vertices.add(p);
			exclusiveNeighbors.retainAll(p.nbrs);
			neighbors.retainAll(p.nbrs);
		}
		neighbors.addAll(verts);
		
		// recompute the reduced bases
		computeReducedBasis();
		computeExclusiveReducedBasis();
	}

	
	/**
	 * Add a vertex to the DSimplex
	 * Updates the common neighbors HashSets
	 * 
	 * @param p
	 */
	public void addVertex(DPoint p) {
		if (p==null) { return; }
		ArrayList<DPoint> ll = new ArrayList<DPoint>();
		ll.add(p);
		addVertices(ll);
	}

	
	/**
	 * Remove some vertices from the DSimplex.
	 * Requires re-computing the common neighbors from scratch.
	 * When we remove vertices from a simplex, unless it removes all
	 * vertices, it actually increases the neighbor set.
	 * 
	 * @param verts
	 */
	public void delVertices(List<DPoint> verts) {
		vertices.removeAll(verts);
		
		if (vertices.size() == 0) {
		    if (reducedBasis != null) { reducedBasis.destroy(); }
		    if (exclusiveReducedBasis != null) { exclusiveReducedBasis.destroy(); }
            reducedBasis = null;
		    exclusiveReducedBasis = null;
		    vsDimension = 0;
		    neighbors.clear();
		    exclusiveNeighbors.clear();
		    return;
		}

        // re-compute the common neighbors, inclusive and exclusive
		exclusiveNeighbors.clear();
		neighbors.clear();
		exclusiveNeighbors.addAll(vertices.iterator().next().nbrs);
		for (DPoint p : vertices) {
		    exclusiveNeighbors.retainAll(p.nbrs);
		    if (exclusiveNeighbors.size()==0) { break; }
		}
		neighbors.addAll(exclusiveNeighbors);
		for (DPoint p : vertices) {
			neighbors.add(p);
		}
		
		// recompute the reduced bases
		computeReducedBasis();
		computeExclusiveReducedBasis();
	}
	
	
	/**
	 * Remove a vertex from the DSimplex.
	 * Requires re-computing the common neighbors from scratch.
	 * 
	 * @param p
	 */
	public void delVertex(DPoint p) {
		if (p==null) { return; }
		ArrayList<DPoint> ll = new ArrayList<DPoint>();
		ll.add(p);
		delVertices(ll);
	}
	
	
	/**
	 * This method contains Mzd objects, so we need a way to to clean them up.
	 */
	public void destroy() {
	    if (reducedBasis != null) {
	        reducedBasis.destroy();
	        reducedBasis = null;
	    }
	    if (exclusiveReducedBasis != null) {
	        exclusiveReducedBasis.destroy();
	        exclusiveReducedBasis = null;
	    }
	    neighbors.clear();
	    exclusiveNeighbors.clear();
	    vsDimension = 0;
	    rank = 0;
	}

	
	/**
	 * toString method.
	 */
	public String toString() {
		String s = "DSimplex: [";
		for (DPoint p : this.vertices) {
			s += " "+p+" ";
		}
		s += "]";
		return s;
	}
	
}
