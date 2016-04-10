/*
 * DComplex.java
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
 * Class to represent a collection of DSimplex.
 * 
 * Don't require hereditaty here.  For our purposes only need to list the "maximal"
 * simplices.
 * This class will support computing U=V(\Sigma), the largest vector space such that
 * \Sigma \subseteq RC(U)
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


public class DComplex {
	
	public int vsDimension = 0;
	public HashSet<DSimplex> simplices = new HashSet<DSimplex>();

	
	/**
	 * Blank constructor
	 * Adds no actual simplices.
	 * 
	 * @param dim
	 */
	public DComplex(int dim) {
		this.vsDimension = dim;
	}
	
	
	/**
	 * Constructor.
	 * Figures out the vsDimension from the simplices passed in.
	 * If you pass in an empty list or empty simplices weird things will happen.
	 * 
	 * @param simpl
	 */
	public DComplex(List<DSimplex> simpl) {
		if (simpl == null) { return; }
		if (simpl.size() == 0) { return; }
		if (simpl.get(0) == null) { return; }
		if (simpl.get(0).vertices.size() == 0) { return; }
		this.vsDimension = DrainageExperiment.vsDimension;

		for (DSimplex ss : simpl) {
			simplices.add(ss);
		}
	}
	
	
	/**
	 * Add a simplex to the complex
	 * 
	 * @param ss
	 */
	public void addSimplex(DSimplex ss) {
		//System.out.println("DComplex.add("+ss+")");
		simplices.add(ss);
	}
	
	

	/**
	 * Compute the largest vector space such that all my DSimplex are in it.
	 * 
	 * An obvious optimization here is that if a simplex has full rank (incl neighbors)
	 * then there is no need to compute the intersection with it; it has no effect.
	 * The only thing you can possibly loos that way is the remnants of the particular basis
	 * the nodes had to start.
	 * 
	 * This is the main reason for creating this class.
	 * @return
	 */
	public Mzd vsCover() {
		if (this.vsDimension == 0) {
			return null;
		}
		
		// initialize the process with the standard basis
		Mzd U = Mzd.identityMatrix(vsDimension);
		
		for (DSimplex ds : simplices) {
			//System.out.println("vsCover(): processing "+ds);
			
			// if this simplex has rank 0 then it's going to kill the whole
			// intersection.  Just return.
			if (ds.rank == 0) {
				//System.out.println("total rank of this simplex is 0 (rank="+ds.rank+")  -  returning");
				U.destroy();
				return new Mzd(0,this.vsDimension);
			}
			
			// if the simplex has full rank then there is no point in computing the intersection.
			// it has no effect.
			if (ds.rank == vsDimension) {
				//continue;
			}
			
			// gather the full inventory covering this simplex
			
			// compute the common neighbors
			// one obvious optimization here is if we kept track of common neighbor
			// sets we had already used, skip processing this loop if:
			// - this common neighbor set was already examined in connection with a different DSimplex
			// - this common neighbor set contains any other common neighbor set we've looked at
			// 
			HashSet<DPoint> neighbors = ds.neighbors;
			
			// collect all the inventories and reduce
			if (neighbors==null || neighbors.size()==0) { return null; }
			int numVectors = 0;
			for (DPoint p : neighbors) {
				numVectors += p.rank;
			}
			Mzd W = new Mzd(numVectors, vsDimension);
			int rowi = 0;
			for (DPoint p : neighbors) {
				Mzd.copyRows(W, rowi, p.reducedBasis, 0, p.rank);
				rowi += p.rank;
			}
			int total_rank = W.echelonize(false);
			
			if (total_rank > 0) {
				Mzd Wr = W.submatrix(null, 0, 0, total_rank, vsDimension);
				W.destroy();
				W = Wr;
			} else {
				//System.out.println("rank of the intersection has been reduced to 0  -  returning");
				//W.print();
				W.destroy();
				U.destroy();
				return new Mzd(0,this.vsDimension);
			}
			
			// compute the running intersection
			Mzd newVsI = Mzd.vsIntersect(U, W);
			U.destroy();
			if (newVsI == null) {
				System.out.println("WARNING: Mzd.vsIntersect() returned null!");
				System.exit(-1);
			}
			if (newVsI.isZero()) {
				return new Mzd(0,this.vsDimension);
			}
			U = newVsI;
		}
		
		return U;
	}
	
	
	
	
	
	
}
