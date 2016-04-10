/*
 * DrainageExperiment.java
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
 * 
 * A base class to handle some common functionality of coverage complex experiments
 * where we try to drain redundancy out of the complex
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


public class DrainageExperiment {

	/*
	 * class data
	 */
	static int vsDimension = 3;  // this always get overwritten
	double radius = 1.0;

	// identity matrix of vsDimension
	// and standard basis
	Mzd V = null;
	Mzd[] basis = null;

	// data structure containing DPoints
	Vector<DPoint> points = null;

	// maintain our own list of simplices that are either 2-simplices or maximal 0 or 1-simplices
	// this will get filled in whcn we compute the coverage complex...?
	DComplex maximalSimplices = null;
	
	// javaplex objects
	ExplicitSimplexStream ripsComplexStream = null;
	ExplicitSimplexStream coverageComplexStream = null;
	ExplicitSimplexStream cycleStream = null;
	ExplicitSimplexStream ripsCycleStream = null;
	ExplicitSimplexStream coverageRipsComplexStream = null;

	// plex-viewer information
	double[][] domainPoints = null;


	/**
	 * Constructor.
	 * Just the basics
	 * 
	 * @param vsDim
	 * @param r
	 */
	public DrainageExperiment(int vsDim, double r) {
		DrainageExperiment.vsDimension = vsDim;
		this.radius = r;

		// seed the m4ri random number generator
		Random generator = new Random();
		Mzd.srandom(generator.nextInt());

		// construct the identity matrix for the full vector space V
		this.V = new Mzd(vsDimension, vsDimension);
		for (int i=0; i<vsDimension; i++) {
			this.V.writeBit(i, i, 1);
		}

		// get the standard basis.  Always handy.
		this.basis = Mzd.standardBasis(vsDimension);
		
		// the DComplex containing "maximal" simplices up to dim 2
		this.maximalSimplices = new DComplex(DrainageExperiment.vsDimension);
	}

	
	/**
	 * Adds up the total size of the inventory of all the points.
	 * The total storage used across the whole network.
	 * 
	 * @return
	 */
	public int totalInventorySize() {
		int totalSize = 0;
		for (DPoint p : this.points) {
			totalSize += p.M.size();
		}
		return totalSize;
	}
	
	
	/**
	 * Note: Must construct the Rips complex before you call this function!!
	 * 
	 * Compute a new inv entory for p such that it is the least number
	 * of vectors necessary to support all of it's maximal or 2-dim cofaces (simplices containing it)
	 * 
	 * This version has the option to fill the inventory too, in the case that it has
	 * an un-spanned coface it raises itself to the full basis and then starts draining.
	 * 
	 * This version drains at most nVectors vectors from the inventory.  It's a little less
	 * efficient this way, because a lot of the work in this method is inital setup of the data
	 * structures, and we don't keep any of that around (can't really).
	 * 
	 * It returns the number of vectors drained.
	 * 
	 * If you pass in nVectors=0 it drains the maximum possible.
	 * 
	 * @param p
	 * @param nVectors
	 * @param fillVertices
	 * @return
	 */
	public int drainVertex(DPoint p, int nVectors) {
		return drainVertex(p, nVectors, false);
	}
	public int drainVertex(DPoint p, int nVectors, boolean fillVertices) {
		System.out.println("----------------\ndrainVertex("+p+")");
	    if (p==null) { return 0; }
	    
	    System.out.println("initial inventory:");
	    for (Mzd v : p.M){
	    	v.print();
	    }
	    
		// an argument of 0 means drain the maximum possible
		if (nVectors == 0) {
			//nVectors = p.M.size();
			nVectors = p.vsDimension;
		}
		
		// if p has no neighbors there's nothing we can do.
		if (p.nbrs.isEmpty()) { return 0; }
		
		// CRUNCH!!!
		// start with the full inventory and reduce it subject to the constraint that is complements U
		Mzd BB = new Mzd(p.M.size(), p.vsDimension);
		for (int i=0; i<p.M.size(); i++) {
			Mzd.copyRow(BB, i, p.M.get(i), 0);
		}
		
		/*
		// figure out how much of the space is spanned to start with
		Mzd tst = new Mzd(2*p.vsDimension, p.vsDimension);
		Mzd.copyRows(tst, 0, U, 0, U.getNrows());
		Mzd.copyRows(tst, p.vsDimension, BB, 0, BB.getNrows());
		int initialRank = tst.echelonize(false);
		tst.destroy();
		*/
		
		// if we're supposed to fill incomplete simplices...
		if (fillVertices) {
			for (DSimplex sigma : p.ripsCofaces) {
				//Mzd tst = new Mzd(2*p.vsDimension, p.vsDimension);
				Mzd U = sigma.exclusiveBasis(p);
				Mzd tst = Mzd.stack(U, BB);
				//System.out.println("U:");
				//U.print();
				//System.out.println("BB:");
				//BB.print();
				//System.out.println("tst:");
				//tst.print();
				//Mzd.copyRows(tst, 0, U, 0, U.getNrows());
				//Mzd.copyRows(tst, p.vsDimension, BB, 0, BB.getNrows());
				int newrank = tst.echelonize(false);
				tst.destroy();
				U.destroy();
				
				// if the rank was reduced bail out
				if (newrank != p.vsDimension) {
					BB.destroy();
					BB = Mzd.identityMatrix(p.vsDimension);
					p.setInventory(BB);
					System.out.println("WARNING: drainVertex() - simplex has unsupported cofaces - filling with standard basis");
					break;
				}
			}
		}
		//System.out.println("done testing for uncovered cofaces.");
		
		int numDrained = 0;
		//Mzd tst = new Mzd(2*p.vsDimension, p.vsDimension);
		Mzd iRow = new Mzd(1, p.vsDimension);
		
		for (int i=0; i<BB.getNrows(); i++) {
			// try adding it to each other vector.
			// move on to next if we succeed
			
			// this is the row we will be trying to add to other things
			// make a back-up copy of it
			Mzd.copyRow(iRow, 0, BB, i);
			
			// if this row is zero don't bother with it
			//System.out.println("iRow(i="+i+"):");
			//iRow.print();
			if (iRow.isZero()) {
				System.out.println("skipping row "+i+" because it is zero");
				continue;
			}
			
			// when i=j we are testing the case of just removing row i
			for (int j=i; j<BB.getNrows(); j++) {
				
				//System.out.println("Test-adding "+i+"+"+j+" ...");
				BB.rowAdd(i, j);
				BB.rowClearOffset(i, 0);

				// test all the possible cofaces to see if this change leaves them less covered
				boolean addFail = false;
				for (DSimplex sigma : p.ripsCofaces) {
					
					//Mzd tst = new Mzd(2*p.vsDimension, p.vsDimension);
					Mzd U = sigma.exclusiveBasis(p);
					Mzd tst = Mzd.stack(U, BB);
					//Mzd.copyRows(tst, 0, U, 0, U.getNrows());
					//Mzd.copyRows(tst, p.vsDimension, BB, 0, BB.getNrows());
					int newrank = tst.echelonize(false);
					tst.destroy();
					U.destroy();
					
					// if the rank was reduced bail out
					if (newrank != p.vsDimension) {
						addFail = true;
						break;
					}
					
				}
				
				// if non of the cofaces were destoyed we have success
				// lock it in by bailing out of this loop
				if (! addFail) {
					numDrained++;
					break;
				}
				
				// otherwise try adding to the next row, but first undo the change
				Mzd.copyRow(BB, i, iRow, 0);
				if (i!=j) {
					BB.rowAdd(i, j);
				}
			}
			if (numDrained >= nVectors) {
				break;
			}
		}
		
		// clean up before letting it go out of scope
		iRow.destroy();
		//tst.destroy();
		
		System.out.println("new inventory:");
		BB.print();
		
		// Note: BB may have several zero rows at the top.
		// That is accounted for in setInventory()
		p.clearInventory();
		p.setInventory(BB);
		
		BB.destroy();
		return numDrained;
	}
	
	
	/**
	 * Does what it says
	 * The neighbor set info is local and is stored inside each DPoint
	 * 
	 * (right now) the neighbor set does not include the point itself...
	 */
	public void computeNeighborSet() {
		double r2 = this.radius*this.radius;

		for (int i=0; i<points.size(); i++) {
			DPoint p1 = points.get(i);
			for (int j=(i+1); j<points.size(); j++) {
				DPoint p2 = points.get(j);
				double d2 = (p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y);
				if (d2 <= r2) {
					p1.nbrsUp.add(p2);
					p1.nbrs.add(p2);
					p2.nbrs.add(p1);
				}
			}
		}
	}


	/**
	 * Method to build the Rips complex of the points
	 * Right now up to dimension 3
	 * Include no filtration information right now
	 * 
	 * Also fills in the DComplex of "maximals"
	 */
	void buildRipsComplex() {
		ripsComplexStream = Plex4.createExplicitSimplexStream();
		
		// add the 0-simplices
		for (DPoint p1 : points) {
			//System.out.println("RIPS: addSimplex [ "+p1.index+" ]");
			ripsComplexStream.addVertex(p1.index, 0);
			
			// see if this vertex is itself maximal
			if (p1.nbrs.size() == 0) {
				DPoint[] verts = { p1 };
				DSimplex ds = new DSimplex(Arrays.asList(verts));
				this.maximalSimplices.addSimplex(ds);
				p1.ripsCofaces.add(ds);
			}
		}
		
		// add the 1-simplices
		for (DPoint p1 : points) {
			for (DPoint p2 : p1.nbrsUp) {
				//System.out.println("RIPS: addSimplex [ "+p1.index+"  "+p2.index+" ]");
				ripsComplexStream.addElement(new int[]{p1.index, p2.index}, 0);
				
				// see if this vertex is itself maximal
				HashSet<DPoint> commonNbrs = new HashSet<DPoint>(p1.nbrs);
				commonNbrs.retainAll(p2.nbrs);
				if (commonNbrs.size() == 0) {
					DPoint[] verts = { p1 , p2 };
					DSimplex ds = new DSimplex(Arrays.asList(verts));
					this.maximalSimplices.addSimplex(ds);
					p1.ripsCofaces.add(ds);
					p2.ripsCofaces.add(ds);
				}
			}
		}
		
		// add the 2-simplices
		for (DPoint p1 : points) {
			for (DPoint p2 : p1.nbrsUp) {
				HashSet<DPoint> nbrs2 = new HashSet<DPoint>(p2.nbrsUp);
				nbrs2.retainAll(p1.nbrsUp);
				for (DPoint p3 : nbrs2) {
					//System.out.println("RIPS: addSimplex [ "+p1.index+"  "+p2.index+"  "+p3.index+" ]");
					ripsComplexStream.addElement(new int[]{p1.index, p2.index, p3.index}, 0);
					
					// we'll consider all 2-simplices to be maximal
					// for the purposes of DComplex
					HashSet<DPoint> commonNbrs = new HashSet<DPoint>(p1.nbrs);
					commonNbrs.retainAll(p2.nbrs);
					commonNbrs.retainAll(p3.nbrs);
					//if (commonNbrs.size() == 0) {
						DPoint[] verts = { p1 , p2 , p3 };
						DSimplex ds = new DSimplex(Arrays.asList(verts));
						this.maximalSimplices.addSimplex(ds);
						p1.ripsCofaces.add(ds);
						p2.ripsCofaces.add(ds);
						p3.ripsCofaces.add(ds);
					//}
				}
			}
		}
		
		// add the 3-simplices
		for (DPoint p1 : points) {
			for (DPoint p2 : p1.nbrsUp) {
				HashSet<DPoint> nbrs2 = new HashSet<DPoint>(p2.nbrsUp);
				nbrs2.retainAll(p1.nbrsUp);
				for (DPoint p3 : nbrs2) {
					HashSet<DPoint> nbrs3 = new HashSet<DPoint>(p3.nbrsUp);
					nbrs3.retainAll(nbrs2);
					for (DPoint p4 : nbrs3) {
						//System.out.println("RIPS: addSimplex [ "+p1.index+"  "+p2.index+"  "+p3.index+"  "+p4.index+" ]");
						ripsComplexStream.addElement(new int[]{p1.index, p2.index, p3.index, p4.index}, 0);
					}
				}
			}
		}
		
		// in case we made a mistake...
		ripsComplexStream.ensureAllFaces();
		ripsComplexStream.finalizeStream();
	}


	/**
	 * Decide if a simplex should be allowed into RC(U)
	 * Namely, do the vertices and their common neighbors together span the
	 * vector space U?
	 * 
	 * @param U
	 * @param vertices
	 * @param nbrs
	 * @return
	 */
	boolean isSimplexCovered(Mzd U, HashSet<DPoint> vertices, HashSet<DPoint> nbrs) {
		//System.out.println("isSimplexCovered()");
		// count up the total vectors amongst the common neighbors
		int total_vectors = 0;
		for (DPoint n : vertices) { total_vectors += n.rank; }
		for (DPoint n : nbrs) { total_vectors += n.rank; }
		//System.out.println("\t total_vectors="+total_vectors);

		// put them all in an m4ri matrix
		// include extra rows to copy the matrix U in as a last step
		Mzd simplex_span = new Mzd(total_vectors+U.getNrows(), vsDimension);
		int row = 0;
		for (DPoint n : vertices) {
			Mzd.copyRows(simplex_span, row, n.reducedBasis, 0, n.rank);
			row += n.rank;
		}
		for (DPoint n : nbrs) {
			Mzd.copyRows(simplex_span, row, n.reducedBasis, 0, n.rank);
			row += n.rank;
		}

		// compute the rank of the nbhd
		//System.out.println("\t simplex_span:");
		//simplex_span.print();
		int rk = simplex_span.echelonize(false);
		//System.out.println("\t simplex_span (reduced):");
		//simplex_span.print();

		// concatenate the basis for U and see if it increases the rank
		Mzd.copyRows(simplex_span, row, U, 0, U.getNrows());
		row += U.getNrows();
		//System.out.println("\t simplex_span after copying in U:");
		//simplex_span.print();

		// if U didn't increase the rank then span(U) is contained in span(nbhd)
		int rk2 = simplex_span.echelonize(false);
		//System.out.println("rk="+rk+"\t rk2="+rk2);

		// before we return, clean up this temp matrix
		simplex_span.destroy();

		if (rk==rk2) {
			//System.out.println("simplex is covered");
			return true;
		}
		//System.out.println("simplex is NOT covered");
		return false;
	}


	/**
	 * Build the coverage complex of the points
	 * There are lots of things here could certainly be done more efficiently:
	 * - start with high-dim simplices and only test their faces if the k-simplex fails
	 * - keep temporary structures for vs spanned and neighbor sets as you go through loops
	 * 
	 * @param U - The vector space to be spanned.  Represented by an m4ri matrix
	 * 			   where the rows (or should it be cols) contain the vectors.
	 */
	void buildCoverageComplex(Mzd U) {
		//System.out.println("buildCoverageComplex\n  U:");
		//U.print();
		//System.out.println("  rank="+U.echelonize(false));

		coverageComplexStream = Plex4.createExplicitSimplexStream();
		HashSet<DPoint> vertices = new HashSet<DPoint>(10);

		// the initial rank of the basis for U
		//int rk = U.echelonize(false);

		// add the 0-simplices
		for (DPoint p1 : points) {
			vertices.add(p1);
			if (isSimplexCovered(U, vertices, p1.nbrs)) {
				//System.out.println("COVERAGE: addSimplex [ "+p1.index+" ]");
				coverageComplexStream.addVertex(p1.index, 0);
			}
			vertices.remove(p1);
		}

		// add the 1-simplices
		for (DPoint p1 : points) {
			vertices.add(p1);
			for (DPoint p2 : p1.nbrsUp) {
				vertices.add(p2);

				// compute the intersection of their neighbor sets
				HashSet<DPoint> common_nbrs = (HashSet<DPoint>)p1.nbrs.clone();
				common_nbrs.retainAll(p2.nbrs);

				if (isSimplexCovered(U, vertices, common_nbrs)) {
					//System.out.println("COVERAGE: addSimplex [ "+p1.index+"  "+p2.index+" ]");
					coverageComplexStream.addElement(new int[]{p1.index, p2.index}, 0);
				}
				vertices.remove(p2);
			}
			vertices.remove(p1);
		}

		// add the 2-simplices
		for (DPoint p1 : points) {
			vertices.add(p1);
			for (DPoint p2 : p1.nbrsUp) {
				vertices.add(p2);
				HashSet<DPoint> nbrs2 = new HashSet<DPoint>(p2.nbrsUp);
				nbrs2.retainAll(p1.nbrsUp);
				for (DPoint p3 : nbrs2) {
					vertices.add(p3);

					// compute the intersection of their neighbor sets
					HashSet<DPoint> common_nbrs = (HashSet<DPoint>)p1.nbrs.clone();
					common_nbrs.retainAll(p2.nbrs);
					common_nbrs.retainAll(p3.nbrs);

					if (isSimplexCovered(U, vertices, common_nbrs)) {
						//System.out.println("COVERAGE: addSimplex [ "+p1.index+"  "+p2.index+"  "+p3.index+" ]");
						coverageComplexStream.addElement(new int[]{p1.index, p2.index, p3.index}, 0);
					}
					vertices.remove(p3);
				}
				vertices.remove(p2);
			}
			vertices.remove(p1);
		}

		// add the 3-simplices
		for (DPoint p1 : points) {
			vertices.add(p1);
			for (DPoint p2 : p1.nbrsUp) {
				vertices.add(p2);
				HashSet<DPoint> nbrs2 = new HashSet<DPoint>(p2.nbrsUp);
				nbrs2.retainAll(p1.nbrsUp);
				for (DPoint p3 : nbrs2) {
					vertices.add(p3);
					HashSet<DPoint> nbrs3 = new HashSet<DPoint>(p3.nbrsUp);
					nbrs3.retainAll(nbrs2);
					for (DPoint p4 : nbrs3) {
						vertices.add(p4);

						// compute the intersection of their neighbor sets
						HashSet<DPoint> common_nbrs = (HashSet<DPoint>)p1.nbrs.clone();
						common_nbrs.retainAll(p2.nbrs);
						common_nbrs.retainAll(p3.nbrs);
						common_nbrs.retainAll(p4.nbrs);

						if (isSimplexCovered(U, vertices, common_nbrs)) {
							//System.out.println("COVERAGE: addSimplex [ "+p1.index+"  "+p2.index+"  "+p3.index+"  "+p4.index+" ]");
							coverageComplexStream.addElement(new int[]{p1.index, p2.index, p3.index, p4.index}, 0);
						}
						vertices.remove(p4);
					}
					vertices.remove(p3);
				}
				vertices.remove(p2);
			}
			vertices.remove(p1);
		}

		coverageComplexStream.ensureAllFaces();
		coverageComplexStream.finalizeStream();
	}


	/**
	 * Build the filtered comples where the coverage complex is the first stage of the
	 * filtration and the Rips complex is the 2nd stage.
	 * 
	 * There are lots of things here could certainly be done more efficiently:
	 * - start with high-dim simplices and only test their faces if the k-simplex fails
	 * - keep temporary structures for vs spanned and neighbor sets as you go through loops
	 * 
	 * @param U - The vector space to be spanned.  Represented by an m4ri matrix
	 * 			   where the rows (or should it be cols) contain the vectors.
	 */
	void buildCoverageRipsComplex(Mzd U) {
		//System.out.println("buildCoverageComplex\n  U:");
		//U.print();
		//System.out.println("  rank="+U.echelonize(false));

		coverageRipsComplexStream = Plex4.createExplicitSimplexStream();
		HashSet<DPoint> vertices = new HashSet<DPoint>(10);

		// the initial rank of the basis for U
		//int rk = U.echelonize(false);

		// add the 0-simplices
		for (DPoint p1 : points) {
			vertices.add(p1);
			if (isSimplexCovered(U, vertices, p1.nbrs)) {
				//System.out.println("COVERAGE_RIPS: addSimplex [ "+p1.index+" ]\t filt=0");
				coverageRipsComplexStream.addVertex(p1.index, 0);
			} else {
				//System.out.println("COVERAGE_RIPS: addSimplex [ "+p1.index+" ]\t filt=1");
				coverageRipsComplexStream.addVertex(p1.index, 1);
			}
			vertices.remove(p1);
		}

		// add the 1-simplices
		for (DPoint p1 : points) {
			vertices.add(p1);
			for (DPoint p2 : p1.nbrsUp) {
				vertices.add(p2);

				// compute the intersection of their neighbor sets
				HashSet<DPoint> common_nbrs = (HashSet<DPoint>)p1.nbrs.clone();
				common_nbrs.retainAll(p2.nbrs);

				if (isSimplexCovered(U, vertices, common_nbrs)) {
					//System.out.println("COVERAGE: addSimplex [ "+p1.index+"  "+p2.index+" ]\t filt=0");
					coverageRipsComplexStream.addElement(new int[]{p1.index, p2.index}, 0);
				} else {
					//System.out.println("COVERAGE: addSimplex [ "+p1.index+"  "+p2.index+" ]\t filt=1");
					coverageRipsComplexStream.addElement(new int[]{p1.index, p2.index}, 1);
				}
				vertices.remove(p2);
			}
			vertices.remove(p1);
		}

		// add the 2-simplices
		for (DPoint p1 : points) {
			vertices.add(p1);
			for (DPoint p2 : p1.nbrsUp) {
				vertices.add(p2);
				HashSet<DPoint> nbrs2 = new HashSet<DPoint>(p2.nbrsUp);
				nbrs2.retainAll(p1.nbrsUp);
				for (DPoint p3 : nbrs2) {
					vertices.add(p3);

					// compute the intersection of their neighbor sets
					HashSet<DPoint> common_nbrs = (HashSet<DPoint>)p1.nbrs.clone();
					common_nbrs.retainAll(p2.nbrs);
					common_nbrs.retainAll(p3.nbrs);

					if (isSimplexCovered(U, vertices, common_nbrs)) {
						//System.out.println("COVERAGE: addSimplex [ "+p1.index+"  "+p2.index+"  "+p3.index+" ]\t filt=0");
						coverageRipsComplexStream.addElement(new int[]{p1.index, p2.index, p3.index}, 0);
					} else {
						//System.out.println("COVERAGE: addSimplex [ "+p1.index+"  "+p2.index+"  "+p3.index+" ]\t filt=1");
						coverageRipsComplexStream.addElement(new int[]{p1.index, p2.index, p3.index}, 1);
					}
					vertices.remove(p3);
				}
				vertices.remove(p2);
			}
			vertices.remove(p1);
		}

		// add the 3-simplices
		for (DPoint p1 : points) {
			vertices.add(p1);
			for (DPoint p2 : p1.nbrsUp) {
				vertices.add(p2);
				HashSet<DPoint> nbrs2 = new HashSet<DPoint>(p2.nbrsUp);
				nbrs2.retainAll(p1.nbrsUp);
				for (DPoint p3 : nbrs2) {
					vertices.add(p3);
					HashSet<DPoint> nbrs3 = new HashSet<DPoint>(p3.nbrsUp);
					nbrs3.retainAll(nbrs2);
					for (DPoint p4 : nbrs3) {
						vertices.add(p4);

						// compute the intersection of their neighbor sets
						HashSet<DPoint> common_nbrs = (HashSet<DPoint>)p1.nbrs.clone();
						common_nbrs.retainAll(p2.nbrs);
						common_nbrs.retainAll(p3.nbrs);
						common_nbrs.retainAll(p4.nbrs);

						if (isSimplexCovered(U, vertices, common_nbrs)) {
							//System.out.println("COVERAGE: addSimplex [ "+p1.index+"  "+p2.index+"  "+p3.index+"  "+p4.index+" ]\t filt=0");
							coverageRipsComplexStream.addElement(new int[]{p1.index, p2.index, p3.index, p4.index}, 0);
						} else {
							//System.out.println("COVERAGE: addSimplex [ "+p1.index+"  "+p2.index+"  "+p3.index+"  "+p4.index+" ]\t filt=1");
							coverageRipsComplexStream.addElement(new int[]{p1.index, p2.index, p3.index, p4.index}, 1);
						}
						vertices.remove(p4);
					}
					vertices.remove(p3);
				}
				vertices.remove(p2);
			}
			vertices.remove(p1);
		}

		coverageRipsComplexStream.ensureAllFaces();
		coverageRipsComplexStream.finalizeStream();
	}



	/**
	 * Computes the persistent homology of the 2-step filtration:
	 *  RC-->R
	 * This allows us to distinguish holes that only exist in the coverage complex
	 * from holes that are also in the Rips complex that we can't hope to cover.
	 * 
	 */
	public Vector<ExplicitSimplexStream> computePersistentHomology() {
		ExplicitSimplexStream stream = coverageRipsComplexStream;
		stream.finalizeStream();

		AbstractPersistenceBasisAlgorithm<Simplex,IntSparseFormalSum<Simplex>> persistence = new IntAbsoluteHomology<Simplex>(ModularIntField.getInstance(2), SimplexComparator.getInstance(), 0, 2);
		System.out.println("Got persistence object...");

		// compute and print the intervals
		BarcodeCollection<Double> intervals = persistence.computeIntervals(stream);
		System.out.println("\nBarcodes for stream: ");
		System.out.println(intervals);

		// compute and print the intervals annotated with a representative cycle
		//AnnotatedBarcodeCollection<Integer,IntSparseFormalSum<Simplex>> circle_intervals2 = persistence.   .computeAnnotatedIntervals(stream);

		//AnnotatedBarcodeCollection<Integer, IntSparseFormalSum<Simplex>> 
		AnnotatedBarcodeCollection<Double, IntSparseFormalSum<Simplex>> annotated_intervals;
		annotated_intervals = persistence.computeAnnotatedIntervals(stream);
		System.out.println("\nAnnotated barcodes for stream: ");
		System.out.println(annotated_intervals);

		AnnotatedBarcodeCollection<Double, IntSparseFormalSum<Simplex>> infinite_intervals = annotated_intervals.getInfiniteIntervals();
		System.out.println("\nJust the infinite intervals for stream: ");
		System.out.println(infinite_intervals);

		List<IntSparseFormalSum<Simplex>> oneCycles = infinite_intervals.getGeneratorsAtDimension(1);

		//System.out.println("\nManually iterate through 1-cycles: ");
		Vector<ExplicitSimplexStream> cycleStreams = new Vector<ExplicitSimplexStream>(oneCycles.size());
		for (IntSparseFormalSum<Simplex> c : oneCycles) {
			//System.out.println(c);

			ExplicitSimplexStream str = new ExplicitSimplexStream();

			for (Simplex ss : c.getSummands()) {
				//System.out.println(ss);
				str.addElement(ss,0);
			}
			cycleStreams.add(str);
		}

		return cycleStreams;
	}


	/**
	 * Compute the homology of the coverage complex
	 * 
	 * Set this.cycleStream to a stream representing the 1-cycles.
	 */
	//public Vector<ExplicitSimplexStream> computeHomology() {
	public void computeHomology() {
		ExplicitSimplexStream stream = coverageComplexStream;
		stream.finalizeStream();

		AbstractPersistenceBasisAlgorithm<Simplex,IntSparseFormalSum<Simplex>> persistence = new IntAbsoluteHomology<Simplex>(ModularIntField.getInstance(2), SimplexComparator.getInstance(), 0, 2);
		System.out.println("Got persistence object...");

		// compute and print the intervals
		BarcodeCollection<Double> intervals = persistence.computeIntervals(stream);
		System.out.println("\nBarcodes for stream: ");
		System.out.println(intervals);

		// compute and print the intervals annotated with a representative cycle
		//AnnotatedBarcodeCollection<Integer,IntSparseFormalSum<Simplex>> circle_intervals2 = persistence.   .computeAnnotatedIntervals(stream);

		//AnnotatedBarcodeCollection<Integer, IntSparseFormalSum<Simplex>> 
		AnnotatedBarcodeCollection<Double, IntSparseFormalSum<Simplex>> annotated_intervals;
		annotated_intervals = persistence.computeAnnotatedIntervals(stream);
		System.out.println("\nAnnotated barcodes for stream: ");
		System.out.println(annotated_intervals);

		AnnotatedBarcodeCollection<Double, IntSparseFormalSum<Simplex>> infinite_intervals = annotated_intervals.getInfiniteIntervals();
		System.out.println("\nJust the infinite intervals for stream: ");
		System.out.println(infinite_intervals);

		List<IntSparseFormalSum<Simplex>> oneCycles = infinite_intervals.getGeneratorsAtDimension(1);

		//System.out.println("\nManually iterate through 1-cycles: ");
		//Vector<ExplicitSimplexStream> cycleStreams = new Vector<ExplicitSimplexStream>(oneCycles.size());
		ExplicitSimplexStream str = new ExplicitSimplexStream();
		for (IntSparseFormalSum<Simplex> c : oneCycles) {
			//System.out.println(c);

			//ExplicitSimplexStream str = new ExplicitSimplexStream();

			for (Simplex ss : c.getSummands()) {
				//System.out.println(ss);
				str.addElement(ss,0);
			}
			//cycleStreams.add(str);
		}

		//return cycleStreams;
		this.cycleStream = str;
	}


	/**
	 * Use plex-viewer to draw the rips and coverage complexes
	 * 
	 * This version uses our specialized drawCoverageComplex() method in plex-viewer
	 * to show the Rips, Coverage, and Cycles in differetn colors.
	 * 
	 */
	public void drawComplex() {
		System.out.println("drawComplex()");

		// first put the points location data into an array of double
		domainPoints = new double[points.size()][];
		double xMean = 0.0;
		double yMean = 0.0;
		for (int i=0; i<points.size(); i++) {
			domainPoints[i] = new double[3];
			domainPoints[i][0] = points.get(i).x;
			domainPoints[i][1] = points.get(i).y;
			domainPoints[i][2] = points.get(i).isElevated() ? 0.5 : 0.0;
			xMean += domainPoints[i][0];
			yMean += domainPoints[i][1];		
		}
		xMean = xMean/points.size();
		yMean = yMean/points.size();

		// adjust the points so they're a little better centered
		for (int i=0; i<points.size(); i++) {
			domainPoints[i][0] -= xMean;
			domainPoints[i][1] -= yMean;
		}

		// pass it to the plex-viewer API along with the stream
		//Api.drawSimplexStream(ripsComplexStream, domainPoints);
		//Api.drawSimplexStream(ss, domainPoints);
		Api.drawCoverageComplex(ripsComplexStream, coverageComplexStream, ripsCycleStream, cycleStream, domainPoints);
	}


	/**
	 * Use plex-viewer to draw the rips and coverage complexes
	 * 
	 * This version takes a filtered complex representing the RC-->R inclusion.
	 */
	public void drawFilteredComplex(ExplicitSimplexStream ss) {
		System.out.println("drawComplex()");
		if (ss==null || coverageComplexStream==null) {
			System.out.println("ERROR: drawComplex() on a null stream!");
			return;
		}

		// first put the points location data into an array of double
		domainPoints = new double[points.size()][];
		double xMean = 0.0;
		double yMean = 0.0;
		for (int i=0; i<points.size(); i++) {
			domainPoints[i] = new double[2];
			domainPoints[i][0] = points.get(i).x;
			domainPoints[i][1] = points.get(i).y;
			xMean += domainPoints[i][0];
			yMean += domainPoints[i][1];		
		}
		xMean = xMean/points.size();
		yMean = yMean/points.size();

		// adjust the points so they're a little better centered
		for (int i=0; i<points.size(); i++) {
			domainPoints[i][0] -= xMean;
			domainPoints[i][1] -= yMean;
		}

		// pass it to the plex-viewer API along with the stream
		Api.drawSimplexStream(ss, domainPoints);
	}



}
