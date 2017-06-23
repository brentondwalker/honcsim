/*
 * CoverageExperiment.java
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
 * CoverageExperiment
 * 
 * A base class to handle some common functionality of coverage complex experiments
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
import edu.stanford.math.primitivelib.autogen.pair.ObjectObjectPair;
import edu.stanford.math.primitivelib.algebraic.impl.*;
import edu.stanford.math.plex_viewer.*;


public class CoverageExperiment {

    /*
     * class data
     */
    int vsDimension = 3;  // this always get overwritten
    double radius = 1.0;

    // identity matrix of vsDimension
    // and standard basis
    Mzd V = null;
    Mzd[] basis = null;
    
    // data structure containing DPoints
    Vector<DPoint> points = null;
    
    // javaplex objects
    ExplicitSimplexStream ripsComplexStream = null;
    ExplicitSimplexStream coverageComplexStream = null;
    ExplicitSimplexStream cycleStream = null;           // the cycles only in coverage complex (finite intervals)
    ExplicitSimplexStream ripsCycleStream = null;       // the cycles also in Rips complex (infinite intervals)
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
    public CoverageExperiment(int vsDim, double r) {
        this.vsDimension = vsDim;
        this.radius = r;

        // seed the m4ri random number generator
        Random generator = new Random();
        Mzd.srandom(generator.nextInt());
        //Mzd.srandom(1);
        
        // construct the identity matrix for the full vector space V
        this.V = new Mzd(vsDimension, vsDimension);
        for (int i=0; i<vsDimension; i++) {
            this.V.writeBit(i, i, 1);
        }

        // get the standard basis.  Always handy.
        this.basis = Mzd.standardBasis(vsDimension);
    }




    /**
     * Does what it says
     * The neighbor set info is local and is stored inside each DPoint
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
     */
    void buildRipsComplex() {
        ripsComplexStream = Plex4.createExplicitSimplexStream();

        // add the 0-simplices
        for (DPoint p1 : points) {
            //System.out.println("RIPS: addSimplex [ "+p1.index+" ]");
            ripsComplexStream.addVertex(p1.index, 0);
        }

        // add the 1-simplices
        for (DPoint p1 : points) {
            for (DPoint p2 : p1.nbrsUp) {
                //System.out.println("RIPS: addSimplex [ "+p1.index+"  "+p2.index+" ]");
                ripsComplexStream.addElement(new int[]{p1.index, p2.index}, 0);
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
        buildCoverageComplex(U, false);
    }


    /**
     * This adds a debatable improvement.
     * 
     * This method optionally excludes 1-simplices (edges) from the coverage
     * complex if they are not the face of a 2-simplex.
     * 
     * In computing coverage complexes we often get 1-simplexes (edges) that are
     * are not bounding any 2-simplex.  We could argue that these do not represent
     * any coverage, since they are 1-dimensional.
     * 
     * Also, in examples like CoverageGrid or FencedCoverageSquare, we often end
     * up with apparent large coverage holes with a web of 1-simplexes running
     * over them, giving a huge number of 1-cycles, when visually, it looks like
     * there should be only one.  Really these smaller holes would give short
     * persistence intervals in a Rips-filtered complex.  But since we're using
     * a persistence computation to distinguish network coverage and data coverage
     * holes in the rips complex --> coverage complex filtration, there isn't a
     * simple way to compute persistent homology along the dual filtrations.
     * There probably is a way, but maybe not a simple way.
     * 
     */
    void buildCoverageComplex(Mzd U, boolean excludeCoverageChords) {
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

        // add the 1-simplices (unless we are excluding coverage chords)
        // if we exclude the 1-simplices, the ones that are faces of a covered
        // 2-simplex will be added at the end anyway.
        if (! excludeCoverageChords) {
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
     * Build the filtered complex where the coverage complex is the first stage of the
     * filtration and the Rips complex is the 2nd stage.
     * 
     * This version uses the (probably inefficient) method of iterating
     * over all simplices in the Rips complex, and checking if each one
     * is also in the coverage complex.
     * 
     * If the Rips or coverage complexes have not already been computed, this
     * method will call the appropriate methods to compute them.
     * 
     * @param U - The vector space to be spanned.  Represented by an m4ri matrix
     * 			   where the rows (or should it be cols) contain the vectors.
     */
    void buildCoverageRipsComplex(Mzd U) {
        buildCoverageRipsComplex(U, false);
    }
    
    
    /**
     * This version includes the debatable enhancement explained above, to
     * optionally exclude 1-simplices that are not faces of a covered 2-simplex
     * from the coverage complex.
     * 
     * @param U
     * @param excludeCoverageChords
     */
    void buildCoverageRipsComplex(Mzd U, boolean excludeCoverageChords) {
        if (ripsComplexStream == null) {
            this.buildRipsComplex();
        }
        
        if (coverageComplexStream == null) {
            this.buildCoverageComplex(U, excludeCoverageChords);
        }
        
        coverageRipsComplexStream = Plex4.createExplicitSimplexStream();

        for (Simplex s : ripsComplexStream) {
            if (! excludeCoverageChords || s.getDimension() != 1) {
                if (coverageComplexStream.containsElement(s)) {
                    coverageRipsComplexStream.addElement(s, 0);
                } else {
                    coverageRipsComplexStream.addElement(s, 1);
                }
            }
        }
        coverageRipsComplexStream.ensureAllFaces();
        coverageRipsComplexStream.finalizeStream();
    }
        
    
    /**
     * Build all the complexes, ripsComplexStream, coverageComplexStream, and
     * the filtered coverageRipsComplexStream.
     * 
     * This version includes the debatable enhancement explained above, to
     * optionally exclude 1-simplices that are not faces of a covered 2-simplex
     * from the coverage complex.
     * 
     * @param U
     */
    void buildAllComplexes(Mzd U) {
        buildAllComplexes(U, false);
    }

    
    /**
     * This version includes the debatable enhancement explained above, to
     * optionally exclude 1-simplices that are not faces of a covered 2-simplex
     * from the coverage complex.
     * 
     * There are lots of things here could certainly be done more efficiently:
     * - start with high-dim simplices and only test their faces if the k-simplex fails
     * - keep temporary structures for vs spanned and neighbor sets as you go through loops

     * @param U
     * @param excludeCoverageChords
     */
    void buildAllComplexes(Mzd U, boolean excludeCoverageChords) {
        //System.out.println("buildCoverageComplex\n  U:");
        //U.print();
        //System.out.println("  rank="+U.echelonize(false));
        
        ripsComplexStream = Plex4.createExplicitSimplexStream();
        coverageComplexStream = Plex4.createExplicitSimplexStream();
        coverageRipsComplexStream = Plex4.createExplicitSimplexStream();
        HashSet<DPoint> vertices = new HashSet<DPoint>(10);

        
        // the initial rank of the basis for U
        //int rk = U.echelonize(false);

        // add the 0-simplices
        for (DPoint p1 : points) {
            vertices.add(p1);
            ripsComplexStream.addVertex(p1.index, 0);
            if (isSimplexCovered(U, vertices, p1.nbrs)) {
                //System.out.println("COVERAGE_RIPS: addSimplex [ "+p1.index+" ]\t filt=0");
                coverageRipsComplexStream.addVertex(p1.index, 0);
                coverageComplexStream.addVertex(p1.index, 0);
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

                ripsComplexStream.addElement(new int[]{p1.index, p2.index}, 0);
                
                // compute the intersection of their neighbor sets
                HashSet<DPoint> common_nbrs = (HashSet<DPoint>)p1.nbrs.clone();
                common_nbrs.retainAll(p2.nbrs);
                
                if (isSimplexCovered(U, vertices, common_nbrs) && (! excludeCoverageChords)) {
                    //System.out.println("COVERAGE: addSimplex [ "+p1.index+"  "+p2.index+" ]\t filt=0");
                    coverageRipsComplexStream.addElement(new int[]{p1.index, p2.index}, 0);
                    coverageComplexStream.addElement(new int[]{p1.index, p2.index}, 0);
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
                    
                    ripsComplexStream.addElement(new int[]{p1.index, p2.index, p3.index}, 0);

                    // compute the intersection of their neighbor sets
                    HashSet<DPoint> common_nbrs = (HashSet<DPoint>)p1.nbrs.clone();
                    common_nbrs.retainAll(p2.nbrs);
                    common_nbrs.retainAll(p3.nbrs);

                    if (isSimplexCovered(U, vertices, common_nbrs)) {
                        //System.out.println("COVERAGE: addSimplex [ "+p1.index+"  "+p2.index+"  "+p3.index+" ]\t filt=0");
                        coverageRipsComplexStream.addElement(new int[]{p1.index, p2.index, p3.index}, 0);
                        coverageComplexStream.addElement(new int[]{p1.index, p2.index, p3.index}, 0);
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
                        
                        ripsComplexStream.addElement(new int[]{p1.index, p2.index, p3.index, p4.index}, 0);

                        // compute the intersection of their neighbor sets
                        HashSet<DPoint> common_nbrs = (HashSet<DPoint>)p1.nbrs.clone();
                        common_nbrs.retainAll(p2.nbrs);
                        common_nbrs.retainAll(p3.nbrs);
                        common_nbrs.retainAll(p4.nbrs);

                        if (isSimplexCovered(U, vertices, common_nbrs)) {
                            //System.out.println("COVERAGE: addSimplex [ "+p1.index+"  "+p2.index+"  "+p3.index+"  "+p4.index+" ]\t filt=0");
                            coverageRipsComplexStream.addElement(new int[]{p1.index, p2.index, p3.index, p4.index}, 0);
                            coverageComplexStream.addElement(new int[]{p1.index, p2.index, p3.index, p4.index}, 0);
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
        
        ripsComplexStream.ensureAllFaces();
        ripsComplexStream.finalizeStream();
        
        coverageComplexStream.ensureAllFaces();
        coverageComplexStream.finalizeStream();
        
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
    //public Vector<ExplicitSimplexStream> computePersistentHomology() {
    public void computePersistentHomology() {
        ExplicitSimplexStream stream = coverageRipsComplexStream;
        stream.finalizeStream();

        AbstractPersistenceBasisAlgorithm<Simplex,IntSparseFormalSum<Simplex>> persistence = new IntAbsoluteHomology<Simplex>(ModularIntField.getInstance(2), SimplexComparator.getInstance(), 1, 2);
        //System.out.println("Got persistence object...");

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

        //List<IntSparseFormalSum<Simplex>> oneCycles = infinite_intervals.getGeneratorsAtDimension(1);
        //List<IntSparseFormalSum<Simplex>> oneCycles = annotated_intervals.getGeneratorsAtDimension(1);
        List<ObjectObjectPair<Interval<Double>, IntSparseFormalSum<Simplex>>> oneCycles = annotated_intervals.getIntervalGeneratorPairsAtDimension(1);

        /*
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

	    return cycleStreams;  */

        //System.out.println("\nManually iterate through 1-cycles: ");
        //Vector<ExplicitSimplexStream> cycleStreams = new Vector<ExplicitSimplexStream>(oneCycles.size());
        ExplicitSimplexStream strInf = new ExplicitSimplexStream();
        ExplicitSimplexStream strFin = new ExplicitSimplexStream();
        if (oneCycles != null) {
            for (ObjectObjectPair<Interval<Double>, IntSparseFormalSum<Simplex>> c : oneCycles) {
                //System.out.println(c);

                if (c.getFirst().isInfinite()) {
                    for (Simplex ss : c.getSecond().getSummands()) {
                        //System.out.println(ss);
                        strInf.addElement(ss,0);
                    }
                } else {
                    for (Simplex ss : c.getSecond().getSummands()) {
                        //System.out.println(ss);
                        strFin.addElement(ss,0);
                    }
                }
            }
        }

        //return cycleStreams;
        this.cycleStream = strFin;
        this.ripsCycleStream = strInf;
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
        //System.out.println("Got persistence object...");

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
        //System.out.println("drawComplex()");
        System.out.println("page-up / page-down to zoom");
        System.out.println("arrows rotate");
        System.out.println("1 - show/hide Rips complex");
        System.out.println("2 - show/hide coverage complex");
        System.out.println("3 - show/hide homology generators");

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
        Api.drawCoverageComplex(ripsComplexStream, coverageComplexStream, cycleStream, ripsCycleStream, domainPoints);
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
