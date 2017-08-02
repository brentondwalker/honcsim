package honcsim;

import static org.junit.Assert.*;

import java.util.ArrayList;

import org.junit.Test;

import m4rjni.Mzd;


public class DSimplexTest {

    /*
     * ======================================
     *        Utility Functions
     * ======================================
     */
    
    /**
     * Makes all the DPoints in the vertex_lists neighbors of each other
     * 
     * @param vertex_lists
     */
    @SafeVarargs
    private final void makeClique(ArrayList<DPoint>... vertex_lists) {
        int total_size = 0;
        for (ArrayList<DPoint> ll : vertex_lists) {
            total_size += ll.size();
        }
        ArrayList<DPoint> biglist = new ArrayList<DPoint>(total_size);
        for (ArrayList<DPoint> ll : vertex_lists) {
            biglist.addAll(ll);
        }
        for (int i=0; i<biglist.size(); i++) {
            for (int j=(i+1); j<biglist.size(); j++) {
                biglist.get(i).addNeighbor(biglist.get(j));
            }
        }
    }
    
    
    /**
     * Creates a DPoint with a single standard basis vector in its inventory.
     * The argument index determines which standard basis vector.
     * 
     * @param dim
     * @param index
     * @return
     */
    private final DPoint unitVertex(int dim, int index) {
        DPoint p = new DPoint(index,index,dim);
        Mzd v = new Mzd(1, dim);
        v.writeBit(0, index, 1);
        p.addInventoryVector(v);
        v.destroy();
        return p;
    }
    
    
    /*
     * ======================================
     *        Actual Testing
     * ======================================
     */

    
    @Test
    public void testDSimplex_Constructor() {
        
        // test the empty constructor
        {
            DSimplex sigma = new DSimplex();
            
            assertNotNull(sigma.vertices);
            assertEquals(0, sigma.vertices.size());
            assertNotNull(sigma.neighbors);
            assertEquals(0, sigma.neighbors.size());
            
            // The reducedBasis matrix does not get allocated until a DPoint vertex
            // is added, which sets the vsDimension
            assertNull(sigma.reducedBasis);
            assertNull(sigma.exclusiveReducedBasis);
            
            sigma.destroy();
        }
        
        // test the constructor with a list of vertices
        // in the case where the simplex is maximal
        {
            int dim = 10;
            int k = 5;
            
            // create some DPoint and put them in a List
            ArrayList<DPoint> vertex_list = new ArrayList<DPoint>(k);
            for (int i=0; i<k; i++) {
                vertex_list.add(i, unitVertex(dim,i));
            }
            
            // make all the vertices neighbors of each other
            makeClique(vertex_list);
            
            // call the constructor
            DSimplex sigma = new DSimplex(vertex_list);
            
            // verify that the key stuff has been created
            assertEquals(dim, sigma.vsDimension);
            assertNotNull(sigma.reducedBasis);
            assertEquals(k, sigma.rank);
            assertNotNull(sigma.exclusiveReducedBasis);
            assertEquals(0, sigma.exclusiveRank);
            
            // clean up
            sigma.destroy();
            for (DPoint p : vertex_list) {
                p.destroy();
            }
            vertex_list.clear();
        }
        
        // test the constructor with a list of vertices
        // in the case where the simplex is *not* maximal
        {
            int dim = 10;
            int k = 5;
            
            // create some DPoint and put them in a List
            ArrayList<DPoint> vertex_list = new ArrayList<DPoint>(k);
            for (int i=0; i<k; i++) {
                vertex_list.add(i, unitVertex(dim,i));
            }
            
            // make all the vertices neighbors of each other
            makeClique(vertex_list);
            
            // create another vertex that is a neighbor to all the points, but
            // not a part of the simplex.  That is, the simplex we create will
            // not be maximal.
            DPoint pp = unitVertex(dim, k+1);
            for (DPoint p : vertex_list) { pp.addNeighbor(p); }
            
            // call the constructor
            DSimplex sigma = new DSimplex(vertex_list);
            
            // verify that the key stuff has been created
            assertEquals(dim, sigma.vsDimension);
            assertNotNull(sigma.reducedBasis);
            assertEquals(k, sigma.rank);
            assertNotNull(sigma.exclusiveReducedBasis);
            assertEquals(1, sigma.exclusiveRank);
            
            // clean up
            sigma.destroy();
            for (DPoint p : vertex_list) {
                p.destroy();
            }
            vertex_list.clear();
            pp.destroy();
        }
    }

    @Test
    public void testDSimplex_addVertices() {
        
        // test adding vertices to an empty simplex
        // use the non-maximal simplex case
        {
            int dim = 10;
            int k = 5;

            ArrayList<DPoint> vertex_list = new ArrayList<DPoint>(k);
            for (int i=0; i<k; i++) {
                vertex_list.add(i, unitVertex(dim,i));
            }

            // make all the vertices neighbors of each other
            makeClique(vertex_list);
            
            // create another vertex that is a neighbor to all the points, but
            // not a part of the simplex.  That is, the simplex we create will
            // not be maximal.
            DPoint pp = unitVertex(dim, k+1);
            for (DPoint p : vertex_list) { pp.addNeighbor(p); }

            // create an empty simplex
            DSimplex sigma = new DSimplex();
            
            // add the vertices
            sigma.addVertices(vertex_list);

            // verify that the key stuff has been created
            assertEquals(dim, sigma.vsDimension);
            assertNotNull(sigma.reducedBasis);
            assertEquals(k, sigma.rank);
            assertNotNull(sigma.exclusiveReducedBasis);
            assertEquals(1, sigma.exclusiveRank);
            
            // clean up
            sigma.destroy();
            for (DPoint p : vertex_list) {
                p.destroy();
            }
            vertex_list.clear();
            pp.destroy();
        }
        
        // test adding vertices to a non-empty simplex
        // use the non-maximal simplex case
        {
            int dim = 20;
            int k = 5;
            
            // create some DPoint and put them in a List
            ArrayList<DPoint> vertex_list1 = new ArrayList<DPoint>(k);
            ArrayList<DPoint> vertex_list2 = new ArrayList<DPoint>(k);

            for (int i=0; i<k; i++) {
                vertex_list1.add(i, unitVertex(dim,i));
                vertex_list2.add(i, unitVertex(dim,k+i));
            }

            // make all the vertices neighbors of each other
            makeClique(vertex_list1, vertex_list2);
            
            // create another vertex that is a neighbor to all the points, but
            // not a part of the simplex.  That is, the simplex we create will
            // not be maximal.
            DPoint pp = unitVertex(dim, 2*k+1);
            for (DPoint p : vertex_list1) { pp.addNeighbor(p); }
            for (DPoint p : vertex_list2) { pp.addNeighbor(p); }
            
            // create a simplex from the first k vertices
            DSimplex sigma = new DSimplex(vertex_list1);
            
            assertEquals(k, sigma.rank);
            assertEquals(k+1, sigma.exclusiveRank);
            
            // add the second list of vertices
            sigma.addVertices(vertex_list2);

            // verify that the key stuff has been created
            assertEquals(dim, sigma.vsDimension);
            assertNotNull(sigma.reducedBasis);
            assertEquals(2*k, sigma.rank);
            assertNotNull(sigma.exclusiveReducedBasis);
            assertEquals(1, sigma.exclusiveRank);
            
            // clean up
            sigma.destroy();
            for (DPoint p : vertex_list1) {
                p.destroy();
            }
            for (DPoint p : vertex_list2) {
                p.destroy();
            }
            vertex_list1.clear();
            vertex_list2.clear();
            pp.destroy();
        }
        
    }
    
    

    
}
