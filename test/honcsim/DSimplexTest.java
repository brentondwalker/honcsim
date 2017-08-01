package honcsim;

import static org.junit.Assert.*;

import java.util.ArrayList;

import org.junit.Test;

import m4rjni.Mzd;


public class DSimplexTest {

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
        {
            int dim = 100;
            int k = 5;
            
            // create some DPoint and put them in a List
            ArrayList<DPoint> vertex_list = new ArrayList<DPoint>(k);
            for (int i=0; i<k; i++) {
                vertex_list.add(i, new DPoint(i,i,dim));
                Mzd v = new Mzd(1, dim);
                v.writeBit(0, i, 1);
                vertex_list.get(i).addInventoryVector(v);
                v.destroy();
            }
            
            // call the constructor
            DSimplex sigma = new DSimplex(vertex_list);
            
            // verify that the key stuff has been created
            assertNotNull(sigma.reducedBasis);
            assertEquals(k, sigma.rank);
            assertNotNull(sigma.exclusiveReducedBasis);
            
        }
    }
    
}
