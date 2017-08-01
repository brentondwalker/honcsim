package honcsim;

import m4rjni.Mzd;
import static org.junit.Assert.*;
import org.junit.Test;


public class DPointTest {

    @Test
    public void testDPoint_Constructor() {
        {
            DPoint p = new DPoint(2.1, 3.4, 100);

            assertNotNull(p.nbrs);
            assertEquals(0, p.nbrs.size());
            assertNotNull(p.nbrsUp);
            assertEquals(0, p.nbrsUp.size());
            assertNotNull(p.ripsCofaces);
            assertEquals(0, p.ripsCofaces.size());
            assertNotNull(p.M);
            assertEquals(0, p.M.size());
            assertNotNull(p.reducedBasis);
            assert(p.reducedBasis.isZero());
            assertEquals(100, p.vsDimension);
            p.destroy();
        }
    }
    
    @Test
    public void testDPoint_addInventoryVector() {
        // test a single unit vector
        {
            int dim = 10;
            DPoint p = new DPoint(2.1, 3.4, dim);
            Mzd v = new Mzd(1, dim);
            v.writeBit(0, 0, 1);
            p.addInventoryVector(v);
            Mzd M1 = new Mzd(dim,dim);
            M1.writeBit(0, 0, 1);
            assertEquals(1, p.rank);
            assertTrue(p.reducedBasis.equals(M1));
            v.destroy();
            M1.destroy();
            p.destroy();
        }
        
        // test a bunch of unit vectors added one at a time
        {
            int dim = 10;
            int k = 5;
            DPoint p = new DPoint(2.1, 3.4, dim);
            for (int i=0; i<k; i++) {
                Mzd v = new Mzd(1, dim);
                v.writeBit(0, i, 1);
                p.addInventoryVector(v);
                v.destroy();
            }
            Mzd Mk = new Mzd(dim,dim);
            for (int i=0; i<k; i++) {
                Mk.writeBit(i, i, 1);
            }
            assertEquals(k, p.rank);
            assertTrue(p.reducedBasis.equals(Mk));
            Mk.destroy();
            p.destroy();
        }
    }
    
    @Test
    public void testDPoint_addInventoryVectors() {
        // test a bunch of unit vectors added together
        {
            int dim = 10;
            int k = 5;
            DPoint p = new DPoint(2.1, 3.4, dim);
            Mzd[] vv = new Mzd[k];
            for (int i=0; i<k; i++) {
                vv[i] = new Mzd(1, dim);
                vv[i].writeBit(0, i, 1);
            }
            p.addInventoryVectors(vv);
            for (int i=0; i<k; i++)
                vv[i].destroy();
            
            Mzd Mk = new Mzd(dim,dim);
            for (int i=0; i<k; i++) {
                Mk.writeBit(i, i, 1);
            }
            assertEquals(k, p.rank);
            assertTrue(p.reducedBasis.equals(Mk));
            Mk.destroy();
            p.destroy();
        }
    }
    
    @Test
    public void testDPoint_addRandomInventoryVector() {
        // test adding a single random vector
        {
            int dim = 10;
            DPoint p = new DPoint(2.1, 3.4, dim);
            p.addRandomInventoryVector();
            assertEquals(1,p.rank);
            p.destroy();
        }
    }
    
    @Test
    public void testDPoint_addRandomInventoryVector_statistics() {
        // test the statistics of adding many random vectors
        // if the rank is (dim-1) then adding another random vector
        // should increase the rank with probability 1/2
        {
            int num_trials = 10000;
            int dim = 10;
            int sum = 0;
            for (int i=0; i<num_trials; i++) {
                DPoint p = new DPoint(2.1, 3.4, dim);
                while (p.rank < (dim-1)) {
                    p.addRandomInventoryVector();
                }
                p.addRandomInventoryVector();
                sum += dim - p.rank;
                p.destroy();
            }
            System.out.println("sum: "+sum);
            double mean = ((double)sum)/((double)num_trials);
            assertTrue(Math.abs(mean - 0.5) < 0.01);
        }
        
        // test the statistics of adding a random vector when the rank is
        // about 1/4 way to completion (25/100).  The probability that new
        // vectors do not increase the rank at this point is very close to
        // zero.  Even in 1000 trials, we do not expect to see it happen.
        {
            int num_trials = 1000;
            int dim = 100;
            int partial_rank = 25;
            int num_additional = 5;
            for (int i=0; i<num_trials; i++) {
                DPoint p = new DPoint(2.1, 3.4, dim);
                while (p.rank < partial_rank) {
                    p.addRandomInventoryVector();
                }
                assertEquals(partial_rank, p.rank);
                for (int j=0; j<num_additional; j++) {
                    p.addRandomInventoryVector();
                }
                assertEquals((partial_rank+num_additional), p.rank);
                p.destroy();
            }
        }
    }
    
    @Test
    public void testDPoint_clearInventory() {
        {
            int dim = 100;
            int partial_rank = 25;
            DPoint p = new DPoint(2.1, 3.4, dim);
            while (p.rank < partial_rank) {
                p.addRandomInventoryVector();
            }
            assertEquals(partial_rank, p.rank);
            assertTrue(p.M.size() >= partial_rank);
            p.clearInventory();
            assertEquals(0, p.rank);
            assertNotNull(p.M);
            assertEquals(0, p.M.size());
            assertNotNull(p.reducedBasis);
            assertTrue(p.reducedBasis.isZero());
            p.destroy();
        }
    }
    
    @Test
    public void testDPoint_addNeighbor() {
        {
            int dim = 10;
            int num_neighbors = 10;
            DPoint[] neighbors = new DPoint[num_neighbors];
            DPoint p = new DPoint(2.1, 3.4, dim);
            for (int i=0; i<num_neighbors; i++) {
                neighbors[i] = new DPoint(i, i, dim);
                p.addNeighbor(neighbors[i]);
            }
            assertEquals(num_neighbors, p.nbrs.size());
            for (int i=0; i<num_neighbors; i++) {
                assertTrue(p.nbrs.contains(neighbors[i]));
                assertTrue(neighbors[i].nbrs.contains(p));
                // for each neighbor pair, exactly one of them is in the other's NbrsUp list.
                assertTrue(p.nbrsUp.contains(neighbors[i]) ^  neighbors[i].nbrsUp.contains(p));
            }
            for (int i=0; i<num_neighbors; i++) {
                neighbors[i].destroy();
            }
            p.destroy();
        }
    }
    
    @Test
    public void testDPoint_removeNeighbor() {
        {
            int dim = 10;
            int num_neighbors = 10;
            DPoint[] neighbors = new DPoint[num_neighbors];
            DPoint p = new DPoint(2.1, 3.4, dim);
            for (int i=0; i<num_neighbors; i++) {
                neighbors[i] = new DPoint(i, i, dim);
                p.addNeighbor(neighbors[i]);
            }
            assertEquals(num_neighbors, p.nbrs.size());
            
            // remove one neighbor and make sure it's gone
            p.removeNeighbor(neighbors[0]);
            assertFalse(p.nbrs.contains(neighbors[0]));
            assertFalse(p.nbrsUp.contains(neighbors[0]) || neighbors[0].nbrsUp.contains(p));
            
            // make sure no other neighbors were messed with
            for (int i=1; i<num_neighbors; i++) {
                assertTrue(p.nbrs.contains(neighbors[i]));
                assertTrue(neighbors[i].nbrs.contains(p));
                // for each neighbor pair, exactly one of them is in the other's NbrsUp list.
                assertTrue(p.nbrsUp.contains(neighbors[i]) ^  neighbors[i].nbrsUp.contains(p));
            }
            
            // clean up
            for (int i=0; i<num_neighbors; i++) {
                neighbors[i].destroy();
            }
            p.destroy();
        }
    }
    
    @Test
    public void testDPoint_clearNeighbors() {
        {
            int dim = 10;
            int num_neighbors = 10;
            DPoint[] neighbors = new DPoint[num_neighbors];
            DPoint p = new DPoint(2.1, 3.4, dim);
            for (int i=0; i<num_neighbors; i++) {
                neighbors[i] = new DPoint(i, i, dim);
                p.addNeighbor(neighbors[i]);
            }
            assertEquals(num_neighbors, p.nbrs.size());
            
            // make everyone a neighbor of everyone else
            for (int i=0; i<num_neighbors; i++) {
                for (int j=(i+1); j<num_neighbors; j++) {
                    neighbors[i].addNeighbor(neighbors[j]);
                }
            }
            
            // clear the neighbors from p
            p.clearNeighbors();

            // verify that p is now alone, and no other DPoint has p as a neighbor
            assertEquals(0, p.nbrs.size());
            assertEquals(0, p.nbrsUp.size());
            for (int i=1; i<num_neighbors; i++) {
                assertFalse(neighbors[0].nbrsUp.contains(p));
            }
            
            // verify that all the other DPoints are still neighbors of each other
            for (int i=0; i<num_neighbors; i++) {
                for (int j=(i+1); j<num_neighbors; j++) {
                    assertTrue(neighbors[i].nbrs.contains(neighbors[j]));
                    assertTrue(neighbors[j].nbrs.contains(neighbors[i]));
                    assertTrue(neighbors[i].nbrsUp.contains(neighbors[j]) ^ neighbors[j].nbrsUp.contains(neighbors[i]));
                }
            }
            
            // clean up
            for (int i=0; i<num_neighbors; i++) {
                neighbors[i].destroy();
            }
            p.destroy();
        }
    }
    
    @Test
    public void testDPoint_destroy() {
        {
            int dim = 10;
            int num_neighbors = 10;
            DPoint[] neighbors = new DPoint[num_neighbors];
            DPoint p = new DPoint(2.1, 3.4, dim);
            for (int i=0; i<num_neighbors; i++) {
                neighbors[i] = new DPoint(i, i, dim);
                p.addNeighbor(neighbors[i]);
            }
            assertEquals(num_neighbors, p.nbrs.size());
            
            // make everyone a neighbor of everyone else
            for (int i=0; i<num_neighbors; i++) {
                for (int j=(i+1); j<num_neighbors; j++) {
                    neighbors[i].addNeighbor(neighbors[j]);
                }
            }
            
            // clear the neighbors from p
            p.destroy();

            // verify that p is now alone, and no other DPoint has p as a neighbor
            assertEquals(0, p.nbrs.size());
            assertEquals(0, p.nbrsUp.size());
            for (int i=1; i<num_neighbors; i++) {
                assertFalse(neighbors[0].nbrsUp.contains(p));
            }
            
            // verify that all the other DPoints are still neighbors of each other
            for (int i=0; i<num_neighbors; i++) {
                for (int j=(i+1); j<num_neighbors; j++) {
                    assertTrue(neighbors[i].nbrs.contains(neighbors[j]));
                    assertTrue(neighbors[j].nbrs.contains(neighbors[i]));
                    assertTrue(neighbors[i].nbrsUp.contains(neighbors[j]) ^ neighbors[j].nbrsUp.contains(neighbors[i]));
                }
            }
            
            // check that all the Mzd objects are gone
            assertNotNull(p.M);
            assertEquals(0, p.M.size());
            assertNull(p.reducedBasis);
            
            // clean up
            for (int i=0; i<num_neighbors; i++) {
                neighbors[i].destroy();
            }
        }
    }
    
}
