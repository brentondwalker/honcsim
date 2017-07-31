package honcsim;

import honcsim.DPoint;
import m4rjni.Mzd;

import static org.junit.Assert.*;

import java.util.Vector;

import static org.hamcrest.CoreMatchers.*;
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
    }

}
