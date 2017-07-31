package honcsim;

import honcsim.DPoint;
import static org.junit.Assert.*;
import org.junit.Test;


public class DPointTest {

    @Test
    public void testDPointConstructor() {
        DPoint p = new DPoint(2.1, 3.4, 100);
        {
            //assertNotEquals(null, p.nbrs);
            assertEquals(0, p.nbrs.size());
        }
        //String s = System.getProperty("java.library.path");
        //System.err.println("java.library.path = "+s);
        //s = System.getProperty("user.dir");
        //System.err.println("user.dir = "+s);
        //{
            //Mzd m = new Mzd(10,10);
            //assertNotEquals(null, m);
            //assertEquals(10, m.getNcols());
            //assertEquals(10, m.getNrows());
            //assert(m.isZero());
            //m.destroy();
        //}
    }
    
}
