/*
 * ThreePointDrainage.java
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
 * Simple experiment with the basis drainage code
 */

import java.util.*;

public class ThreePointDrainage extends DrainageExperiment {
	
	static final int vsDimension = 2;

	public ThreePointDrainage() {
		super(vsDimension, 1.0);
		
		points = new Vector<DPoint>(3);
		points.add(new DPoint(-0.9, 0.0, vsDimension));
		points.add(new DPoint(0.0, 0.0, vsDimension));
		points.add(new DPoint(0.9, 0.0, vsDimension));

		// give everything a full basis to start
		for (DPoint p : points) {
			p.addInventoryVectors(this.basis);
		}
		
		this.computeNeighborSet();
		this.buildRipsComplex();
		this.buildCoverageComplex(this.V);
	}
	
	
	/* **************************************
	 *
	 * Main Routine
	 *
	 * **************************************
	 */
	public static void main(String[] args) {
		
		ThreePointDrainage tpd = new ThreePointDrainage();
		
		tpd.drainVertex(tpd.points.get(0),0);
		tpd.drainVertex(tpd.points.get(1),0);
		tpd.drainVertex(tpd.points.get(2),0);
		
	}
	
}
