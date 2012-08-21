/* $Revision$ $Author$ $Date$
 *
 * Copyright (C) 2011-2012       Syed Asad Rahman <asad@ebi.ac.uk>
 *           
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package graph.atom.model;

import java.util.Collections;
import java.util.Map;

/**
 *
 * @author Syed Asad Rahman (2012) 
 * @cdk.keyword fingerprint 
 * @cdk.keyword similarity 
 * @cdk.module standard
 * @cdk.githash
 */
public class ShortestPathContainer {

    private final AtomGraph atomContainerGraph;
    private final AtomVertex start;
    private final Map<AtomVertex, Integer> path;
    private final Map<AtomVertex, AtomVertex> ancestors;

    /**
     *
     * @param path
     * @param ancestors
     * @param atomContainerGraph
     * @param start
     */
    public ShortestPathContainer(
            Map<AtomVertex, Integer> path,
            Map<AtomVertex, AtomVertex> ancestors,
            AtomGraph atomContainerGraph,
            AtomVertex start) {
        super();
        this.path = path;
        this.ancestors = ancestors;
        this.atomContainerGraph = atomContainerGraph;
        this.start = start;
    }

    /**
     * @return the atomContainerGraph
     */
    public AtomGraph getAtomContainerGraph() {
        return atomContainerGraph;
    }

    /**
     * @return the start
     */
    public AtomVertex getStart() {
        return start;
    }

    /**
     * @return the path
     */
    public Map<AtomVertex, Integer> getPath() {
        return Collections.unmodifiableMap(path);
    }

    /**
     * @return the ancestors
     */
    public Map<AtomVertex, AtomVertex> getAncestors() {
        return Collections.unmodifiableMap(ancestors);
    }
}
