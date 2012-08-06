/* $Revision$ $Author$ $Date$
 *
 * Copyright (C) 2011-2012       Syed Asad Rahman <asad@ebi.ac.uk>
 *           
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute itAtom and/or
 * modify itAtom under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that itAtom will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package graph.model;

import graph.algorithm.helper.VertexCanonicalisation;
import java.io.Serializable;
import java.util.*;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 *
 * @author Syed Asad Rahman (2012) 
 * @cdk.keyword fingerprint 
 * @cdk.keyword similarity 
 * @cdk.module standard
 * @cdk.githash
 */
public class AtomContainerGraph extends ExampleGraphContainers implements Serializable {

    private static final long serialVersionUID = 88898677862424221L;
    private Set<AtomVertex> atomVertexSet;
    private AdjacencyMatrix adjacencyMatrix;
    private Map<IAtom, AtomVertex> vertexLookupMap;

    public AtomContainerGraph() {
        this.atomVertexSet = new HashSet<AtomVertex>();
        this.vertexLookupMap = new HashMap<IAtom, AtomVertex>();
    }

    /**
     * AtomContainer to be converted as graph
     *
     * @param container
     * @param weighted
     */
    public AtomContainerGraph(IAtomContainer container, boolean weighted) {
        this();
        adjacencyMatrix = new AdjacencyMatrix(container.getAtomCount());
        setAtomContainer(container, weighted);
    }

    private void setAtomContainer(IAtomContainer container, boolean weighted) {
        int i = 1;
        /*
         * Canonicalized vertexs as SP may return only one path out of k-sp Hence this path has to be same for similar
         * molecules
         */
        Collection<IAtom> canonicalizeAtoms = new VertexCanonicalisation().canonicalizeAtoms(container);
        for (Iterator<IAtom> it = canonicalizeAtoms.iterator(); it.hasNext();) {
            IAtom atom = it.next();
            AtomVertex atomVertex = new AtomVertex(atom, i);
            i += 1;
            addVertex(atomVertex, true);
        }


        for (Iterator<IBond> itBond = container.bonds().iterator(); itBond.hasNext();) {
            IBond bond = itBond.next();
            AtomVertex v1 = getVertexLookupMap().get(bond.getAtom(0));
            AtomVertex v2 = getVertexLookupMap().get(bond.getAtom(1));
            if (weighted) {
                addEdge(v1, new Edge(v2, (bond.getOrder().ordinal() + 1)));
                addEdge(v2, new Edge(v1, (bond.getOrder().ordinal() + 1)));
            } else {
                addEdge(v1, new Edge(v2, 1));
                addEdge(v2, new Edge(v1, 1));
            }

        }
    }

    /**
     *
     * @param atomVertex
     * @param intoMatrix
     */
    private void addVertex(AtomVertex atomVertex, boolean intoMatrix) {
        atomVertexSet.add(atomVertex);
        getVertexLookupMap().put(atomVertex.getAtom(), atomVertex);
        if (intoMatrix) {
            if (adjacencyMatrix == null) {
                throw new NullPointerException(
                        "Adjacency Matrix is not initialized!");
            }
            adjacencyMatrix.addVertex(atomVertex);
        }
    }

    /**
     *
     * @param atomVertex
     * @return
     */
    public Map<AtomVertex, Integer> getAdjacentVertices(AtomVertex atomVertex) {
        return adjacencyMatrix.getAdjacentVertices(atomVertex);
    }

    /**
     * Add edge to the vertex
     *
     * @param v
     * @param l
     */
    private void addEdge(AtomVertex v, Edge... l) {
        for (Edge e : l) {
            adjacencyMatrix.addEdge(v, e.getSinkVertex(), e.getWeight());
        }
    }

    public Map<IAtom, AtomVertex> getVertexLookupMap() {
        return Collections.synchronizedMap(this.vertexLookupMap);
    }

    /**
     *
     * @return
     */
    public Collection<AtomVertex> getVertexSet() {
        return Collections.unmodifiableSet(atomVertexSet);
    }
}
