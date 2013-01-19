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
package fingerprints.model;

import java.io.Serializable;
import java.util.*;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;

/**
 *
 * @author Syed Asad Rahman (2012) 
 * @cdk.keyword fingerprint 
 * @cdk.keyword similarity 
 * @cdk.module standard 
 * @cdk.githash
 */
public class AtomGraph extends ExampleGraphContainers implements Serializable {

    private static final long serialVersionUID = 88898677862424221L;
    private final Set<AtomVertex> atomVertexSet;
    private AdjacencyAtomMatrix adjacencyMatrix;
    private final Map<IAtom, AtomVertex> vertexLookupMap;
    private final Collection<IAtom> skipVertices;

    public AtomGraph() {
        this.atomVertexSet = new HashSet<AtomVertex>();
        this.vertexLookupMap = new HashMap<IAtom, AtomVertex>();
        this.skipVertices = new HashSet<IAtom>();
    }

    /**
     * AtomContainer to be converted as graph
     *
     * @param container
     * @param weighted
     */
    public AtomGraph(IAtomContainer container, boolean weighted) {
        this();
        adjacencyMatrix = new AdjacencyAtomMatrix(container.getAtomCount());
        setAtomContainer(container, weighted);
    }

    /**
     * AtomContainer to be converted as graph
     *
     * @param container
     * @param skipVertices
     * @param weighted
     */
    public AtomGraph(IAtomContainer container, IAtomContainer skipVertices, boolean weighted) {
        this();
        adjacencyMatrix = new AdjacencyAtomMatrix(container.getAtomCount());
        for (IAtom atom : skipVertices.atoms()) {
            this.skipVertices.add(atom);
        }
        setAtomContainer(container, weighted);
    }

    /**
     * AtomContainer to be converted as graph
     *
     * @param container
     * @param skipVertices
     * @param weighted
     */
    public AtomGraph(IAtomContainer container, IAtomContainerSet skipVertices, boolean weighted) {
        this();
        adjacencyMatrix = new AdjacencyAtomMatrix(container.getAtomCount());
        for (IAtomContainer ac : skipVertices.atomContainers()) {
            for (IAtom atom : ac.atoms()) {
                this.skipVertices.add(atom);
            }
        }
        setAtomContainer(container, weighted);
    }

    private void setAtomContainer(IAtomContainer container, boolean weighted) {
        int i = 1;
        /*
         * Canonicalized vertexs as SP may return only one path out of k-sp Hence this path has to be same for similar
         * molecules
         */
        Collection<IAtom> canonicalizeAtoms = new SimpleAtomCanonicalizer().canonicalizeAtoms(container);
        for (Iterator<IAtom> it = canonicalizeAtoms.iterator(); it.hasNext();) {
            IAtom atom = it.next();
            if (skipVertices.contains(atom)) {
                continue;
            }
            AtomVertex atomVertex = new AtomVertex(atom, i);
            addVertex(atomVertex, true);
            i += 1;
        }


        for (Iterator<IBond> itBond = container.bonds().iterator(); itBond.hasNext();) {
            IBond bond = itBond.next();
            if (skipVertices.contains(bond.getAtom(0)) || skipVertices.contains(bond.getAtom(1))) {
                continue;
            }
            AtomVertex v1 = getVertexLookupMap().get(bond.getAtom(0));
            AtomVertex v2 = getVertexLookupMap().get(bond.getAtom(1));

            if (weighted) {
                addEdge(v1, new AtomEdge(v2, (bond.getOrder().ordinal() + 1)));
                addEdge(v2, new AtomEdge(v1, (bond.getOrder().ordinal() + 1)));
            } else {
                addEdge(v1, new AtomEdge(v2, 1));
                addEdge(v2, new AtomEdge(v1, 1));
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
        //System.out.println("S: "+ atomVertex + ", N: " + adjacencyMatrix.getAdjacentVertices(atomVertex).keySet());
        return adjacencyMatrix.getAdjacentVertices(atomVertex);
    }

    /**
     * Add edge to the vertex
     *
     * @param v
     * @param l
     */
    private void addEdge(AtomVertex v, AtomEdge... l) {
        for (AtomEdge e : l) {
            adjacencyMatrix.addEdge(v, e.getSinkVertex(), e.getWeight());
        }
    }

    public Map<IAtom, AtomVertex> getVertexLookupMap() {
        return this.vertexLookupMap;
    }

    /**
     *
     * @return
     */
    public Collection<AtomVertex> getVertexSet() {
        return atomVertexSet;
    }
}
