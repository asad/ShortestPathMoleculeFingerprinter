/* $Revision$ $Author$ $Date$
 *
 * Copyright (C) 2012   Syed Asad Rahman <asad@ebi.ac.uk>
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
package fingerprints;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import org._3pq.jgrapht.Edge;
import org._3pq.jgrapht.Graph;
import org._3pq.jgrapht.graph.SimpleGraph;
import org._3pq.jgrapht.traverse.BreadthFirstIterator;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 * @author Syed Asad Rahman (2012)
 * @cdk.keyword fingerprint 
 * @cdk.keyword similarity 
 * @cdk.module standard 
 * @cdk.githash
 *
 */
public class PathFinder {

    /**
     * Creates a undirected, unweighted molecule graph for use with jgrapht. Bond orders are not respected.
     *
     * @param molecule the specified molecule
     * @return a graph representing the molecule
     */
    static public SimpleGraph getMoleculeGraph(IAtomContainer molecule) {
        SimpleGraph graph = new SimpleGraph();
        for (IAtom atom : molecule.atoms()) {
            graph.addVertex(atom);
        }

        for (IBond bond : molecule.bonds()) {
            graph.addEdge(bond.getAtom(0), bond.getAtom(1));
            graph.addEdge(bond.getAtom(1), bond.getAtom(0));
        }
        return graph;
    }

    /**
     * Finds BFS Shortest Path between Source and Sink atoms
     *
     * @param graph
     * @param startVertex
     * @param endVertex
     * @return
     */
    public static LinkedList<Edge> findPathBetween(Graph graph, Object startVertex, Object endVertex) {
        BFSIterator iter =
                new BFSIterator(graph, startVertex);

        while (iter.hasNext()) {
            Object vertex = iter.next();
            if (vertex.equals(endVertex)) {
                return createPath(iter, endVertex);
            }
        }

        return null;
    }

    /**
     * Finds BFS Shortest Path between Source and all sink atoms
     *
     * @param graph
     * @param startVertex
     * @return
     */
    public static List<LinkedList<Edge>> findPaths(Graph graph, Object startVertex) {
        List<LinkedList<Edge>> paths = new ArrayList<LinkedList<Edge>>();
        BFSIterator iter =
                new BFSIterator(graph, startVertex);

        while (iter.hasNext()) {
            Object sinkVertex = iter.next();
            paths.add(createPath(iter, sinkVertex));
        }

        return paths;
    }

    private static LinkedList<Edge> createPath(BFSIterator iter, Object endVertex) {
        LinkedList<Edge> path = new LinkedList<Edge>();

        while (true) {
            Edge edge = iter.getSpanningTreeEdge(endVertex);
            if (edge == null) {
                break;
            }
            path.add(edge);
            endVertex = edge.oppositeVertex(endVertex);
        }
        Collections.reverse(path);
        return path;
    }

    private static class BFSIterator extends BreadthFirstIterator {

        public BFSIterator(Graph g, Object startVertex) {
            super(g, startVertex);
        }

        @Override
        protected void encounterVertex(Object vertex, Edge edge) {
            super.encounterVertex(vertex, edge);
            putSeenData(vertex, edge);
        }

        public Edge getSpanningTreeEdge(Object vertex) {
            return (Edge) getSeenData(vertex);
        }
    }
}
