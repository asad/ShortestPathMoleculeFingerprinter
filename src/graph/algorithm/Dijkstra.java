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
package graph.algorithm;

import graph.algorithm.helper.Printer;
import graph.model.AtomContainerGraph;
import graph.model.AtomVertex;
import graph.model.ShortestPathContainer;
import java.util.*;
import java.util.Map.Entry;

/**
 *
 * @author Syed Asad Rahman (2012) 
 * @cdk.keyword fingerprint 
 * @cdk.keyword similarity 
 * @cdk.module standard
 * @cdk.githash
 */

public class Dijkstra extends Printer {

    private ShortestPathContainer shortestPaths;
    private long computingTime;
    private final AtomContainerGraph graph;
    private final AtomVertex source;

    /**
     * Initialize with graph atom container and the source vertex
     *
     * @param graph
     * @param start
     */
    public Dijkstra(AtomContainerGraph graph, AtomVertex start) {
        this.graph = graph;
        this.source = start;
        long startTime = System.currentTimeMillis();
        this.shortestPaths = findShortestPaths();
        long endTime = System.currentTimeMillis();
        this.computingTime = endTime - startTime;
    }

    private Stack<AtomVertex> reconstructShortestPath(
            ShortestPathContainer container,
            AtomVertex destination) {
        Stack<AtomVertex> pathStack = new Stack<AtomVertex>();
        AtomVertex current = destination;
        while (current != null) {
            pathStack.push(current);
            current = container.getAncestors().get(current);
        }
        return pathStack;
    }

    /**
     *
     * @param graph
     * @param source
     * @param sink
     * @param weight
     * @param path
     * @param ancestors
     */
    private void updateDistance(
            AtomVertex source,
            AtomVertex sink,
            Integer weight,
            Map<AtomVertex, Integer> path,
            Map<AtomVertex, AtomVertex> ancestors) {
        int summedWeight = path.get(source) + weight;
        if (summedWeight < path.get(sink)) {
            path.put(sink, summedWeight);
            ancestors.put(sink, source);
        }
    }

    private ShortestPathContainer findShortestPaths() {
        Map<AtomVertex, Integer> path = new HashMap<AtomVertex, Integer>();
        Map<AtomVertex, AtomVertex> ancestors = new HashMap<AtomVertex, AtomVertex>();
        Collection<AtomVertex> vertices = new LinkedList<AtomVertex>();
        initializeCost(graph, source, path, ancestors, vertices);

        /*
         * Dynamic recursive call for the main Dijkstra
         */
        while (!vertices.isEmpty()) {
            AtomVertex u = findCheapestPath(vertices, path);
            vertices.remove(u);
            for (Iterator<Entry<AtomVertex, Integer>> it =
                    graph.getAdjacentVertices(u).entrySet().iterator(); it.hasNext();) {
                Entry<AtomVertex, Integer> v = it.next();
                if (vertices.contains(v.getKey())) {
                    updateDistance(u, v.getKey(), v.getValue(), path,
                            ancestors);
                }
            }
        }
        return new ShortestPathContainer(path, ancestors, graph, source);
    }

    /**
     *
     * @param graph
     * @param start
     * @param pathMap
     * @param ancestors
     * @param vertices
     */
    private void initializeCost(
            AtomContainerGraph graph,
            AtomVertex start,
            Map<AtomVertex, Integer> pathMap,
            Map<AtomVertex, AtomVertex> ancestors,
            Collection<AtomVertex> vertices) {
        // initialize the matrix with infinity (Integer.MAX_VALUE) cost
        for (AtomVertex v : graph.getVertexSet()) {
            pathMap.put(v, Integer.valueOf(Integer.MAX_VALUE));
            ancestors.put(v, v);
            vertices.add(v.clone());
        }
        // set the distance from start to itself as zero
        pathMap.put(start, Integer.valueOf(0));
        // set the ancestors for start to null (root)
        ancestors.put(start, null);
    }

    /**
     *
     * @param q
     * @param path
     * @return
     */
    private AtomVertex findCheapestPath(
            Collection<AtomVertex> q,
            Map<AtomVertex, Integer> path) {

        AtomVertex lowest = null;
        int currentLowest = -1;

        for (AtomVertex v : q) {
            if (lowest == null) {
                lowest = v;
                currentLowest = path.get(lowest);
            } else {
                Integer u = path.get(v);
                if (u.intValue() < currentLowest) {
                    lowest = v;
                    currentLowest = u;
                }
            }
        }

        return lowest;
    }

    /**
     * @return the shortestPaths
     */
    public ShortestPathContainer getShortestPaths() {
        return shortestPaths;
    }

    /**
     * @return the Dijkstra computing Time
     */
    public long getComputingTime() {
        return computingTime;
    }

    /**
     * Return all source sink shortest paths
     *
     * @return
     */
    public List<Stack<AtomVertex>> getAllSinkShorestPath() {
        List<Stack<AtomVertex>> paths = new ArrayList<Stack<AtomVertex>>();
        for (AtomVertex destination : graph.getVertexSet()) {
            Stack<AtomVertex> reconstructedShortestPath = reconstructShortestPath(shortestPaths, destination);
            paths.add(reconstructedShortestPath);
        }
        return paths;
    }

    /**
     * Return source sink shortest paths
     *
     * @param sink
     * @return
     */
    public Stack<AtomVertex> getSinkShorestPath(AtomVertex sink) {
        Stack<AtomVertex> reconstructedShortestPath;
        reconstructedShortestPath = reconstructShortestPath(shortestPaths, sink);
        return reconstructedShortestPath;
    }
}
