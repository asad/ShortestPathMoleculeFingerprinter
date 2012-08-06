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
package graph.model;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author Syed Asad Rahman (2012) 
 * @cdk.keyword fingerprint 
 * @cdk.keyword similarity 
 * @cdk.module standard
 * @cdk.githash
 */
public class AdjacencyMatrix implements Serializable {

    private static final long serialVersionUID = 7688786252424151L;
    private Map<Integer, AtomVertex> lookupMap;
    private Map<AtomVertex, Integer> atomIndexMap;
    private int currentIndex;
    private Integer[][] matrix;

    public AdjacencyMatrix() {
    }

    public AdjacencyMatrix(int vertices) {
        super();
        this.matrix = new Integer[vertices][vertices];
        this.atomIndexMap = new HashMap<AtomVertex, Integer>();
        this.lookupMap = new HashMap<Integer, AtomVertex>();
        initializeMatrix();
    }

    /*
     * initialize the matrix without element as null
     */
    private void initializeMatrix() {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix.length; j++) {
                matrix[i][j] = null;
            }
        }
    }
    
    /**
     * Add an unweighted edge
     * @param source
     * @param sink
     */
    public void addEdge(AtomVertex source, AtomVertex sink) {
        addEdge(source, sink, 0);
    }
    
    /**
     * Add a weighted edge
     * @param source
     * @param sink
     * @param weight
     */
    public void addEdge(AtomVertex source, AtomVertex sink, int weight) {
        Integer i = atomIndexMap.get(source);
        Integer j = atomIndexMap.get(sink);

        if (i == null || j == null || i > matrix.length - 1
                || j > matrix.length - 1) {
            return;
        }

        this.matrix[i][j] = Integer.valueOf(weight);
        this.lookupMap.put(j, sink);
    }

    /*
     * Get the adjacent vertices of the query vertex
     *
     *
     * @param queryVertex atom vertex @return
     */
    public Map<AtomVertex, Integer> getAdjacentVertices(AtomVertex queryVertex) {
        Integer i = atomIndexMap.get(queryVertex);
        return i == null ? null : getAdjacentVertices(i);
    }

    private Map<AtomVertex, Integer> getAdjacentVertices(Integer i) {
        Map<AtomVertex, Integer> adjacentAtomVertexMap = new HashMap<AtomVertex, Integer>();
        for (int j = 0; j < matrix.length; j++) {
            if (matrix[i][j] != null) {
                adjacentAtomVertexMap.put(lookupMap.get(j), matrix[i][j]);
            }
        }
        return adjacentAtomVertexMap;
    }

    /**
     * Add a vertex to the adjacency matrix
     *
     * @param atomVertex
     * @return previous value associated with the key 
     */
    public Integer addVertex(AtomVertex atomVertex) {
        Integer value = this.atomIndexMap.put(atomVertex, currentIndex);
        currentIndex += 1;
        return value;
    }
}
