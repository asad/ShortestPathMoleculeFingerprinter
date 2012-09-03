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

import fingerprints.model.AtomGraph;
import fingerprints.model.AtomVertex;
import java.util.*;
import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.Bond;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

/**
 * @author Syed Asad Rahman (2012) 
 * @cdk.keyword fingerprint 
 * @cdk.keyword similarity 
 * @cdk.module standard 
 * @cdk.githash
 *
 */
public class PathFinder {
    
    private static ILoggingTool logger =
            LoggingToolFactory.createLoggingTool(PathFinder.class);
    private final AtomGraph graph;
    private final AtomVertex startNode;
    private final boolean heuristic;
    private int depth;
    private boolean allPath;
    
    public PathFinder(AtomGraph graph, AtomVertex startNode) {
        this.graph = graph;
        this.startNode = startNode;
        this.heuristic = false;
    }
    
    PathFinder(AtomGraph moleculeGraph, AtomVertex startNode, boolean heuristic) {
        this.graph = moleculeGraph;
        this.startNode = startNode;
        this.heuristic = heuristic;
    }

    /*
     * recursive call to build the path
     */
    protected Stack<AtomVertex> constructPath(fingerprints.model.Path node) {
        Stack<AtomVertex> path = new Stack<AtomVertex>();
        int localDepth = 1;
        while (node.getParent() != null) {
            path.push(node.getNode());
            node = node.getParent();
            localDepth++;
        }
        path.push(startNode);
        return path;
    }

    /**
     *
     * @return
     */
    public Collection<Stack<AtomVertex>> getSinkKShorestPath() {
        
        this.depth = 0;
        this.allPath = false;
        
        List<Stack<AtomVertex>> paths = new ArrayList<Stack<AtomVertex>>();
        // list of visited nodes
        LinkedList<fingerprints.model.Path> closedList = new LinkedList<fingerprints.model.Path>();
        // list of nodes to visit (sorted)
        LinkedList<fingerprints.model.Path> openList = new LinkedList<fingerprints.model.Path>();
        openList.add(new fingerprints.model.Path(startNode));
        while (!openList.isEmpty()) {
            fingerprints.model.Path currentPath = openList.removeFirst();
            logger.debug("Path found");
            //Has covered all shortest path at the given depth
            if (!allPath) {
                // path found!
                Stack<AtomVertex> constructedPath = constructPath(currentPath);
                if (!constructedPath.isEmpty()) {
                    paths.add(constructedPath);
                }
                logger.debug("New Path found " + constructedPath);
                logger.debug("Depth " + this.depth);
            }
            
            closedList.add(currentPath);
            // add neighbors to the open list
            Iterator<AtomVertex> i = graph.getAdjacentVertices(currentPath.getNode()).keySet().iterator();
            while (i.hasNext()) {
                AtomVertex neighborNode = i.next();
                fingerprints.model.Path neighborPath = new fingerprints.model.Path(neighborNode);
                if (heuristic && !closedList.contains(neighborPath) && !openList.contains(neighborPath)) {
                    neighborPath.setPathParent(currentPath);
                    openList.add(neighborPath);
                } else if (!heuristic && !closedList.contains(neighborPath)) {
                    neighborPath.setPathParent(currentPath);
                    openList.add(neighborPath);
                }
            }
        }
        logger.debug("\nDepth: " + this.depth);
        /*
         * No paths found
         */
        return Collections.unmodifiableList(paths);
    }
    
    public static void main(String[] args) throws InvalidSmilesException {
        IAtomContainer atomContainer = new AtomContainer();
        IAtom atom1 = new Atom("Cl");
        IAtom atom2 = new Atom("C");
        IAtom atom3 = new Atom("O");
        IAtom atom4 = new Atom("N");
        IAtom atom5 = new Atom("N");
        IAtom atom6 = new Atom("S");
        IAtom atom7 = new Atom("Br");
        
        IBond bond1 = new Bond(atom1, atom2, IBond.Order.SINGLE);
        IBond bond2 = new Bond(atom2, atom3, IBond.Order.SINGLE);
        IBond bond3 = new Bond(atom2, atom4, IBond.Order.SINGLE);
        IBond bond4 = new Bond(atom2, atom5, IBond.Order.SINGLE);
        IBond bond5 = new Bond(atom3, atom6, IBond.Order.SINGLE);
        IBond bond6 = new Bond(atom1, atom3, IBond.Order.SINGLE);
        IBond bond7 = new Bond(atom5, atom7, IBond.Order.SINGLE);
        IBond bond8 = new Bond(atom4, atom7, IBond.Order.SINGLE);
        
        atomContainer.addAtom(atom1);
        atomContainer.addAtom(atom2);
        atomContainer.addAtom(atom3);
        atomContainer.addAtom(atom4);
        atomContainer.addAtom(atom5);
        atomContainer.addAtom(atom6);
        atomContainer.addAtom(atom7);
        
        
        atomContainer.addBond(bond1);
        atomContainer.addBond(bond2);
        atomContainer.addBond(bond3);
        atomContainer.addBond(bond4);
        atomContainer.addBond(bond5);
        atomContainer.addBond(bond6);
        atomContainer.addBond(bond7);
        atomContainer.addBond(bond8);
        
        AtomGraph g = new AtomGraph(atomContainer, true);
        PathFinder bfsksP = new PathFinder(g, g.getVertexLookupMap().get(atom1));
        Collection<Stack<AtomVertex>> kShortestPaths = bfsksP.getSinkKShorestPath();
        
        
        System.out.println("BFS S");
        for (Iterator<Stack<AtomVertex>> it = kShortestPaths.iterator(); it.hasNext();) {
            Stack<AtomVertex> paths = it.next();
            while (!paths.isEmpty()) {
                System.out.print(paths.pop() + ", ");
            }
            System.out.println();
        }
        System.out.println();
    }
}
