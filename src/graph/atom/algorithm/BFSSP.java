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
package graph.atom.algorithm;

import graph.atom.model.AtomGraph;
import graph.atom.model.AtomVertex;
import graph.atom.model.Path;
import java.util.*;
import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.Bond;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

/**
 * This is a k-shortest path as mentioned in my PhD for Pathway Hunter Tool (PHT). This version is for unweighted graph
 * but the PHT implementation is a weighted one.
 *
 * @author Syed Asad Rahman (2012) @cdk.keyword fingerprint @cdk.keyword similarity @cdk.module standard @cdk.githash
 */
public class BFSSP {

    private static ILoggingTool logger =
            LoggingToolFactory.createLoggingTool(BFSSP.class);
    private final AtomGraph graph;
    private final AtomVertex startNode;
    private int depth;
    private boolean firstPath;
    private boolean allPath;

    public BFSSP(AtomGraph graph, AtomVertex startNode) {
        this.graph = graph;
        this.startNode = startNode;
    }

    /*
     * recursive call to build the path
     */
    protected Stack<AtomVertex> constructPath(Path node) {
        Stack<AtomVertex> path = new Stack<AtomVertex>();
        int localDepth = 1;
        while (node.getParent() != null) {
            path.push(node.getNode());
            node = node.getParent();
            localDepth++;
        }
        path.push(startNode);
        /*
         * Do not report path longer than first shortest path
         */
        if (firstPath && localDepth != depth) {
            return new Stack<AtomVertex>();
        }
        if (!firstPath) {
            this.firstPath = true;
            this.depth = localDepth;
        }
        if (firstPath && localDepth > depth) {
            this.allPath = true;
        }
        return path;
    }

    /**
     *
     * @param goalNode
     * @return
     */
    public Collection<Stack<AtomVertex>> getSinkKShorestPath(AtomVertex goalNode) {

        this.depth = 0;
        this.firstPath = false;
        this.allPath = false;

        List<Stack<AtomVertex>> paths = new ArrayList<Stack<AtomVertex>>();
        // list of visited nodes
        LinkedList<Path> closedList = new LinkedList<Path>();
        // list of nodes to visit (sorted)
        LinkedList<Path> openList = new LinkedList<Path>();
        openList.add(new Path(startNode));
        while (!openList.isEmpty()) {
            Path currentPath = openList.removeFirst();
            //System.out.println("\nVisiting " + currentPath.getNode().getAtom().getSymbol());

            if (currentPath.getNode() == goalNode) {
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
            } else {
                closedList.add(currentPath);
                //                System.out.println("\nCurrent: " + currentPath.getNode());
                // add neighbors to the open list
                Set<AtomVertex> container = graph.getAdjacentVertices(currentPath.getNode()).keySet();
                Iterator<AtomVertex> i = container.iterator();
//                System.out.println("NeighborPath:");
                while (i.hasNext()) {
                    AtomVertex neighborNode = i.next();
                    Path neighborPath = new Path(neighborNode);
                    if (!closedList.contains(neighborPath) && !openList.contains(neighborPath)) {
                        neighborPath.setPathParent(currentPath);
//                        System.out.print(" " + neighborPath.getNode() + ":");
                        openList.add(neighborPath);
                    }
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
        BFSSP bfsksP = new BFSSP(g, g.getVertexLookupMap().get(atom1));
        Collection<Stack<AtomVertex>> kShortestPaths = bfsksP.getSinkKShorestPath(g.getVertexLookupMap().get(atom7));

//        GraphOutputWriter.GraphMLWriter(pg.getGraph(), new File("example.gml"));

        System.out.println("BFS S");
        for (Iterator<Stack<AtomVertex>> it = kShortestPaths.iterator(); it.hasNext();) {
            Stack<AtomVertex> paths = it.next();
            while (!paths.isEmpty()) {
                System.out.print(paths.pop() + ", ");
            }
            System.out.println();
        }
        Collection<List<IAtom>> shortestPath = PathTools.getPathsOfLengthUpto(atomContainer, atom1, 4);
        for (List<IAtom> l : shortestPath) {
            System.out.println("CDK shortestPath: ");
            for (IAtom a : l) {
                System.out.print(a.getSymbol() + ", ");
            }
            System.out.println();
        }
    }
}
