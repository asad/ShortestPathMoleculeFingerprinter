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

import graph.model.AtomContainerGraph;
import graph.model.AtomVertex;
import graph.model.Path;
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
 * This is a k-shortest path as mentioned in my PhD for Pathway Hunter Tool (PHT). This version is for unweighted graph
 * but the PHT implementation is a weighted one.
 *
 * @author Syed Asad Rahman (2012) @cdk.keyword fingerprint @cdk.keyword similarity @cdk.module standard @cdk.githash
 */
public class BFSKSP {

    private static ILoggingTool logger =
            LoggingToolFactory.createLoggingTool(BFSKSP.class);
    private final AtomContainerGraph graph;
    private final AtomVertex startNode;
    private int depth;
    private boolean firstPath;
    private boolean allPath;

    public BFSKSP(AtomContainerGraph graph, AtomVertex startNode) {
        this.graph = graph;
        this.startNode = startNode;
    }

    /*
     * recursive call to build the path
     */
    protected List<AtomVertex> constructPath(Path node) {
        LinkedList<AtomVertex> path = new LinkedList<AtomVertex>();
        int localDepth = 1;
        while (node.getParent() != null) {
            path.addFirst(node.getNode());
            node = node.getParent();
            localDepth++;
        }
        path.addFirst(startNode);

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
    public Collection<List<AtomVertex>> getKShortestPaths(AtomVertex goalNode) {

        this.depth = 0;
        this.firstPath = false;
        this.allPath = false;

        List<List<AtomVertex>> paths = new ArrayList<List<AtomVertex>>();
        // list of visited nodes
        LinkedList<Path> closedList = new LinkedList<Path>();
        // list of nodes to visit (sorted)
        LinkedList<Path> openList = new LinkedList<Path>();
        openList.add(new Path(startNode));
        while (!openList.isEmpty()) {
            Path currentPath = openList.removeFirst();
            logger.debug("\nVisiting " + currentPath.getNode().getAtom().getSymbol());

            if (currentPath.getNode() == goalNode) {
                logger.debug("Path found");
                //Has covered all shortest path at the given depth
                if (!allPath) {
                    // path found!
                    List<AtomVertex> constructedPath = constructPath(currentPath);
                    paths.add(constructedPath);
                    logger.debug("New Path found " + constructedPath);
                }
            } else {
                closedList.add(currentPath);
                // add neighbors to the open list
                Iterator<AtomVertex> i = graph.getAdjacentVertices(currentPath.getNode()).keySet().iterator();
                logger.debug("added  neighborPath:");
                while (i.hasNext()) {
                    AtomVertex neighborNode = i.next();
                    Path neighborPath = new Path(neighborNode);
                    if (!closedList.contains(neighborPath)
                            && !openList.contains(neighborPath)) {
                        neighborPath.setPathParent(currentPath);
                        logger.debug(" " + neighborPath.getNode().getAtom().getSymbol() + ", ");
                        openList.add(neighborPath);
                    }
                }
            }
        }
        logger.debug("\nDepth: " + this.depth);
        /*
         * No paths found
         */
        return Collections.unmodifiableCollection(paths);
    }

    public static void main(String[] args) throws InvalidSmilesException {
        IAtomContainer atomContainer = new AtomContainer();
        IAtom atom1 = new Atom("C");
        IAtom atom2 = new Atom("N");
        IAtom atom3 = new Atom("O");
        IAtom atom4 = new Atom("S");
        IAtom atom5 = new Atom("C");

        IBond bond1 = new Bond(atom1, atom2, IBond.Order.SINGLE);
        IBond bond2 = new Bond(atom1, atom3, IBond.Order.SINGLE);
        IBond bond3 = new Bond(atom1, atom4, IBond.Order.SINGLE);
        IBond bond4 = new Bond(atom2, atom5, IBond.Order.SINGLE);
        IBond bond5 = new Bond(atom3, atom5, IBond.Order.SINGLE);

        atomContainer.addAtom(atom1);
        atomContainer.addAtom(atom2);
        atomContainer.addAtom(atom3);
        atomContainer.addAtom(atom4);
        atomContainer.addAtom(atom5);


        atomContainer.addBond(bond1);
        atomContainer.addBond(bond2);
        atomContainer.addBond(bond3);
        atomContainer.addBond(bond4);
        atomContainer.addBond(bond5);

        AtomContainerGraph g = new AtomContainerGraph(atomContainer, true);
        BFSKSP bfsksP = new BFSKSP(g, g.getVertexLookupMap().get(atom1));
        Collection<List<AtomVertex>> kShortestPaths = bfsksP.getKShortestPaths(g.getVertexLookupMap().get(atom5));

        for (List<AtomVertex> path : kShortestPaths) {
            System.out.println(path);
        }
    }
}
