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
package graph.algorithm.helper;

import graph.model.AtomVertex;
import graph.model.ShortestPathContainer;
import java.util.Map;
import java.util.Stack;

/**
 *
 * @author Asad
 */
public class Printer {

    /**
     *
     * @param pathContainer
     * @param withCosts
     */
    public void printShortestPaths(ShortestPathContainer pathContainer, boolean withCosts) {
        for (Map.Entry<AtomVertex, AtomVertex> entry : pathContainer.getAncestors().entrySet()) {
            if (withCosts) {
                System.out.println("sink " + entry.getKey().toString() + " cost: "
                        + pathContainer.getPath().get(entry.getKey()));
            } else {
                System.out.println("\t" + entry.getKey());
            }
        }
    }

    /**
     *
     * @param pathStack
     */
    public void printShortestPath(Stack<AtomVertex> pathStack) {
        while (!pathStack.isEmpty()) {
            System.out.println("\t" + pathStack.pop());
        }
    }
}
