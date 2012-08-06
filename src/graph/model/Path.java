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

import graph.model.AtomVertex;
import java.io.Serializable;

/**
 *
 * @author Syed Asad Rahman (2012) @cdk.keyword fingerprint @cdk.keyword similarity @cdk.module standard @cdk.githash
 */
public class Path implements Cloneable, Serializable, Comparable<Path> {

    private static final long serialVersionUID = 136373272772822828L;

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Path other = (Path) obj;
        if (this.currentNode != other.getNode() && (this.currentNode == null || !this.currentNode.equals(other.getNode()))) {
            return false;
        }
        if (this.pathParent != other.getParent() && (this.pathParent == null || !this.pathParent.equals(other.getParent()))) {
            return false;
        }
        return true;
    }

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 43 * hash + (this.currentNode != null ? this.currentNode.hashCode() : 0);
        hash = 43 * hash + (this.pathParent != null ? this.pathParent.hashCode() : 0);
        return hash;
    }

    public AtomVertex getNode() {
        return currentNode;
    }

    @Override
    public String toString() {
        return "Path{" + "currentNode=" + currentNode + ", pathParent=" + pathParent + '}';
    }

    public Path(AtomVertex currentNode) {
        this.currentNode = currentNode;
        this.pathParent = null;
    }
    private final AtomVertex currentNode;
    private Path pathParent;

    public Path getParent() {
        return pathParent;
    }

    public void setPathParent(Path pathParent) {
        this.pathParent = pathParent;
    }

    @Override
    public int compareTo(Path t) {
        if (this.pathParent.getNode() == null || t.getParent().getNode() == null) {
            return 10 * this.currentNode.compareTo(t.getNode());
        } else if (this.currentNode.compareTo(t.getNode()) == 0
                && this.pathParent.getNode().compareTo(t.getParent().getNode()) == 0) {
            return 0;
        } else if (this.currentNode.compareTo(t.getNode()) == 0
                && this.pathParent.getNode().compareTo(t.getParent().getNode()) != 0) {
            return this.pathParent.getNode().compareTo(t.getParent().getNode());
        } else if (this.currentNode.compareTo(t.getNode()) != 0
                && this.pathParent.getNode().compareTo(t.getParent().getNode()) == 0) {
            return this.currentNode.compareTo(t.getNode());
        } else {
            return 100;
        }
    }
}
