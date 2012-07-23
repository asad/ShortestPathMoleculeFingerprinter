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
import org.openscience.cdk.interfaces.IAtom;

/**
 *
 * @author Asad
 */
class Vertex implements
        Cloneable, Serializable, Comparable<AtomVertex> {

    private static final long serialVersionUID = 786786786110723560L;
    private int vertexId;

    public Vertex() {
        super();
    }

    public void setVertexId(int vertexId) {
        this.vertexId = vertexId;
    }

    public int getVertexId() {
        return vertexId;
    }

    @Override
    public int compareTo(AtomVertex t) {
        if (this.getVertexId() > t.getVertexId()) {
            return 1;
        } else if (this.getVertexId() < t.getVertexId()) {
            return -1;
        } else {
            return 0;
        }
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Vertex other = (Vertex) obj;
        if (this.vertexId != other.getVertexId()) {
            return false;
        }
        return true;
    }

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 37 * hash + this.getVertexId();
        return hash;
    }
}

public class AtomVertex extends Vertex {

    private static final long serialVersionUID = 136767675678688282L;
    private IAtom atom;

    public AtomVertex() {
    }

    @Override
    public int compareTo(AtomVertex o) {
        return this.getAtom().getSymbol().compareTo(o.getAtom().getSymbol());
    }

    public AtomVertex(IAtom name) {
        super();
        this.setAtom(name);
    }

    public AtomVertex(IAtom name, int id) {
        super();
        this.setAtom(name);
        this.setVertexId(id);
    }

    @Override
    public AtomVertex clone() {
        return new AtomVertex(this.getAtom(), getVertexId());
    }

    @Override
    public String toString() {
        return "Vertex {" + getAtom().getSymbol() + ", " + getVertexId() + "}";
    }

    public void setAtom(IAtom name) {
        this.atom = name;
    }

    public IAtom getAtom() {
        return atom;
    }
}
