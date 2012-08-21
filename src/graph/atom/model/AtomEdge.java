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

package graph.atom.model;
/**
 *
 * @author Syed Asad Rahman (2012) 
 * @cdk.keyword fingerprint 
 * @cdk.keyword similarity 
 * @cdk.module standard
 * @cdk.githash
 */
public class AtomEdge {

    private AtomVertex toVertex;
    private int weight;

    public AtomEdge() {
    }

    public AtomEdge(AtomVertex sinkVertex) {
        super();
        this.toVertex = sinkVertex;
    }

    public AtomEdge(AtomVertex to, int weight) {
        super();
        this.toVertex = to;
        this.weight = weight;
    }

    @Override
    public String toString() {
        return "Edge [sink=" + getSinkVertex() + ", weight=" + weight + "]";
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result
                + ((getSinkVertex() == null) ? 0 : getSinkVertex().hashCode());
        result = prime * result + weight;
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        AtomEdge other = (AtomEdge) obj;
        if (getSinkVertex() == null) {
            if (other.getSinkVertex() != null) {
                return false;
            }
        } else if (!toVertex.equals(other.getSinkVertex())) {
            return false;
        }
        if (weight != other.getWeight()) {
            return false;
        }
        return true;
    }

    public int getWeight() {
        return weight;
    }

    /**
     * @return the toVertex
     */
    public AtomVertex getSinkVertex() {
        return toVertex;
    }
}
