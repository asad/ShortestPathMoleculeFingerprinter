/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package graph.algorithm.helper;

import graph.model.AtomVertex;
import java.util.*;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 *
 * @author Asad
 */
public class VertexCanonicalisation {

    public VertexCanonicalisation() {
    }

    /**
     * @param inputVertexSet the input Vertex Set
     * @return canonicalized VertexList
     */
    public Collection<AtomVertex> canonicalizeAtomVertex(Collection<AtomVertex> inputVertexSet) {
        List<AtomVertex> canonicalizedVertexList = new LinkedList<AtomVertex>(inputVertexSet);
        Collections.sort(canonicalizedVertexList, new VertexComparator());
        return Collections.synchronizedCollection(canonicalizedVertexList);
    }

    /**
     * @param atomSet the atomSet to set
     * @return canonicalized atoms
     */
    public Collection<IAtom> canonicalizeAtoms(Collection<IAtom> atomSet) {
        List<IAtom> canonicalizedVertexList = new LinkedList<IAtom>(atomSet);
        Collections.sort(canonicalizedVertexList, new AtomComparator());
        return Collections.synchronizedCollection(canonicalizedVertexList);
    }

    /**
     * @param container the container
     * @return canonicalized atoms
     */
    public Collection<IAtom> canonicalizeAtoms(IAtomContainer container) {

        List<IAtom> canonicalizedVertexList = new LinkedList<IAtom>();
        int i = 0;
        for (Iterator<IAtom> it = container.atoms().iterator(); it.hasNext();) {
            IAtom atom = it.next();
            canonicalizedVertexList.add(i, atom);
            i++;
        }
        Collections.sort(canonicalizedVertexList, new AtomComparator());
        return Collections.synchronizedCollection(canonicalizedVertexList);
    }
}

class VertexComparator implements Comparator<AtomVertex> {

    @Override
    public int compare(AtomVertex o1, AtomVertex o2) {
        if (!(o1 instanceof AtomVertex) || !(o2 instanceof AtomVertex)) {
            throw new ClassCastException();
        }

        if (o1.getAtom().getSymbol().compareToIgnoreCase(o2.getAtom().getSymbol()) == 0) {
            if ((o1.getAtom().getHybridization() != null
                    && o2.getAtom().getHybridization() != null)) {
                return o1.getAtom().getHybridization().compareTo(o2.getAtom().getHybridization());
            }
            return 0;
        } else {
            if ((o1.getAtom().getHybridization() != null
                    && o2.getAtom().getHybridization() != null)) {
                return 10 * o1.getAtom().getHybridization().compareTo(o2.getAtom().getHybridization());
            }
        }
        return -100;
    }
}

class AtomComparator implements Comparator<IAtom> {

    @Override
    public int compare(IAtom o1, IAtom o2) {
        if (!(o1 instanceof IAtom) || !(o2 instanceof IAtom)) {
            throw new ClassCastException();
        }
        if (o1.getSymbol().compareToIgnoreCase(o2.getSymbol()) == 0) {
            if ((o1.getHybridization() != null
                    && o2.getHybridization() != null)) {
                return o1.getHybridization().compareTo(o2.getHybridization());
            }
            return 0;
        } else {
            if ((o1.getHybridization() != null
                    && o2.getHybridization() != null)) {
                return 100 * o1.getHybridization().compareTo(o2.getHybridization());
            }
            return 10 * o1.getSymbol().compareToIgnoreCase(o2.getSymbol());
        }
    }
}
