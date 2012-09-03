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
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.NoSuchAtomException;
import org.openscience.cdk.graph.SpanningTree;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.tools.periodictable.PeriodicTable;

/**
 *
 * @author Syed Asad Rahman (2012) 
 * @cdk.keyword fingerprint 
 * @cdk.keyword similarity 
 * @cdk.module standard 
 * @cdk.githash
 *
 */
final class ShortestPathWalker {

    private static final long serialVersionUID = 0x3b728f46;
    private final IAtomContainer atomContainer;
    private final Set<String> cleanPath;
    private final List<String> pseudoAtoms;
    private int pseduoAtomCounter;
    private final Set<String> allPaths;
    private final IRingSet basicRings;
    private final Map<IAtom, Map<IAtom, IBond>> cache;

    /**
     *
     * @param atomContainer
     * @throws CloneNotSupportedException
     * @throws CDKException
     */
    public ShortestPathWalker(IAtomContainer atomContainer) throws NoSuchAtomException {
        SpanningTree spanningTree = new SpanningTree(atomContainer);
        this.basicRings = spanningTree.getBasicRings();
        this.cleanPath = new HashSet<String>();
        this.atomContainer = atomContainer;
        this.pseudoAtoms = new ArrayList<String>();
        this.pseduoAtomCounter = 0;
        this.allPaths = new HashSet<String>();
        this.cache = new HashMap<IAtom, Map<IAtom, IBond>>();
        findPaths();
    }

    /**
     * @return the cleanPath
     */
    public Set<String> getPaths() {
        return Collections.unmodifiableSet(cleanPath);
    }

    /**
     * @return the cleanPath
     */
    public int getPathCount() {
        return cleanPath.size();
    }

    private void findPaths() {
        pseudoAtoms.clear();
        getAtomPaths();
        getRingAtomPaths();
        for (String s1 : allPaths) {
            if (s1.equals("")) {
                continue;
            }
            if (!cleanPath.contains(s1)) {
//                System.out.println(s1);
                cleanPath.add(s1);
            }
        }
    }

    private void setAtom(IAtom atomCurrent, StringBuilder sb) {
        if (atomCurrent instanceof IPseudoAtom) {
            if (!pseudoAtoms.contains(atomCurrent.getSymbol())) {
                pseudoAtoms.add(pseduoAtomCounter, atomCurrent.getSymbol());
                pseduoAtomCounter += 1;
            }
            sb.append((char) (PeriodicTable.getElementCount()
                    + pseudoAtoms.indexOf(atomCurrent.getSymbol()) + 1));
        } else {
            Integer atnum = PeriodicTable.getAtomicNumber(atomCurrent.getSymbol());
            if (atnum != null) {
                sb.append(toAtomPattern(atomCurrent));
            } else {
                sb.append((char) PeriodicTable.getElementCount() + 1);
            }
        }
    }

    private String toAtomPattern(IAtom atom) {
        return atom.getSymbol();
    }

    /**
     * Gets the bondSymbol attribute of the HashedFingerprinter class
     *
     * @param bond Description of the Parameter
     * @return The bondSymbol value
     */
    private char getBondSymbol(IBond bond) {
        if (isSP2Bond(bond)) {
            return '@';
        } else {
            switch (bond.getOrder()) {
                case SINGLE:
                    return '1';
                case DOUBLE:
                    return '2';
                case TRIPLE:
                    return '3';
                case QUADRUPLE:
                    return '4';
                default:
                    return '5';
            }
        }
    }

    /**
     * Returns true if the bond binds two atoms, and both atoms are SP2 in a ring system.
     */
    private boolean isSP2Bond(IBond bond) {
        if (bond.getAtomCount() == 2
                && bond.getAtom(0).getHybridization() == IAtomType.Hybridization.SP2
                && bond.getAtom(1).getHybridization() == IAtomType.Hybridization.SP2
                && isRingBond(bond)) {
            return true;
        }
        return false;
    }

    /**
     * Returns true if the bond binds two atoms, and both atoms are SP2 in a ring system.
     */
    private boolean isRingBond(IBond bond) {
        if (bond.getAtomCount() == 2
                && basicRings.contains(bond.getAtom(0))
                && basicRings.contains(bond.getAtom(1))) {
            return true;
        }
        return false;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (String path : cleanPath) {
            sb.append(path).append("->");
        }
        return sb.toString();
    }
    /*
     * This module generates shortest path between atoms using BFS
     */

    private void getAtomPaths() {

        /*
         * Canonicalisation of atoms for reporting unique pathsToDestination with consistency
         */
        AtomGraph moleculeGraph = new AtomGraph(atomContainer, basicRings, false);

        for (IAtom sourceAtom : atomContainer.atoms()) {
            if (basicRings.contains(sourceAtom)) {
                continue;
            }
            StringBuilder sb = new StringBuilder();
            setAtom(sourceAtom, sb);
            allPaths.add(sb.toString().trim());
            PathFinder pathFinder = new PathFinder(moleculeGraph, moleculeGraph.getVertexLookupMap().get(sourceAtom), true);
            Collection<Stack<AtomVertex>> shortestPaths = pathFinder.getSinkKShorestPath();

            for (Stack<AtomVertex> shortestPath : shortestPaths) {
                sb = new StringBuilder();
                if (shortestPath == null || shortestPath.isEmpty() || shortestPath.size() > 8) {
                    continue;
                }
                IAtom atomCurrent = shortestPath.pop().getAtom();
                while (!shortestPath.empty()) {
                    final IAtom atomNext = shortestPath.pop().getAtom();
                    Map<IAtom, IBond> m = cache.get(atomCurrent);
                    final IBond[] b = {m != null ? m.get(atomNext) : null};
                    if (b[0] == null) {
                        b[0] = atomContainer.getBond(atomCurrent, atomNext);
                        cache.put(atomCurrent,
                                new HashMap<IAtom, IBond>() {

                                    {
                                        put(atomNext, b[0]);
                                    }
                                    private static final long serialVersionUID = 0xb3a7a32449fL;
                                });
                    }

                    setAtom(atomCurrent, sb);
                    sb.append(getBondSymbol(b[0]));
                    setAtom(atomNext, sb);
                    atomCurrent = atomNext;

                }
                allPaths.add(sb.toString().trim());
            }
        }
    }

    /*
     * This module generates shortest path between atoms using BFS
     */
    private void getRingAtomPaths() {

        /*
         * Canonicalisation of atoms for reporting unique pathsToDestination with consistency
         */

        for (IAtomContainer ac : basicRings.atomContainers()) {
            AtomGraph moleculeGraph = new AtomGraph(ac, false);
            for (IAtom sourceAtom : ac.atoms()) {
                StringBuilder sb = new StringBuilder();
                setAtom(sourceAtom, sb);
                allPaths.add(sb.toString().trim());
                PathFinder pathFinder = new PathFinder(moleculeGraph, moleculeGraph.getVertexLookupMap().get(sourceAtom), true);
                Collection<Stack<AtomVertex>> shortestPaths = pathFinder.getSinkKShorestPath();
                for (Stack<AtomVertex> shortestPath : shortestPaths) {
                    sb = new StringBuilder();
                    if (shortestPath == null || shortestPath.isEmpty()) {
                        continue;
                    }
                    IAtom atomCurrent = shortestPath.pop().getAtom();

                    while (!shortestPath.empty()) {
                        final IAtom atomNext = shortestPath.pop().getAtom();
                        Map<IAtom, IBond> m = cache.get(atomCurrent);
                        final IBond[] b = {m != null ? m.get(atomNext) : null};
                        if (b[0] == null) {
                            b[0] = ac.getBond(atomCurrent, atomNext);
                            cache.put(atomCurrent,
                                    new HashMap<IAtom, IBond>() {

                                        {
                                            put(atomNext, b[0]);
                                        }
                                        private static final long serialVersionUID = 0xb3a7a32449fL;
                                    });
                        }

                        setAtom(atomCurrent, sb);
                        sb.append(getBondSymbol(b[0]));
                        setAtom(atomNext, sb);
                        atomCurrent = atomNext;
                    }
                    allPaths.add(sb.toString().trim());
                }
            }
        }
    }
}
