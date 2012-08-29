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
package fingerprints;

import fingerprints.helper.ShortestPathWalker;
import fingerprints.helper.RandomNumber;
import java.io.Serializable;
import java.util.*;
import java.util.logging.Level;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.*;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.interfaces.IAtomType.Hybridization;
import org.openscience.cdk.ringsearch.SSSRFinder;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.RingSetManipulator;
import org.openscience.cdk.tools.periodictable.PeriodicTable;

/**
 * Generates a fingerprint for a given {@link IAtomContainer}. Fingerprints are one-dimensional bit arrays, where bits
 * are set according to a the occurrence of a particular structural feature (See for example the Daylight inc. theory
 * manual for more information). Fingerprints allow for a fast screening step to exclude candidates for a substructure
 * search in a database. They are also a means for determining the similarity of chemical structures.
 *
 * A fingerprint is generated for an AtomContainer with this code:
 * <pre>
 *   AtomContainer molecule = new AtomContainer();
 *   IFingerprinter fingerprinter = new ShortestPathFingerprinter();
 *   IBitFingerprint fingerprint = fingerprinter.getFingerprint(molecule);
 *   fingerprint.fingerprintLength(); // returns 1024 by default
 *   fingerprint.length(); // returns the highest set bit
 * </pre> <p>
 *
 * <p>The FingerPrinter calculates fingerprint based on the Shortest Paths between two atoms. It also takes care of the
 * ring system and charges, if desired.
 *
 * <p>The FingerPrinter assumes that hydrogens are explicitly given! Furthermore, if pseudo atoms or atoms with
 * malformed symbols are present, their atomic number is taken as one more than the last element currently supported in {@link PeriodicTable}.
 *
 * <p>Unlike the {@link Fingerprinter}, this fingerprinter does not take into account aromaticity. Instead, it takes
 * into account SP2
 * {@link Hybridization}.
 *
 * @author Syed Asad Rahman (2012) @cdk.keyword fingerprint @cdk.keyword similarity @cdk.module standard @cdk.githash
 */
public class ShortestPathFingerprinter extends RandomNumber implements IFingerprinter, Serializable {

    /**
     * The default length of created fingerprints.
     */
    public final static int DEFAULT_SIZE = 1024;
    private static final long serialVersionUID = 7867864332244557861L;
    /**
     * The default length of created fingerprints.
     */
    private int fingerprintLength;
    private static ILoggingTool logger =
            LoggingToolFactory.createLoggingTool(ShortestPathFingerprinter.class);

    /**
     * Creates a fingerprint generator of length
     * <code>DEFAULT_SIZE</code>
     */
    public ShortestPathFingerprinter() {
        this(DEFAULT_SIZE);
    }

    /**
     * Constructs a fingerprint generator that creates fingerprints of the given fingerprintLength, using a generation
     * algorithm with shortest paths.
     *
     * @param fingerprintLength The desired fingerprintLength of the fingerprint
     */
    public ShortestPathFingerprinter(int fingerprintLength) {
        this.fingerprintLength = fingerprintLength;
    }

    /**
     * Generates a fingerprint of the default fingerprintLength for the given AtomContainer.
     *
     * @param atomContainer The AtomContainer for which a Fingerprint is generated
     * @exception CDKException if there is a timeout in ring or aromaticity perception
     * @return A {@link BitSet} representing the fingerprint
     */
    @Override
    public IBitFingerprint getBitFingerprint(
            IAtomContainer atomContainer)
            throws CDKException {
        logger.debug("Entering Fingerprinter");
        logger.debug("Starting Aromaticity Detection");
        long before = System.currentTimeMillis();
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(atomContainer);
        CDKHueckelAromaticityDetector.detectAromaticity(atomContainer);
        long after = System.currentTimeMillis();
        logger.debug("time for aromaticity calculation: "
                + (after - before) + " milliseconds");
        logger.debug("Finished Aromaticity Detection");
        BitSet bitSet = new BitSet(fingerprintLength);
        if (!ConnectivityChecker.isConnected(atomContainer)) {
            IAtomContainerSet partitionedMolecules = ConnectivityChecker.partitionIntoMolecules(atomContainer);
            for (IAtomContainer container : partitionedMolecules.atomContainers()) {
                addUniquePath(container, bitSet);
            }
        } else {
            addUniquePath(atomContainer, bitSet);
        }
        return new BitSetFingerprint(bitSet);
    }

    private void addUniquePath(IAtomContainer container, BitSet bitSet) {
        try {
            Integer[] hashes = findPaths(container);
            for (Integer hash : hashes) {
                int position = (int) generateMersenneTwisterRandomNumber(fingerprintLength, hash.intValue());
                bitSet.set(position);
            }
        } catch (CloneNotSupportedException ex) {
            logger.error(Level.SEVERE, null, ex);
        } catch (CDKException ex) {
            logger.error(Level.SEVERE, null, ex);
        }
    }

    /**
     * {@inheritDoc}
     *
     * @param atomContainer
     * @return
     * @throws CDKException
     */
    @Override
    public Map<String, Integer> getRawFingerprint(IAtomContainer atomContainer) throws CDKException {
        Map<String, Integer> uniquePaths = new TreeMap<String, Integer>();
        if (!ConnectivityChecker.isConnected(atomContainer)) {
            IAtomContainerSet partitionedMolecules = ConnectivityChecker.partitionIntoMolecules(atomContainer);
            for (IAtomContainer container : partitionedMolecules.atomContainers()) {
                addUniquePath(container, uniquePaths);
            }
        } else {
            addUniquePath(atomContainer, uniquePaths);
        }
        return uniquePaths;
    }

    private void addUniquePath(IAtomContainer atomContainer, Map<String, Integer> uniquePaths) {
        Integer[] hashes;
        try {
            hashes = findPaths(atomContainer);
            for (Integer hash : hashes) {
                int position = (int) generateMersenneTwisterRandomNumber(fingerprintLength, hash.intValue());
                uniquePaths.put(new Integer(position).toString(), hash);
            }
        } catch (CloneNotSupportedException ex) {
            logger.error(Level.SEVERE, null, ex);
        } catch (CDKException ex) {
            logger.error(Level.SEVERE, null, ex);
        }
    }

    /**
     * Get all paths of lengths 0 to the specified length.
     *
     * This method will find all paths upto length N starting from each atom in the molecule and return the unique set
     * of such paths.
     *
     * @param container The molecule to search
     * @return A map of path strings, keyed on themselves
     * @throws CloneNotSupportedException
     * @throws CDKException
     */
    protected Integer[] findPaths(IAtomContainer container) throws CloneNotSupportedException, CDKException {

        ShortestPathWalker walker = new ShortestPathWalker(container);
        // convert paths to hashes
        List<Integer> paths = new ArrayList<Integer>();
        int patternIndex = 0;

        for (String s : walker.getPaths()) {
            int toHashCode = s.hashCode();
            paths.add(patternIndex, toHashCode);
            patternIndex++;
        }
        SSSRFinder finder = new SSSRFinder(container);
        IRingSet sssr = finder.findEssentialRings();
        RingSetManipulator.sort(sssr);
        int ringCounter = sssr.getAtomContainerCount();
        for (Iterator<IAtomContainer> it = sssr.atomContainers().iterator(); it.hasNext();) {
            IAtomContainer ring = it.next();
            int toHashCode = String.valueOf(ringCounter * ring.getAtomCount()).hashCode();
            paths.add(patternIndex, toHashCode);
            patternIndex++;
            ringCounter--;
        }

        List<String> l = new ArrayList<String>();
        for (Iterator<IAtom> it = container.atoms().iterator(); it.hasNext();) {
            IAtom atom = it.next();
            int charge = atom.getFormalCharge() == null ? 0 : atom.getFormalCharge().intValue();
            if (charge != 0) {
                l.add(atom.getSymbol().concat(String.valueOf(charge)));
            }
        }
        Collections.sort(l);
        int toHashCode = l.hashCode();
        paths.add(patternIndex, toHashCode);
        patternIndex++;

        l = new ArrayList<String>();
        /*
         * atom stereo parity
         */
        for (Iterator<IAtom> it = container.atoms().iterator(); it.hasNext();) {
            IAtom atom = it.next();
            int st = atom.getStereoParity() == null ? 0 : atom.getStereoParity().intValue();
            if (st != 0) {
                l.add(atom.getSymbol().concat(String.valueOf(st)));
            }
        }
        Collections.sort(l);
        toHashCode = l.hashCode();
        paths.add(patternIndex, toHashCode);
        patternIndex++;

        if (container.getSingleElectronCount() > 0) {
            StringBuilder radicalInformation = new StringBuilder();
            radicalInformation.append("RAD: ".concat(String.valueOf(container.getSingleElectronCount())));
            paths.add(patternIndex, radicalInformation.toString().hashCode());
            patternIndex++;
        }
        if (container.getLonePairCount() > 0) {
            StringBuilder lpInformation = new StringBuilder();
            lpInformation.append("LP: ".concat(String.valueOf(container.getLonePairCount())));
            paths.add(patternIndex, lpInformation.toString().hashCode());
            patternIndex++;
        }

        return paths.toArray(new Integer[paths.size()]);
    }

    @Override
    public int getSize() {
        return fingerprintLength;
    }

    @Override
    public ICountFingerprint getCountFingerprint(IAtomContainer iac) throws CDKException {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
