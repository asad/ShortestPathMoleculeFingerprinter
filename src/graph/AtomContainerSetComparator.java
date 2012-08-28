/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package graph;

import java.io.Serializable;
import java.util.Comparator;
import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.config.IsotopeFactory;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 *
 * @author Asad
 *
 * AtomContainerSetComparator Comparator
 */
public class AtomContainerSetComparator implements Comparator<IAtomContainer>, Serializable {

    private static final long serialVersionUID = 74747883636787861L;
    /**
     * Configure LoggingTool
     */
    private ILoggingTool logger =
            LoggingToolFactory.createLoggingTool(AtomContainerSetComparator.class);

    /**
     * Creates a new instance of AtomContainerComparator
     */
    public AtomContainerSetComparator() {
    }

    /*
     * <p>Compares two IAtomContainers for order with the following criteria with decreasing priority:</p> <ul>
     * <li>Compare atom count <li>Compare molecular weight (heavy atoms only) <li>Compare bond count <li>Compare sum of
     * bond orders (heavy atoms only) </ul> <p>If no difference can be found with the above criteria, the
     * IAtomContainers are considered equal.</p> <p>Returns a negative integer, zero, or a positive integer as the first
     * argument is less than, equal to, or greater than the second.</p> <p>This method is null safe.</p>
     *
     * @param o1 the first IAtomContainer @param o2 the second IAtomContainer @return a negative integer, zero, or a
     * positive integer as the first argument is less than, equal to, or greater than the second.
     */
    /**
     *
     * @param o1
     * @param o2
     * @return
     */
    @TestMethod("testCompare_Object_Object")
    @Override
    public int compare(IAtomContainer o1, IAtomContainer o2) {
        // Check for nulls
        if (o1 == null && o2 == null) {
            return 0;
        }
        if (o1 == null) {
            return -1;
        }
        if (o2 == null) {
            return 1;
        }

        // Check for correct instances
        if (!(o1 instanceof IAtomContainer) && !(o2 instanceof IAtomContainer)) {
            return 0;
        }
        if (!(o1 instanceof IAtomContainer)) {
            return -1;
        }
        if (!(o2 instanceof IAtomContainer)) {
            return 1;
        }

        // Check for correct instances
        if (!(o1 instanceof IAtomContainer) && !(o2 instanceof IAtomContainer)) {
            return 0;
        }
        if (!(o1 instanceof IAtomContainer)) {
            return -1;
        }
        if (!(o2 instanceof IAtomContainer)) {
            return 1;
        }

        IAtomContainer atomContainer1 = o1;
        IAtomContainer atomContainer2 = o2;

        // 1. Compare atom count
        if (atomContainer1.getAtomCount() > atomContainer2.getAtomCount()) {
            return -1;
        } else if (atomContainer1.getAtomCount() < atomContainer2.getAtomCount()) {
            return 1;
        } else {
            // 2. Atom count equal, compare molecular weight (heavy atoms only)
            double mw1;
            double mw2;
            try {
                mw1 = getMolecularWeight(atomContainer1);
                mw2 = getMolecularWeight(atomContainer2);
            } catch (CDKException e) {
                logger.warn("Exception in molecular mass calculation.");
                return 0;
            }
            if (mw1 > mw2) {
                return -1;
            } else if (mw1 < mw2) {
                return 1;
            } else {
                // 3. Molecular weight equal, compare bond count
                if (atomContainer1.getBondCount() > atomContainer2.getBondCount()) {
                    return -1;
                } else if (atomContainer1.getBondCount() < atomContainer2.getBondCount()) {
                    return 1;
                } else {
                    // 4. Bond count equal, compare sum of bond orders (heavy atoms only)
                    double bondOrderSum1 = AtomContainerManipulator.getSingleBondEquivalentSum(atomContainer1);
                    double bondOrderSum2 = AtomContainerManipulator.getSingleBondEquivalentSum(atomContainer2);
                    if (bondOrderSum1 > bondOrderSum2) {
                        return -1;
                    } else if (bondOrderSum1 < bondOrderSum2) {
                        return 1;
                    }
                }

            }
        }
        // AtomContainers are equal in terms of this comparator
        return 0;
    }

    /**
     * Returns the molecular weight (exact mass) of the major isotopes of all heavy atoms of the given IAtomContainer.
     *
     * @param atomContainer an IAtomContainer to calculate the molecular weight for
     * @throws org.openscience.cdk.exception.CDKException if an error occurs with the IsotopeFactory
     * @return the molecular weight (exact mass) of the major isotopes of all heavy atoms of the given IAtomContainer
     */
    private double getMolecularWeight(IAtomContainer atomContainer) throws CDKException {
        double mw = 0.0;
        try {
            for (IAtom atom : atomContainer.atoms()) {
                if (!atom.getSymbol().equals("H") && !atom.getSymbol().equals("R")) {
                    try {
                        mw += IsotopeFactory.getInstance(atomContainer.getBuilder()).getMajorIsotope(atom.getSymbol()).getExactMass();
                    } catch (Exception e) {
                        logger.warn("Molecular weight calculation failes for atom " + atom.getSymbol());
                    }
                } else if (atom.getSymbol().equals("R")) {
                    mw += IsotopeFactory.getInstance(atomContainer.getBuilder()).getMajorIsotope("C").getExactMass();
                }
            }
        } catch (Exception e) {
        }
        return mw;
    }
}
