/* Copyright (C) 1997-2007,2011  Egon Willighagen <egonw@users.sf.net>
 * 
 * Contact: cdk-devel@lists.sourceforge.net
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All I ask is that proper credit is given for my work, which includes
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

import java.math.BigInteger;
import java.util.BitSet;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.Atom;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.FingerprinterTool;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.graph.AtomContainerAtomPermutor;
import org.openscience.cdk.graph.AtomContainerBondPermutor;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.templates.MoleculeFactory;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

/**
 * @cdk.module test-standard
 */
public class FingerprinterTest {

    boolean standAlone = false;
    private static ILoggingTool logger =
            LoggingToolFactory.createLoggingTool(FingerprinterTest.class);

    public IFingerprinter getBitFingerprinter() {
        return new ShortestPathFingerprinter();
    }

    @Test
    public void testRegression() throws Exception {
        IAtomContainer mol1 = MoleculeFactory.makeIndole();
        IAtomContainer mol2 = MoleculeFactory.makePyrrole();
        ShortestPathFingerprinter fingerprinter = new ShortestPathFingerprinter();
        IBitFingerprint bs1 = fingerprinter.getBitFingerprint(mol1);
        Assert.assertEquals("Seems the fingerprint code has changed. This will cause a number of other tests to fail too!", 18, bs1.cardinality());
        IBitFingerprint bs2 = fingerprinter.getBitFingerprint(mol2);
        Assert.assertEquals("Seems the fingerprint code has changed. This will cause a number of other tests to fail too!", 9, bs2.cardinality());
    }

    @Test
    public void testGetSize() throws java.lang.Exception {
        IFingerprinter fingerprinter = new ShortestPathFingerprinter(512);
        Assert.assertNotNull(fingerprinter);
        Assert.assertEquals(512, fingerprinter.getSize());
    }

    @Test
    public void testgetBitFingerprint_IAtomContainer() throws java.lang.Exception {
        ShortestPathFingerprinter fingerprinter = new ShortestPathFingerprinter();

        IAtomContainer mol = MoleculeFactory.makeIndole();
        IBitFingerprint bs = fingerprinter.getBitFingerprint(mol);
        Assert.assertNotNull(bs);
        Assert.assertEquals(fingerprinter.getSize(), bs.size());
    }

    @Test
    public void testFingerprinter() throws java.lang.Exception {
        ShortestPathFingerprinter fingerprinter = new ShortestPathFingerprinter();
        Assert.assertNotNull(fingerprinter);

        IAtomContainer mol = MoleculeFactory.makeIndole();
        BitSet bs = fingerprinter.getBitFingerprint(mol).asBitSet();
        IAtomContainer frag1 = MoleculeFactory.makePyrrole();
        BitSet bs1 = fingerprinter.getBitFingerprint(frag1).asBitSet();
        Assert.assertTrue(FingerprinterTool.isSubset(bs, bs1));
    }

    @Test
    public void testFingerprinter_int() throws java.lang.Exception {
        ShortestPathFingerprinter fingerprinter = new ShortestPathFingerprinter(512);
        Assert.assertNotNull(fingerprinter);

        IAtomContainer mol = MoleculeFactory.makeIndole();
        BitSet bs = fingerprinter.getBitFingerprint(mol).asBitSet();
        IAtomContainer frag1 = MoleculeFactory.makePyrrole();
        BitSet bs1 = fingerprinter.getBitFingerprint(frag1).asBitSet();
        Assert.assertTrue(FingerprinterTool.isSubset(bs, bs1));
    }

    @Test
    public void testFingerprinter_int_int() throws java.lang.Exception {
        ShortestPathFingerprinter fingerprinter = new ShortestPathFingerprinter(1024);
        Assert.assertNotNull(fingerprinter);

        IAtomContainer mol = MoleculeFactory.makeIndole();
        BitSet bs = fingerprinter.getBitFingerprint(mol).asBitSet();
        IAtomContainer frag1 = MoleculeFactory.makePyrrole();
        BitSet bs1 = fingerprinter.getBitFingerprint(frag1).asBitSet();
        Assert.assertTrue(FingerprinterTool.isSubset(bs, bs1));
    }

    @Test
    public void testFingerprinterBitSetSize() throws Exception {
        ShortestPathFingerprinter fingerprinter = new ShortestPathFingerprinter(1024);
        Assert.assertNotNull(fingerprinter);
        IAtomContainer mol = MoleculeFactory.makeIndole();
        BitSet bs = fingerprinter.getBitFingerprint(mol).asBitSet();
        Assert.assertEquals(1009, bs.length()); // highest set bit
        Assert.assertEquals(1024, bs.size()); // actual bit set size
    }

    /**
     * @cdk.bug 2819557
     *
     * @throws org.openscience.cdk.exception.CDKException
     */
    @Test
    public void testBug2819557() throws CDKException {
        IAtomContainer butane = makeButane();
        IAtomContainer propylAmine = makePropylAmine();

        ShortestPathFingerprinter fp = new ShortestPathFingerprinter();
        BitSet b1 = fp.getBitFingerprint(butane).asBitSet();
        BitSet b2 = fp.getBitFingerprint(propylAmine).asBitSet();

        System.out.println("f1:" + b1);
        System.out.println("f2:" + b2);
        
        Assert.assertFalse(FingerprinterTool.isSubset(b2, b1));
        Assert.assertFalse("butane should not be a substructure of propylamine", FingerprinterTool.isSubset(b2, b1));
    }

    @Test
    public void testBondPermutation() throws CDKException {
        IAtomContainer pamine = makePropylAmine();
        ShortestPathFingerprinter fp = new ShortestPathFingerprinter();
        IBitFingerprint bs1 = fp.getBitFingerprint(pamine);

        AtomContainerBondPermutor acp = new AtomContainerBondPermutor(pamine);
        while (acp.hasNext()) {
            IAtomContainer container = acp.next();
            IBitFingerprint bs2 = fp.getBitFingerprint(container);
            Assert.assertTrue(bs1.equals(bs2));
        }
    }

    @Test
    public void testAtomPermutation() throws CDKException {
        IAtomContainer pamine = makePropylAmine();
        ShortestPathFingerprinter fp = new ShortestPathFingerprinter();
        IBitFingerprint bs1 = fp.getBitFingerprint(pamine);

        AtomContainerAtomPermutor acp = new AtomContainerAtomPermutor(pamine);
        while (acp.hasNext()) {
            IAtomContainer container = acp.next();
            IBitFingerprint bs2 = fp.getBitFingerprint(container);
            Assert.assertTrue(bs1.equals(bs2));
        }
    }

    @Test
    public void testBondPermutation2() throws CDKException {
        IAtomContainer pamine = MoleculeFactory.makeCyclopentane();
        ShortestPathFingerprinter fp = new ShortestPathFingerprinter();
        IBitFingerprint bs1 = fp.getBitFingerprint(pamine);

        AtomContainerBondPermutor acp = new AtomContainerBondPermutor(pamine);
        while (acp.hasNext()) {
            IAtomContainer container = acp.next();
            IBitFingerprint bs2 = fp.getBitFingerprint(container);
            Assert.assertTrue(bs1.equals(bs2));
        }
    }

    @Test
    public void testAtomPermutation2() throws CDKException {
        IAtomContainer pamine = MoleculeFactory.makeCyclopentane();
        ShortestPathFingerprinter fp = new ShortestPathFingerprinter();
        IBitFingerprint bs1 = fp.getBitFingerprint(pamine);

        AtomContainerAtomPermutor acp = new AtomContainerAtomPermutor(pamine);
        while (acp.hasNext()) {
            IAtomContainer container = acp.next();
            IBitFingerprint bs2 = fp.getBitFingerprint(container);
            Assert.assertTrue(bs1.equals(bs2));
        }
    }

    public static IAtomContainer makeFragment1() {
        IAtomContainer mol = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
        mol.addAtom(new Atom("C")); // 0
        mol.addAtom(new Atom("C")); // 1
        mol.addAtom(new Atom("C")); // 2
        mol.addAtom(new Atom("C")); // 3
        mol.addAtom(new Atom("C")); // 4
        mol.addAtom(new Atom("C")); // 5
        mol.addAtom(new Atom("C")); // 6

        mol.addBond(0, 1, IBond.Order.SINGLE); // 1
        mol.addBond(0, 2, IBond.Order.SINGLE); // 2
        mol.addBond(0, 3, IBond.Order.SINGLE); // 3
        mol.addBond(0, 4, IBond.Order.SINGLE); // 4
        mol.addBond(3, 5, IBond.Order.SINGLE); // 5
        mol.addBond(5, 6, IBond.Order.DOUBLE); // 6
        return mol;
    }

    public static IAtomContainer makeFragment4() {
        IAtomContainer mol = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
        mol.addAtom(new Atom("C")); // 0
        mol.addAtom(new Atom("C")); // 1

        mol.addBond(0, 1, IBond.Order.SINGLE); // 1
        return mol;
    }

    public static IAtomContainer makeFragment2() {
        IAtomContainer mol = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
        mol.addAtom(new Atom("C")); // 0
        mol.addAtom(new Atom("C")); // 1
        mol.addAtom(new Atom("C")); // 2
        mol.addAtom(new Atom("S")); // 3
        mol.addAtom(new Atom("O")); // 4
        mol.addAtom(new Atom("C")); // 5
        mol.addAtom(new Atom("C")); // 6

        mol.addBond(0, 1, IBond.Order.DOUBLE); // 1
        mol.addBond(0, 2, IBond.Order.SINGLE); // 2
        mol.addBond(0, 3, IBond.Order.SINGLE); // 3
        mol.addBond(0, 4, IBond.Order.SINGLE); // 4
        mol.addBond(3, 5, IBond.Order.SINGLE); // 5
        mol.addBond(5, 6, IBond.Order.DOUBLE); // 6
        mol.addBond(5, 6, IBond.Order.DOUBLE); // 7
        return mol;
    }

    public static IAtomContainer makeFragment3() {
        IAtomContainer mol = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
        mol.addAtom(new Atom("C")); // 0
        mol.addAtom(new Atom("C")); // 1
        mol.addAtom(new Atom("C")); // 2
        mol.addAtom(new Atom("C")); // 3
        mol.addAtom(new Atom("C")); // 4
        mol.addAtom(new Atom("C")); // 5
        mol.addAtom(new Atom("C")); // 6

        mol.addBond(0, 1, IBond.Order.SINGLE); // 1
        mol.addBond(0, 2, IBond.Order.SINGLE); // 2
        mol.addBond(0, 3, IBond.Order.SINGLE); // 3
        mol.addBond(0, 4, IBond.Order.SINGLE); // 4
        mol.addBond(3, 5, IBond.Order.DOUBLE); // 5
        mol.addBond(5, 6, IBond.Order.SINGLE); // 6
        return mol;
    }

    public static IAtomContainer makeButane() {
        IAtomContainer mol = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
        mol.addAtom(new Atom("C")); // 0
        mol.addAtom(new Atom("C")); // 1
        mol.addAtom(new Atom("C")); // 2
        mol.addAtom(new Atom("C")); // 3

        mol.addBond(0, 1, IBond.Order.SINGLE); // 1
        mol.addBond(1, 2, IBond.Order.SINGLE); // 2
        mol.addBond(2, 3, IBond.Order.SINGLE); // 3

        return mol;
    }

    public static IAtomContainer makePropylAmine() {
        IAtomContainer mol = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
        mol.addAtom(new Atom("C")); // 0
        mol.addAtom(new Atom("C")); // 1
        mol.addAtom(new Atom("C")); // 2
        mol.addAtom(new Atom("N")); // 3

        mol.addBond(0, 1, IBond.Order.SINGLE); // 1
        mol.addBond(1, 2, IBond.Order.SINGLE); // 2
        mol.addBond(2, 3, IBond.Order.SINGLE); // 3

        return mol;
    }

    public static void main(String[] args) throws Exception {
        BigInteger bi = new BigInteger("0");
        bi = bi.add(BigInteger.valueOf((long) Math.pow(2, 63)));
        System.err.println(bi.toString());
        bi = bi.add(BigInteger.valueOf((long) Math.pow(2, 0)));
        System.err.println(bi.toString());
        FingerprinterTest fpt = new FingerprinterTest();
        fpt.standAlone = true;
        //fpt.testFingerprinter();
        //fpt.testFingerprinterArguments();
        //fpt.testBug706786();
        //fpt.testBug771485();
        //fpt.testBug853254();
        //fpt.testBug931608();
//        fpt.testBug934819();
    }
}
