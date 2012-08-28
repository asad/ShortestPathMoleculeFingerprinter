/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package fingerprints;

import java.util.BitSet;
import java.util.Map;
import java.util.TreeMap;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.fingerprint.FingerprinterTool;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.similarity.Tanimoto;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 *
 * @author Asad
 */
public class RingSearchSPTestCases {

    @Test
    public void testGenerateFingerprintNaphthalene() throws InvalidSmilesException, Exception {

        String smiles = "C1=CC2=CC=CC=C2C=C1";
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        System.out.println("Atom count " + molecule.getAtomCount());
        ShortestPathFingerprinter fingerprint = new ShortestPathFingerprinter(1024);
        fingerprint.setRespectRingMatches(true);
        fingerprint.setRespectFormalCharges(true);
        BitSet fingerprint1;
        fingerprint1 = fingerprint.getBitFingerprint(molecule).asBitSet();
        System.out.println("Naphthalene fp " + fingerprint1.toString());
    }

    @Test
    public void testGenerateFingerprintAnthracene() throws InvalidSmilesException, Exception {

        String smiles = "C1=CC2=CC3=CC=CC=C3C=C2C=C1";
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        System.out.println("Atom count " + molecule.getAtomCount());
        ShortestPathFingerprinter fingerprint = new ShortestPathFingerprinter(1024);
        fingerprint.setRespectRingMatches(true);
        fingerprint.setRespectFormalCharges(true);
        BitSet fingerprint1;
        fingerprint1 = fingerprint.getBitFingerprint(molecule).asBitSet();
        System.out.println("Anthracene fp " + fingerprint1.toString());
    }

    @Test
    public void testGenerateFingerprintMultiphtalene() throws InvalidSmilesException, Exception {

        String smiles = "C1=CC2=CC=C3C4=CC5=CC6=CC=CC=C6C=C5C=C4C=CC3=C2C=C1";
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
        System.out.println("Atom count " + molecule.getAtomCount());
        ShortestPathFingerprinter fingerprint = new ShortestPathFingerprinter(1024);
        fingerprint.setRespectRingMatches(true);
        fingerprint.setRespectFormalCharges(true);
        BitSet fingerprint1;
        fingerprint1 = fingerprint.getBitFingerprint(molecule).asBitSet();
        System.out.println(" Multiphtalene fp " + fingerprint1.toString());
    }

    @Test
    public void testGenerateFingerprintSimilarity() throws InvalidSmilesException, Exception {

        ShortestPathFingerprinter fingerprint = new ShortestPathFingerprinter(1024);
//        fingerprint.setRespectRingMatches(true);
//        fingerprint.setRespectFormalCharges(true);
        Map<String, String> mols = new TreeMap<String, String>();
        mols.put("Multiphtalene", "C1=CC2=CC=C3C4=CC5=CC6=CC=CC=C6C=C5C=C4C=CC3=C2C=C1");
        mols.put("Anthracene", "C1=CC2=CC3=CC=CC=C3C=C2C=C1");
        mols.put("Naphthalene", "C1=CC2=CC=CC=C2C=C1");

        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());

        System.out.println("Query\tTarget\tSimilarity\tSubgraph");

        for (String key1 : mols.keySet()) {
            IAtomContainer molecule1 = smilesParser.parseSmiles(mols.get(key1));
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule1);
            BitSet fingerprint1 = fingerprint.getBitFingerprint(molecule1).asBitSet();
            for (String key2 : mols.keySet()) {
                IAtomContainer molecule2 = smilesParser.parseSmiles(mols.get(key2));
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule2);
                BitSet fingerprint2 = fingerprint.getBitFingerprint(molecule2).asBitSet();
                int flag = (FingerprinterTool.isSubset(fingerprint1, fingerprint2)) ? 1 : 0;
                float calculate = Tanimoto.calculate(fingerprint1, fingerprint2) + 0.0f;
                String format = String.format("%s\t%s\t%f\t%d", key1, key2, calculate, flag);
                System.out.println(format);

            }
        }
    }
}
