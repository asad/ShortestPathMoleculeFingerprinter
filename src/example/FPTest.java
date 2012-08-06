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

package example;


import fingerprints.HashedSPFingerprinter;
import java.util.BitSet;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.inchi.InChIToStructure;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

/**
 *
 * @author Syed Asad Rahman (2012) 
 * @cdk.keyword fingerprint 
 * @cdk.keyword similarity 
 * @cdk.module standard
 * @cdk.githash
 */
public class FPTest {

    /**
     * @param args the command line arguments
     * @throws CDKException  
     */
    public static void main(String[] args) throws CDKException {
        int n = 7; //don't run with n > 6, unless on a supercomputer
        String[] problemStructures = {
            "InChI=1S/C23H30/c24-15(25-24)4-2-1-3(2)6(1,2,4,15)5(2,4,15,36(15)32(15)45(4,5,15,36)43(4,5,15)28(4)15)7(1,2,3,4,6,15,27(2)41(1,2,7)35(1,2)7)12(1,2,3,6)10-11(12)13(1,3,10,12,30-31(13)46(3,12,13,30)42(1,3,12)13)19(3,10,11,12)9(6)8(4,5,6)14(4,5,6,9,15)17(8,9,19)16(9,10,11,19)18(8,9,14,17,19,47(8,9,14,17)48(8,9,14)37(8,14)38(8,14)48)20(3,9,10,11,12,13,16,17,19)21(10,11,13,19)22(10,11,19,20,26-39(10,21)22,33(21)34(21)22,40(21)44(11,21,22)49(10,11,21,22)40)23(9,10,11,16,17,18,19,20,21)29(16)52(16,17,18,23)50(16,17,18,23)51(16,17,18,23)52/h28H",
            "InChI=1S/C4Cl12/c5-1-2(5)3(1)4(1,2)7(3)8(3,4)12(1,2,3,4)10(1,2,3,4)6(1,2,5,14(1,2,3,4,5,10,12)16(1,2,3,4,7,8,10)12)13(1,2,3,4,5)9(1,2,3,4,5)11(1,2,3,4,7,13)15(1,2,3,4,7,8,9)13",
            "InChI=1S/C30H28N6O6S4/c37-27-23(54(27)55(23)27)25(27,60(23,27)61(23,27)37)11-12(25)24(11,43-44-25)14-19(9-6-2-3(6,48-2,51(2)6)8(2,6,9)17(6,9,19,52(8)9)28(9,14,19,24,63(8,9,17)19)34(14,17,19,24,57(14,19)28)36(11,12,14,24,28,66(14,24,28)34)40(11,12,24)33(11,12,23,25,27,59(23,25)27)41(11,12,25,36)40)21-15-4-1(47-4)5(4,15,50(1)4)16(4,15,21)29(15,21)22(21)26(29)10-7-18(10)13-20(18,42(13,49-13)58(13)64(13,20,42)65(13,20,42)58)31(7,10,18,38(7)18,53(18)20)46(18,20)45(10,26)39(7,10,26)32(7,10,22,26)35(21,22,26,29,56(22)32)30(15,16,21,22,26,29,62(5,15,16)29,67(15,16,26)29)68(21,22,29)35/h1-2,14,20,22,37H",
            "InChI=1S/C76H68/c77-10(5-15-17(5)27(15,107(15,17)81-17)21-22-25(15,21,27,99(21)27)26(21,22)16(5)18(5)28(16,22,26,100(22)26)108(16,18)82-18)19(77)20(10,78-10)29(19)30(19,20)41(29)31-32(41)46-40-48(32,46,86-106(46)48)58(29,30,31,32,41,102(32)48,110(29,30)41)57(29,30,31,32,41,109(29,30)41)47(31,101(31)57)39-45(31,47,105(47)85-47)53(39)65(39,45,125(39,53,89-53)123(39,45,47)53)66(40,46)54(40,46,90-126(40,54,66)124(40,46,48)54)70(65,66)61-35-36(61,62(35,61,70)69(53,61,65,66,70)72(61,62,65,66,70,129(61,62,69,93-61)133(61,65,69,72)121(53,65)69,130(61,62,70,94-62)134(62,66,70,72)122(54,66)70)76(39,40,45,46,53,54,65,66,69,70)73(31,32,39,40,41,45,46,47,48,57)58)42(35)49-33-37-43(33,49)51(33,49,87-117(33,51)113(33,43)51)59(35,36,42,49,103(49)51,111(35,36)42)50(42,49)34-38-44(34,50)52(34,50,88-118(34,52)114(34,44)52)60(35,36,42,49,50,59,104(50)52,112(35,36)42)75(33,34,42,43,44,49,50,51,52,59)74(33,34,37,38,43,44)63(37,43,115(37,43)83-37)64(38,44,74,116(38,44)84-38)67(37,63,74)55-2-1(4-6-8(11(4,6)79-97(6)11)13-14-9-7(4,12(4,9)80-98(7)12)24(9,13,14,96(9)14)23(6,8,13,14)95(8)13)3(2,55)56(2,55,67)68(38,55,63,64,67,74)71(55,56,63,64,67,74,127(55,56,67,91-55)131(55,63,67,71)119(37,63)67)128(55,56,68,92-56)132(56,64,68,71)120(38,64)68/h8-9,13-14,19-22,27-28H",
            "InChI=1S/C12H37O3/c16-5-1-3-6(1,5)9(1,3,5,16,21-39(9)22-9,10(1,3,5,6,17-5,29(1)30(1)10,41(1,5)18-5,26-48(1,3,6,9)10)43(6)37(6)44(6,10)43)13(3)4-2-7(33(5)34(5)7)8(2,4,13)11(2,4,7,13,23-40(11)24-11,35(7)36(7)11,14(3,4,6,9,13,27(3)6)15(3,4,8,9,11,13)28(4)8)12(2,4,7,8,19-7,31(2)32(2)12,42(2,7)20-7,25-47(2,4,8)12)45(8)38(8)46(8,12)45/h1-4H",
            "InChI=1S/C18H44/c19-45-10-6-7(10,33(6)34(6)7,35(6)36(6)7)11(6,10)14(10,45)3-5(14,31(3)32(3)5)8(3,10,14)15(3,5,14,22-5,44(3,5,8)14,17(6,7,8,10,11,14,45)18(6,7,8,10,11,14,23-8,24-8,38(10)17,48(6,7,10,11,17)49(6,7,10,11,17)18)39(11)50(8,10,11,17,18)52(8,10,11,14,17,18)42(11,17)25-11)4(14,41(14,15)45)9(3,5,15)1-2(9,27(1)28(1)2,29(1)30(1)2)12(1,3,5,9,15)13(1,2,4,9,15,20-4,21-4)16(1,2,3,4,5,9,12,14,15,37(9)13,46(1,2,9,12)13)43(12,26-12)51(4,9,12,13,15,16)47(4,9,12,13,16)40(12)13/h1-8,14,19H",
            "InChI=1S/C32.C18H13PS2.C10H15.C2Cl2.2CF3O3S.2Ru/c1-6-9(1)18(1,6)2-7-8(2)15(7)3-5-4-16-10-11(16)12(10)19-13-14(19)17(6,13,19)21(13,14,19)20(9,18)22(15)23(16)24(3,4,5,22)27(3,5,15,22,23)26(2,7,8,18,20,22)25(1,6,9,18,20,21,32(7,8,15,20,22,24,26)27)31(13,14,17,19,21,23)29(10,11,12,16,23,24)28(4,5,16,22,23,24,27)30(10,11,12,19,21,23,29)31;20-17-12-6-4-10-15(17)19(14-8-2-1-3-9-14)16-11-5-7-13-18(16)21;1-6-7(2)9(4)10(5)8(6)3;3-1-2(3)4-1;2*2-1(3,4)8(5,6)7;;/h;1-13H;1-5H3;;;;;",
            "InChI=1S/C8H16Cl16/c25-1-3-5(1)7(1,3,25,29(3)5)11(1)9(1,19(1,3,5,7,11)17(3,5,7)21(1,3,5,7,19,33(3,5)17,37(5,17)39(3,5,17)21)23(1,3,5,7,9,11,17,19,27(1)19)31(1,7)25)15(3,5,7,11)13(3,5,7,35(3,7)15)14-4-2-6(4,14)8(2,4,14,30(4)6)12(2)10(2,16(4,6,8,12,14)36(4,8)14)20(2,4,6,8,12)18(4,6,8)22(2,4,6,8,20,34(4,6)18,38(6,18)40(4,6,18)22)24(2,4,6,8,10,12,18,20,28(2)20)32(2,8)26(2)8"
        };

        for (String inchi : problemStructures) {

            InChIGeneratorFactory f = InChIGeneratorFactory.getInstance();
            InChIToStructure c = f.getInChIToStructure(inchi, SilentChemObjectBuilder.getInstance());
            //Assert.assertEquals(INCHI_RET.OKAY, c.getReturnStatus());
            IAtomContainer mol = c.getAtomContainer();

            HashedSPFingerprinter fp = new HashedSPFingerprinter(1024);

// here add code to configure atom types  - and make sure to set atomic numbers as well, 
//because perceiveAtomType() method doesn't set them, and they are null after JNIInchi 


            long now = System.currentTimeMillis();
            BitSet bitset = fp.getBitFingerprint(mol).asBitSet();
            System.out.println(String.format(
                    "Atoms\t%d"
                    + "\tBonds\t%d"
                    + "\tFP Density\t%d"
                    + "\tElapsed time\t%d"
                    + "\tBitset\t%s",
                    mol.getAtomCount(), mol.getBondCount(),
                    bitset.cardinality(), (System.currentTimeMillis() - now), bitset));
            /*
             * depth 1 : 188 ms depth 2 : 291 ms depth 3 : 525 ms depth 4 : 1474 ms depth 5 : 9147 ms depth 6 : 71634 ms
             */

            System.out.println();
        }
    }
}
