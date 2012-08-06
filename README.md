ShortestPathMoleculeFingerprinter
=================================

This is a hashed fingerprint based on the shortest path. The idea is to skip the exhaustive path based finger printers where the runtime is exponential.

Nina's example for CDK Hashed Fingerprinter failure.


            "InChI=1S/C23H30/c24-15(25-24)4-2-1-3(2)6(1,2,4,15)5(2,4,15,36(15)32(15)45(4,5,15,36)43(4,5,15)28(4)15)7(1,2,3,4,6,15,27(2)41(1,2,7)35(1,2)7)12(1,2,3,6)10-11(12)13(1,3,10,12,30-31(13)46(3,12,13,30)42(1,3,12)13)19(3,10,11,12)9(6)8(4,5,6)14(4,5,6,9,15)17(8,9,19)16(9,10,11,19)18(8,9,14,17,19,47(8,9,14,17)48(8,9,14)37(8,14)38(8,14)48)20(3,9,10,11,12,13,16,17,19)21(10,11,13,19)22(10,11,19,20,26-39(10,21)22,33(21)34(21)22,40(21)44(11,21,22)49(10,11,21,22)40)23(9,10,11,16,17,18,19,20,21)29(16)52(16,17,18,23)50(16,17,18,23)51(16,17,18,23)52/h28H",
            "InChI=1S/C4Cl12/c5-1-2(5)3(1)4(1,2)7(3)8(3,4)12(1,2,3,4)10(1,2,3,4)6(1,2,5,14(1,2,3,4,5,10,12)16(1,2,3,4,7,8,10)12)13(1,2,3,4,5)9(1,2,3,4,5)11(1,2,3,4,7,13)15(1,2,3,4,7,8,9)13",
            "InChI=1S/C30H28N6O6S4/c37-27-23(54(27)55(23)27)25(27,60(23,27)61(23,27)37)11-12(25)24(11,43-44-25)14-19(9-6-2-3(6,48-2,51(2)6)8(2,6,9)17(6,9,19,52(8)9)28(9,14,19,24,63(8,9,17)19)34(14,17,19,24,57(14,19)28)36(11,12,14,24,28,66(14,24,28)34)40(11,12,24)33(11,12,23,25,27,59(23,25)27)41(11,12,25,36)40)21-15-4-1(47-4)5(4,15,50(1)4)16(4,15,21)29(15,21)22(21)26(29)10-7-18(10)13-20(18,42(13,49-13)58(13)64(13,20,42)65(13,20,42)58)31(7,10,18,38(7)18,53(18)20)46(18,20)45(10,26)39(7,10,26)32(7,10,22,26)35(21,22,26,29,56(22)32)30(15,16,21,22,26,29,62(5,15,16)29,67(15,16,26)29)68(21,22,29)35/h1-2,14,20,22,37H",
            "InChI=1S/C76H68/c77-10(5-15-17(5)27(15,107(15,17)81-17)21-22-25(15,21,27,99(21)27)26(21,22)16(5)18(5)28(16,22,26,100(22)26)108(16,18)82-18)19(77)20(10,78-10)29(19)30(19,20)41(29)31-32(41)46-40-48(32,46,86-106(46)48)58(29,30,31,32,41,102(32)48,110(29,30)41)57(29,30,31,32,41,109(29,30)41)47(31,101(31)57)39-45(31,47,105(47)85-47)53(39)65(39,45,125(39,53,89-53)123(39,45,47)53)66(40,46)54(40,46,90-126(40,54,66)124(40,46,48)54)70(65,66)61-35-36(61,62(35,61,70)69(53,61,65,66,70)72(61,62,65,66,70,129(61,62,69,93-61)133(61,65,69,72)121(53,65)69,130(61,62,70,94-62)134(62,66,70,72)122(54,66)70)76(39,40,45,46,53,54,65,66,69,70)73(31,32,39,40,41,45,46,47,48,57)58)42(35)49-33-37-43(33,49)51(33,49,87-117(33,51)113(33,43)51)59(35,36,42,49,103(49)51,111(35,36)42)50(42,49)34-38-44(34,50)52(34,50,88-118(34,52)114(34,44)52)60(35,36,42,49,50,59,104(50)52,112(35,36)42)75(33,34,42,43,44,49,50,51,52,59)74(33,34,37,38,43,44)63(37,43,115(37,43)83-37)64(38,44,74,116(38,44)84-38)67(37,63,74)55-2-1(4-6-8(11(4,6)79-97(6)11)13-14-9-7(4,12(4,9)80-98(7)12)24(9,13,14,96(9)14)23(6,8,13,14)95(8)13)3(2,55)56(2,55,67)68(38,55,63,64,67,74)71(55,56,63,64,67,74,127(55,56,67,91-55)131(55,63,67,71)119(37,63)67)128(55,56,68,92-56)132(56,64,68,71)120(38,64)68/h8-9,13-14,19-22,27-28H",
            "InChI=1S/C12H37O3/c16-5-1-3-6(1,5)9(1,3,5,16,21-39(9)22-9,10(1,3,5,6,17-5,29(1)30(1)10,41(1,5)18-5,26-48(1,3,6,9)10)43(6)37(6)44(6,10)43)13(3)4-2-7(33(5)34(5)7)8(2,4,13)11(2,4,7,13,23-40(11)24-11,35(7)36(7)11,14(3,4,6,9,13,27(3)6)15(3,4,8,9,11,13)28(4)8)12(2,4,7,8,19-7,31(2)32(2)12,42(2,7)20-7,25-47(2,4,8)12)45(8)38(8)46(8,12)45/h1-4H",
            "InChI=1S/C18H44/c19-45-10-6-7(10,33(6)34(6)7,35(6)36(6)7)11(6,10)14(10,45)3-5(14,31(3)32(3)5)8(3,10,14)15(3,5,14,22-5,44(3,5,8)14,17(6,7,8,10,11,14,45)18(6,7,8,10,11,14,23-8,24-8,38(10)17,48(6,7,10,11,17)49(6,7,10,11,17)18)39(11)50(8,10,11,17,18)52(8,10,11,14,17,18)42(11,17)25-11)4(14,41(14,15)45)9(3,5,15)1-2(9,27(1)28(1)2,29(1)30(1)2)12(1,3,5,9,15)13(1,2,4,9,15,20-4,21-4)16(1,2,3,4,5,9,12,14,15,37(9)13,46(1,2,9,12)13)43(12,26-12)51(4,9,12,13,15,16)47(4,9,12,13,16)40(12)13/h1-8,14,19H",
            "InChI=1S/C32.C18H13PS2.C10H15.C2Cl2.2CF3O3S.2Ru/c1-6-9(1)18(1,6)2-7-8(2)15(7)3-5-4-16-10-11(16)12(10)19-13-14(19)17(6,13,19)21(13,14,19)20(9,18)22(15)23(16)24(3,4,5,22)27(3,5,15,22,23)26(2,7,8,18,20,22)25(1,6,9,18,20,21,32(7,8,15,20,22,24,26)27)31(13,14,17,19,21,23)29(10,11,12,16,23,24)28(4,5,16,22,23,24,27)30(10,11,12,19,21,23,29)31;20-17-12-6-4-10-15(17)19(14-8-2-1-3-9-14)16-11-5-7-13-18(16)21;1-6-7(2)9(4)10(5)8(6)3;3-1-2(3)4-1;2*2-1(3,4)8(5,6)7;;/h;1-13H;1-5H3;;;;;",
            "InChI=1S/C8H16Cl16/c25-1-3-5(1)7(1,3,25,29(3)5)11(1)9(1,19(1,3,5,7,11)17(3,5,7)21(1,3,5,7,19,33(3,5)17,37(5,17)39(3,5,17)21)23(1,3,5,7,9,11,17,19,27(1)19)31(1,7)25)15(3,5,7,11)13(3,5,7,35(3,7)15)14-4-2-6(4,14)8(2,4,14,30(4)6)12(2)10(2,16(4,6,8,12,14)36(4,8)14)20(2,4,6,8,12)18(4,6,8)22(2,4,6,8,20,34(4,6)18,38(6,18)40(4,6,18)22)24(2,4,6,8,10,12,18,20,28(2)20)32(2,8)26(2)8"
     
     I have resolved this by Shortest path based new fingerprinter.
     
     Here is the output.
     

Atoms       52	Bonds	194	FP Density	17	Elapsed time	713	Bitset	{18, 86, 98, 201, 262, 307, 337, 372, 379, 396, 492, 577, 587, 664, 671, 726, 769}

Atoms	16	Bonds	72	FP Density	10	Elapsed time	35	Bitset	{1, 105, 151, 160, 247, 328, 443, 550, 577, 587}

Atoms	68	Bonds	197	FP Density	339	Elapsed time	317	Bitset	{4, 5, 7, 11, 18, 19, 22, 29, 32, 33, 34, 35, 36, 37, 44, 47, 49, 51, 52, 55, 57, 59, 64, 66, 67, 68, 86, 92, 94, 96, 98, 100, 102, 103, 108, 110, 116, 119, 121, 122, 127, 128, 129, 135, 137, 139, 140, 144, 145, 148, 158, 163, 167, 168, 170, 171, 176, 180, 182, 183, 187, 192, 198, 201, 207, 212, 216, 220, 222, 223, 225, 230, 231, 234, 237, 238, 243, 244, 255, 256, 261, 262, 265, 269, 273, 274, 278, 282, 283, 286, 293, 296, 297, 299, 301, 302, 303, 305, 307, 311, 315, 317, 318, 321, 325, 326, 329, 333, 334, 335, 337, 338, 339, 341, 343, 344, 345, 349, 354, 355, 358, 372, 373, 374, 375, 377, 379, 381, 383, 384, 387, 389, 394, 396, 397, 400, 404, 407, 410, 411, 416, 418, 420, 421, 432, 437, 438, 439, 440, 441, 444, 445, 448, 454, 460, 466, 467, 468, 472, 477, 491, 492, 493, 495, 497, 509, 511, 512, 513, 514, 516, 518, 520, 522, 533, 536, 544, 545, 553, 554, 555, 563, 566, 570, 572, 577, 581, 584, 587, 589, 590, 591, 592, 594, 599, 600, 601, 614, 618, 623, 625, 627, 630, 631, 639, 646, 653, 654, 655, 657, 658, 659, 660, 662, 663, 664, 666, 669, 670, 671, 672, 677, 679, 681, 687, 690, 691, 694, 697, 700, 703, 710, 711, 713, 715, 716, 720, 726, 727, 728, 733, 736, 741, 745, 746, 747, 752, 753, 759, 767, 769, 773, 777, 787, 788, 799, 800, 801, 804, 806, 809, 811, 812, 817, 819, 821, 822, 824, 826, 829, 832, 834, 836, 839, 840, 842, 843, 849, 855, 856, 860, 863, 864, 872, 875, 876, 882, 883, 884, 885, 887, 891, 892, 894, 900, 904, 905, 907, 919, 924, 930, 932, 934, 935, 937, 940, 948, 952, 956, 957, 958, 959, 961, 963, 964, 965, 967, 968, 969, 970, 972, 974, 978, 979, 985, 990, 991, 994, 997, 1001, 1004, 1005, 1009, 1010, 1011, 1012, 1013, 1019, 1020}

Atoms	134	Bonds	402	FP Density	285	Elapsed time	773	Bitset	{6, 9, 10, 14, 16, 18, 19, 20, 21, 22, 24, 29, 32, 33, 36, 41, 44, 48, 49, 57, 59, 63, 72, 73, 78, 83, 85, 86, 87, 91, 92, 98, 99, 101, 110, 116, 118, 121, 122, 128, 131, 134, 136, 144, 161, 162, 175, 178, 182, 183, 190, 193, 199, 200, 201, 202, 205, 208, 209, 210, 214, 215, 229, 232, 233, 236, 237, 239, 244, 246, 247, 254, 256, 258, 261, 262, 272, 273, 277, 278, 280, 285, 288, 290, 294, 296, 298, 301, 302, 303, 304, 307, 310, 317, 321, 322, 323, 328, 337, 341, 343, 347, 348, 353, 354, 355, 358, 360, 369, 372, 378, 379, 380, 386, 387, 389, 390, 393, 396, 403, 404, 410, 412, 416, 421, 425, 432, 438, 440, 442, 444, 454, 465, 466, 468, 476, 478, 483, 486, 492, 493, 495, 496, 499, 502, 511, 513, 514, 527, 534, 540, 541, 544, 546, 547, 548, 555, 557, 558, 562, 565, 569, 572, 573, 576, 577, 578, 582, 584, 587, 592, 594, 599, 601, 614, 617, 621, 630, 631, 638, 640, 641, 645, 648, 649, 651, 655, 658, 664, 666, 669, 671, 677, 680, 690, 695, 696, 698, 699, 706, 710, 712, 715, 716, 723, 726, 729, 733, 734, 741, 743, 745, 754, 760, 761, 762, 763, 765, 769, 776, 782, 787, 791, 796, 808, 816, 820, 825, 826, 827, 830, 832, 839, 843, 844, 845, 849, 850, 854, 861, 864, 875, 877, 879, 882, 884, 892, 894, 895, 898, 903, 905, 906, 907, 911, 913, 916, 917, 918, 923, 930, 933, 937, 938, 943, 948, 951, 952, 955, 959, 961, 964, 966, 973, 976, 978, 983, 988, 989, 990, 993, 997, 1003, 1009, 1014}

Atoms	48	Bonds	127	FP Density	42	Elapsed time	46	Bitset	{44, 86, 142, 145, 238, 256, 259, 273, 305, 307, 337, 372, 390, 396, 420, 437, 444, 486, 492, 558, 577, 587, 636, 657, 664, 666, 674, 719, 720, 745, 760, 769, 775, 790, 795, 846, 925, 963, 964, 967, 983, 1009}

Atoms	52	Bonds	177	FP Density	24	Elapsed time	40	Bitset	{18, 44, 86, 98, 201, 202, 262, 307, 337, 372, 379, 396, 444, 492, 540, 558, 573, 577, 587, 664, 671, 726, 754, 769}

Atoms	85	Bonds	157	FP Density	82	Elapsed time	38	Bitset	{2, 8, 9, 17, 59, 98, 108, 117, 120, 137, 146, 148, 171, 184, 210, 239, 247, 255, 258, 271, 274, 280, 297, 299, 325, 371, 375, 379, 380, 387, 415, 416, 418, 433, 434, 450, 456, 470, 472, 492, 495, 511, 525, 531, 544, 550, 555, 558, 577, 579, 580, 587, 600, 621, 651, 657, 658, 664, 687, 708, 735, 743, 746, 752, 757, 781, 782, 791, 793, 799, 804, 816, 832, 834, 844, 854, 864, 941, 969, 976, 978, 979}

Atoms	40	Bonds	145	FP Density	71	Elapsed time	101	Bitset	{1, 6, 8, 21, 26, 28, 40, 46, 51, 69, 73, 86, 94, 97, 105, 120, 160, 182, 209, 247, 255, 267, 289, 299, 305, 307, 309, 313, 328, 337, 362, 368, 377, 393, 396, 416, 427, 438, 443, 461, 492, 506, 524, 550, 567, 573, 577, 579, 587, 594, 623, 656, 670, 678, 686, 688, 691, 751, 769, 784, 808, 822, 825, 843, 889, 902, 912, 975, 981, 982, 995}