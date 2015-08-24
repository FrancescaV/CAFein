//JACOBIAN R = 6X6 BRUTE FORCE. DO NOT ERASE
//
//double a00 = gsl_matrix_get(p -> matA_struct,0,0),
//a01 = gsl_matrix_get(p -> matA_struct,0,1),
//a02 = gsl_matrix_get(p -> matA_struct,0,2),
//a03 = gsl_matrix_get(p -> matA_struct,0,3),
//a04 = gsl_matrix_get(p -> matA_struct,0,4),
//a05 = gsl_matrix_get(p -> matA_struct,0,5),
//
//a10 = gsl_matrix_get(p -> matA_struct,1,0),
//a11 = gsl_matrix_get(p -> matA_struct,1,1),
//a12 = gsl_matrix_get(p -> matA_struct,1,2),
//a13 = gsl_matrix_get(p -> matA_struct,1,3),
//a14 = gsl_matrix_get(p -> matA_struct,1,4),
//a15 = gsl_matrix_get(p -> matA_struct,1,5),
//
//a20 = gsl_matrix_get(p -> matA_struct,2,0),
//a21 = gsl_matrix_get(p -> matA_struct,2,1),
//a22 = gsl_matrix_get(p -> matA_struct,2,2),
//a23 = gsl_matrix_get(p -> matA_struct,2,3),
//a24 = gsl_matrix_get(p -> matA_struct,2,4),
//a25 = gsl_matrix_get(p -> matA_struct,2,5),
//
//a30 = gsl_matrix_get(p -> matA_struct,3,0),
//a31 = gsl_matrix_get(p -> matA_struct,3,1),
//a32 = gsl_matrix_get(p -> matA_struct,3,2),
//a33 = gsl_matrix_get(p -> matA_struct,3,3),
//a34 = gsl_matrix_get(p -> matA_struct,3,4),
//a35 = gsl_matrix_get(p -> matA_struct,3,5),
//
//a40 = gsl_matrix_get(p -> matA_struct,4,0),
//a41 = gsl_matrix_get(p -> matA_struct,4,1),
//a42 = gsl_matrix_get(p -> matA_struct,4,2),
//a43 = gsl_matrix_get(p -> matA_struct,4,3),
//a44 = gsl_matrix_get(p -> matA_struct,4,4),
//a45 = gsl_matrix_get(p -> matA_struct,4,5),
//
//a50 = gsl_matrix_get(p -> matA_struct,5,0),
//a51 = gsl_matrix_get(p -> matA_struct,5,1),
//a52 = gsl_matrix_get(p -> matA_struct,5,2),
//a53 = gsl_matrix_get(p -> matA_struct,5,3),
//a54 = gsl_matrix_get(p -> matA_struct,5,4),
//a55 = gsl_matrix_get(p -> matA_struct,5,5),
//
//c00 = gsl_matrix_get(p -> matC_struct,0,0),
//c01 = gsl_matrix_get(p -> matC_struct,0,1),
//c02 = gsl_matrix_get(p -> matC_struct,0,2),
//c03 = gsl_matrix_get(p -> matC_struct,0,3),
//c04 = gsl_matrix_get(p -> matC_struct,0,4),
//c05 = gsl_matrix_get(p -> matC_struct,0,5),
//c10 = gsl_matrix_get(p -> matC_struct,1,0),
//c11 = gsl_matrix_get(p -> matC_struct,1,1),
//c12 = gsl_matrix_get(p -> matC_struct,1,2),
//c13 = gsl_matrix_get(p -> matC_struct,1,3),
//c14 = gsl_matrix_get(p -> matC_struct,1,4),
//c15 = gsl_matrix_get(p -> matC_struct,1,5),
//c20 = gsl_matrix_get(p -> matC_struct,2,0),
//c21 = gsl_matrix_get(p -> matC_struct,2,1),
//c22 = gsl_matrix_get(p -> matC_struct,2,2),
//c23 = gsl_matrix_get(p -> matC_struct,2,3),
//c24 = gsl_matrix_get(p -> matC_struct,2,4),
//c25 = gsl_matrix_get(p -> matC_struct,2,5),
//c30 = gsl_matrix_get(p -> matC_struct,3,0),
//c31 = gsl_matrix_get(p -> matC_struct,3,1),
//c32 = gsl_matrix_get(p -> matC_struct,3,2),
//c33 = gsl_matrix_get(p -> matC_struct,3,3),
//c34 = gsl_matrix_get(p -> matC_struct,3,4),
//c35 = gsl_matrix_get(p -> matC_struct,3,5),
//c40 = gsl_matrix_get(p -> matC_struct,4,0),
//c41 = gsl_matrix_get(p -> matC_struct,4,1),
//c42 = gsl_matrix_get(p -> matC_struct,4,2),
//c43 = gsl_matrix_get(p -> matC_struct,4,3),
//c44 = gsl_matrix_get(p -> matC_struct,4,4),
//c45 = gsl_matrix_get(p -> matC_struct,4,5),
//c50 = gsl_matrix_get(p -> matC_struct,5,0),
//c51 = gsl_matrix_get(p -> matC_struct,5,1),
//c52 = gsl_matrix_get(p -> matC_struct,5,2),
//c53 = gsl_matrix_get(p -> matC_struct,5,3),
//c54 = gsl_matrix_get(p -> matC_struct,5,4),
//c55 = gsl_matrix_get(p -> matC_struct,5,5),
//
//
//d00 = gsl_matrix_get(p -> matD_struct,0,0),
//d01 = gsl_matrix_get(p -> matD_struct,0,1),
//d02 = gsl_matrix_get(p -> matD_struct,0,2),
//d03 = gsl_matrix_get(p -> matD_struct,0,3),
//d04 = gsl_matrix_get(p -> matD_struct,0,4),
//d05 = gsl_matrix_get(p -> matD_struct,0,5),
//d10 = gsl_matrix_get(p -> matD_struct,1,0),
//d11 = gsl_matrix_get(p -> matD_struct,1,1),
//d12 = gsl_matrix_get(p -> matD_struct,1,2),
//d13 = gsl_matrix_get(p -> matD_struct,1,3),
//d14 = gsl_matrix_get(p -> matD_struct,1,4),
//d15 = gsl_matrix_get(p -> matD_struct,1,5),
//d20 = gsl_matrix_get(p -> matD_struct,2,0),
//d21 = gsl_matrix_get(p -> matD_struct,2,1),
//d22 = gsl_matrix_get(p -> matD_struct,2,2),
//d23 = gsl_matrix_get(p -> matD_struct,2,3),
//d24 = gsl_matrix_get(p -> matD_struct,2,4),
//d25 = gsl_matrix_get(p -> matD_struct,2,5),
//d30 = gsl_matrix_get(p -> matD_struct,3,0),
//d31 = gsl_matrix_get(p -> matD_struct,3,1),
//d32 = gsl_matrix_get(p -> matD_struct,3,2),
//d33 = gsl_matrix_get(p -> matD_struct,3,3),
//d34 = gsl_matrix_get(p -> matD_struct,3,4),
//d35 = gsl_matrix_get(p -> matD_struct,3,5),
//d40 = gsl_matrix_get(p -> matD_struct,4,0),
//d41 = gsl_matrix_get(p -> matD_struct,4,1),
//d42 = gsl_matrix_get(p -> matD_struct,4,2),
//d43 = gsl_matrix_get(p -> matD_struct,4,3),
//d44 = gsl_matrix_get(p -> matD_struct,4,4),
//d45 = gsl_matrix_get(p -> matD_struct,4,5),
//d50 = gsl_matrix_get(p -> matD_struct,5,0),
//d51 = gsl_matrix_get(p -> matD_struct,5,1),
//d52 = gsl_matrix_get(p -> matD_struct,5,2),
//d53 = gsl_matrix_get(p -> matD_struct,5,3),
//d54 = gsl_matrix_get(p -> matD_struct,5,4),
//d55 = gsl_matrix_get(p -> matD_struct,5,5),
//
//r00 = y[0],
//r01 = y[1],
//r02 = y[2],
//r03 = y[3],
//r04 = y[4],
//r05 = y[5],
//r10 = y[6],
//r11 = y[7],
//r12 = y[8],
//r13 = y[9],
//r14 = y[10],
//r15 = y[11],
//r20 = y[12],
//r21 = y[13],
//r22 = y[14],
//r23 = y[15],
//r24 = y[16],
//r25 = y[17],
//r30 = y[18],
//r31 = y[19],
//r32 = y[20],
//r33 = y[21],
//r34 = y[22],
//r35 = y[23],
//r40 = y[24],
//r41 = y[25],
//r42 = y[26],
//r43 = y[27],
//r44 = y[28],
//r45 = y[29],
//r50 = y[30],
//r51 = y[31],
//r52 = y[32],
//r53 = y[33],
//r54 = y[34],
//r55 = y[35],
//v0 = y[36],
//v1 = y[37],
//v2 = y[38],
//v3 = y[39],
//v4 = y[40],
//v5 = y[41];
//
//gsl_matrix_set(m, 0, 0, a00 - d00 - 2*c00*r00 - c10*r01 - c20*r02 - c30*r03 - c40*r04 - c50*r05 - c01*r10 - c02*r20 - c03*r30 - c04*r40 - c05*r50);
//gsl_matrix_set(m, 0, 1, -d10 - c10*r00 - c11*r10 - c12*r20 - c13*r30 - c14*r40 - c15*r50);
//gsl_matrix_set(m, 0, 2,-d20 - c20*r00 - c21*r10 - c22*r20 - c23*r30 - c24*r40 - c25*r50);
//gsl_matrix_set(m, 0, 3,-d30 - c30*r00 - c31*r10 - c32*r20 - c33*r30 - c34*r40 - c35*r50);
//gsl_matrix_set(m, 0, 4,-d40 - c40*r00 - c41*r10 - c42*r20 - c43*r30 - c44*r40 - c45*r50);
//gsl_matrix_set(m, 0, 5,-d50 - c50*r00 - c51*r10 - c52*r20 - c53*r30 - c54*r40 - c55*r50);
//gsl_matrix_set(m, 0, 6, a01 - c01*r00 - c11*r01 - c21*r02 - c31*r03 - c41*r04 - c51*r05);
//gsl_matrix_set(m, 0, 12, a02 - c02*r00 - c12*r01 - c22*r02 - c32*r03 - c42*r04 - c52*r05);
//gsl_matrix_set(m, 0, 18, a03 - c03*r00 - c13*r01 - c23*r02 - c33*r03 - c43*r04 - c53*r05);
//gsl_matrix_set(m, 0, 24, a04 - c04*r00 - c14*r01 - c24*r02 - c34*r03 - c44*r04 - c54*r05);
//gsl_matrix_set(m, 0, 30, a05 - c05*r00 - c15*r01 - c25*r02 - c35*r03 - c45*r04 - c55*r05);
//
//gsl_matrix_set(m, 1, 0, -d01 - c00*r01 - c01*r11 - c02*r21 - c03*r31 - c04*r41 - c05*r51);
//gsl_matrix_set(m, 1, 1, a00 - d11 - c00*r00 - 2*c10*r01 - c20*r02 - c30*r03 - c40*r04 - c50*r05 - c11*r11 - c12*r21 - c13*r31 - c14*r41 - c15*r51);
//gsl_matrix_set(m, 1, 2, -d21 - c20*r01 - c21*r11 - c22*r21 - c23*r31 - c24*r41 - c25*r51);
//gsl_matrix_set(m, 1, 3, -d31 - c30*r01 - c31*r11 - c32*r21 - c33*r31 - c34*r41 - c35*r51);
//gsl_matrix_set(m, 1, 4, -d41 - c40*r01 - c41*r11 - c42*r21 - c43*r31 - c44*r41 - c45*r51);
//gsl_matrix_set(m, 1, 5, -d51 - c50*r01 - c51*r11 - c52*r21 - c53*r31 - c54*r41 - c55*r51);
//
//gsl_matrix_set(m, 1, 7, a01 - c01*r00 - c11*r01 - c21*r02 - c31*r03 - c41*r04 - c51*r05);
//gsl_matrix_set(m, 1, 13, a02 - c02*r00 - c12*r01 - c22*r02 - c32*r03 - c42*r04 - c52*r05);
//gsl_matrix_set(m, 1, 19, a03 - c03*r00 - c13*r01 - c23*r02 - c33*r03 - c43*r04 - c53*r05);
//gsl_matrix_set(m, 1, 25, a04 - c04*r00 - c14*r01 - c24*r02 - c34*r03 - c44*r04 - c54*r05);
//gsl_matrix_set(m, 1, 31, a05 - c05*r00 - c15*r01 - c25*r02 - c35*r03 - c45*r04 - c55*r05);
//
//
//gsl_matrix_set(m, 2, 0, -d02 - c00*r02 - c01*r12 - c02*r22 - c03*r32 - c04*r42 - c05*r52);
//gsl_matrix_set(m, 2, 1, -d12 - c10*r02 - c11*r12 - c12*r22 - c13*r32 - c14*r42 - c15*r52);
//gsl_matrix_set(m, 2, 2, a00 - d22 - c00*r00 - c10*r01 - 2*c20*r02 - c30*r03 - c40*r04 - c50*r05 - c21*r12 - c22*r22 - c23*r32 - c24*r42 - c25*r52);
//gsl_matrix_set(m, 2, 3, -d32 - c30*r02 - c31*r12 - c32*r22 - c33*r32 - c34*r42 - c35*r52);
//gsl_matrix_set(m, 2, 4, -d42 - c40*r02 - c41*r12 - c42*r22 - c43*r32 - c44*r42 - c45*r52);
//gsl_matrix_set(m, 2, 5, -d52 - c50*r02 - c51*r12 - c52*r22 - c53*r32 - c54*r42 - c55*r52);
//
//gsl_matrix_set(m, 2, 8, a01 - c01*r00 - c11*r01 - c21*r02 - c31*r03 - c41*r04 - c51*r05);
//gsl_matrix_set(m, 2, 14, a02 - c02*r00 - c12*r01 - c22*r02 - c32*r03 - c42*r04 - c52*r05);
//gsl_matrix_set(m, 2, 20, a03 - c03*r00 - c13*r01 - c23*r02 - c33*r03 - c43*r04 - c53*r05);
//gsl_matrix_set(m, 2, 26, a04 - c04*r00 - c14*r01 - c24*r02 - c34*r03 - c44*r04 - c54*r05);
//gsl_matrix_set(m, 2, 32, a05 - c05*r00 - c15*r01 - c25*r02 - c35*r03 - c45*r04 - c55*r05);
//
//
//gsl_matrix_set(m, 3, 0, -d03 - c00*r03 - c01*r13 - c02*r23 - c03*r33 - c04*r43 - c05*r53);
//gsl_matrix_set(m, 3, 1, -d13 - c10*r03 - c11*r13 - c12*r23 - c13*r33 - c14*r43 - c15*r53);
//gsl_matrix_set(m, 3, 2, -d23 - c20*r03 - c21*r13 - c22*r23 - c23*r33 - c24*r43 - c25*r53);
//gsl_matrix_set(m, 3, 3, a00 - d33 - c00*r00 - c10*r01 - c20*r02 - 2*c30*r03 - c40*r04 - c50*r05 - c31*r13 - c32*r23 - c33*r33 - c34*r43 - c35*r53);
//gsl_matrix_set(m, 3, 4, -d43 - c40*r03 - c41*r13 - c42*r23 - c43*r33 - c44*r43 - c45*r53);
//gsl_matrix_set(m, 3, 5, -d53 - c50*r03 - c51*r13 - c52*r23 - c53*r33 - c54*r43 - c55*r53);
//
//gsl_matrix_set(m, 3, 9, a01 - c01*r00 - c11*r01 - c21*r02 - c31*r03 - c41*r04 - c51*r05);
//gsl_matrix_set(m, 3, 15, a02 - c02*r00 - c12*r01 - c22*r02 - c32*r03 - c42*r04 - c52*r05);
//gsl_matrix_set(m, 3, 21, a03 - c03*r00 - c13*r01 - c23*r02 - c33*r03 - c43*r04 - c53*r05);
//gsl_matrix_set(m, 3, 27, a04 - c04*r00 - c14*r01 - c24*r02 - c34*r03 - c44*r04 - c54*r05);
//gsl_matrix_set(m, 3, 33, a05 - c05*r00 - c15*r01 - c25*r02 - c35*r03 - c45*r04 - c55*r05);
//
//
//gsl_matrix_set(m, 4, 0, -d04 - c00*r04 - c01*r14 - c02*r24 - c03*r34 - c04*r44 - c05*r54);
//gsl_matrix_set(m, 4, 1, -d14 - c10*r04 - c11*r14 - c12*r24 - c13*r34 - c14*r44 - c15*r54);
//gsl_matrix_set(m, 4, 2, -d24 - c20*r04 - c21*r14 - c22*r24 - c23*r34 - c24*r44 - c25*r54);
//gsl_matrix_set(m, 4, 3, -d34 - c30*r04 - c31*r14 - c32*r24 - c33*r34 - c34*r44 - c35*r54);
//gsl_matrix_set(m, 4, 4, a00 - d44 - c00*r00 - c10*r01 - c20*r02 - c30*r03 - 2*c40*r04 - c50*r05 - c41*r14 - c42*r24 - c43*r34 - c44*r44 - c45*r54);
//gsl_matrix_set(m, 4, 5, -d54 - c50*r04 - c51*r14 - c52*r24 - c53*r34 - c54*r44 - c55*r54);
//
//gsl_matrix_set(m, 4, 10, a01 - c01*r00 - c11*r01 - c21*r02 - c31*r03 - c41*r04 - c51*r05);
//gsl_matrix_set(m, 4, 16, a02 - c02*r00 - c12*r01 - c22*r02 - c32*r03 - c42*r04 - c52*r05);
//gsl_matrix_set(m, 4, 22, a03 - c03*r00 - c13*r01 - c23*r02 - c33*r03 - c43*r04 - c53*r05);
//gsl_matrix_set(m, 4, 28, a04 - c04*r00 - c14*r01 - c24*r02 - c34*r03 - c44*r04 - c54*r05);
//gsl_matrix_set(m, 4, 34, a05 - c05*r00 - c15*r01 - c25*r02 - c35*r03 - c45*r04 - c55*r05);
//
//gsl_matrix_set(m, 5, 0, -d05 - c00*r05 - c01*r15 - c02*r25 - c03*r35 - c04*r45 - c05*r55);
//gsl_matrix_set(m, 5, 1, -d15 - c10*r05 - c11*r15 - c12*r25 - c13*r35 - c14*r45 - c15*r55);
//gsl_matrix_set(m, 5, 2, -d25 - c20*r05 - c21*r15 - c22*r25 - c23*r35 - c24*r45 - c25*r55);
//gsl_matrix_set(m, 5, 3, -d35 - c30*r05 - c31*r15 - c32*r25 - c33*r35 - c34*r45 - c35*r55);
//gsl_matrix_set(m, 5, 4, -d45 - c40*r05 - c41*r15 - c42*r25 - c43*r35 - c44*r45 - c45*r55);
//gsl_matrix_set(m, 5, 5, a00 - d55 - c00*r00 - c10*r01 - c20*r02 - c30*r03 - c40*r04 - 2*c50*r05 - c51*r15 - c52*r25 - c53*r35 - c54*r45 - c55*r55);
//
//gsl_matrix_set(m, 5, 11, a01 - c01*r00 - c11*r01 - c21*r02 - c31*r03 - c41*r04 - c51*r05);
//gsl_matrix_set(m, 5, 17, a02 - c02*r00 - c12*r01 - c22*r02 - c32*r03 - c42*r04 - c52*r05);
//gsl_matrix_set(m, 5, 23, a03 - c03*r00 - c13*r01 - c23*r02 - c33*r03 - c43*r04 - c53*r05);
//gsl_matrix_set(m, 5, 29, a04 - c04*r00 - c14*r01 - c24*r02 - c34*r03 - c44*r04 - c54*r05);
//gsl_matrix_set(m, 5, 35, a05 - c05*r00 - c15*r01 - c25*r02 - c35*r03 - c45*r04 - c55*r05);
//
//
//gsl_matrix_set(m, 6, 0, a10 - c00*r10 - c10*r11 - c20*r12 - c30*r13 - c40*r14 - c50*r15);
//gsl_matrix_set(m, 6, 6, a11 - d00 - c00*r00 - 2*c01*r10 - c11*r11 - c21*r12 - c31*r13 - c41*r14 - c51*r15 - c02*r20 - c03*r30 - c04*r40 - c05*r50);
//gsl_matrix_set(m, 6, 7, -d10 - c10*r00 - c11*r10 - c12*r20 - c13*r30 - c14*r40 - c15*r50);
//gsl_matrix_set(m, 6, 8, -d20 - c20*r00 - c21*r10 - c22*r20 - c23*r30 - c24*r40 - c25*r50);
//gsl_matrix_set(m, 6, 9, -d30 - c30*r00 - c31*r10 - c32*r20 - c33*r30 - c34*r40 - c35*r50);
//gsl_matrix_set(m, 6, 10, -d40 - c40*r00 - c41*r10 - c42*r20 - c43*r30 - c44*r40 - c45*r50);
//gsl_matrix_set(m, 6, 11, -d50 - c50*r00 - c51*r10 - c52*r20 - c53*r30 - c54*r40 - c55*r50);
//
//gsl_matrix_set(m, 6, 12, a12 - c02*r10 - c12*r11 - c22*r12 - c32*r13 - c42*r14 - c52*r15);
//gsl_matrix_set(m, 6, 18, a13 - c03*r10 - c13*r11 - c23*r12 - c33*r13 - c43*r14 - c53*r15);
//gsl_matrix_set(m, 6, 24, a14 - c04*r10 - c14*r11 - c24*r12 - c34*r13 - c44*r14 - c54*r15);
//gsl_matrix_set(m, 6, 30, a15 - c05*r10 - c15*r11 - c25*r12 - c35*r13 - c45*r14 - c55*r15);
//
//gsl_matrix_set(m, 7, 1, a10 - c00*r10 - c10*r11 - c20*r12 - c30*r13 - c40*r14 - c50*r15);
//gsl_matrix_set(m, 7, 6, -d01 - c00*r01 - c01*r11 - c02*r21 - c03*r31 - c04*r41 - c05*r51);
//gsl_matrix_set(m, 7, 7, a11 - d11 - c10*r01 - c01*r10 - 2*c11*r11 - c21*r12 - c31*r13 - c41*r14 - c51*r15 - c12*r21 - c13*r31 - c14*r41 - c15*r51);
//gsl_matrix_set(m, 7, 8, -d21 - c20*r01 - c21*r11 - c22*r21 - c23*r31 - c24*r41 - c25*r51);
//gsl_matrix_set(m, 7, 9, -d31 - c30*r01 - c31*r11 - c32*r21 - c33*r31 - c34*r41 - c35*r51);
//gsl_matrix_set(m, 7, 10, -d41 - c40*r01 - c41*r11 - c42*r21 - c43*r31 - c44*r41 - c45*r51);
//gsl_matrix_set(m, 7, 11, -d51 - c50*r01 - c51*r11 - c52*r21 - c53*r31 - c54*r41 - c55*r51);
//gsl_matrix_set(m, 7, 13, a12 - c02*r10 - c12*r11 - c22*r12 - c32*r13 - c42*r14 - c52*r15);
//gsl_matrix_set(m, 7, 19, a13 - c03*r10 - c13*r11 - c23*r12 - c33*r13 - c43*r14 - c53*r15);
//gsl_matrix_set(m, 7, 25, a14 - c04*r10 - c14*r11 - c24*r12 - c34*r13 - c44*r14 - c54*r15);
//gsl_matrix_set(m, 7, 31, a15 - c05*r10 - c15*r11 - c25*r12 - c35*r13 - c45*r14 - c55*r15);
//
//
//
//gsl_matrix_set(m, 8, 2, a10 - c00*r10 - c10*r11 - c20*r12 - c30*r13 - c40*r14 - c50*r15);
//gsl_matrix_set(m, 8, 6, -d02 - c00*r02 - c01*r12 - c02*r22 - c03*r32 - c04*r42 - c05*r52);
//gsl_matrix_set(m, 8, 7, -d12 - c10*r02 - c11*r12 - c12*r22 - c13*r32 - c14*r42 - c15*r52);
//gsl_matrix_set(m, 8, 8, a11 - d22 - c20*r02 - c01*r10 - c11*r11 - 2*c21*r12 - c31*r13 - c41*r14 - c51*r15 - c22*r22 - c23*r32 - c24*r42 - c25*r52);
//gsl_matrix_set(m, 8, 9, -d32 - c30*r02 - c31*r12 - c32*r22 - c33*r32 - c34*r42 - c35*r52);
//gsl_matrix_set(m, 8, 10, -d42 - c40*r02 - c41*r12 - c42*r22 - c43*r32 - c44*r42 - c45*r52);
//gsl_matrix_set(m, 8, 11, -d52 - c50*r02 - c51*r12 - c52*r22 - c53*r32 - c54*r42 - c55*r52);
//gsl_matrix_set(m, 8, 14, a12 - c02*r10 - c12*r11 - c22*r12 - c32*r13 - c42*r14 - c52*r15);
//gsl_matrix_set(m, 8, 20, a13 - c03*r10 - c13*r11 - c23*r12 - c33*r13 - c43*r14 - c53*r15);
//gsl_matrix_set(m, 8, 26, a14 - c04*r10 - c14*r11 - c24*r12 - c34*r13 - c44*r14 - c54*r15);
//gsl_matrix_set(m, 8, 32, a15 - c05*r10 - c15*r11 - c25*r12 - c35*r13 - c45*r14 - c55*r15);
//
//
//gsl_matrix_set(m, 9, 3, a10 - c00*r10 - c10*r11 - c20*r12 - c30*r13 - c40*r14 - c50*r15);
//gsl_matrix_set(m, 9, 6, -d03 - c00*r03 - c01*r13 - c02*r23 - c03*r33 - c04*r43 - c05*r53);
//gsl_matrix_set(m, 9, 7, -d13 - c10*r03 - c11*r13 - c12*r23 - c13*r33 - c14*r43 - c15*r53);
//gsl_matrix_set(m, 9, 8, -d23 - c20*r03 - c21*r13 - c22*r23 - c23*r33 - c24*r43 - c25*r53);
//gsl_matrix_set(m, 9, 9, a11 - d33 - c30*r03 - c01*r10 - c11*r11 - c21*r12 - 2*c31*r13 - c41*r14 - c51*r15 - c32*r23 - c33*r33 - c34*r43 - c35*r53);
//gsl_matrix_set(m, 9, 10, -d43 - c40*r03 - c41*r13 - c42*r23 - c43*r33 - c44*r43 - c45*r53);
//gsl_matrix_set(m, 9, 11, -d53 - c50*r03 - c51*r13 - c52*r23 - c53*r33 - c54*r43 - c55*r53);
//gsl_matrix_set(m, 9, 15, a12 - c02*r10 - c12*r11 - c22*r12 - c32*r13 - c42*r14 - c52*r15);
//gsl_matrix_set(m, 9, 21, a13 - c03*r10 - c13*r11 - c23*r12 - c33*r13 - c43*r14 - c53*r15);
//gsl_matrix_set(m, 9, 27, a14 - c04*r10 - c14*r11 - c24*r12 - c34*r13 - c44*r14 - c54*r15);
//gsl_matrix_set(m, 9, 33, a15 - c05*r10 - c15*r11 - c25*r12 - c35*r13 - c45*r14 - c55*r15);
//
//
//gsl_matrix_set(m, 10, 4, a10 - c00*r10 - c10*r11 - c20*r12 - c30*r13 - c40*r14 - c50*r15);
//gsl_matrix_set(m, 10, 6, -d04 - c00*r04 - c01*r14 - c02*r24 - c03*r34 - c04*r44 - c05*r54);
//gsl_matrix_set(m, 10, 7, -d14 - c10*r04 - c11*r14 - c12*r24 - c13*r34 - c14*r44 - c15*r54);
//gsl_matrix_set(m, 10, 8, -d24 - c20*r04 - c21*r14 - c22*r24 - c23*r34 - c24*r44 - c25*r54);
//gsl_matrix_set(m, 10, 9, -d34 - c30*r04 - c31*r14 - c32*r24 - c33*r34 - c34*r44 - c35*r54);
//gsl_matrix_set(m, 10, 10, a11 - d44 - c40*r04 - c01*r10 - c11*r11 - c21*r12 - c31*r13 - 2*c41*r14 - c51*r15 - c42*r24 - c43*r34 - c44*r44 - c45*r54);
//gsl_matrix_set(m, 10, 11, -d54 - c50*r04 - c51*r14 - c52*r24 - c53*r34 - c54*r44 - c55*r54);
//gsl_matrix_set(m, 10, 16, a12 - c02*r10 - c12*r11 - c22*r12 - c32*r13 - c42*r14 - c52*r15);
//gsl_matrix_set(m, 10, 22, a13 - c03*r10 - c13*r11 - c23*r12 - c33*r13 - c43*r14 - c53*r15);
//gsl_matrix_set(m, 10, 28, a14 - c04*r10 - c14*r11 - c24*r12 - c34*r13 - c44*r14 - c54*r15);
//gsl_matrix_set(m, 10, 34, a15 - c05*r10 - c15*r11 - c25*r12 - c35*r13 - c45*r14 - c55*r15);
//
//
//
//
//gsl_matrix_set(m, 11, 5, a10 - c00*r10 - c10*r11 - c20*r12 - c30*r13 - c40*r14 - c50*r15);
//
//
//gsl_matrix_set(m, 11, 6, -d05 - c00*r05 - c01*r15 - c02*r25 - c03*r35 - c04*r45 - c05*r55);
//gsl_matrix_set(m, 11, 7, -d15 - c10*r05 - c11*r15 - c12*r25 - c13*r35 - c14*r45 - c15*r55);
//gsl_matrix_set(m, 11, 8, -d25 - c20*r05 - c21*r15 - c22*r25 - c23*r35 - c24*r45 - c25*r55);
//gsl_matrix_set(m, 11, 9, -d35 - c30*r05 - c31*r15 - c32*r25 - c33*r35 - c34*r45 - c35*r55);
//gsl_matrix_set(m, 11, 10, -d45 - c40*r05 - c41*r15 - c42*r25 - c43*r35 - c44*r45 - c45*r55);
//gsl_matrix_set(m, 11, 11, a11 - d55 - c50*r05 - c01*r10 - c11*r11 - c21*r12 - c31*r13 - c41*r14 - 2*c51*r15 - c52*r25 - c53*r35 - c54*r45 - c55*r55);
//gsl_matrix_set(m, 11, 17, a12 - c02*r10 - c12*r11 - c22*r12 - c32*r13 - c42*r14 - c52*r15);
//gsl_matrix_set(m, 11, 23, a13 - c03*r10 - c13*r11 - c23*r12 - c33*r13 - c43*r14 - c53*r15);
//gsl_matrix_set(m, 11, 29, a14 - c04*r10 - c14*r11 - c24*r12 - c34*r13 - c44*r14 - c54*r15);
//gsl_matrix_set(m, 11, 35, a15 - c05*r10 - c15*r11 - c25*r12 - c35*r13 - c45*r14 - c55*r15);
//
//
//gsl_matrix_set(m, 12, 0, a20 - c00*r20 - c10*r21 - c20*r22 - c30*r23 - c40*r24 - c50*r25);
//gsl_matrix_set(m, 12, 6, a21 - c01*r20 - c11*r21 - c21*r22 - c31*r23 - c41*r24 - c51*r25);
//gsl_matrix_set(m, 12, 12, a22 - d00 - c00*r00 - c01*r10 - 2*c02*r20 - c12*r21 - c22*r22 - c32*r23 - c42*r24 - c52*r25 - c03*r30 - c04*r40 - c05*r50);
//gsl_matrix_set(m, 12, 13, -d10 - c10*r00 - c11*r10 - c12*r20 - c13*r30 - c14*r40 - c15*r50);
//gsl_matrix_set(m, 12, 14, -d20 - c20*r00 - c21*r10 - c22*r20 - c23*r30 - c24*r40 - c25*r50);
//gsl_matrix_set(m, 12, 15, -d30 - c30*r00 - c31*r10 - c32*r20 - c33*r30 - c34*r40 - c35*r50);
//gsl_matrix_set(m, 12, 16, -d40 - c40*r00 - c41*r10 - c42*r20 - c43*r30 - c44*r40 - c45*r50);
//gsl_matrix_set(m, 12, 17, -d50 - c50*r00 - c51*r10 - c52*r20 - c53*r30 - c54*r40 - c55*r50);
//
//gsl_matrix_set(m, 12, 18, a23 - c03*r20 - c13*r21 - c23*r22 - c33*r23 - c43*r24 - c53*r25);
//gsl_matrix_set(m, 12, 24, a24 - c04*r20 - c14*r21 - c24*r22 - c34*r23 - c44*r24 - c54*r25);
//gsl_matrix_set(m, 12, 30, a25 - c05*r20 - c15*r21 - c25*r22 - c35*r23 - c45*r24 - c55*r25);
//
//
//
//gsl_matrix_set(m, 13, 1, a20 - c00*r20 - c10*r21 - c20*r22 - c30*r23 - c40*r24 - c50*r25);
//gsl_matrix_set(m, 13, 7, a21 - c01*r20 - c11*r21 - c21*r22 - c31*r23 - c41*r24 - c51*r25);
//gsl_matrix_set(m, 13, 12, -d01 - c00*r01 - c01*r11 - c02*r21 - c03*r31 - c04*r41 - c05*r51);
//gsl_matrix_set(m, 13, 13, a22 - d11 - c10*r01 - c11*r11 - c02*r20 - 2*c12*r21 - c22*r22 - c32*r23 - c42*r24 - c52*r25 - c13*r31 - c14*r41 - c15*r51);
//gsl_matrix_set(m, 13, 14, -d21 - c20*r01 - c21*r11 - c22*r21 - c23*r31 - c24*r41 - c25*r51);
//gsl_matrix_set(m, 13, 15, -d31 - c30*r01 - c31*r11 - c32*r21 - c33*r31 - c34*r41 - c35*r51);
//gsl_matrix_set(m, 13, 16, -d41 - c40*r01 - c41*r11 - c42*r21 - c43*r31 - c44*r41 - c45*r51);
//gsl_matrix_set(m, 13, 17, -d51 - c50*r01 - c51*r11 - c52*r21 - c53*r31 - c54*r41 - c55*r51);
//gsl_matrix_set(m, 13, 19, a23 - c03*r20 - c13*r21 - c23*r22 - c33*r23 - c43*r24 - c53*r25);
//gsl_matrix_set(m, 13, 25, a24 - c04*r20 - c14*r21 - c24*r22 - c34*r23 - c44*r24 - c54*r25);
//gsl_matrix_set(m, 13, 31, a25 - c05*r20 - c15*r21 - c25*r22 - c35*r23 - c45*r24 - c55*r25);
//
//
//
//gsl_matrix_set(m, 14, 0, 0.0);
//gsl_matrix_set(m, 14, 1, 0.0);
//gsl_matrix_set(m, 14, 2, a20 - c00*r20 - c10*r21 - c20*r22 - c30*r23 - c40*r24 - c50*r25);
//gsl_matrix_set(m, 14, 3, 0.0);
//gsl_matrix_set(m, 14, 4, 0.0);
//gsl_matrix_set(m, 14, 5, 0.0);
//
//gsl_matrix_set(m, 14, 6, 0.0);
//gsl_matrix_set(m, 14, 7, 0.0);
//gsl_matrix_set(m, 14, 8, a21 - c01*r20 - c11*r21 - c21*r22 - c31*r23 - c41*r24 - c51*r25);
//gsl_matrix_set(m, 14, 9, 0.0);
//gsl_matrix_set(m, 14, 10, 0.0);
//gsl_matrix_set(m, 14, 11, 0.0);
//
//gsl_matrix_set(m, 14, 12, -d02 - c00*r02 - c01*r12 - c02*r22 - c03*r32 - c04*r42 - c05*r52);
//gsl_matrix_set(m, 14, 13, -d12 - c10*r02 - c11*r12 - c12*r22 - c13*r32 - c14*r42 - c15*r52);
//gsl_matrix_set(m, 14, 14, a22 - d22 - c20*r02 - c21*r12 - c02*r20 - c12*r21 - 2*c22*r22 - c32*r23 - c42*r24 - c52*r25 - c23*r32 - c24*r42 - c25*r52);
//gsl_matrix_set(m, 14, 15, -d32 - c30*r02 - c31*r12 - c32*r22 - c33*r32 - c34*r42 - c35*r52);
//gsl_matrix_set(m, 14, 16, -d42 - c40*r02 - c41*r12 - c42*r22 - c43*r32 - c44*r42 - c45*r52);
//gsl_matrix_set(m, 14, 17, -d52 - c50*r02 - c51*r12 - c52*r22 - c53*r32 - c54*r42 - c55*r52);
//
//gsl_matrix_set(m, 14, 18, 0.0);
//gsl_matrix_set(m, 14, 19, 0.0);
//gsl_matrix_set(m, 14, 20, a23 - c03*r20 - c13*r21 - c23*r22 - c33*r23 - c43*r24 - c53*r25);
//gsl_matrix_set(m, 14, 21, 0.0);
//gsl_matrix_set(m, 14, 22, 0.0);
//gsl_matrix_set(m, 14, 23, 0.0);
//
//gsl_matrix_set(m, 14, 24, 0.0);
//gsl_matrix_set(m, 14, 25, 0.0);
//gsl_matrix_set(m, 14, 26, a24 - c04*r20 - c14*r21 - c24*r22 - c34*r23 - c44*r24 - c54*r25);
//gsl_matrix_set(m, 14, 27, 0.0);
//gsl_matrix_set(m, 14, 28, 0.0);
//gsl_matrix_set(m, 14, 29, 0.0);
//
//gsl_matrix_set(m, 14, 30, 0.0);
//gsl_matrix_set(m, 14, 31, 0.0);
//gsl_matrix_set(m, 14, 32, a25 - c05*r20 - c15*r21 - c25*r22 - c35*r23 - c45*r24 - c55*r25);
//gsl_matrix_set(m, 14, 33, 0.0);
//gsl_matrix_set(m, 14, 34, 0.0);
//gsl_matrix_set(m, 14, 35, 0.0);
//
//
//
//
//gsl_matrix_set(m, 15, 0, 0.0);
//gsl_matrix_set(m, 15, 1, 0.0);
//gsl_matrix_set(m, 15, 2, 0.0);
//gsl_matrix_set(m, 15, 3, a20 - c00*r20 - c10*r21 - c20*r22 - c30*r23 - c40*r24 - c50*r25);
//gsl_matrix_set(m, 15, 4, 0.0);
//gsl_matrix_set(m, 15, 5, 0.0);
//
//gsl_matrix_set(m, 15, 6, 0.0);
//gsl_matrix_set(m, 15, 7, 0.0);
//gsl_matrix_set(m, 15, 8, 0.0);
//gsl_matrix_set(m, 15, 9, a21 - c01*r20 - c11*r21 - c21*r22 - c31*r23 - c41*r24 - c51*r25);
//gsl_matrix_set(m, 15, 10, 0.0);
//gsl_matrix_set(m, 15, 11, 0.0);
//
//gsl_matrix_set(m, 15, 12, -d03 - c00*r03 - c01*r13 - c02*r23 - c03*r33 - c04*r43 - c05*r53);
//gsl_matrix_set(m, 15, 13, -d13 - c10*r03 - c11*r13 - c12*r23 - c13*r33 - c14*r43 - c15*r53);
//gsl_matrix_set(m, 15, 14, -d23 - c20*r03 - c21*r13 - c22*r23 - c23*r33 - c24*r43 - c25*r53);
//gsl_matrix_set(m, 15, 15, a22 - d33 - c30*r03 - c31*r13 - c02*r20 - c12*r21 - c22*r22 - 2*c32*r23 - c42*r24 - c52*r25 - c33*r33 - c34*r43 - c35*r53);
//gsl_matrix_set(m, 15, 16, -d43 - c40*r03 - c41*r13 - c42*r23 - c43*r33 - c44*r43 - c45*r53);
//gsl_matrix_set(m, 15, 17, -d53 - c50*r03 - c51*r13 - c52*r23 - c53*r33 - c54*r43 - c55*r53);
//
//gsl_matrix_set(m, 15, 18, 0.0);
//gsl_matrix_set(m, 15, 19, 0.0);
//gsl_matrix_set(m, 15, 20, 0.0);
//gsl_matrix_set(m, 15, 21, a23 - c03*r20 - c13*r21 - c23*r22 - c33*r23 - c43*r24 - c53*r25);
//gsl_matrix_set(m, 15, 22, 0.0);
//gsl_matrix_set(m, 15, 23, 0.0);
//
//gsl_matrix_set(m, 15, 24, 0.0);
//gsl_matrix_set(m, 15, 25, 0.0);
//gsl_matrix_set(m, 15, 26, 0.0);
//gsl_matrix_set(m, 15, 27, a24 - c04*r20 - c14*r21 - c24*r22 - c34*r23 - c44*r24 - c54*r25);
//gsl_matrix_set(m, 15, 28, 0.0);
//gsl_matrix_set(m, 15, 29, 0.0);
//
//gsl_matrix_set(m, 15, 30, 0.0);
//gsl_matrix_set(m, 15, 31, 0.0);
//gsl_matrix_set(m, 15, 32, 0.0);
//gsl_matrix_set(m, 15, 33, a25 - c05*r20 - c15*r21 - c25*r22 - c35*r23 - c45*r24 - c55*r25);
//gsl_matrix_set(m, 15, 34, 0.0);
//gsl_matrix_set(m, 15, 35, 0.0);
//
//
//
//gsl_matrix_set(m, 16, 0, 0.0);
//gsl_matrix_set(m, 16, 1, 0.0);
//gsl_matrix_set(m, 16, 2, 0.0);
//gsl_matrix_set(m, 16, 3, 0.0);
//gsl_matrix_set(m, 16, 4, a20 - c00*r20 - c10*r21 - c20*r22 - c30*r23 - c40*r24 - c50*r25);
//gsl_matrix_set(m, 16, 5, 0.0);
//
//gsl_matrix_set(m, 16, 6, 0.0);
//gsl_matrix_set(m, 16, 7, 0.0);
//gsl_matrix_set(m, 16, 8, 0.0);
//gsl_matrix_set(m, 16, 9, 0.0);
//gsl_matrix_set(m, 16, 10, a21 - c01*r20 - c11*r21 - c21*r22 - c31*r23 - c41*r24 - c51*r25);
//gsl_matrix_set(m, 16, 11, 0.0);
//
//gsl_matrix_set(m, 16, 12, -d04 - c00*r04 - c01*r14 - c02*r24 - c03*r34 - c04*r44 - c05*r54);
//gsl_matrix_set(m, 16, 13, -d14 - c10*r04 - c11*r14 - c12*r24 - c13*r34 - c14*r44 - c15*r54);
//gsl_matrix_set(m, 16, 14, -d24 - c20*r04 - c21*r14 - c22*r24 - c23*r34 - c24*r44 - c25*r54);
//gsl_matrix_set(m, 16, 15, -d34 - c30*r04 - c31*r14 - c32*r24 - c33*r34 - c34*r44 - c35*r54);
//gsl_matrix_set(m, 16, 16, a22 - d44 - c40*r04 - c41*r14 - c02*r20 - c12*r21 - c22*r22 - c32*r23 - 2*c42*r24 - c52*r25 - c43*r34 - c44*r44 - c45*r54);
//gsl_matrix_set(m, 16, 17, -d54 - c50*r04 - c51*r14 - c52*r24 - c53*r34 - c54*r44 - c55*r54);
//
//gsl_matrix_set(m, 16, 18, 0.0);
//gsl_matrix_set(m, 16, 19, 0.0);
//gsl_matrix_set(m, 16, 20, 0.0);
//gsl_matrix_set(m, 16, 21, 0.0);
//gsl_matrix_set(m, 16, 22, a23 - c03*r20 - c13*r21 - c23*r22 - c33*r23 - c43*r24 - c53*r25);
//gsl_matrix_set(m, 16, 23, 0.0);
//
//gsl_matrix_set(m, 16, 24, 0.0);
//gsl_matrix_set(m, 16, 25, 0.0);
//gsl_matrix_set(m, 16, 26, 0.0);
//gsl_matrix_set(m, 16, 27, 0.0);
//gsl_matrix_set(m, 16, 28, a24 - c04*r20 - c14*r21 - c24*r22 - c34*r23 - c44*r24 - c54*r25);
//gsl_matrix_set(m, 16, 29, 0.0);
//
//gsl_matrix_set(m, 16, 30, 0.0);
//gsl_matrix_set(m, 16, 31, 0.0);
//gsl_matrix_set(m, 16, 32, 0.0);
//gsl_matrix_set(m, 16, 33, 0.0);
//gsl_matrix_set(m, 16, 34, a25 - c05*r20 - c15*r21 - c25*r22 - c35*r23 - c45*r24 - c55*r25);
//gsl_matrix_set(m, 16, 35, 0.0);
//
//
//gsl_matrix_set(m, 17, 0, 0.0);
//gsl_matrix_set(m, 17, 1, 0.0);
//gsl_matrix_set(m, 17, 2, 0.0);
//gsl_matrix_set(m, 17, 3, 0.0);
//gsl_matrix_set(m, 17, 4, 0.0);
//gsl_matrix_set(m, 17, 5, a20 - c00*r20 - c10*r21 - c20*r22 - c30*r23 - c40*r24 - c50*r25);
//
//gsl_matrix_set(m, 17, 6, 0.0);
//gsl_matrix_set(m, 17, 7, 0.0);
//gsl_matrix_set(m, 17, 8, 0.0);
//gsl_matrix_set(m, 17, 9, 0.0);
//gsl_matrix_set(m, 17, 10, 0.0);
//gsl_matrix_set(m, 17, 11, a21 - c01*r20 - c11*r21 - c21*r22 - c31*r23 - c41*r24 - c51*r25);
//
//gsl_matrix_set(m, 17, 12, -d05 - c00*r05 - c01*r15 - c02*r25 - c03*r35 - c04*r45 - c05*r55);
//gsl_matrix_set(m, 17, 13, -d15 - c10*r05 - c11*r15 - c12*r25 - c13*r35 - c14*r45 - c15*r55);
//gsl_matrix_set(m, 17, 14, -d25 - c20*r05 - c21*r15 - c22*r25 - c23*r35 - c24*r45 - c25*r55);
//gsl_matrix_set(m, 17, 15, -d35 - c30*r05 - c31*r15 - c32*r25 - c33*r35 - c34*r45 - c35*r55);
//gsl_matrix_set(m, 17, 16, -d45 - c40*r05 - c41*r15 - c42*r25 - c43*r35 - c44*r45 - c45*r55);
//gsl_matrix_set(m, 17, 17, a22 - d55 - c50*r05 - c51*r15 - c02*r20 - c12*r21 - c22*r22 - c32*r23 - c42*r24 - 2*c52*r25 - c53*r35 - c54*r45 - c55*r55);
//
//gsl_matrix_set(m, 17, 18, 0.0);
//gsl_matrix_set(m, 17, 19, 0.0);
//gsl_matrix_set(m, 17, 20, 0.0);
//gsl_matrix_set(m, 17, 21, 0.0);
//gsl_matrix_set(m, 17, 22, 0.0);
//gsl_matrix_set(m, 17, 23, a23 - c03*r20 - c13*r21 - c23*r22 - c33*r23 - c43*r24 - c53*r25);
//
//gsl_matrix_set(m, 17, 24, 0.0);
//gsl_matrix_set(m, 17, 25, 0.0);
//gsl_matrix_set(m, 17, 26, 0.0);
//gsl_matrix_set(m, 17, 27, 0.0);
//gsl_matrix_set(m, 17, 28, 0.0);
//gsl_matrix_set(m, 17, 29, a24 - c04*r20 - c14*r21 - c24*r22 - c34*r23 - c44*r24 - c54*r25);
//
//gsl_matrix_set(m, 17, 30, 0.0);
//gsl_matrix_set(m, 17, 31, 0.0);
//gsl_matrix_set(m, 17, 32, 0.0);
//gsl_matrix_set(m, 17, 33, 0.0);
//gsl_matrix_set(m, 17, 34, 0.0);
//gsl_matrix_set(m, 17, 35, a25 - c05*r20 - c15*r21 - c25*r22 - c35*r23 - c45*r24 - c55*r25);
//
//
//
//
//gsl_matrix_set(m, 18, 0, a30 - c00*r30 - c10*r31 - c20*r32 - c30*r33 - c40*r34 - c50*r35);
//gsl_matrix_set(m, 18, 1, 0.0);
//gsl_matrix_set(m, 18, 2, 0.0);
//gsl_matrix_set(m, 18, 3, 0.0);
//gsl_matrix_set(m, 18, 4, 0.0);
//gsl_matrix_set(m, 18, 5, 0.0);
//
//gsl_matrix_set(m, 18, 6, a31 - c01*r30 - c11*r31 - c21*r32 - c31*r33 - c41*r34 - c51*r35);
//gsl_matrix_set(m, 18, 7, 0.0);
//gsl_matrix_set(m, 18, 8, 0.0);
//gsl_matrix_set(m, 18, 9, 0.0);
//gsl_matrix_set(m, 18, 10, 0.0);
//gsl_matrix_set(m, 18, 11, 0.0);
//
//gsl_matrix_set(m, 18, 12, a32 - c02*r30 - c12*r31 - c22*r32 - c32*r33 - c42*r34 - c52*r35);
//gsl_matrix_set(m, 18, 13, 0.0);
//gsl_matrix_set(m, 18, 14, 0.0);
//gsl_matrix_set(m, 18, 15, 0.0);
//gsl_matrix_set(m, 18, 16, 0.0);
//gsl_matrix_set(m, 18, 17, 0.0);
//
//gsl_matrix_set(m, 18, 18, a33 - d00 - c00*r00 - c01*r10 - c02*r20 - 2*c03*r30 - c13*r31 - c23*r32 - c33*r33 - c43*r34 - c53*r35 - c04*r40 - c05*r50);
//gsl_matrix_set(m, 18, 19, -d10 - c10*r00 - c11*r10 - c12*r20 - c13*r30 - c14*r40 - c15*r50);
//gsl_matrix_set(m, 18, 20, -d20 - c20*r00 - c21*r10 - c22*r20 - c23*r30 - c24*r40 - c25*r50);
//gsl_matrix_set(m, 18, 21, -d30 - c30*r00 - c31*r10 - c32*r20 - c33*r30 - c34*r40 - c35*r50);
//gsl_matrix_set(m, 18, 22, -d40 - c40*r00 - c41*r10 - c42*r20 - c43*r30 - c44*r40 - c45*r50);
//gsl_matrix_set(m, 18, 23, -d50 - c50*r00 - c51*r10 - c52*r20 - c53*r30 - c54*r40 - c55*r50);
//
//gsl_matrix_set(m, 18, 24, a34 - c04*r30 - c14*r31 - c24*r32 - c34*r33 - c44*r34 - c54*r35);
//gsl_matrix_set(m, 18, 25, 0.0);
//gsl_matrix_set(m, 18, 26, 0.0);
//gsl_matrix_set(m, 18, 27, 0.0);
//gsl_matrix_set(m, 18, 28, 0.0);
//gsl_matrix_set(m, 18, 29, 0.0);
//
//gsl_matrix_set(m, 18, 30, a35 - c05*r30 - c15*r31 - c25*r32 - c35*r33 - c45*r34 - c55*r35);
//gsl_matrix_set(m, 18, 31, 0.0);
//gsl_matrix_set(m, 18, 32, 0.0);
//gsl_matrix_set(m, 18, 33, 0.0);
//gsl_matrix_set(m, 18, 34, 0.0);
//gsl_matrix_set(m, 18, 35, 0.0);
//
//
//
//
//gsl_matrix_set(m, 19, 0, 0.0);
//gsl_matrix_set(m, 19, 1, a30 - c00*r30 - c10*r31 - c20*r32 - c30*r33 - c40*r34 - c50*r35);
//gsl_matrix_set(m, 19, 2, 0.0);
//gsl_matrix_set(m, 19, 3, 0.0);
//gsl_matrix_set(m, 19, 4, 0.0);
//gsl_matrix_set(m, 19, 5, 0.0);
//
//gsl_matrix_set(m, 19, 6, 0.0);
//gsl_matrix_set(m, 19, 7, a31 - c01*r30 - c11*r31 - c21*r32 - c31*r33 - c41*r34 - c51*r35);
//gsl_matrix_set(m, 19, 8, 0.0);
//gsl_matrix_set(m, 19, 9, 0.0);
//gsl_matrix_set(m, 19, 10, 0.0);
//gsl_matrix_set(m, 19, 11, 0.0);
//
//gsl_matrix_set(m, 19, 12, 0.0);
//gsl_matrix_set(m, 19, 13, a32 - c02*r30 - c12*r31 - c22*r32 - c32*r33 - c42*r34 - c52*r35);
//gsl_matrix_set(m, 19, 14, 0.0);
//gsl_matrix_set(m, 19, 15, 0.0);
//gsl_matrix_set(m, 19, 16, 0.0);
//gsl_matrix_set(m, 19, 17, 0.0);
//
//gsl_matrix_set(m, 19, 18, -d01 - c00*r01 - c01*r11 - c02*r21 - c03*r31 - c04*r41 - c05*r51);
//gsl_matrix_set(m, 19, 19, a33 - d11 - c10*r01 - c11*r11 - c12*r21 - c03*r30 - 2*c13*r31 - c23*r32 - c33*r33 - c43*r34 - c53*r35 - c14*r41 - c15*r51);
//gsl_matrix_set(m, 19, 20, -d21 - c20*r01 - c21*r11 - c22*r21 - c23*r31 - c24*r41 - c25*r51);
//gsl_matrix_set(m, 19, 21, -d31 - c30*r01 - c31*r11 - c32*r21 - c33*r31 - c34*r41 - c35*r51);
//gsl_matrix_set(m, 19, 22, -d41 - c40*r01 - c41*r11 - c42*r21 - c43*r31 - c44*r41 - c45*r51);
//gsl_matrix_set(m, 19, 23, -d51 - c50*r01 - c51*r11 - c52*r21 - c53*r31 - c54*r41 - c55*r51);
//
//gsl_matrix_set(m, 19, 24, 0.0);
//gsl_matrix_set(m, 19, 25, a34 - c04*r30 - c14*r31 - c24*r32 - c34*r33 - c44*r34 - c54*r35);
//gsl_matrix_set(m, 19, 26, 0.0);
//gsl_matrix_set(m, 19, 27, 0.0);
//gsl_matrix_set(m, 19, 28, 0.0);
//gsl_matrix_set(m, 19, 29, 0.0);
//
//gsl_matrix_set(m, 19, 30, 0.0);
//gsl_matrix_set(m, 19, 31, a35 - c05*r30 - c15*r31 - c25*r32 - c35*r33 - c45*r34 - c55*r35);
//gsl_matrix_set(m, 19, 32, 0.0);
//gsl_matrix_set(m, 19, 33, 0.0);
//gsl_matrix_set(m, 19, 34, 0.0);
//gsl_matrix_set(m, 19, 35, 0.0);
//
//
//
//gsl_matrix_set(m, 20, 0, 0.0);
//gsl_matrix_set(m, 20, 1, 0.0);
//gsl_matrix_set(m, 20, 2, a30 - c00*r30 - c10*r31 - c20*r32 - c30*r33 - c40*r34 - c50*r35);
//gsl_matrix_set(m, 20, 3, 0.0);
//gsl_matrix_set(m, 20, 4, 0.0);
//gsl_matrix_set(m, 20, 5, 0.0);
//
//gsl_matrix_set(m, 20, 6, 0.0);
//gsl_matrix_set(m, 20, 7, 0.0);
//gsl_matrix_set(m, 20, 8, a31 - c01*r30 - c11*r31 - c21*r32 - c31*r33 - c41*r34 - c51*r35);
//gsl_matrix_set(m, 20, 9, 0.0);
//gsl_matrix_set(m, 20, 10, 0.0);
//gsl_matrix_set(m, 20, 11, 0.0);
//
//gsl_matrix_set(m, 20, 12, 0.0);
//gsl_matrix_set(m, 20, 13, 0.0);
//gsl_matrix_set(m, 20, 14, a32 - c02*r30 - c12*r31 - c22*r32 - c32*r33 - c42*r34 - c52*r35);
//gsl_matrix_set(m, 20, 15, 0.0);
//gsl_matrix_set(m, 20, 16, 0.0);
//gsl_matrix_set(m, 20, 17, 0.0);
//
//gsl_matrix_set(m, 20, 18, -d02 - c00*r02 - c01*r12 - c02*r22 - c03*r32 - c04*r42 - c05*r52);
//gsl_matrix_set(m, 20, 19, -d12 - c10*r02 - c11*r12 - c12*r22 - c13*r32 - c14*r42 - c15*r52);
//gsl_matrix_set(m, 20, 20, a33 - d22 - c20*r02 - c21*r12 - c22*r22 - c03*r30 - c13*r31 - 2*c23*r32 - c33*r33 - c43*r34 - c53*r35 - c24*r42 - c25*r52);
//gsl_matrix_set(m, 20, 21, -d32 - c30*r02 - c31*r12 - c32*r22 - c33*r32 - c34*r42 - c35*r52);
//gsl_matrix_set(m, 20, 22, -d42 - c40*r02 - c41*r12 - c42*r22 - c43*r32 - c44*r42 - c45*r52);
//gsl_matrix_set(m, 20, 23, -d52 - c50*r02 - c51*r12 - c52*r22 - c53*r32 - c54*r42 - c55*r52);
//
//gsl_matrix_set(m, 20, 24, 0.0);
//gsl_matrix_set(m, 20, 25, 0.0);
//gsl_matrix_set(m, 20, 26, a34 - c04*r30 - c14*r31 - c24*r32 - c34*r33 - c44*r34 - c54*r35);
//gsl_matrix_set(m, 20, 27, 0.0);
//gsl_matrix_set(m, 20, 28, 0.0);
//gsl_matrix_set(m, 20, 29, 0.0);
//
//gsl_matrix_set(m, 20, 30, 0.0);
//gsl_matrix_set(m, 20, 31, 0.0);
//gsl_matrix_set(m, 20, 32, a35 - c05*r30 - c15*r31 - c25*r32 - c35*r33 - c45*r34 - c55*r35);
//gsl_matrix_set(m, 20, 33, 0.0);
//gsl_matrix_set(m, 20, 34, 0.0);
//gsl_matrix_set(m, 20, 35, 0.0);
//
//
//gsl_matrix_set(m, 21, 0, 0.0);
//gsl_matrix_set(m, 21, 1, 0.0);
//gsl_matrix_set(m, 21, 2, 0.0);
//gsl_matrix_set(m, 21, 3, a30 - c00*r30 - c10*r31 - c20*r32 - c30*r33 - c40*r34 - c50*r35);
//gsl_matrix_set(m, 21, 4, 0.0);
//gsl_matrix_set(m, 21, 5, 0.0);
//
//gsl_matrix_set(m, 21, 6, 0.0);
//gsl_matrix_set(m, 21, 7, 0.0);
//gsl_matrix_set(m, 21, 8, 0.0);
//gsl_matrix_set(m, 21, 9, a31 - c01*r30 - c11*r31 - c21*r32 - c31*r33 - c41*r34 - c51*r35);
//gsl_matrix_set(m, 21, 10, 0.0);
//gsl_matrix_set(m, 21, 11, 0.0);
//
//gsl_matrix_set(m, 21, 12, 0.0);
//gsl_matrix_set(m, 21, 13, 0.0);
//gsl_matrix_set(m, 21, 14, 0.0);
//gsl_matrix_set(m, 21, 15, a32 - c02*r30 - c12*r31 - c22*r32 - c32*r33 - c42*r34 - c52*r35);
//gsl_matrix_set(m, 21, 16, 0.0);
//gsl_matrix_set(m, 21, 17, 0.0);
//
//gsl_matrix_set(m, 21, 18, -d03 - c00*r03 - c01*r13 - c02*r23 - c03*r33 - c04*r43 - c05*r53);
//gsl_matrix_set(m, 21, 19, -d13 - c10*r03 - c11*r13 - c12*r23 - c13*r33 - c14*r43 - c15*r53);
//gsl_matrix_set(m, 21, 20, -d23 - c20*r03 - c21*r13 - c22*r23 - c23*r33 - c24*r43 - c25*r53);
//gsl_matrix_set(m, 21, 21, a33 - d33 - c30*r03 - c31*r13 - c32*r23 - c03*r30 - c13*r31 - c23*r32 - 2*c33*r33 - c43*r34 - c53*r35 - c34*r43 - c35*r53);
//gsl_matrix_set(m, 21, 22, -d43 - c40*r03 - c41*r13 - c42*r23 - c43*r33 - c44*r43 - c45*r53);
//gsl_matrix_set(m, 21, 23, -d53 - c50*r03 - c51*r13 - c52*r23 - c53*r33 - c54*r43 - c55*r53);
//
//gsl_matrix_set(m, 21, 24, 0.0);
//gsl_matrix_set(m, 21, 25, 0.0);
//gsl_matrix_set(m, 21, 26, 0.0);
//gsl_matrix_set(m, 21, 27, a34 - c04*r30 - c14*r31 - c24*r32 - c34*r33 - c44*r34 - c54*r35);
//gsl_matrix_set(m, 21, 28, 0.0);
//gsl_matrix_set(m, 21, 29, 0.0);
//
//gsl_matrix_set(m, 21, 30, 0.0);
//gsl_matrix_set(m, 21, 31, 0.0);
//gsl_matrix_set(m, 21, 32, 0.0);
//gsl_matrix_set(m, 21, 33, a35 - c05*r30 - c15*r31 - c25*r32 - c35*r33 - c45*r34 - c55*r35);
//gsl_matrix_set(m, 21, 34, 0.0);
//gsl_matrix_set(m, 21, 35, 0.0);
//
//
//
//gsl_matrix_set(m, 22, 0, 0.0);
//gsl_matrix_set(m, 22, 1, 0.0);
//gsl_matrix_set(m, 22, 2, 0.0);
//gsl_matrix_set(m, 22, 3, 0.0);
//gsl_matrix_set(m, 22, 4, a30 - c00*r30 - c10*r31 - c20*r32 - c30*r33 - c40*r34 - c50*r35);
//gsl_matrix_set(m, 22, 5, 0.0);
//
//gsl_matrix_set(m, 22, 6, 0.0);
//gsl_matrix_set(m, 22, 7, 0.0);
//gsl_matrix_set(m, 22, 8, 0.0);
//gsl_matrix_set(m, 22, 9, 0.0);
//gsl_matrix_set(m, 22, 10, a31 - c01*r30 - c11*r31 - c21*r32 - c31*r33 - c41*r34 - c51*r35);
//gsl_matrix_set(m, 22, 11, 0.0);
//
//gsl_matrix_set(m, 22, 12, 0.0);
//gsl_matrix_set(m, 22, 13, 0.0);
//gsl_matrix_set(m, 22, 14, 0.0);
//gsl_matrix_set(m, 22, 15, 0.0);
//gsl_matrix_set(m, 22, 16, a32 - c02*r30 - c12*r31 - c22*r32 - c32*r33 - c42*r34 - c52*r35);
//gsl_matrix_set(m, 22, 17, 0.0);
//
//gsl_matrix_set(m, 22, 18, -d04 - c00*r04 - c01*r14 - c02*r24 - c03*r34 - c04*r44 - c05*r54);
//gsl_matrix_set(m, 22, 19, -d14 - c10*r04 - c11*r14 - c12*r24 - c13*r34 - c14*r44 - c15*r54);
//gsl_matrix_set(m, 22, 20, -d24 - c20*r04 - c21*r14 - c22*r24 - c23*r34 - c24*r44 - c25*r54);
//gsl_matrix_set(m, 22, 21, -d34 - c30*r04 - c31*r14 - c32*r24 - c33*r34 - c34*r44 - c35*r54);
//gsl_matrix_set(m, 22, 22, a33 - d44 - c40*r04 - c41*r14 - c42*r24 - c03*r30 - c13*r31 - c23*r32 - c33*r33 - 2*c43*r34 - c53*r35 - c44*r44 - c45*r54);
//gsl_matrix_set(m, 22, 23, -d54 - c50*r04 - c51*r14 - c52*r24 - c53*r34 - c54*r44 - c55*r54);
//
//gsl_matrix_set(m, 22, 24, 0.0);
//gsl_matrix_set(m, 22, 25, 0.0);
//gsl_matrix_set(m, 22, 26, 0.0);
//gsl_matrix_set(m, 22, 27, 0.0);
//gsl_matrix_set(m, 22, 28, a34 - c04*r30 - c14*r31 - c24*r32 - c34*r33 - c44*r34 - c54*r35);
//gsl_matrix_set(m, 22, 29, 0.0);
//
//gsl_matrix_set(m, 22, 30, 0.0);
//gsl_matrix_set(m, 22, 31, 0.0);
//gsl_matrix_set(m, 22, 32, 0.0);
//gsl_matrix_set(m, 22, 33, 0.0);
//gsl_matrix_set(m, 22, 34, a35 - c05*r30 - c15*r31 - c25*r32 - c35*r33 - c45*r34 - c55*r35);
//gsl_matrix_set(m, 22, 35, 0.0);
//
//
//gsl_matrix_set(m, 23, 0, 0.0);
//gsl_matrix_set(m, 23, 1, 0.0);
//gsl_matrix_set(m, 23, 2, 0.0);
//gsl_matrix_set(m, 23, 3, 0.0);
//gsl_matrix_set(m, 23, 4, 0.0);
//gsl_matrix_set(m, 23, 5, a30 - c00*r30 - c10*r31 - c20*r32 - c30*r33 - c40*r34 - c50*r35);
//
//gsl_matrix_set(m, 23, 6, 0.0);
//gsl_matrix_set(m, 23, 7, 0.0);
//gsl_matrix_set(m, 23, 8, 0.0);
//gsl_matrix_set(m, 23, 9, 0.0);
//gsl_matrix_set(m, 23, 10, 0.0);
//gsl_matrix_set(m, 23, 11, a31 - c01*r30 - c11*r31 - c21*r32 - c31*r33 - c41*r34 - c51*r35);
//
//gsl_matrix_set(m, 23, 12, 0.0);
//gsl_matrix_set(m, 23, 13, 0.0);
//gsl_matrix_set(m, 23, 14, 0.0);
//gsl_matrix_set(m, 23, 15, 0.0);
//gsl_matrix_set(m, 23, 16, 0.0);
//gsl_matrix_set(m, 23, 17, a32 - c02*r30 - c12*r31 - c22*r32 - c32*r33 - c42*r34 - c52*r35);
//
//gsl_matrix_set(m, 23, 18, -d05 - c00*r05 - c01*r15 - c02*r25 - c03*r35 - c04*r45 - c05*r55);
//gsl_matrix_set(m, 23, 19, -d15 - c10*r05 - c11*r15 - c12*r25 - c13*r35 - c14*r45 - c15*r55);
//gsl_matrix_set(m, 23, 20, -d25 - c20*r05 - c21*r15 - c22*r25 - c23*r35 - c24*r45 - c25*r55);
//gsl_matrix_set(m, 23, 21, -d35 - c30*r05 - c31*r15 - c32*r25 - c33*r35 - c34*r45 - c35*r55);
//gsl_matrix_set(m, 23, 22, -d45 - c40*r05 - c41*r15 - c42*r25 - c43*r35 - c44*r45 - c45*r55);
//gsl_matrix_set(m, 23, 23, a33 - d55 - c50*r05 - c51*r15 - c52*r25 - c03*r30 - c13*r31 - c23*r32 - c33*r33 - c43*r34 - 2*c53*r35 - c54*r45 - c55*r55);
//
//gsl_matrix_set(m, 23, 24, 0.0);
//gsl_matrix_set(m, 23, 25, 0.0);
//gsl_matrix_set(m, 23, 26, 0.0);
//gsl_matrix_set(m, 23, 27, 0.0);
//gsl_matrix_set(m, 23, 28, 0.0);
//gsl_matrix_set(m, 23, 29, a34 - c04*r30 - c14*r31 - c24*r32 - c34*r33 - c44*r34 - c54*r35);
//
//gsl_matrix_set(m, 23, 30, 0.0);
//gsl_matrix_set(m, 23, 31, 0.0);
//gsl_matrix_set(m, 23, 32, 0.0);
//gsl_matrix_set(m, 23, 33, 0.0);
//gsl_matrix_set(m, 23, 34, 0.0);
//gsl_matrix_set(m, 23, 35, a35 - c05*r30 - c15*r31 - c25*r32 - c35*r33 - c45*r34 - c55*r35);
//
//
//
//gsl_matrix_set(m, 24, 0, a40 - c00*r40 - c10*r41 - c20*r42 - c30*r43 - c40*r44 - c50*r45);
//gsl_matrix_set(m, 24, 1, 0.0);
//gsl_matrix_set(m, 24, 2, 0.0);
//gsl_matrix_set(m, 24, 3, 0.0);
//gsl_matrix_set(m, 24, 4, 0.0);
//gsl_matrix_set(m, 24, 5, 0.0);
//
//gsl_matrix_set(m, 24, 6, a41 - c01*r40 - c11*r41 - c21*r42 - c31*r43 - c41*r44 - c51*r45);
//gsl_matrix_set(m, 24, 7, 0.0);
//gsl_matrix_set(m, 24, 8, 0.0);
//gsl_matrix_set(m, 24, 9, 0.0);
//gsl_matrix_set(m, 24, 10, 0.0);
//gsl_matrix_set(m, 24, 11, 0.0);
//
//gsl_matrix_set(m, 24, 12, a42 - c02*r40 - c12*r41 - c22*r42 - c32*r43 - c42*r44 - c52*r45);
//gsl_matrix_set(m, 24, 13, 0.0);
//gsl_matrix_set(m, 24, 14, 0.0);
//gsl_matrix_set(m, 24, 15, 0.0);
//gsl_matrix_set(m, 24, 16, 0.0);
//gsl_matrix_set(m, 24, 17, 0.0);
//
//gsl_matrix_set(m, 24, 18, a43 - c03*r40 - c13*r41 - c23*r42 - c33*r43 - c43*r44 - c53*r45);
//gsl_matrix_set(m, 24, 19, 0.0);
//gsl_matrix_set(m, 24, 20, 0.0);
//gsl_matrix_set(m, 24, 21, 0.0);
//gsl_matrix_set(m, 24, 22, 0.0);
//gsl_matrix_set(m, 24, 23, 0.0);
//
//gsl_matrix_set(m, 24, 24, a44 - d00 - c00*r00 - c01*r10 - c02*r20 - c03*r30 - 2*c04*r40 - c14*r41 - c24*r42 - c34*r43 - c44*r44 - c54*r45 - c05*r50);
//gsl_matrix_set(m, 24, 25, -d10 - c10*r00 - c11*r10 - c12*r20 - c13*r30 - c14*r40 - c15*r50);
//gsl_matrix_set(m, 24, 26, -d20 - c20*r00 - c21*r10 - c22*r20 - c23*r30 - c24*r40 - c25*r50);
//gsl_matrix_set(m, 24, 27, -d30 - c30*r00 - c31*r10 - c32*r20 - c33*r30 - c34*r40 - c35*r50);
//gsl_matrix_set(m, 24, 28, -d40 - c40*r00 - c41*r10 - c42*r20 - c43*r30 - c44*r40 - c45*r50);
//gsl_matrix_set(m, 24, 29, -d50 - c50*r00 - c51*r10 - c52*r20 - c53*r30 - c54*r40 - c55*r50);
//
//gsl_matrix_set(m, 24, 30, a45 - c05*r40 - c15*r41 - c25*r42 - c35*r43 - c45*r44 - c55*r45);
//gsl_matrix_set(m, 24, 31, 0.0);
//gsl_matrix_set(m, 24, 32, 0.0);
//gsl_matrix_set(m, 24, 33, 0.0);
//gsl_matrix_set(m, 24, 34, 0.0);
//gsl_matrix_set(m, 24, 35, 0.0);
//
//
//gsl_matrix_set(m, 25, 0, 0.0);
//gsl_matrix_set(m, 25, 1, a40 - c00*r40 - c10*r41 - c20*r42 - c30*r43 - c40*r44 - c50*r45);
//gsl_matrix_set(m, 25, 2, 0.0);
//gsl_matrix_set(m, 25, 3, 0.0);
//gsl_matrix_set(m, 25, 4, 0.0);
//gsl_matrix_set(m, 25, 5, 0.0);
//
//gsl_matrix_set(m, 25, 6, 0.0);
//gsl_matrix_set(m, 25, 7, a41 - c01*r40 - c11*r41 - c21*r42 - c31*r43 - c41*r44 - c51*r45);
//gsl_matrix_set(m, 25, 8, 0.0);
//gsl_matrix_set(m, 25, 9, 0.0);
//gsl_matrix_set(m, 25, 10, 0.0);
//gsl_matrix_set(m, 25, 11, 0.0);
//
//gsl_matrix_set(m, 25, 12, 0.0);
//gsl_matrix_set(m, 25, 13, a42 - c02*r40 - c12*r41 - c22*r42 - c32*r43 - c42*r44 - c52*r45);
//gsl_matrix_set(m, 25, 14, 0.0);
//gsl_matrix_set(m, 25, 15, 0.0);
//gsl_matrix_set(m, 25, 16, 0.0);
//gsl_matrix_set(m, 25, 17, 0.0);
//
//gsl_matrix_set(m, 25, 18, 0.0);
//gsl_matrix_set(m, 25, 19, a43 - c03*r40 - c13*r41 - c23*r42 - c33*r43 - c43*r44 - c53*r45);
//gsl_matrix_set(m, 25, 20, 0.0);
//gsl_matrix_set(m, 25, 21, 0.0);
//gsl_matrix_set(m, 25, 22, 0.0);
//gsl_matrix_set(m, 25, 23, 0.0);
//
//gsl_matrix_set(m, 25, 24, -d01 - c00*r01 - c01*r11 - c02*r21 - c03*r31 - c04*r41 - c05*r51);
//gsl_matrix_set(m, 25, 25, a44 - d11 - c10*r01 - c11*r11 - c12*r21 - c13*r31 - c04*r40 - 2*c14*r41 - c24*r42 - c34*r43 - c44*r44 - c54*r45 - c15*r51);
//gsl_matrix_set(m, 25, 26, -d21 - c20*r01 - c21*r11 - c22*r21 - c23*r31 - c24*r41 - c25*r51);
//gsl_matrix_set(m, 25, 27, -d31 - c30*r01 - c31*r11 - c32*r21 - c33*r31 - c34*r41 - c35*r51);
//gsl_matrix_set(m, 25, 28, -d41 - c40*r01 - c41*r11 - c42*r21 - c43*r31 - c44*r41 - c45*r51);
//gsl_matrix_set(m, 25, 29, -d51 - c50*r01 - c51*r11 - c52*r21 - c53*r31 - c54*r41 - c55*r51);
//
//gsl_matrix_set(m, 25, 30, 0.0);
//gsl_matrix_set(m, 25, 31, a45 - c05*r40 - c15*r41 - c25*r42 - c35*r43 - c45*r44 - c55*r45);
//gsl_matrix_set(m, 25, 32, 0.0);
//gsl_matrix_set(m, 25, 33, 0.0);
//gsl_matrix_set(m, 25, 34, 0.0);
//gsl_matrix_set(m, 25, 35, 0.0);
//
//
//
//
//gsl_matrix_set(m, 26, 0, 0.0);
//gsl_matrix_set(m, 26, 1, 0.0);
//gsl_matrix_set(m, 26, 2, a40 - c00*r40 - c10*r41 - c20*r42 - c30*r43 - c40*r44 - c50*r45);
//gsl_matrix_set(m, 26, 3, 0.0);
//gsl_matrix_set(m, 26, 4, 0.0);
//gsl_matrix_set(m, 26, 5, 0.0);
//
//gsl_matrix_set(m, 26, 6, 0.0);
//gsl_matrix_set(m, 26, 7, 0.0);
//gsl_matrix_set(m, 26, 8, a41 - c01*r40 - c11*r41 - c21*r42 - c31*r43 - c41*r44 - c51*r45);
//gsl_matrix_set(m, 26, 9, 0.0);
//gsl_matrix_set(m, 26, 10, 0.0);
//gsl_matrix_set(m, 26, 11, 0.0);
//
//gsl_matrix_set(m, 26, 12, 0.0);
//gsl_matrix_set(m, 26, 13, 0.0);
//gsl_matrix_set(m, 26, 14, a42 - c02*r40 - c12*r41 - c22*r42 - c32*r43 - c42*r44 - c52*r45);
//gsl_matrix_set(m, 26, 15, 0.0);
//gsl_matrix_set(m, 26, 16, 0.0);
//gsl_matrix_set(m, 26, 17, 0.0);
//
//gsl_matrix_set(m, 26, 18, 0.0);
//gsl_matrix_set(m, 26, 19, 0.0);
//gsl_matrix_set(m, 26, 20, a43 - c03*r40 - c13*r41 - c23*r42 - c33*r43 - c43*r44 - c53*r45);
//gsl_matrix_set(m, 26, 21, 0.0);
//gsl_matrix_set(m, 26, 22, 0.0);
//gsl_matrix_set(m, 26, 23, 0.0);
//
//gsl_matrix_set(m, 26, 24, -d02 - c00*r02 - c01*r12 - c02*r22 - c03*r32 - c04*r42 - c05*r52);
//gsl_matrix_set(m, 26, 25, -d12 - c10*r02 - c11*r12 - c12*r22 - c13*r32 - c14*r42 - c15*r52);
//gsl_matrix_set(m, 26, 26, a44 - d22 - c20*r02 - c21*r12 - c22*r22 - c23*r32 - c04*r40 - c14*r41 - 2*c24*r42 - c34*r43 - c44*r44 - c54*r45 - c25*r52);
//gsl_matrix_set(m, 26, 27, -d32 - c30*r02 - c31*r12 - c32*r22 - c33*r32 - c34*r42 - c35*r52);
//gsl_matrix_set(m, 26, 28, -d42 - c40*r02 - c41*r12 - c42*r22 - c43*r32 - c44*r42 - c45*r52);
//gsl_matrix_set(m, 26, 29, -d52 - c50*r02 - c51*r12 - c52*r22 - c53*r32 - c54*r42 - c55*r52);
//
//gsl_matrix_set(m, 26, 30, 0.0);
//gsl_matrix_set(m, 26, 31, 0.0);
//gsl_matrix_set(m, 26, 32, a45 - c05*r40 - c15*r41 - c25*r42 - c35*r43 - c45*r44 - c55*r45);
//gsl_matrix_set(m, 26, 33, 0.0);
//gsl_matrix_set(m, 26, 34, 0.0);
//gsl_matrix_set(m, 26, 35, 0.0);
//
//
//
//
//
//
//gsl_matrix_set(m, 27, 0, 0.0);
//gsl_matrix_set(m, 27, 1, 0.0);
//gsl_matrix_set(m, 27, 2, 0.0);
//gsl_matrix_set(m, 27, 3, a40 - c00*r40 - c10*r41 - c20*r42 - c30*r43 - c40*r44 - c50*r45);
//gsl_matrix_set(m, 27, 4, 0.0);
//gsl_matrix_set(m, 27, 5, 0.0);
//
//gsl_matrix_set(m, 27, 6, 0.0);
//gsl_matrix_set(m, 27, 7, 0.0);
//gsl_matrix_set(m, 27, 8, 0.0);
//gsl_matrix_set(m, 27, 9, a41 - c01*r40 - c11*r41 - c21*r42 - c31*r43 - c41*r44 - c51*r45);
//gsl_matrix_set(m, 27, 10, 0.0);
//gsl_matrix_set(m, 27, 11, 0.0);
//
//gsl_matrix_set(m, 27, 12, 0.0);
//gsl_matrix_set(m, 27, 13, 0.0);
//gsl_matrix_set(m, 27, 14, 0.0);
//gsl_matrix_set(m, 27, 15, a42 - c02*r40 - c12*r41 - c22*r42 - c32*r43 - c42*r44 - c52*r45);
//gsl_matrix_set(m, 27, 16, 0.0);
//gsl_matrix_set(m, 27, 17, 0.0);
//
//gsl_matrix_set(m, 27, 18, 0.0);
//gsl_matrix_set(m, 27, 19, 0.0);
//gsl_matrix_set(m, 27, 20, 0.0);
//gsl_matrix_set(m, 27, 21, a43 - c03*r40 - c13*r41 - c23*r42 - c33*r43 - c43*r44 - c53*r45);
//gsl_matrix_set(m, 27, 22, 0.0);
//gsl_matrix_set(m, 27, 23, 0.0);
//
//gsl_matrix_set(m, 27, 24, -d03 - c00*r03 - c01*r13 - c02*r23 - c03*r33 - c04*r43 - c05*r53);
//gsl_matrix_set(m, 27, 25, -d13 - c10*r03 - c11*r13 - c12*r23 - c13*r33 - c14*r43 - c15*r53);
//gsl_matrix_set(m, 27, 26, -d23 - c20*r03 - c21*r13 - c22*r23 - c23*r33 - c24*r43 - c25*r53);
//gsl_matrix_set(m, 27, 27, a44 - d33 - c30*r03 - c31*r13 - c32*r23 - c33*r33 - c04*r40 - c14*r41 - c24*r42 - 2*c34*r43 - c44*r44 - c54*r45 - c35*r53);
//gsl_matrix_set(m, 27, 28, -d43 - c40*r03 - c41*r13 - c42*r23 - c43*r33 - c44*r43 - c45*r53);
//gsl_matrix_set(m, 27, 29, -d53 - c50*r03 - c51*r13 - c52*r23 - c53*r33 - c54*r43 - c55*r53);
//
//gsl_matrix_set(m, 27, 30, 0.0);
//gsl_matrix_set(m, 27, 31, 0.0);
//gsl_matrix_set(m, 27, 32, 0.0);
//gsl_matrix_set(m, 27, 33, a45 - c05*r40 - c15*r41 - c25*r42 - c35*r43 - c45*r44 - c55*r45);
//gsl_matrix_set(m, 27, 34, 0.0);
//gsl_matrix_set(m, 27, 35, 0.0);
//
//
//
//
//
//
//gsl_matrix_set(m, 28, 0, 0.0);
//gsl_matrix_set(m, 28, 1, 0.0);
//gsl_matrix_set(m, 28, 2, 0.0);
//gsl_matrix_set(m, 28, 3, 0.0);
//gsl_matrix_set(m, 28, 4, a40 - c00*r40 - c10*r41 - c20*r42 - c30*r43 - c40*r44 - c50*r45);
//gsl_matrix_set(m, 28, 5, 0.0);
//
//gsl_matrix_set(m, 28, 6, 0.0);
//gsl_matrix_set(m, 28, 7, 0.0);
//gsl_matrix_set(m, 28, 8, 0.0);
//gsl_matrix_set(m, 28, 9, 0.0);
//gsl_matrix_set(m, 28, 10, a41 - c01*r40 - c11*r41 - c21*r42 - c31*r43 - c41*r44 - c51*r45);
//gsl_matrix_set(m, 28, 11, 0.0);
//
//gsl_matrix_set(m, 28, 12, 0.0);
//gsl_matrix_set(m, 28, 13, 0.0);
//gsl_matrix_set(m, 28, 14, 0.0);
//gsl_matrix_set(m, 28, 15, 0.0);
//gsl_matrix_set(m, 28, 16, a42 - c02*r40 - c12*r41 - c22*r42 - c32*r43 - c42*r44 - c52*r45);
//gsl_matrix_set(m, 28, 17, 0.0);
//
//gsl_matrix_set(m, 28, 18, 0.0);
//gsl_matrix_set(m, 28, 19, 0.0);
//gsl_matrix_set(m, 28, 20, 0.0);
//gsl_matrix_set(m, 28, 21, 0.0);
//gsl_matrix_set(m, 28, 22, a43 - c03*r40 - c13*r41 - c23*r42 - c33*r43 - c43*r44 - c53*r45);
//gsl_matrix_set(m, 28, 23, 0.0);
//
//gsl_matrix_set(m, 28, 24, -d04 - c00*r04 - c01*r14 - c02*r24 - c03*r34 - c04*r44 - c05*r54);
//gsl_matrix_set(m, 28, 25, -d14 - c10*r04 - c11*r14 - c12*r24 - c13*r34 - c14*r44 - c15*r54);
//gsl_matrix_set(m, 28, 26, -d24 - c20*r04 - c21*r14 - c22*r24 - c23*r34 - c24*r44 - c25*r54);
//gsl_matrix_set(m, 28, 27, -d34 - c30*r04 - c31*r14 - c32*r24 - c33*r34 - c34*r44 - c35*r54);
//gsl_matrix_set(m, 28, 28, a44 - d44 - c40*r04 - c41*r14 - c42*r24 - c43*r34 - c04*r40 - c14*r41 - c24*r42 - c34*r43 - 2*c44*r44 - c54*r45 - c45*r54);
//gsl_matrix_set(m, 28, 29, -d54 - c50*r04 - c51*r14 - c52*r24 - c53*r34 - c54*r44 - c55*r54);
//
//gsl_matrix_set(m, 28, 30, 0.0);
//gsl_matrix_set(m, 28, 31, 0.0);
//gsl_matrix_set(m, 28, 32, 0.0);
//gsl_matrix_set(m, 28, 33, 0.0);
//gsl_matrix_set(m, 28, 34, a45 - c05*r40 - c15*r41 - c25*r42 - c35*r43 - c45*r44 - c55*r45);
//gsl_matrix_set(m, 28, 35, 0.0);
//
//
//
//
//
//
//gsl_matrix_set(m, 29, 0, 0.0);
//gsl_matrix_set(m, 29, 1, 0.0);
//gsl_matrix_set(m, 29, 2, 0.0);
//gsl_matrix_set(m, 29, 3, 0.0);
//gsl_matrix_set(m, 29, 4, 0.0);
//gsl_matrix_set(m, 29, 5, a40 - c00*r40 - c10*r41 - c20*r42 - c30*r43 - c40*r44 - c50*r45);
//
//gsl_matrix_set(m, 29, 6, 0.0);
//gsl_matrix_set(m, 29, 7, 0.0);
//gsl_matrix_set(m, 29, 8, 0.0);
//gsl_matrix_set(m, 29, 9, 0.0);
//gsl_matrix_set(m, 29, 10, 0.0);
//gsl_matrix_set(m, 29, 11, a41 - c01*r40 - c11*r41 - c21*r42 - c31*r43 - c41*r44 - c51*r45);
//
//gsl_matrix_set(m, 29, 12, 0.0);
//gsl_matrix_set(m, 29, 13, 0.0);
//gsl_matrix_set(m, 29, 14, 0.0);
//gsl_matrix_set(m, 29, 15, 0.0);
//gsl_matrix_set(m, 29, 16, 0.0);
//gsl_matrix_set(m, 29, 17, a42 - c02*r40 - c12*r41 - c22*r42 - c32*r43 - c42*r44 - c52*r45);
//
//gsl_matrix_set(m, 29, 18, 0.0);
//gsl_matrix_set(m, 29, 19, 0.0);
//gsl_matrix_set(m, 29, 20, 0.0);
//gsl_matrix_set(m, 29, 21, 0.0);
//gsl_matrix_set(m, 29, 22, 0.0);
//gsl_matrix_set(m, 29, 23, a43 - c03*r40 - c13*r41 - c23*r42 - c33*r43 - c43*r44 - c53*r45);
//
//gsl_matrix_set(m, 29, 24, -d05 - c00*r05 - c01*r15 - c02*r25 - c03*r35 - c04*r45 - c05*r55);
//gsl_matrix_set(m, 29, 25, -d15 - c10*r05 - c11*r15 - c12*r25 - c13*r35 - c14*r45 - c15*r55);
//gsl_matrix_set(m, 29, 26, -d25 - c20*r05 - c21*r15 - c22*r25 - c23*r35 - c24*r45 - c25*r55);
//gsl_matrix_set(m, 29, 27, -d35 - c30*r05 - c31*r15 - c32*r25 - c33*r35 - c34*r45 - c35*r55);
//gsl_matrix_set(m, 29, 28, -d45 - c40*r05 - c41*r15 - c42*r25 - c43*r35 - c44*r45 - c45*r55);
//gsl_matrix_set(m, 29, 29, a44 - d55 - c50*r05 - c51*r15 - c52*r25 - c53*r35 - c04*r40 - c14*r41 - c24*r42 - c34*r43 - c44*r44 - 2*c54*r45 - c55*r55);
//
//gsl_matrix_set(m, 29, 30, 0.0);
//gsl_matrix_set(m, 29, 31, 0.0);
//gsl_matrix_set(m, 29, 32, 0.0);
//gsl_matrix_set(m, 29, 33, 0.0);
//gsl_matrix_set(m, 29, 34, 0.0);
//gsl_matrix_set(m, 29, 35, a45 - c05*r40 - c15*r41 - c25*r42 - c35*r43 - c45*r44 - c55*r45);
//
//
//
//
//
//
//gsl_matrix_set(m, 30, 0, a50 - c00*r50 - c10*r51 - c20*r52 - c30*r53 - c40*r54 - c50*r55);
//gsl_matrix_set(m, 30, 1, 0.0);
//gsl_matrix_set(m, 30, 2, 0.0);
//gsl_matrix_set(m, 30, 3, 0.0);
//gsl_matrix_set(m, 30, 4, 0.0);
//gsl_matrix_set(m, 30, 5, 0.0);
//
//gsl_matrix_set(m, 30, 6, a51 - c01*r50 - c11*r51 - c21*r52 - c31*r53 - c41*r54 - c51*r55);
//gsl_matrix_set(m, 30, 7, 0.0);
//gsl_matrix_set(m, 30, 8, 0.0);
//gsl_matrix_set(m, 30, 9, 0.0);
//gsl_matrix_set(m, 30, 10, 0.0);
//gsl_matrix_set(m, 30, 11, 0.0);
//
//gsl_matrix_set(m, 30, 12, a52 - c02*r50 - c12*r51 - c22*r52 - c32*r53 - c42*r54 - c52*r55);
//gsl_matrix_set(m, 30, 13, 0.0);
//gsl_matrix_set(m, 30, 14, 0.0);
//gsl_matrix_set(m, 30, 15, 0.0);
//gsl_matrix_set(m, 30, 16, 0.0);
//gsl_matrix_set(m, 30, 17, 0.0);
//
//gsl_matrix_set(m, 30, 18, a53 - c03*r50 - c13*r51 - c23*r52 - c33*r53 - c43*r54 - c53*r55);
//gsl_matrix_set(m, 30, 19, 0.0);
//gsl_matrix_set(m, 30, 20, 0.0);
//gsl_matrix_set(m, 30, 21, 0.0);
//gsl_matrix_set(m, 30, 22, 0.0);
//gsl_matrix_set(m, 30, 23, 0.0);
//
//gsl_matrix_set(m, 30, 24, a54 - c04*r50 - c14*r51 - c24*r52 - c34*r53 - c44*r54 - c54*r55);
//gsl_matrix_set(m, 30, 25, 0.0);
//gsl_matrix_set(m, 30, 26, 0.0);
//gsl_matrix_set(m, 30, 27, 0.0);
//gsl_matrix_set(m, 30, 28, 0.0);
//gsl_matrix_set(m, 30, 29, 0.0);
//
//gsl_matrix_set(m, 30, 30, a55 - d00 - c00*r00 - c01*r10 - c02*r20 - c03*r30 - c04*r40 - 2*c05*r50 - c15*r51 - c25*r52 - c35*r53 - c45*r54 - c55*r55);
//gsl_matrix_set(m, 30, 31, -d10 - c10*r00 - c11*r10 - c12*r20 - c13*r30 - c14*r40 - c15*r50);
//gsl_matrix_set(m, 30, 32, -d20 - c20*r00 - c21*r10 - c22*r20 - c23*r30 - c24*r40 - c25*r50);
//gsl_matrix_set(m, 30, 33, -d30 - c30*r00 - c31*r10 - c32*r20 - c33*r30 - c34*r40 - c35*r50);
//gsl_matrix_set(m, 30, 34, -d40 - c40*r00 - c41*r10 - c42*r20 - c43*r30 - c44*r40 - c45*r50);
//gsl_matrix_set(m, 30, 35, -d50 - c50*r00 - c51*r10 - c52*r20 - c53*r30 - c54*r40 - c55*r50);
//
//
//
//
//
//gsl_matrix_set(m, 31, 0, 0.0);
//gsl_matrix_set(m, 31, 1, a50 - c00*r50 - c10*r51 - c20*r52 - c30*r53 - c40*r54 - c50*r55);
//gsl_matrix_set(m, 31, 2, 0.0);
//gsl_matrix_set(m, 31, 3, 0.0);
//gsl_matrix_set(m, 31, 4, 0.0);
//gsl_matrix_set(m, 31, 5, 0.0);
//
//gsl_matrix_set(m, 31, 6, 0.0);
//gsl_matrix_set(m, 31, 7, a51 - c01*r50 - c11*r51 - c21*r52 - c31*r53 - c41*r54 - c51*r55);
//gsl_matrix_set(m, 31, 8, 0.0);
//gsl_matrix_set(m, 31, 9, 0.0);
//gsl_matrix_set(m, 31, 10, 0.0);
//gsl_matrix_set(m, 31, 11, 0.0);
//
//gsl_matrix_set(m, 31, 12, 0.0);
//gsl_matrix_set(m, 31, 13, a52 - c02*r50 - c12*r51 - c22*r52 - c32*r53 - c42*r54 - c52*r55);
//gsl_matrix_set(m, 31, 14, 0.0);
//gsl_matrix_set(m, 31, 15, 0.0);
//gsl_matrix_set(m, 31, 16, 0.0);
//gsl_matrix_set(m, 31, 17, 0.0);
//
//gsl_matrix_set(m, 31, 18, 0.0);
//gsl_matrix_set(m, 31, 19, a53 - c03*r50 - c13*r51 - c23*r52 - c33*r53 - c43*r54 - c53*r55);
//gsl_matrix_set(m, 31, 20, 0.0);
//gsl_matrix_set(m, 31, 21, 0.0);
//gsl_matrix_set(m, 31, 22, 0.0);
//gsl_matrix_set(m, 31, 23, 0.0);
//
//gsl_matrix_set(m, 31, 24, 0.0);
//gsl_matrix_set(m, 31, 25, a54 - c04*r50 - c14*r51 - c24*r52 - c34*r53 - c44*r54 - c54*r55);
//gsl_matrix_set(m, 31, 26, 0.0);
//gsl_matrix_set(m, 31, 27, 0.0);
//gsl_matrix_set(m, 31, 28, 0.0);
//gsl_matrix_set(m, 31, 29, 0.0);
//
//gsl_matrix_set(m, 31, 30, -d01 - c00*r01 - c01*r11 - c02*r21 - c03*r31 - c04*r41 - c05*r51);
//gsl_matrix_set(m, 31, 31, a55 - d11 - c10*r01 - c11*r11 - c12*r21 - c13*r31 - c14*r41 - c05*r50 - 2*c15*r51 - c25*r52 - c35*r53 - c45*r54 - c55*r55);
//gsl_matrix_set(m, 31, 32, -d21 - c20*r01 - c21*r11 - c22*r21 - c23*r31 - c24*r41 - c25*r51);
//gsl_matrix_set(m, 31, 33, -d31 - c30*r01 - c31*r11 - c32*r21 - c33*r31 - c34*r41 - c35*r51);
//gsl_matrix_set(m, 31, 34, -d41 - c40*r01 - c41*r11 - c42*r21 - c43*r31 - c44*r41 - c45*r51);
//gsl_matrix_set(m, 31, 35, -d51 - c50*r01 - c51*r11 - c52*r21 - c53*r31 - c54*r41 - c55*r51);
//
//
//
//
//gsl_matrix_set(m, 32, 0, 0.0);
//gsl_matrix_set(m, 32, 1, 0.0);
//gsl_matrix_set(m, 32, 2, a50 - c00*r50 - c10*r51 - c20*r52 - c30*r53 - c40*r54 - c50*r55);
//gsl_matrix_set(m, 32, 3, 0.0);
//gsl_matrix_set(m, 32, 4, 0.0);
//gsl_matrix_set(m, 32, 5, 0.0);
//
//gsl_matrix_set(m, 32, 6, 0.0);
//gsl_matrix_set(m, 32, 7, 0.0);
//gsl_matrix_set(m, 32, 8, a51 - c01*r50 - c11*r51 - c21*r52 - c31*r53 - c41*r54 - c51*r55);
//gsl_matrix_set(m, 32, 9, 0.0);
//gsl_matrix_set(m, 32, 10, 0.0);
//gsl_matrix_set(m, 32, 11, 0.0);
//
//gsl_matrix_set(m, 32, 12, 0.0);
//gsl_matrix_set(m, 32, 13, 0.0);
//gsl_matrix_set(m, 32, 14, a52 - c02*r50 - c12*r51 - c22*r52 - c32*r53 - c42*r54 - c52*r55);
//gsl_matrix_set(m, 32, 15, 0.0);
//gsl_matrix_set(m, 32, 16, 0.0);
//gsl_matrix_set(m, 32, 17, 0.0);
//
//gsl_matrix_set(m, 32, 18, 0.0);
//gsl_matrix_set(m, 32, 19, 0.0);
//gsl_matrix_set(m, 32, 20, a53 - c03*r50 - c13*r51 - c23*r52 - c33*r53 - c43*r54 - c53*r55);
//gsl_matrix_set(m, 32, 21, 0.0);
//gsl_matrix_set(m, 32, 22, 0.0);
//gsl_matrix_set(m, 32, 23, 0.0);
//
//gsl_matrix_set(m, 32, 24, 0.0);
//gsl_matrix_set(m, 32, 25, 0.0);
//gsl_matrix_set(m, 32, 26, a54 - c04*r50 - c14*r51 - c24*r52 - c34*r53 - c44*r54 - c54*r55);
//gsl_matrix_set(m, 32, 27, 0.0);
//gsl_matrix_set(m, 32, 28, 0.0);
//gsl_matrix_set(m, 32, 29, 0.0);
//
//gsl_matrix_set(m, 32, 30, -d02 - c00*r02 - c01*r12 - c02*r22 - c03*r32 - c04*r42 - c05*r52);
//gsl_matrix_set(m, 32, 31, -d12 - c10*r02 - c11*r12 - c12*r22 - c13*r32 - c14*r42 - c15*r52);
//gsl_matrix_set(m, 32, 32, a55 - d22 - c20*r02 - c21*r12 - c22*r22 - c23*r32 - c24*r42 - c05*r50 - c15*r51 - 2*c25*r52 - c35*r53 - c45*r54 - c55*r55);
//gsl_matrix_set(m, 32, 33, -d32 - c30*r02 - c31*r12 - c32*r22 - c33*r32 - c34*r42 - c35*r52);
//gsl_matrix_set(m, 32, 34, -d42 - c40*r02 - c41*r12 - c42*r22 - c43*r32 - c44*r42 - c45*r52);
//gsl_matrix_set(m, 32, 35, -d52 - c50*r02 - c51*r12 - c52*r22 - c53*r32 - c54*r42 - c55*r52);
//
//
//
//gsl_matrix_set(m, 33, 0, 0.0);
//gsl_matrix_set(m, 33, 1, 0.0);
//gsl_matrix_set(m, 33, 2, 0.0);
//gsl_matrix_set(m, 33, 3, a50 - c00*r50 - c10*r51 - c20*r52 - c30*r53 - c40*r54 - c50*r55);
//gsl_matrix_set(m, 33, 4, 0.0);
//gsl_matrix_set(m, 33, 5, 0.0);
//
//gsl_matrix_set(m, 33, 6, 0.0);
//gsl_matrix_set(m, 33, 7, 0.0);
//gsl_matrix_set(m, 33, 8, 0.0);
//gsl_matrix_set(m, 33, 9, a51 - c01*r50 - c11*r51 - c21*r52 - c31*r53 - c41*r54 - c51*r55);
//gsl_matrix_set(m, 33, 10, 0.0);
//gsl_matrix_set(m, 33, 11, 0.0);
//
//gsl_matrix_set(m, 33, 12, 0.0);
//gsl_matrix_set(m, 33, 13, 0.0);
//gsl_matrix_set(m, 33, 14, 0.0);
//gsl_matrix_set(m, 33, 15, a52 - c02*r50 - c12*r51 - c22*r52 - c32*r53 - c42*r54 - c52*r55);
//gsl_matrix_set(m, 33, 16, 0.0);
//gsl_matrix_set(m, 33, 17, 0.0);
//
//gsl_matrix_set(m, 33, 18, 0.0);
//gsl_matrix_set(m, 33, 19, 0.0);
//gsl_matrix_set(m, 33, 20, 0.0);
//gsl_matrix_set(m, 33, 21, a53 - c03*r50 - c13*r51 - c23*r52 - c33*r53 - c43*r54 - c53*r55);
//gsl_matrix_set(m, 33, 22, 0.0);
//gsl_matrix_set(m, 33, 23, 0.0);
//
//gsl_matrix_set(m, 33, 24, 0.0);
//gsl_matrix_set(m, 33, 25, 0.0);
//gsl_matrix_set(m, 33, 26, 0.0);
//gsl_matrix_set(m, 33, 27, a54 - c04*r50 - c14*r51 - c24*r52 - c34*r53 - c44*r54 - c54*r55);
//gsl_matrix_set(m, 33, 28, 0.0);
//gsl_matrix_set(m, 33, 29, 0.0);
//
//gsl_matrix_set(m, 33, 30, -d03 - c00*r03 - c01*r13 - c02*r23 - c03*r33 - c04*r43 - c05*r53);
//gsl_matrix_set(m, 33, 31, -d13 - c10*r03 - c11*r13 - c12*r23 - c13*r33 - c14*r43 - c15*r53);
//gsl_matrix_set(m, 33, 32, -d23 - c20*r03 - c21*r13 - c22*r23 - c23*r33 - c24*r43 - c25*r53);
//gsl_matrix_set(m, 33, 33, a55 - d33 - c30*r03 - c31*r13 - c32*r23 - c33*r33 - c34*r43 - c05*r50 - c15*r51 - c25*r52 - 2*c35*r53 - c45*r54 - c55*r55);
//gsl_matrix_set(m, 33, 34, -d43 - c40*r03 - c41*r13 - c42*r23 - c43*r33 - c44*r43 - c45*r53);
//gsl_matrix_set(m, 33, 35, -d53 - c50*r03 - c51*r13 - c52*r23 - c53*r33 - c54*r43 - c55*r53);
//
//
//
//
//gsl_matrix_set(m, 34, 0, 0.0);
//gsl_matrix_set(m, 34, 1, 0.0);
//gsl_matrix_set(m, 34, 2, 0.0);
//gsl_matrix_set(m, 34, 3, 0.0);
//gsl_matrix_set(m, 34, 4, a50 - c00*r50 - c10*r51 - c20*r52 - c30*r53 - c40*r54 - c50*r55);
//gsl_matrix_set(m, 34, 5, 0.0);
//
//gsl_matrix_set(m, 34, 6, 0.0);
//gsl_matrix_set(m, 34, 7, 0.0);
//gsl_matrix_set(m, 34, 8, 0.0);
//gsl_matrix_set(m, 34, 9, 0.0);
//gsl_matrix_set(m, 34, 10, a51 - c01*r50 - c11*r51 - c21*r52 - c31*r53 - c41*r54 - c51*r55);
//gsl_matrix_set(m, 34, 11, 0.0);
//
//gsl_matrix_set(m, 34, 12, 0.0);
//gsl_matrix_set(m, 34, 13, 0.0);
//gsl_matrix_set(m, 34, 14, 0.0);
//gsl_matrix_set(m, 34, 15, 0.0);
//gsl_matrix_set(m, 34, 16, a52 - c02*r50 - c12*r51 - c22*r52 - c32*r53 - c42*r54 - c52*r55);
//gsl_matrix_set(m, 34, 17, 0.0);
//
//gsl_matrix_set(m, 34, 18, 0.0);
//gsl_matrix_set(m, 34, 19, 0.0);
//gsl_matrix_set(m, 34, 20, 0.0);
//gsl_matrix_set(m, 34, 21, 0.0);
//gsl_matrix_set(m, 34, 22, a53 - c03*r50 - c13*r51 - c23*r52 - c33*r53 - c43*r54 - c53*r55);
//gsl_matrix_set(m, 34, 23, 0.0);
//
//gsl_matrix_set(m, 34, 24, 0.0);
//gsl_matrix_set(m, 34, 25, 0.0);
//gsl_matrix_set(m, 34, 26, 0.0);
//gsl_matrix_set(m, 34, 27, 0.0);
//gsl_matrix_set(m, 34, 28, a54 - c04*r50 - c14*r51 - c24*r52 - c34*r53 - c44*r54 - c54*r55);
//gsl_matrix_set(m, 34, 29, 0.0);
//
//
//gsl_matrix_set(m, 34, 30, -d04 - c00*r04 - c01*r14 - c02*r24 - c03*r34 - c04*r44 - c05*r54);
//gsl_matrix_set(m, 34, 31, -d14 - c10*r04 - c11*r14 - c12*r24 - c13*r34 - c14*r44 - c15*r54);
//gsl_matrix_set(m, 34, 32, -d24 - c20*r04 - c21*r14 - c22*r24 - c23*r34 - c24*r44 - c25*r54);
//gsl_matrix_set(m, 34, 33, -d34 - c30*r04 - c31*r14 - c32*r24 - c33*r34 - c34*r44 - c35*r54);
//gsl_matrix_set(m, 34, 34, a55 - d44 - c40*r04 - c41*r14 - c42*r24 - c43*r34 - c44*r44 - c05*r50 - c15*r51 - c25*r52 - c35*r53 - 2*c45*r54 - c55*r55);
//gsl_matrix_set(m, 34, 35, -d54 - c50*r04 - c51*r14 - c52*r24 - c53*r34 - c54*r44 - c55*r54);
//
//
//
//
//gsl_matrix_set(m, 35, 0, 0.0);
//gsl_matrix_set(m, 35, 1, 0.0);
//gsl_matrix_set(m, 35, 2, 0.0);
//gsl_matrix_set(m, 35, 3, 0.0);
//gsl_matrix_set(m, 35, 4, 0.0);
//gsl_matrix_set(m, 35, 5, a50 - c00*r50 - c10*r51 - c20*r52 - c30*r53 - c40*r54 - c50*r55);
//
//gsl_matrix_set(m, 35, 6, 0.0);
//gsl_matrix_set(m, 35, 7, 0.0);
//gsl_matrix_set(m, 35, 8, 0.0);
//gsl_matrix_set(m, 35, 9, 0.0);
//gsl_matrix_set(m, 35, 10, 0.0);
//gsl_matrix_set(m, 35, 11, a51 - c01*r50 - c11*r51 - c21*r52 - c31*r53 - c41*r54 - c51*r55);
//
//gsl_matrix_set(m, 35, 12, 0.0);
//gsl_matrix_set(m, 35, 13, 0.0);
//gsl_matrix_set(m, 35, 14, 0.0);
//gsl_matrix_set(m, 35, 15, 0.0);
//gsl_matrix_set(m, 35, 16, 0.0);
//gsl_matrix_set(m, 35, 17, a52 - c02*r50 - c12*r51 - c22*r52 - c32*r53 - c42*r54 - c52*r55);
//
//gsl_matrix_set(m, 35, 18, 0.0);
//gsl_matrix_set(m, 35, 19, 0.0);
//gsl_matrix_set(m, 35, 20, 0.0);
//gsl_matrix_set(m, 35, 21, 0.0);
//gsl_matrix_set(m, 35, 22, 0.0);
//gsl_matrix_set(m, 35, 23, a53 - c03*r50 - c13*r51 - c23*r52 - c33*r53 - c43*r54 - c53*r55);
//
//gsl_matrix_set(m, 35, 24, 0.0);
//gsl_matrix_set(m, 35, 25, 0.0);
//gsl_matrix_set(m, 35, 26, 0.0);
//gsl_matrix_set(m, 35, 27, 0.0);
//gsl_matrix_set(m, 35, 28, 0.0);
//gsl_matrix_set(m, 35, 29, a54 - c04*r50 - c14*r51 - c24*r52 - c34*r53 - c44*r54 - c54*r55);
//
//gsl_matrix_set(m, 35, 30, -d05 - c00*r05 - c01*r15 - c02*r25 - c03*r35 - c04*r45 - c05*r55);
//gsl_matrix_set(m, 35, 31, -d15 - c10*r05 - c11*r15 - c12*r25 - c13*r35 - c14*r45 - c15*r55);
//gsl_matrix_set(m, 35, 32, -d25 - c20*r05 - c21*r15 - c22*r25 - c23*r35 - c24*r45 - c25*r55);
//gsl_matrix_set(m, 35, 33, -d35 - c30*r05 - c31*r15 - c32*r25 - c33*r35 - c34*r45 - c35*r55);
//gsl_matrix_set(m, 35, 34, -d45 - c40*r05 - c41*r15 - c42*r25 - c43*r35 - c44*r45 - c45*r55);
//gsl_matrix_set(m, 35, 35, a55 - d55 - c50*r05 - c51*r15 - c52*r25 - c53*r35 - c54*r45 - c05*r50 - c15*r51 - c25*r52 - c35*r53 - c45*r54 - 2*c55*r55);
//
//
//
//gsl_matrix_set(m, 36, 0, c00*v0);
//gsl_matrix_set(m, 36, 1, c00*v1);
//gsl_matrix_set(m, 36, 2, c00*v2);
//gsl_matrix_set(m, 36, 3, c00*v3);
//gsl_matrix_set(m, 36, 4, c00*v4);
//gsl_matrix_set(m, 36, 5, c00*v5);
//gsl_matrix_set(m, 36, 6, c01*v0);
//gsl_matrix_set(m, 36, 7, c01*v1);
//gsl_matrix_set(m, 36, 8, c01*v2);
//gsl_matrix_set(m, 36, 9, c01*v3);
//gsl_matrix_set(m, 36, 10, c01*v4);
//gsl_matrix_set(m, 36, 11, c01*v5);
//gsl_matrix_set(m, 36, 12, c02*v0);
//gsl_matrix_set(m, 36, 13, c02*v1);
//gsl_matrix_set(m, 36, 14, c02*v2);
//gsl_matrix_set(m, 36, 15, c02*v3);
//gsl_matrix_set(m, 36, 16, c02*v4);
//gsl_matrix_set(m, 36, 17, c02*v5);
//gsl_matrix_set(m, 36, 18, c03*v0);
//gsl_matrix_set(m, 36, 19, c03*v1);
//gsl_matrix_set(m, 36, 20, c03*v2);
//gsl_matrix_set(m, 36, 21, c03*v3);
//gsl_matrix_set(m, 36, 22, c03*v4);
//gsl_matrix_set(m, 36, 23, c03*v5);
//gsl_matrix_set(m, 36, 24, c04*v0);
//gsl_matrix_set(m, 36, 25, c04*v1);
//gsl_matrix_set(m, 36, 26, c04*v2);
//gsl_matrix_set(m, 36, 27, c04*v3);
//gsl_matrix_set(m, 36, 28, c04*v4);
//gsl_matrix_set(m, 36, 29, c04*v5);
//gsl_matrix_set(m, 36, 30, c05*v0);
//gsl_matrix_set(m, 36, 31, c05*v1);
//gsl_matrix_set(m, 36, 32, c05*v2);
//gsl_matrix_set(m, 36, 33, c05*v3);
//gsl_matrix_set(m, 36, 34, c05*v4);
//gsl_matrix_set(m, 36, 35, c05*v5);
//gsl_matrix_set(m, 36, 36, d00 + c00*r00 + c01*r10 + c02*r20 + c03*r30 + c04*r40 + c05*r50);
//gsl_matrix_set(m, 36, 37, d01 + c00*r01 + c01*r11 + c02*r21 + c03*r31 + c04*r41 + c05*r51);
//gsl_matrix_set(m, 36, 38, d02 + c00*r02 + c01*r12 + c02*r22 + c03*r32 + c04*r42 + c05*r52);
//gsl_matrix_set(m, 36, 39, d03 + c00*r03 + c01*r13 + c02*r23 + c03*r33 + c04*r43 + c05*r53);
//gsl_matrix_set(m, 36, 40, d04 + c00*r04 + c01*r14 + c02*r24 + c03*r34 + c04*r44 + c05*r54);
//gsl_matrix_set(m, 36, 41, d05 + c00*r05 + c01*r15 + c02*r25 + c03*r35 + c04*r45 + c05*r55);
//
//
//
//gsl_matrix_set(m, 37, 0, c10*v0);
//gsl_matrix_set(m, 37, 1, c10*v1);
//gsl_matrix_set(m, 37, 2, c10*v2);
//gsl_matrix_set(m, 37, 3, c10*v3);
//gsl_matrix_set(m, 37, 4, c10*v4);
//gsl_matrix_set(m, 37, 5, c10*v5);
//gsl_matrix_set(m, 37, 6, c11*v0);
//gsl_matrix_set(m, 37, 7, c11*v1);
//gsl_matrix_set(m, 37, 8, c11*v2);
//gsl_matrix_set(m, 37, 9, c11*v3);
//gsl_matrix_set(m, 37, 10, c11*v4);
//gsl_matrix_set(m, 37, 11, c11*v5);
//gsl_matrix_set(m, 37, 12, c12*v0);
//gsl_matrix_set(m, 37, 13, c12*v1);
//gsl_matrix_set(m, 37, 14, c12*v2);
//gsl_matrix_set(m, 37, 15, c12*v3);
//gsl_matrix_set(m, 37, 16, c12*v4);
//gsl_matrix_set(m, 37, 17, c12*v5);
//gsl_matrix_set(m, 37, 18, c13*v0);
//gsl_matrix_set(m, 37, 19, c13*v1);
//gsl_matrix_set(m, 37, 20, c13*v2);
//gsl_matrix_set(m, 37, 21, c13*v3);
//gsl_matrix_set(m, 37, 22, c13*v4);
//gsl_matrix_set(m, 37, 23, c13*v5);
//gsl_matrix_set(m, 37, 24, c14*v0);
//gsl_matrix_set(m, 37, 25, c14*v1);
//gsl_matrix_set(m, 37, 26, c14*v2);
//gsl_matrix_set(m, 37, 27, c14*v3);
//gsl_matrix_set(m, 37, 28, c14*v4);
//gsl_matrix_set(m, 37, 29, c14*v5);
//gsl_matrix_set(m, 37, 30, c15*v0);
//gsl_matrix_set(m, 37, 31, c15*v1);
//gsl_matrix_set(m, 37, 32, c15*v2);
//gsl_matrix_set(m, 37, 33, c15*v3);
//gsl_matrix_set(m, 37, 34, c15*v4);
//gsl_matrix_set(m, 37, 35, c15*v5);
//gsl_matrix_set(m, 37, 36, d10 + c10*r00 + c11*r10 + c12*r20 + c13*r30 + c14*r40 + c15*r50);
//gsl_matrix_set(m, 37, 37, d11 + c10*r01 + c11*r11 + c12*r21 + c13*r31 + c14*r41 + c15*r51);
//gsl_matrix_set(m, 37, 38, d12 + c10*r02 + c11*r12 + c12*r22 + c13*r32 + c14*r42 + c15*r52);
//gsl_matrix_set(m, 37, 39, d13 + c10*r03 + c11*r13 + c12*r23 + c13*r33 + c14*r43 + c15*r53);
//gsl_matrix_set(m, 37, 40, d14 + c10*r04 + c11*r14 + c12*r24 + c13*r34 + c14*r44 + c15*r54);
//gsl_matrix_set(m, 37, 41, d15 + c10*r05 + c11*r15 + c12*r25 + c13*r35 + c14*r45 + c15*r55);
//
//
//
//
//gsl_matrix_set(m, 38, 0, c20*v0);
//gsl_matrix_set(m, 38, 1, c20*v1);
//gsl_matrix_set(m, 38, 2, c20*v2);
//gsl_matrix_set(m, 38, 3, c20*v3);
//gsl_matrix_set(m, 38, 4, c20*v4);
//gsl_matrix_set(m, 38, 5, c20*v5);
//gsl_matrix_set(m, 38, 6, c21*v0);
//gsl_matrix_set(m, 38, 7, c21*v1);
//gsl_matrix_set(m, 38, 8, c21*v2);
//gsl_matrix_set(m, 38, 9, c21*v3);
//gsl_matrix_set(m, 38, 10, c21*v4);
//gsl_matrix_set(m, 38, 11, c21*v5);
//gsl_matrix_set(m, 38, 12, c22*v0);
//gsl_matrix_set(m, 38, 13, c22*v1);
//gsl_matrix_set(m, 38, 14, c22*v2);
//gsl_matrix_set(m, 38, 15, c22*v3);
//gsl_matrix_set(m, 38, 16, c22*v4);
//gsl_matrix_set(m, 38, 17, c22*v5);
//gsl_matrix_set(m, 38, 18, c23*v0);
//gsl_matrix_set(m, 38, 19, c23*v1);
//gsl_matrix_set(m, 38, 20, c23*v2);
//gsl_matrix_set(m, 38, 21, c23*v3);
//gsl_matrix_set(m, 38, 22, c23*v4);
//gsl_matrix_set(m, 38, 23, c23*v5);
//gsl_matrix_set(m, 38, 24, c24*v0);
//gsl_matrix_set(m, 38, 25, c24*v1);
//gsl_matrix_set(m, 38, 26, c24*v2);
//gsl_matrix_set(m, 38, 27, c24*v3);
//gsl_matrix_set(m, 38, 28, c24*v4);
//gsl_matrix_set(m, 38, 29, c24*v5);
//gsl_matrix_set(m, 38, 30, c25*v0);
//gsl_matrix_set(m, 38, 31, c25*v1);
//gsl_matrix_set(m, 38, 32, c25*v2);
//gsl_matrix_set(m, 38, 33, c25*v3);
//gsl_matrix_set(m, 38, 34, c25*v4);
//gsl_matrix_set(m, 38, 35, c25*v5);
//gsl_matrix_set(m, 38, 36, d20 + c20*r00 + c21*r10 + c22*r20 + c23*r30 + c24*r40 + c25*r50);
//gsl_matrix_set(m, 38, 37, d21 + c20*r01 + c21*r11 + c22*r21 + c23*r31 + c24*r41 + c25*r51);
//gsl_matrix_set(m, 38, 38, d22 + c20*r02 + c21*r12 + c22*r22 + c23*r32 + c24*r42 + c25*r52);
//gsl_matrix_set(m, 38, 39, d23 + c20*r03 + c21*r13 + c22*r23 + c23*r33 + c24*r43 + c25*r53);
//gsl_matrix_set(m, 38, 40, d24 + c20*r04 + c21*r14 + c22*r24 + c23*r34 + c24*r44 + c25*r54);
//gsl_matrix_set(m, 38, 41, d25 + c20*r05 + c21*r15 + c22*r25 + c23*r35 + c24*r45 + c25*r55);
//
//
//
//
//
//
//
//gsl_matrix_set(m, 39, 0, c30*v0);
//gsl_matrix_set(m, 39, 1, c30*v1);
//gsl_matrix_set(m, 39, 2, c30*v2);
//gsl_matrix_set(m, 39, 3, c30*v3);
//gsl_matrix_set(m, 39, 4, c30*v4);
//gsl_matrix_set(m, 39, 5, c30*v5);
//gsl_matrix_set(m, 39, 6, c31*v0);
//gsl_matrix_set(m, 39, 7, c31*v1);
//gsl_matrix_set(m, 39, 8, c31*v2);
//gsl_matrix_set(m, 39, 9, c31*v3);
//gsl_matrix_set(m, 39, 10, c31*v4);
//gsl_matrix_set(m, 39, 11, c31*v5);
//gsl_matrix_set(m, 39, 12, c32*v0);
//gsl_matrix_set(m, 39, 13, c32*v1);
//gsl_matrix_set(m, 39, 14, c32*v2);
//gsl_matrix_set(m, 39, 15, c32*v3);
//gsl_matrix_set(m, 39, 16, c32*v4);
//gsl_matrix_set(m, 39, 17, c32*v5);
//gsl_matrix_set(m, 39, 18, c33*v0);
//gsl_matrix_set(m, 39, 19, c33*v1);
//gsl_matrix_set(m, 39, 20, c33*v2);
//gsl_matrix_set(m, 39, 21, c33*v3);
//gsl_matrix_set(m, 39, 22, c33*v4);
//gsl_matrix_set(m, 39, 23, c33*v5);
//gsl_matrix_set(m, 39, 24, c34*v0);
//gsl_matrix_set(m, 39, 25, c34*v1);
//gsl_matrix_set(m, 39, 26, c34*v2);
//gsl_matrix_set(m, 39, 27, c34*v3);
//gsl_matrix_set(m, 39, 28, c34*v4);
//gsl_matrix_set(m, 39, 29, c34*v5);
//gsl_matrix_set(m, 39, 30, c35*v0);
//gsl_matrix_set(m, 39, 31, c35*v1);
//gsl_matrix_set(m, 39, 32, c35*v2);
//gsl_matrix_set(m, 39, 33, c35*v3);
//gsl_matrix_set(m, 39, 34, c35*v4);
//gsl_matrix_set(m, 39, 35, c35*v5);
//gsl_matrix_set(m, 39, 36, d30 + c30*r00 + c31*r10 + c32*r20 + c33*r30 + c34*r40 + c35*r50);
//gsl_matrix_set(m, 39, 37, d31 + c30*r01 + c31*r11 + c32*r21 + c33*r31 + c34*r41 + c35*r51);
//gsl_matrix_set(m, 39, 38, d32 + c30*r02 + c31*r12 + c32*r22 + c33*r32 + c34*r42 + c35*r52);
//gsl_matrix_set(m, 39, 39, d33 + c30*r03 + c31*r13 + c32*r23 + c33*r33 + c34*r43 + c35*r53);
//gsl_matrix_set(m, 39, 40, d34 + c30*r04 + c31*r14 + c32*r24 + c33*r34 + c34*r44 + c35*r54);
//gsl_matrix_set(m, 39, 41, d35 + c30*r05 + c31*r15 + c32*r25 + c33*r35 + c34*r45 + c35*r55);
//
//
//
//
//
//
//gsl_matrix_set(m, 40, 0, c40*v0);
//gsl_matrix_set(m, 40, 1, c40*v1);
//gsl_matrix_set(m, 40, 2, c40*v2);
//gsl_matrix_set(m, 40, 3, c40*v3);
//gsl_matrix_set(m, 40, 4, c40*v4);
//gsl_matrix_set(m, 40, 5, c40*v5);
//gsl_matrix_set(m, 40, 6, c41*v0);
//gsl_matrix_set(m, 40, 7, c41*v1);
//gsl_matrix_set(m, 40, 8, c41*v2);
//gsl_matrix_set(m, 40, 9, c41*v3);
//gsl_matrix_set(m, 40, 10, c41*v4);
//gsl_matrix_set(m, 40, 11, c41*v5);
//gsl_matrix_set(m, 40, 12, c42*v0);
//gsl_matrix_set(m, 40, 13, c42*v1);
//gsl_matrix_set(m, 40, 14, c42*v2);
//gsl_matrix_set(m, 40, 15, c42*v3);
//gsl_matrix_set(m, 40, 16, c42*v4);
//gsl_matrix_set(m, 40, 17, c42*v5);
//gsl_matrix_set(m, 40, 18, c43*v0);
//gsl_matrix_set(m, 40, 19, c43*v1);
//gsl_matrix_set(m, 40, 20, c43*v2);
//gsl_matrix_set(m, 40, 21, c43*v3);
//gsl_matrix_set(m, 40, 22, c43*v4);
//gsl_matrix_set(m, 40, 23, c43*v5);
//gsl_matrix_set(m, 40, 24, c44*v0);
//gsl_matrix_set(m, 40, 25, c44*v1);
//gsl_matrix_set(m, 40, 26, c44*v2);
//gsl_matrix_set(m, 40, 27, c44*v3);
//gsl_matrix_set(m, 40, 28, c44*v4);
//gsl_matrix_set(m, 40, 29, c44*v5);
//gsl_matrix_set(m, 40, 30, c45*v0);
//gsl_matrix_set(m, 40, 31, c45*v1);
//gsl_matrix_set(m, 40, 32, c45*v2);
//gsl_matrix_set(m, 40, 33, c45*v3);
//gsl_matrix_set(m, 40, 34, c45*v4);
//gsl_matrix_set(m, 40, 35, c45*v5);
//gsl_matrix_set(m, 40, 36, d40 + c40*r00 + c41*r10 + c42*r20 + c43*r30 + c44*r40 + c45*r50);
//gsl_matrix_set(m, 40, 37, d41 + c40*r01 + c41*r11 + c42*r21 + c43*r31 + c44*r41 + c45*r51);
//gsl_matrix_set(m, 40, 38, d42 + c40*r02 + c41*r12 + c42*r22 + c43*r32 + c44*r42 + c45*r52);
//gsl_matrix_set(m, 40, 39, d43 + c40*r03 + c41*r13 + c42*r23 + c43*r33 + c44*r43 + c45*r53);
//gsl_matrix_set(m, 40, 40, d44 + c40*r04 + c41*r14 + c42*r24 + c43*r34 + c44*r44 + c45*r54);
//gsl_matrix_set(m, 40, 41, d45 + c40*r05 + c41*r15 + c42*r25 + c43*r35 + c44*r45 + c45*r55);
//
//gsl_matrix_set(m, 41, 0, c50*v0);
//gsl_matrix_set(m, 41, 1, c50*v1);
//gsl_matrix_set(m, 41, 2, c50*v2);
//gsl_matrix_set(m, 41, 3, c50*v3);
//gsl_matrix_set(m, 41, 4, c50*v4);
//gsl_matrix_set(m, 41, 5, c50*v5);
//gsl_matrix_set(m, 41, 6, c51*v0);
//gsl_matrix_set(m, 41, 7, c51*v1);
//gsl_matrix_set(m, 41, 8, c51*v2);
//gsl_matrix_set(m, 41, 9, c51*v3);
//gsl_matrix_set(m, 41, 10, c51*v4);
//gsl_matrix_set(m, 41, 11, c51*v5);
//gsl_matrix_set(m, 41, 12, c52*v0);
//gsl_matrix_set(m, 41, 13, c52*v1);
//gsl_matrix_set(m, 41, 14, c52*v2);
//gsl_matrix_set(m, 41, 15, c52*v3);
//gsl_matrix_set(m, 41, 16, c52*v4);
//gsl_matrix_set(m, 41, 17, c52*v5);
//gsl_matrix_set(m, 41, 18, c53*v0);
//gsl_matrix_set(m, 41, 19, c53*v1);
//gsl_matrix_set(m, 41, 20, c53*v2);
//gsl_matrix_set(m, 41, 21, c53*v3);
//gsl_matrix_set(m, 41, 22, c53*v4);
//gsl_matrix_set(m, 41, 23, c53*v5);
//gsl_matrix_set(m, 41, 24, c54*v0);
//gsl_matrix_set(m, 41, 25, c54*v1);
//gsl_matrix_set(m, 41, 26, c54*v2);
//gsl_matrix_set(m, 41, 27, c54*v3);
//gsl_matrix_set(m, 41, 28, c54*v4);
//gsl_matrix_set(m, 41, 29, c54*v5);
//gsl_matrix_set(m, 41, 30, c55*v0);
//gsl_matrix_set(m, 41, 31, c55*v1);
//gsl_matrix_set(m, 41, 32, c55*v2);
//gsl_matrix_set(m, 41, 33, c55*v3);
//gsl_matrix_set(m, 41, 34, c55*v4);
//gsl_matrix_set(m, 41, 35, c55*v5);
//gsl_matrix_set(m, 41, 36, d50 + c50*r00 + c51*r10 + c52*r20 + c53*r30 + c54*r40 + c55*r50);
//gsl_matrix_set(m, 41, 37, d51 + c50*r01 + c51*r11 + c52*r21 + c53*r31 + c54*r41 + c55*r51);
//gsl_matrix_set(m, 41, 38, d52 + c50*r02 + c51*r12 + c52*r22 + c53*r32 + c54*r42 + c55*r52);
//gsl_matrix_set(m, 41, 39, d53 + c50*r03 + c51*r13 + c52*r23 + c53*r33 + c54*r43 + c55*r53);
//gsl_matrix_set(m, 41, 40, d54 + c50*r04 + c51*r14 + c52*r24 + c53*r34 + c54*r44 + c55*r54);
//gsl_matrix_set(m, 41, 41, d55 + c50*r05 + c51*r15 + c52*r25 + c53*r35 + c54*r45 + c55*r55);
//


