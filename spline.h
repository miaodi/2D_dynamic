#ifndef _spline_h
#define _spline_h
#endif


namespace spline {
    void DersBasisFuns(int i, double u, int p, double* U, int n, double** & ders);
    int Findspan(int m, int p, double* U, double u);
    double DersOneBasisFUn(int p, int m, double* U, int i, double u, int n);
    static double x3[] = { -0.77459666924148340427791, 0.00000000000000000000000, 0.77459666924148340427791	};
    static double w3[] = { 0.55555555555555546920488, 0.88888888888888895056795, 0.55555555555555546920488	};

    static double x4[] = { -0.86113631159405257253781, - 0.33998104358485625731134, 0.33998104358485625731134, 0.86113631159405257253781	};
    static double w4[] = { 0.34785484513745379420158,	0.65214515486254631682073, 0.65214515486254631682073, 0.34785484513745379420158	};

    static double x5[] = { -0.90617984593866396370032, - 0.53846931010568310771447, 0.00000000000000000000000, 0.53846931010568310771447, 0.90617984593866396370032	};
    static double w5[] = { 0.23692688505618964001087, 0.47862867049936608232485, 0.56888888888888855532855, 0.47862867049936608232485, 0.23692688505618964001087	};

    static double x6[] = { -0.93246951420315205005807, - 0.66120938646626448154109, - 0.23861918608319693246855, 0.23861918608319693246855, 0.66120938646626448154109, 0.93246951420315205005807	};
    static double w6[] = { 0.17132449237916899664746, 0.36076157304813905035701, 0.46791393457269181421765, 0.46791393457269181421765, 0.36076157304813905035701, 0.17132449237916899664746	};

    static double x7[] = { -0.94910791234275848626822, - 0.74153118559939446008400, - 0.40584515137739718415588, 0.00000000000000000000000, 0.40584515137739718415588, 0.74153118559939446008400, 0.94910791234275848626822	};
    static double w7[] = { 0.12948496616887075760793, 0.27970539148927631156738, 0.38183005050511836797611, 0.41795918367346879263025, 0.38183005050511836797611, 0.27970539148927631156738,	0.12948496616887075760793 };

    static double x8[] = { -0.96028985649753628717207,
                           - 0.79666647741362672796583,
                           - 0.52553240991632899081765,
                           - 0.18343464249564980783624,
                           0.18343464249564980783624,
                           0.52553240991632899081765,
                           0.79666647741362672796583,
                           0.96028985649753628717207
    };
    static double w8[] = { 0.10122853629037641132182,
                           0.22238103445337426000705,
                           0.31370664587788738009166,
                           0.36268378337836187919052,
                           0.36268378337836187919052,
                           0.31370664587788738009166,
                           0.22238103445337426000705,
                           0.10122853629037641132182
    };

    static double x9[] = { -0.96816023950762608585308,
                           - 0.83603110732663576953883,
                           - 0.61337143270059035771169,
                           - 0.32425342340380891581475,
                           0.00000000000000000000000,
                           0.32425342340380891581475,
                           0.61337143270059035771169,
                           0.83603110732663576953883,
                           0.96816023950762608585308
    };
    static double w9[] = { 0.08127438836157502288771,
                           0.18064816069485722938026,
                           0.26061069640293543780984,
                           0.31234707704000252981302,
                           0.33023935500125939368488,
                           0.31234707704000252981302,
                           0.26061069640293543780984,
                           0.18064816069485722938026,
                           0.08127438836157502288771
    };
}