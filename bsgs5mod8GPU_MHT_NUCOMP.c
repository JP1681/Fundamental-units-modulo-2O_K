long dwr(long a, long b){
    return (a - (((a % b) + b) % b)) / b;
};
long floorintsqrt(long D){
    long floorRootD = floor(sqrt((float)D));
    if (floorRootD*floorRootD>D){
        return floorRootD - 1;
    } else if ((floorRootD+1)*(floorRootD+1)<D){
        return floorRootD + 1;
    } else {
        return floorRootD;
    }
}
void extgcd(long a, long b, long *o_s, long *o_t, long *o_r){
    long old_r = a;
    long r = b;
    long old_s = 1;
    long s = 0;
    long old_t = 0;
    long t = 1;
    while (r!=0){
        long q = dwr(old_r,r);
        long tempr = r;
        r = old_r - q*r;
        old_r = tempr;
        long temps = s;
        s = old_s - q*s;
        old_s = temps;
        long tempt = t;
        t = old_t - q*t;
        old_t = tempt;
    };
    *o_s = old_s;
    *o_t = old_t;
    *o_r = old_r;
}

long mod(long a, long b){
    return ((a%b)+b)%b;
}
void NUCOMP(long D, long Q1, long P1, long Q2, long P2, long *newQprime, long *newPprime, long *A, long *B, long *C, int *check){
    long floorRootD = floorintsqrt(D);
    int r = 2;
    if (Q1<=Q2){
        long tempQ = Q2;
        Q2 = Q1;
        Q1 = tempQ;
        long tempP = P2;
        P2 = P1;
        P1 = tempP;
    };
    long Q1onr = dwr(Q1,r);
    long Q2onr = dwr(Q2,r);
    long s;
    long X;
    long G;
    extgcd(Q1onr,Q2onr,&s,&X,&G);
    X = mod(X,Q1onr);
    long Y;
    long Z;
    long S;
    extgcd(dwr(P1+P2,r),G,&Y,&Z,&S);
    long R2 = dwr(D-P2*P2,Q2);
    long Q1onS = dwr(Q1,S);
    long U = mod((X * Z * (P1-P2) + Y * R2),Q1onS);
    long Rm1 = U;
    long R0 = Q1onS;
    long Cm1 = -1;
    long C0 = 0;
    long newQ;
    long newP;
    long NewQprime;
    long NewPprime;
    long Bm1;
    long B0;
    long q;
    long newC;
    long newR;
    long M1;
    long M2;
    long newB;
    int i = -1;
    long Q2onrS = dwr(Q2onr,S);
    long f = floorintsqrt(4*floorintsqrt(D));
    //*check = f;
    //return;
    if (R0<=f){
        newQ = dwr(Q1*Q2,r*S*S);
        newP = mod(P2+U*Q2onrS,newQ);
    } else {
        while (R0>f){
            i++;
            if (i>100){
                *check = -996;
                return;
            }
            q = dwr(Rm1,R0);
            newC = Cm1 - q*C0;
            Cm1 = C0;
            C0 = newC;
            newR = Rm1 - q*R0;
            Rm1 = R0;
            R0 = newR;
            
        };
        M1 = dwr(Q2onrS*R0+(P1-P2)*C0,Q1onS);
        M2 = dwr((P1+P2)*R0+r*S*R2*C0,Q1onS);
        int power;
        if (mod(i+1,2)==0){
            power = 1;
        } else {
            power = -1;
        }
        newQ = power*(R0*M1-C0*M2);
        newP = dwr(Q2onrS*R0+newQ*Cm1,C0)-P2;
    };
    NewQprime = abs(newQ);
    long k = dwr(floorRootD-newP,NewQprime);
    NewPprime = k * NewQprime + newP;
    Bm1 = abs(Cm1);
    B0 = abs(C0);
    int l = 0;
    while (true){
        if (NewPprime + floorRootD >= NewQprime){
            *newQprime = NewQprime;
            *newPprime = mod(NewPprime,NewQprime);
            *A = S*(newQ*Bm1+newP*B0);
            *B = -1*S*B0;
            *C = newQ;
            return;
        }
        l++;
        if (l>300){
            *check = -995;
            return;
        }
        q = dwr(newP+floorRootD,NewQprime);
        newP = q * NewQprime - newP;
        newQ = dwr(D-newP*newP,NewQprime);
        NewQprime = abs(newQ);
        k = dwr(floorRootD-newP,NewQprime);
        NewPprime = k * NewQprime + newP;
        newB = q * B0 + Bm1;
        Bm1 = mod(B0,4);
        B0 = mod(newB,4);
    }
}

struct QPTtuple {
    long Q;
    long P;
    int theta;
};

long gcd(long a, long b){
    long old_r = a;
    long r = b;
    while (r!=0){
        long q = dwr(old_r,r);
        long tempr = r;
        r = old_r - q*r;
        old_r = tempr;
    }
    return old_r;
}
int Ccoset5M8(long x, long y, long z, bool s){
    if (!s){
        long f = gcd(x,gcd(y,z));
        if (z<0){
            f = -1*f;
        };
        x = dwr(x,f);
        y = dwr(y,f);
        z = dwr(z,f);
    };
    if (z==1){
        x = 2*x;
        y = 2*y;
    };
    if (mod(x,4)==mod(y,4)){
        if (mod(y,2)==0){
            return -1;
        } else {
            return 1;
        }
    } else {
        if (mod(y,2)==0){
            return 0;
        } else {
            return 2;
        }
    };
}
int CFFcoset5M8(long x, long y, long z){
    if (z<0){
        x=-x;y=-y,z=-z;
    };
    if (mod(z,4)==0){
        return -1;
    };
    if (mod(z,2)==0){
        return Ccoset5M8(x,y,2,true);
    }
    else {
        return Ccoset5M8(x,y,1,true);
    };
}
int reduceQPfull(long D, long Q, long P){
    long floorRootD = floorintsqrt(D);
    int sD = 2;
    long q0 = dwr(P+floorRootD,Q);
    int i = 1;
    int Q0 = Q;
    long Bm1 = 0;
    long B0 = 0;
    long B1 = 1;
    long newB;
    while (i==1 || Q!=sD){
        P = q0*Q-P;
        Q = dwr(D-P*P,Q);
        q0 = dwr(P+floorRootD,Q);
        long tempB0 = B0;
        long tempB1 = B1;
        newB = mod(q0*B1+B0,4);
        Bm1 = tempB0;
        B0 = tempB1;
        B1 = newB;
        i++;
    }
    long G1 = mod(P*B0+Q*Bm1,4);
    int theta = Ccoset5M8(G1,B0,Q0,false);
    return theta;
}
void reduceQPbaby(long D, long Q, long P, struct QPTtuple CQPdict[], struct QPTtuple spareDict[], int* spareDictPos,  long correctIdeal[], int* cReducedTheta, bool* endEarlyFlag, int* BS){
    float rootD = sqrt((float)D);
    long floorRootD = floorintsqrt(D);
    long q0 = dwr(P+floorRootD,Q);
    int i = 1;
    bool flag = false;
    int c = 0;
    float theta = 0;
    long Q0 = Q;
    long RBm2 = 0;
    long RBm1 = 0;
    long RB0 = 0;
    long RB1 = 1;
    long RG0 = mod(Q,4);
    long RG1 = mod(Q*q0-P,4);
    long lastP, lastQ;
    int cyclicTheta;
    int sdp = *spareDictPos;
    float fourthRootD = pow((float)D,(float)1/(float)4);

    unsigned int mask = (1 << $nBits) - 1;
    unsigned int qmod = mask & Q;
    if (CQPdict[qmod].Q==0){
        CQPdict[qmod].Q = Q;
        CQPdict[qmod].P = mod(P,Q);
        CQPdict[qmod].theta = 0;
    } else {
        spareDict[sdp].Q = Q;
        spareDict[sdp].P = mod(P,Q);
        spareDict[sdp].theta = 0;
        (sdp)++;
    }
    while (i==1 || c<3){
        lastP = P;
        P = q0*Q-P;
        lastQ = Q;
        Q = dwr(D-P*P,Q);
        q0 = dwr(P+floorRootD,Q);
        RG0 = RG1;
        RG1 = mod(P*RB1+Q*RB0,4);
        RBm2 = RBm1;
        RBm1 = RB0;
        RB0 = RB1;
        RB1 = mod(q0*RB0+RBm1,4);
        if (P==lastP){
            *endEarlyFlag = true;
            *cReducedTheta = Ccoset5M8(RG1*RBm1+RG0*RBm2,RBm1*(RB0+RBm2),Q0,false);
            break;
        };
        if (Q==lastQ){
            *endEarlyFlag = true;
            *cReducedTheta = Ccoset5M8(RG1*RB0+RG0*RBm1,RB0*RB0+RBm1*RBm1,Q0,false);
            break;
        };
        cyclicTheta = Ccoset5M8(RG1,RB0,Q0,false);
        i++;

        unsigned int mask = (1 << $nBits) - 1;
        unsigned int qmod = mask & Q;
        if (CQPdict[qmod].Q==0){
            CQPdict[qmod].Q = Q;
            CQPdict[qmod].P = mod(P,Q);
            CQPdict[qmod].theta = cyclicTheta;
        } else {
            spareDict[sdp].Q = Q;
            spareDict[sdp].P = mod(P,Q);
            spareDict[sdp].theta = cyclicTheta;
            (sdp)++;
        }
        theta = theta + native_log((float)(P+rootD)/(float)lastQ);
        if (theta>fourthRootD){
            flag = true;
        };
        if (flag){
            c++;
        };
        if (c==1){
            correctIdeal[0] = Q;
            correctIdeal[1] = mod(P,Q);
            *cReducedTheta = cyclicTheta;
        };
    }
    *BS = i;
    *spareDictPos = sdp;
}
int checkContains(long Q, long P, struct QPTtuple cqpdict[], struct QPTtuple spareDict[], int spareDictPos){
    unsigned int mask = (1 << $nBits) - 1;
    unsigned int qmod = mask & Q;
    long qval = cqpdict[qmod].Q;
    if (qval==0){
        return -1;
    } else if (qval==Q){
        if (cqpdict[qmod].P == P){
            return cqpdict[qmod].theta;
        };
    };
    for (int i = 0; i<spareDictPos; ++i){
        if (spareDict[i].Q==Q){
            if (spareDict[i].P==P){
                return spareDict[i].theta;
            }
        }
    }
    return -1;
}
void find_FUnit(long D, int *FU, int *BS, int *GS, int *SDS){
    long P = 1;
    long Q = 2;
    const int N = $mainDictSize;
    const int spareN = $spareDictSize;
    struct QPTtuple Cdict[N];
    struct QPTtuple spareDict[spareN];
    int spareDictPos = 0;
    for (int i=0; i<N; i++){
        Cdict[i].Q = 0;
    };
    long firstQP[2];
    int Ctheta;
    bool eef = false;
    
    reduceQPbaby(D,Q,P,Cdict,spareDict,&spareDictPos,firstQP,&Ctheta,&eef,BS);
    *SDS = spareDictPos;
    if (spareDictPos>=spareN){
        *FU = -999;
        return;
    };
    if (eef){
        *FU = Ctheta;
        return;
    };
    int firstCTheta = Ctheta;
    long firstQ, firstP;
    firstQ = firstQP[0];
    firstP = firstQP[1];
    long newQ = firstQ;
    long newP = firstP;
    int i = 0;
    int cc = -1;
    int Cd;
    long S;
    int c;
    int check=0;
    long testA;
    long testB;
    long testC;
    while (i==0 || cc == -1){
        NUCOMP(D,newQ,newP,firstQ,firstP,&newQ,&newP,&testA,&testB,&testC,&check);
        if (check<0){
            *FU = check;
            return;
        };
        Cd = CFFcoset5M8(testA,testB,testC);
        Ctheta = mod(firstCTheta+Ctheta-Cd,3);
        i++;
        cc = checkContains(newQ,mod(newP,newQ),Cdict,spareDict,spareDictPos);
        if (i>$giantStepMax){
            *FU = -998;
            *GS = i;
            return;
        };
    };
    int Cfunit = mod(Ctheta - cc,3);
    *FU = Cfunit;
    *GS = i;
    return;    
};
__kernel void sum_array(__global const long *array, __global int *result, __global int *babystep, __global int *giantstep, __global int *sparedictsize) {
    int gid = get_global_id(0);
    long D = array[gid];
    int FU;
    int BS=0;
    int GS=0;
    int SDS=0;
    find_FUnit(D,&FU,&BS,&GS,&SDS);
    result[gid] = FU;
    babystep[gid] = BS;
    giantstep[gid] = GS;
    sparedictsize[gid] = SDS;
};