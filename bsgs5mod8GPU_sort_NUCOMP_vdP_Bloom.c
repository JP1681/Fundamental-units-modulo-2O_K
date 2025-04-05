// NOTE: throughout this code, it is assumed that we are working in the real quadratic number field with discriminant
// a positive integer D that is squarefree, and equal to 5 modulo 8. We are also reducing the fundamental unit modulo 2O_K.
// Other values of D might not function correctly.
long dwr(long a, long b){
    // Floor division and modulo are implemented differently in C to in Python.
    // This function simulates the Python operation.
    return (a - (((a % b) + b) % b)) / b;
};
long floorintsqrt(long D){
    // If D is close to a square number, sometimes the sqrt function will return a value that is off by 1.
    // This function checks whether the value is wrong and fixes it if it is.
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
    // This is (roughly) the standard extended Euclidean algorithm.
    // The code is modified from the pseudocode in the Wikipedia article, see the following link.
    // https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Pseudocode
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
    if (old_r < 0){
        old_r = old_r * -1;
        old_s = old_s * -1;
        old_t = old_t * -1;
    }
    *o_s = old_s;
    *o_t = old_t;
    *o_r = old_r;
}
long mod(long a, long b){
    // Floor division and modulo are implemented differently in C to in Python.
    // This function simulates the Python operation.
    return ((a%b)+b)%b;
}
void abcsolve(long a, long b, long c, long *x1, long *x2, long *x3, long *g){
    // Solves the equation a*x1 + b*x2 + c*x3 = g, where g is gcd(a,b,c).
    long n,m,G1,x,Y,G2;
    extgcd(a,b,&n,&m,&G1);
    extgcd(G1,c,&x,&Y,&G2);
    *x1 = x*n;
    *x2 = x*m;
    *x3 = Y;
    *g = G2;
}
void idealProduct(long D, long Q1, long P1, long Q2, long P2, long *S, long *Q, long *P){
    // Computes the ideal product (Q1,P1) * (Q2,P2) without reduction, in the real quadratic field of discriminant D.
    // This algorithm is taken from Section 5.4 of Solving the Pell Equation, by Jacobson and Williams. MR2466979.
    long r = 2;
    long Q1onr = dwr(Q1,r);
    long Q2onr = dwr(Q2,r);
    long P12onr = dwr(P1+P2,r);
    long V,W,Y;
    abcsolve(Q1onr,Q2onr,P12onr,&V,&W,&Y,&(*S));
    long Q1onS = dwr(Q1,*S);
    long Q2onrS = dwr(Q2onr,*S);
    *Q = Q1onS*Q2onrS;
    long R2 = dwr(D-P2*P2,Q2);
    long U = (W*(P1-P2)+Y*R2);
    *P = mod(P2+U*Q2onrS,*Q);
}
void NUDUPLJvdP(long D, long u, long v, long w, long *u3, long *v3, long *a, long *b, long *c){
    // Computes a near-reduced double, (u3,v3,w3) of a binary quadratic form, (u,v,w) with discriminant D. Also computes the
    // factor induced by reduction, \frac{1}{\gamma}, where \gamma is \frac{a+b\sqrt{D}}{c}.
    // This algorithm is the NUDUPL algorithm described in Computational aspects of NUCOMP, by 
    // Jacobson and van der Poorten, but computes \gamma as in the algorithm from A note on NUCOMP, by van der Poorten.
    long L,x,y,G,Ax,By,Dy,Bx,bx,by,z,q,t,ax,ay,dx,U3,V3,w3,Q1,dy;
    L = floorintsqrt(floorintsqrt(dwr(D,4)));
    //Step 1
    extgcd(u,v,&x,&y,&G);
    Ax = G;
    By = dwr(u,G);
    Dy = dwr(v,G);
    //Step 2
    Bx = mod(y*w,By);
    //Step 3
    bx = Bx;
    by = By;
    //Step 3a
    x = 1;
    y = 0;
    z = 0;
    while (abs(by)>L & bx!=0){
        //Step 3c
        q = dwr(by,bx);
        t = mod(by,bx);
        by = bx;
        bx = t;
        t = y - q*x;
        y = x;
        x = t;
        z = z + 1;
    }
    //Step 3b
    if (mod(z,2)!=0){
        by = by * -1;
        y = y * -1;
    }
    ax = G*x;
    ay = G*y;
    //Step 4
    if (z==0){
        dx = dwr(bx * Dy - w,By);
        U3 = by*by;
        w3 = bx*bx;
        V3 = v - (bx + by)*(bx + by) + U3 + w3;
        w3 = w3 - G*dx;
    } else {
        dx = dwr(bx * Dy - w * x,By);
        Q1 = dx * y;
        //*a = dy;
        //return;
        dy = Q1 + Dy;
        V3 = G * (dy + Q1);
        dy = dwr(dy,x);
        U3 = by*by;
        w3 = bx*bx;
        V3 = V3 - (bx+by)*(bx+by) + U3 + w3;
        U3 = U3 - ay * dy;
        w3 = w3 - ax * dx;
    }
    t = 2*U3;
    *u3 = U3;
    *v3 = V3;
    *a = G*(x*t+y*V3);
    *b = G*y;
    *c = t;
    return;
}
void NUCOMPJvdP(long D,long u1, long v1, long w1, long u2, long v2, long w2, long *u3, long *v3, long *A, long *B, long *C){
    // Computes a near-reduced composite, (u3,v3,w3) of two binary quadratic forms, (u1,v1,w1) and (u2,v2,w2), both with discriminant D.
    // Also computes the factor induced by reduction, \frac{1}{\gamma}, where \gamma is \frac{a+b\sqrt{D}}{c}.
    // This algorithm is the NUCOMP algorithm described in Computational aspects of NUCOMP, by 
    // Jacobson and van der Poorten, but computes \gamma as in the algorithm from A note on NUCOMP, by van der Poorten.
    long Ax,Bx,By,Cy,Dy,x,y,z,G,H,l,bx,by,q,t,ax,ay,Q1,Q2,Q3,Q4,cx,dx,temp,L,s,m,F,w3,cy,dy,U3,V3,b,c;
    if (w1<w2){
        temp = u1;
        u1 = u2;
        u2 = temp;
        temp = v1;
        v1 = v2;
        v2 = temp;
        temp = w1;
        w1 = w2;
        w2 = temp;
    }
    L = floorintsqrt(floorintsqrt(dwr(D,4)));
    //Step 1
    s = dwr(v1+v2,2);
    m = v2 - s;
    //Step 2
    extgcd(u2,u1,&b,&c,&F);

    if (mod(s,F)==0){
        G = F;
        Ax = G;
        Bx = m * b;
        By = dwr(u1,G);
        Cy = dwr(u2,G);
        Dy = dwr(s,G);
    } else { //Step 3
        extgcd(F,s,&x,&y,&G);
        H = dwr(F,G);
        By = dwr(u1,G);
        Cy = dwr(u2,G);
        Dy = dwr(s,G);
        //Step 4
        l = mod(y*(b*w1+c*w2),H);
        Bx = dwr(b*m+l*By,H);
    }
    //Step 5
    bx = mod(Bx,By);
    by = By;
    //Step 5a
    x = 1;
    y = 0;
    z = 0;
    while (abs(by)>L & bx!=0){
        //Step 5c
        q = dwr(by,bx);
        t = mod(by,bx);
        by = bx;
        bx = t;
        t = y - q*x;
        y = x;
        x = t;
        z = z + 1;
    }
    //Step 5b
    if (mod(z,2)!=0){
        by = by * -1;
        y = y * -1;
    }
    ax = G*x;
    ay = G*y;
    //Step 6
    if (z==0){
        Q1 = Cy * bx;
        cx = dwr(Q1-m,By);
        dx = dwr(bx*Dy-w2,By);
        U3 = by*Cy;
        w3 = bx*cx - G*dx;
        V3 = v2 - 2*Q1;
    } else {
        cx = dwr(Cy * bx - m * x,By);
        Q1 = by * cx;
        Q2 = Q1 + m;
        dx = dwr(Dy * bx - w2 * x,By);
        Q3 = y * dx;
        Q4 = Q3 + Dy;
        dy = dwr(Q4,x);
        if (bx != 0){
            cy = dwr(Q2,bx);
        } else {
            cy = dwr(cx * dy - w1,dx);
        }
        U3 = by * cy - ay * dy;
        w3 = bx * cx - ax * dx;
        V3 = G * (Q3+Q4) - Q1 - Q2;
    }
    *u3 = U3;
    *v3 = V3;
    t = 2 * U3;
    *A = G*(x*t+y*V3);
    *B = G*y;
    *C = t;
    return;
}
void NUCOMPchoose(long D, long u1, long v1, long u2, long v2, long *u, long *v, long *a, long *b, long *c){
    // This algorithm takes two ideals, (u1,v1) and (u2,v2). If the norm of either ideal is small (i.e. <=50), it is preferable to
    // use the ideal product algorithm to get a composite. Otherwise, the ideals are converted to their equivalent
    // binary quadratic forms, and passed to NUDUPL (if they are the same) or NUCOMP (if they are different). The resulting composite
    // binary quadratic form is converted back to an ideal and returned, along with the factor induced by reduction during NUCOMP/NUDUPL.
    long w1,w2,u3,v3,A,B,C,U;
    u1 = abs(u1);
    u2 = abs(u2);
    v1 = mod(v1,u1);
    v2 = mod(v2,u2);
    w1 = dwr(v1*v1-D,2*u1);
    if (u1<=50 | u2<= 50){
        long S,Q,P;
        idealProduct(D,u1,v1,u2,v2,&S,&Q,&P);
        *u = Q;
        *v = P;
        *a = 1;
        *b = 0;
        *c = 1;
        return;
    }
    if (u1==u2 & v1==v2){
        NUDUPLJvdP(D,dwr(u1,2),-1*v1,w1,&u3,&v3,&A,&B,&C);
    } else {
        w2 = dwr(v2*v2-D,2*u2);
        NUCOMPJvdP(D,dwr(u1,2),-1*v1,w1,dwr(u2,2),-1*v2,w2,&u3,&v3,&A,&B,&C);
    }
    U = abs(2*u3);
    *u = U;
    *v = mod(-1*v3,U);
    *a = A;
    *b = B;
    *c = C;
    return;
}
void NUCOMP(long D, long Q1, long P1, long Q2, long P2, long *newQprime, long *newPprime, long *A, long *B, long *C, int *check){
    // This algorithm does not seem to work correctly for large values of D.
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
    // Each ideal is stored as a struct, containing the values Q and P such that the ideal is [\frac{Q}{2}, \frac{P+\sqrt{D}}{2}].
    // For brevity we notate this as (Q,P). theta is the ideal's generator, reduced mod 2O_K.
    long Q;
    long P;
    int theta;
};
long gcd(long a, long b){
    // The standard Euclidean algorithm for gcd.
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
    // This determines which coset of 2O_K the values \frac{x+y\sqrt{D}}{z} is in.
    // If it is in the same coset as 0, return -1 (this should never happen).
    // If in the same coset as 1, return 0. If in the same coset as \frac{1+\sqrt{D}}{2}, return 1.
    // If in the same coset as \frac{3+\sqrt{D}}{2}, return 2.
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
    // Does the same as Ccoset5M8, but is extended to the localisation of O_K at the prime 2.
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
void reduceQPred(long D,long Q,long P,long *q, long *p, int *ctheta){
    // Given an ideal (Q,P), computes a reduced ideal (q,p) equivalent to it. ctheta is the ratio of their generators, reduced mod 2O_K.
    // This algorithm is modified from the reduction algorithm described in Section 5.2 and Section 3.1 of Solving the Pell Equation,
    // by Jacobson and Williams. MR2466979.
    long floorRootD = floorintsqrt(D);
    long floorRootD32 = floorintsqrt(dwr(D*9,4));
    int sD = 2;
    long q0 = dwr(P+floorRootD,Q);
    int i = 1;
    bool flag = false;
    int c = 0;
    int Ctheta = 0;
    int Cpsi;
    while (i==1 | c<1){
        if (Q<floorRootD32 & 0<Q){
            flag = true;
        }
        if (flag){
            c = c + 1;
        }
        P = q0 * Q - P;
        Cpsi = CFFcoset5M8(P,1,Q);
        Q = dwr(D-P*P,Q);
        Ctheta = mod(Ctheta+Cpsi,3);
        q0 = dwr(P+floorRootD,Q);
        i = i + 1;
    }
    *q = Q;
    *p = P;
    *ctheta = Ctheta;
}
int reduceQPfull(long D, long Q, long P){
    // This function is not used in the code but is helpful for checking. The function reduceQPfull(D,2,1)
    // will give the reduced fundamental unit, using continued fractions only (no infrastructure). Adapted from 
    // the algorithm on page 59 of Solving the Pell Equation, by Jacobson and Williams. MR2466979.
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

// The next three functions are hashes of an ideal (i.e. of the values Q and P). These functions were written by GPT-4.
uint hash1(long Q, long P) {
    return (uint)((Q ^ P) & ($BLOOM_SIZE - 1));
}

uint hash2(long Q, long P) {
    return (uint)(((Q * 31) + P) & ($BLOOM_SIZE - 1));
}

uint hash3(long Q, long P) {
    return (uint)((Q + P * 17) & ($BLOOM_SIZE - 1));
}
void bloom_insert(uchar* bloom, long Q, long P) {
    // Inserts an ideal (Q,P) into the Bloom filter. This function was written by GPT-4.
    uint h1 = hash1(Q, P);
    uint h2 = hash2(Q, P);
    uint h3 = hash3(Q, P);

    bloom[h1 / 8] |= (1 << (h1 % 8));
    bloom[h2 / 8] |= (1 << (h2 % 8));
    bloom[h3 / 8] |= (1 << (h3 % 8));
}
int bloom_maybe_contains(uchar* bloom, long Q, long P) {
    // This function checks (probabilistically) whether a given ideal is in the Bloom filter. If the answer is no,
    // then the ideal is definitely not there. If the answer is yes, the list of ideals must be checked to confirm this.
    // This function was written by GPT-4.
    uint h1 = hash1(Q, P);
    uint h2 = hash2(Q, P);
    uint h3 = hash3(Q, P);

    return ((bloom[h1 / 8] & (1 << (h1 % 8))) &&
            (bloom[h2 / 8] & (1 << (h2 % 8))) &&
            (bloom[h3 / 8] & (1 << (h3 % 8))));
}

void reduceQPbaby(long D, long Q, long P, struct QPTtuple CQPdict[], uchar bloom[], long correctIdeal[], int* cReducedTheta, bool* endEarlyFlag, int* BS){
    // This function constructs the baby step list of ideals as described in Section 7.4 of 
    // Solving the Pell Equation, by Jacobson and Williams. MR2466979. It will return this list, along with a Bloom filter for the ideals,
    // the starting ideal \mathfrak{b}_1 (correctIdeal) and its generator (cReducedTheta). If the fundamental unit is found without
    // proceeding to the baby step stage (see page 176 of Solving the Pell Equation), then we return with endEarlyFlag marked as true,
    // along with the reduced unit (cReducedTheta). The number of baby steps is also returned (BS).
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
    float fourthRootD = pow((float)D,(float)1/(float)4);
    int dictIndex = 0;
    struct QPTtuple new;
    new.Q = Q; new.P=mod(P,Q); new.theta=0;
    bloom_insert(bloom,Q,mod(P,Q));
    CQPdict[dictIndex] = new;
    dictIndex++;
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
        new.Q = Q; new.P = mod(P,Q); new.theta = cyclicTheta;
        bloom_insert(bloom,Q,mod(P,Q));
        CQPdict[dictIndex] = new;
        dictIndex++;
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
    *BS = dictIndex;
}

int is_less(struct QPTtuple a, struct QPTtuple b) {
    // This function determines whether an ideal is 'less' than another (for sorting purposes). We sort by Q, then by P.
    // This function was written by GPT-4.
    if (a.Q < b.Q) return 1;
    if (a.Q > b.Q) return 0;
    // Qs are equal, compare P
    return a.P < b.P;
}

void quicksort(struct QPTtuple *data, int length) {
    // This function sorts the list of ideals using quicksort. Without a Bloom filter, this is much faster than a simple search on the
    // unsorted list of ideals. But with a Bloom filter, most of the requests to search the baby step list are quickly dismissed, and
    // so the cost of sorting the list outweighs the benefit of using binary search. At present, this function is not used.
    // This function was written by GPT-4.
    int stack[600];
    int top = -1;

    stack[++top] = 0;
    stack[++top] = length - 1;

    while (top >= 0) {
        int high = stack[top--];
        int low = stack[top--];

        struct QPTtuple pivot = data[high];
        int i = low - 1;

        for (int j = low; j < high; j++) {
            if (is_less(data[j], pivot)) {
                i++;
                // Swap data[i] and data[j]
                struct QPTtuple temp = data[i];
                data[i] = data[j];
                data[j] = temp;
            }
        }

        // Put pivot in correct place
        struct QPTtuple temp = data[i + 1];
        data[i + 1] = data[high];
        data[high] = temp;

        int p = i + 1;

        if (p - 1 > low) {
            stack[++top] = low;
            stack[++top] = p - 1;
        }

        if (p + 1 < high) {
            stack[++top] = p + 1;
            stack[++top] = high;
        }
    }
}
int checkContains(long Q, long P, struct QPTtuple *data, int length) {
    // Checks whether a sorted list of ideals contains the ideal (Q,P), using binary search.
    // This functions was written by GPT-4 and is presently not used (see the comments for quicksort).
    int low = 0;
    int high = length - 1;

    while (low <= high) {
        int mid = (low + high) / 2;
        long midQ = data[mid].Q;
        long midP = data[mid].P;

        if (midQ == Q && midP == P) {
            return data[mid].theta;
        }

        if (midQ < Q || (midQ == Q && midP < P)) {
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }

    return -1; // Not found
}

int checkContainsUnsorted(long Q, long P, struct QPTtuple *data, int length) {
    // Checks whether an unsorted list of ideals contains the ideal (Q,P), using linear search.
    // This function was written by GPT-4.
    for (int i = 0; i < length; i++) {
        if (data[i].Q == Q && data[i].P == P) {
            return data[i].theta;
        }
    }
    return -1;
}
int fastCheck(long Q, long P, struct QPTtuple* data, int length, uchar* bloom) {
    // Checks (definitively) whether the list of ideals contains the ideal (Q,P), by first consulting the Bloom filter,
    // and checking the list manually if necessary. This function was written by GPT-4.
    if (!bloom_maybe_contains(bloom, Q, P)) {
        return -1;  // definitely not present
    }
    return checkContainsUnsorted(Q, P, data, length);  // maybe present
}
void find_FUnit(long D, int *FU, int *BS, int *GS){
    // Modified from the infrastructure algorithm described in Section 7.4 of Solving the Pell Equation,
    // by Jacobson and Williams. MR2466979.
    long P = 1;
    long Q = 2;
    const int N = $mainDictSize;
    const int bloom_bytes = $BLOOM_SIZE / 8;
    uchar bloom[bloom_bytes];
    for (int i = 0; i<bloom_bytes; i++){
        bloom[i] = 0;
    };
    struct QPTtuple Cdict[N];
    for (int i=0; i<N; i++){
        Cdict[i].Q = 0;
    };
    long firstQP[2];
    int Ctheta;
    bool eef = false;
    
    reduceQPbaby(D,Q,P,Cdict,bloom,firstQP,&Ctheta,&eef,BS);
    //quicksort(Cdict,*BS);
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
    int Cd,Cd2;
    long S;
    int c;
    int check=0;
    long testA;
    long testB;
    long testC;
    while (i==0 || cc == -1){
        NUCOMPchoose(D,newQ,newP,firstQ,firstP,&newQ,&newP,&testA,&testB,&testC);
        reduceQPred(D,newQ,newP,&newQ,&newP,&Cd2);
    
        Cd = CFFcoset5M8(testA,testB,testC);
        Ctheta = mod(firstCTheta+Ctheta-Cd+Cd2,3);
        i++;
        cc = fastCheck(newQ,mod(newP,newQ),Cdict,*BS,bloom);
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
__kernel void sum_array(__global const long *array, __global int *result, __global int *babystep, __global int *giantstep) {
    // This is the function called from Python. It just runs the find_FUnit functions and writes the output to the buffer arrays.
    // This function was modified from code given by GPT-4.
    int gid = get_global_id(0);
    long D = array[gid];
    int FU;
    int BS=0;
    int GS=0;
    find_FUnit(D,&FU,&BS,&GS);
    result[gid] = FU;
    babystep[gid] = BS;
    giantstep[gid] = GS;
};