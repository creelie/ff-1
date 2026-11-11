"""Re-check every finite computation in
"Isomorphic Invariant Orders Without Hermite Equivalence" (D. Bhattacharjee).

Requires cypari2 (pip install cypari2). Run: python3 verify_resultants.py
"""
import math
import cypari2

pari = cypari2.Pari()
x = pari('x')

def poly(c):
    k = len(c) - 1
    return sum(int(c[i]) * x**(k - i) for i in range(k + 1))

def hres(F, G):
    # resultant of binary forms with nonzero leading coefficients
    return int(pari.polresultant(poly(F), poly(G)))

H7 = [9233654208009611348491360869664999, -158785974126310386817980985635657596,
      351336977243533453505499571533318677, 2956841106545254256188116418780082341,
      -535022898010655232998946727047265289, 228203582350855930536981792812157305,
      -62503182984543134439438862552002751, 20218560847845476602195651556634169]

ROWS = [
    (40,  (1, 0, -10),  (2, 4, -3),    {'H4': [-213, 254, 1679, -1135, 67],
                                        'H6': [-1, -2, 11, 19, -9, 9, -7]}),
    (229, (1, 1, -57),  (3, -5, -17),  {'H3': [14, -33, -845, 2678]}),
    (401, (1, 1, -100), (2, 23, 16),   {'H5': [526505, 10598024, 47546124, -52212324,
                                               69923173, -93789359]}),
    (577, (1, 1, -144), (2, -11, -57), {'H7': H7}),
]

ok = True
def check(cond, msg):
    global ok
    print(('ok   ' if cond else 'FAIL ') + msg)
    ok &= bool(cond)

def disc2(q):
    return q[1]**2 - 4*q[0]*q[2]

def cycle(d, start):
    # cycle of reduced indefinite forms; step (A,B,C) -> (C,B',C')
    s = math.isqrt(d)
    def red_b(B, C):
        m = 2*abs(C)
        Bp = (-B) % m
        while Bp <= s - m:   # want sqrt(d)-2|C| < B' < sqrt(d)
            Bp += m
        while Bp > s:
            Bp -= m
        return Bp
    out, f = [], start
    while True:
        out.append(f)
        A, B, C = f
        Bp = red_b(B, C)
        Cp = (Bp*Bp - d) // (4*C)
        assert Bp*Bp - 4*C*Cp == d
        f = (C, Bp, Cp)
        if f == start or len(out) > 100:
            return out

for d, q0, q1, Hs in ROWS:
    check(disc2(q0) == d and disc2(q1) == d, f'd={d}: disc q0 = disc q1 = {d}')
    for name, H in Hs.items():
        r0, r1 = hres(list(q0), H), hres(list(q1), H)
        P = poly(H)
        check(abs(r0) == 1 and abs(r1) == 1, f'd={d} {name}: Res(q0)={r0}, Res(q1)={r1}')
        check(H[0] != 0 and pari.poldisc(P) != 0 and pari.content(P) == 1
              and pari.polisirreducible(P), f'd={d} {name}: leading coeff, D(H)!=0, primitive, irreducible')
    K = pari.bnfinit(pari(f't^2-{d}'), 1)
    a, b = q1[0], q1[1]
    I = pari.idealhnf(K, a, pari(f'({-b}+t)/2'))   # a*I_{q1}
    cls = pari.bnfisprincipal(K, I, 0)
    check(any(c != 0 for c in cls), f'd={d}: a*I_q1 not principal (class {cls}, Cl={K.bnf_get_cyc()})')

# principal cycles of Lemma 6.2
for d, start in [(229, (1, 15, -1)), (401, (1, 19, -10)), (577, (1, 23, -12))]:
    cyc = cycle(d, start)
    firsts = sorted({f[0] for f in cyc})
    print(f'     d={d} cycle: {cyc}')
    check(not ({2, -2, 3, -3} & set(firsts)), f'd={d}: first coefficients {firsts} avoid +-2, +-3')
check(all((t*t) % 5 not in (2, 3) for t in range(5)) is True, 'x^2-10y^2 = +-2, +-3 impossible mod 5')

# Theorem B
H = [1, 0, -8, -6, -1]
q1s = [3, 8, 2]
check(disc2(q1s) == 40, 'q1* = 3x^2+8xy+2y^2 has discriminant 40')
check(hres([1, 0, -10], H) == 1 and hres(q1s, H) == 1, 'Res(q0,H) = Res(q1*,H) = 1')
P = poly([1, 0, -10]) * poly(H)
Pp = poly(q1s) * poly(H)
check(P == pari('x^6-18*x^4-6*x^3+79*x^2+60*x+10'), 'P as stated')
check(Pp == pari('3*x^6+8*x^5-22*x^4-82*x^3-67*x^2-20*x-2'), "P' as stated")
check(pari.poldisc(P) == pari.poldisc(Pp) != 0, f'D(P) = D(P\') = {pari.poldisc(P)}')
K = pari.bnfinit(pari('t^2-40'), 1)
check(any(c != 0 for c in pari.bnfisprincipal(K, pari.idealhnf(K, 3, pari('-4+t/2')), 0)),
      '3*I_{q1*} = (3, -4+sqrt10) not principal')
H3 = ROWS[1][3]['H3']
Q = poly([1, 1, -57]) * poly(H3)
Qp = poly([3, -5, -17]) * poly(H3)
check(Q == pari('14*x^5-19*x^4-1676*x^3+3714*x^2+50843*x-152646'), 'Q as stated')
check(Qp == pari('42*x^5-169*x^4-2608*x^3+12820*x^2+975*x-45526'), "Q' as stated")
check(pari.poldisc(Q) == pari.poldisc(Qp) != 0, f"D(Q) = D(Q') = {pari.poldisc(Q)}")
print('ALL CHECKS PASSED' if ok else 'SOME CHECK FAILED')
