import sys, math, mpmath, numpy as N

def HB(z1, z2):
   u = z1/z2
   return (u**4 - 3*u**3 - u**2 + 3*u + 1) *\
          (u**8 + 4*u**7 + 7*u**6 + 2*u**5 + 15*u**4 - 2*u**3 + 7*u**2 - 4*u + 1)

def B(z1, z2):
   u = z1/z2
   return -1 - u - 7*(u**2 - u**3 + u**5 + u**6) + u**7 - u**8

def D(z1, z2):
   u = z1/z2
   return -1 + u * (2 + u * (5 + u*u * (5 + u * (-2 - u))))

def f(z1, z2):
   u = (z1/z2)**5
   return z1/z2 * (-1 + u * (11 + u))

def H(z1, z2):
   u = (z1/z2)**5
   return -1 + u * (-228 + u * (-494 + u * (228 - u)))

def T(z1, z2):
   u = (z1/z2)**5
   return 1 + u * (-522 + u * (-10005 + u * u * (-10005 + u * (522 + u))))

def f1f2(a, b, c):
   return a**4 - b**3 + a*b*c

def H1H2(a, b, c):
   return c**4 + 40*a**2*b*c**2 - 192*a**5*c - 120*a*b**3*c + 640*a**4*b**2 - 144*b**5

def T1T2(a, b, c):
   return c**6 + 60*a**2*b*c**4 + 576*a**5*c**3 - 180*a*b**3*c**3 + 648*b**5*c**2 - 2760*a**4*b**2*c**2 +\
          7200*a**7*b*c - 1728*a**10 + 9360*a**3*b**4*c - 2080*a**6*b**3 - 16200*a**2*b**6

def nablasq(a, b, c):
   return 108*a**5*c - 135*a**4*b**2 + 90*a**2*b*c**2 - 320*a*b**3*c + 256*b**5 + c**4

def p(a, b, c):
   return (12**3*f1f2(a, b, c)**5 + H1H2(a, b, c)**3/(12**3) - T1T2(a, b, c)**2/(12**3)) / 2

def q(a, b, c):
   return (-8*a**5*c - 40*a**4*b**2 + 10*a**2*b*c**2 + 45*a*b**3*c - 81*b**5 - c**4) *\
          (64*a**10 + 40*a**7*b*c - 160*a**6*b**3 + a**5*c**3 - 5*a**4*b**2*c**2 +\
                                             5*a**3*b**4*c - 25*a**2*b**6 - b**5*c**2) / 2

def r(a, b, c):
   return (a**2*c**5 - a*b**2*c**4 + 53*a**4*b*c**3 + 64*a**7*c**2 - 7*b**4*c**3 - 225*a**3*b**3*c**2 -\
          12*a**6*b**2*c + 216*a**9*b + 717*a**2*b**5*c - 464*a**5*b**4 - 720*a*b**7) / 2

def s(a, b, c):
   return (-a**2*c**3 + 3*a*b**2*c**2 - 9*b**4*c - 4*a**4*b*c - 8*a**7 - 80*a**3*b**3) / 2

def H1cubf2fiv(a, b, c, nabla):
   return p(a, b, c) + nabla * q(a, b, c)

def M1f2(a, b, c, nabla):
   return (11*a**3*b + 2*b**2*c - a*c**2 - a*nabla) / 2

def N1f1sqT2(a, b, c, nabla):
   return r(a, b, c) + nabla * s(a, b, c)

def Z(a, b, c, nabla):
   return H1cubf2fiv(a, b, c, nabla) / (1728 * f1f2(a, b, c)**5)

def I(z1, z2):
   return H(z1, z2)**3 / (1728*f(z1, z2)**5)

def KZ(a, b, c, nabla):
   m = (11*a**3*b + 2*b**2*c - a*c**2 + a*nabla) / (24*(a**4 - b**3 + a*b*c))
   return (48*a*m**2 - 12*b*m - c)**3/(64*a**2*(12*(a*c-b**2)*m - b*c))

def y(a, b, c, nabla):
   eps = N.exp(2*N.pi*1j / 5)
   l1, l2 = sFunction(Z(a, b, c, nabla)), 1.0
   ys = []
   for nu in range(5):
      ys.append(f(l1, l2) / HB(l1, l2) *\
                M1f2(a, b, c, nabla) / f1f2(a, b, c) +
                D(l1, l2) * T(l1, l2) / (HB(l1, l2) * f(l1, l2)**2) *\
                N1f1sqT2(a, b, c, nabla) / T1T2(a, b, c))
      l1 *= eps
   return ys

def sFunction(Z):
   mpmath.mp.dps = 20
   return N.exp(-N.log(1728*Z) / 5)*\
          mpmath.hyp2f1(31.0/60, 11.0/60, 6.0/5, 1.0/Z) /\
          mpmath.hyp2f1(19.0/60,  -1.0/60, 4.0/5,  1.0/Z)

def randomQuintic():
   while True:
      a, b, c = tuple(N.random.uniform(-1.0, 1.0, 3))
      nabla = N.exp(N.log(nablasq(a, b, c)+0j) / 2.0)
      ZZ = Z(a, b, c, nabla)
      if ZZ.imag > 0.0 or (ZZ.imag == 0 and ZZ.real <= 1):
         continue
      break
   return (a, b, c, nabla, ZZ)

def quintic(z, a, b, c):
   return z**5 + 5*a*z**2 + 5*b*z + c

(a, b, c, nabla, ZZ) = randomQuintic()

print 'Solving quintic: y^5 + %f y^2 + %f y + %f = 0' % (5*a, 5*b, c)
print 'Z:    ', ZZ
print 'Klein:', KZ(a, b, c, -nabla), KZ(a, b, c, nabla)
print 'I:    ', I(sFunction(ZZ), 1.0)
print '\n'.join(['root=%.6f + %.6fi (error=%.10f)' %\
                 (root.real, root.imag, abs(quintic(root, a, b, c)))
                  for root in y(a, b, c, nabla)])
