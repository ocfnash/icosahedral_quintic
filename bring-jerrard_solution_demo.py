import mpmath, numpy as N

def nabla(cc):
    return N.exp(N.log(cc**4 + 256 + 0j)/2)

def Z(cc):
    return (2*1728 + cc**4*(207 + cc**4) - cc*cc*(81 + cc**4)*nabla(cc))/(2*1728)

def z(ZZ):
    mpmath.mp.dps = 20
    return N.exp(-N.log(1728*ZZ) / 5) * mpmath.hyp2f1(31.0/60, 11.0/60, 6.0/5, 1.0/ZZ) /\
                                        mpmath.hyp2f1(19.0/60, -1.0/60, 4.0/5, 1.0/ZZ)

def HB(zz):
    return (zz**4 - 3*zz**3 - zz**2 + 3*zz + 1) *\
           (zz**8 + 4*zz**7 + 7*zz**6 + 2*zz**5 + 15*zz**4 - 2*zz**3 + 7*zz**2 - 4*zz + 1)

def D(zz):
    return -1 + 2*zz + 5*zz**2 + 5*zz**4 - 2*zz**5 - zz**6

def f(zz):
    u = zz**5
    return zz*(-1 + u * (11 + u))

def T(zz):
    u = zz**5
    return 1 + u * (-522 + u * (-10005 + u * u * (-10005 + u * (522 + u))))

def y(cc):
    ys = []
    eps = N.exp(2*N.pi*1j / 5)
    zz = z(Z(cc))
    for n in range(5):
        ys.append(-cc * f(zz) / HB(zz) -\
                  (7*cc*cc + 9*nabla(cc))/(2*cc**5 + 2*648*cc)*D(zz)*T(zz)/(HB(zz)*f(zz)**2))
        zz *= eps
    return [(yy, abs(yy**5 + 5*yy + cc)) for yy in ys]

for (z, err) in y(-2.8234): # Randomly chosen value for quintic (see xkcd #221).
    print '%.6f + %.6fi   (%.10f)' % (z.real, z.imag, err)
