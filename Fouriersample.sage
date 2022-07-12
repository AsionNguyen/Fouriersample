import cmath
def ring(q, alpha):
    return q/alpha^2, alpha^2*q

def is_in_ring(t, ring):
    norm = t[0]*t[0].conjugate()+t[1]*t[1].conjugate()
    if norm > ring[0] and norm < ring[1]:
        return True
    else:
        return False

def Is_in_ring(x, y, d, ring):
    for i in range(d):
        if is_in_ring([x[i],y[i]], ring) == False:
            return (i, False)
    return True
    
def sample_annulus(alpha, q):
    a, b = ring(q, alpha)
    T_r = RealDistribution('uniform', [sqrt(a), sqrt(b)])
    T_theta = RealDistribution('uniform', [0, 2*pi])
    
    r = T_r.get_random_element()
    theta = T_theta.get_random_element()
    x = r*cos(theta)
    y = r*sin(theta)
    
    mu = T_theta.get_random_element()
    f = abs(x)*cmath.exp(1j*mu)
    mu = T_theta.get_random_element()
    g = abs(y)*cmath.exp(1j*mu)
    
    return f,g

def sample_Fourier(x,y,mu):
    f = abs(x)*exp(1j*mu)
    g = abs(y)*exp(1j*mu)
    
    return f,g





'''Settings
phi_n : nth cyclotomic polynome
d     : deg(phi_n)=#{k|(k,n)=1}
'''
n = 150
d = euler_phi(n)
q = 12289
alpha = 1.15
a,b = ring(q,alpha)
omega = [cmath.exp(2j*pi*k/n) for k in range(n) if gcd(k,n)==1]
A = matrix(CC, d, d, lambda i,j: omega[i]^j)
A_inv = ~A

'''Choose Fourier representations'''
F=[]; G=[]
dd = d/2
for _ in range(dd):
    f,g = sample_annulus(alpha,q)
    F.append(f)
    G.append(g)

F2 = [F[len(F)-1-i].conjugate() for i in range(len(F))]
G2 = [G[len(G)-1-i].conjugate() for i in range(len(F))]

FF = matrix(CC, d, 1, F+F2)
GG = matrix(CC, d, 1, G+G2)

'''Integral Fourier inverse'''
#theoretically, f,g are real but they are complex in this computation
f = A_inv*FF           
g = A_inv*GG

#rounding f to ff (closest integral vector)
ff = [round(f.numpy()[i][0].real) for i in range(d)]
gg = [round(g.numpy()[i][0].real) for i in range(d)]

#correct if Fourier representations of ff and gg do not stay in the annulus
def correct(ff,gg):
    f2 = matrix(CC, d, 1, ff)
    g2 = matrix(CC, d, 1, gg)
    FF2 = A*f2
    GG2 = A*g2

    while Is_in_ring(FF2,GG2,d,[a,b])!= True:
        i = Is_in_ring(FF2,GG2,d,[a,b])[0]
        if ff[i]*ff[i].conjugate()+gg[i]*gg[i].conjugate()<a:
            ff[i]+=1
        if ff[i]*ff[i].conjugate()+gg[i]*gg[i].conjugate()>b:
            ff[i]-=1
        f2 = matrix(CC, d, 1, ff)
        g2 = matrix(CC, d, 1, gg)
        FF2 = A*f2
        GG2 = A*g2
        
    return ff,gg

correct(ff,gg)