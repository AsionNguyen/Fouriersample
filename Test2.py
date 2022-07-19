#!/usr/bin/env python
# coding: utf-8

# In[2]:


def ring(q, alpha):
    return q/alpha^2, alpha^2*q


# In[8]:


center = 2/sqrt(24)*sqrt(d)
error = sqrt(2)/12*alpha*sqrt(q*n)
new_ring = [a-center+error, b-center-error]; new_ring


# In[40]:


def norm(x):
    return x*x.conjugate()

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

def sample_annulus(ring):
    a = ring[0]
    b = ring[1]
    T_r = RealDistribution('uniform', [sqrt(a), sqrt(b)])
    T_theta = RealDistribution('uniform', [0, 2*pi])
    
    r = T_r.get_random_element()
    theta = T_theta.get_random_element()
    x = r*cos(theta)
    y = r*sin(theta)
    
    mu = T_theta.get_random_element()
    f = abs(x)*exp(1j*mu)
    mu = T_theta.get_random_element()
    g = abs(y)*exp(1j*mu)
    
    return f,g


# In[60]:


n = 100
d = euler_phi(n)
q = 12289
alpha = 1.15
a,b = ring(q,alpha)
omega = [exp(2j*pi*k/n) for k in range(n) if gcd(k,n)==1]
A = matrix(CC, d, d, lambda i,j: omega[i]^j)
A_inv = ~A

'''Choose Fourier representations'''
F=[]; G=[]
dd = d/2
for _ in range(dd):
    f,g = sample_annulus(new_ring)
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
f2 = matrix(CC, d, 1, ff)
g2 = matrix(CC, d, 1, gg)
FF2 = A*f2
GG2 = A*g2


# In[61]:


Is_in_ring(FF2,GG2,d,[a,b])


# In[69]:


num_exp =20000
num_suc = 0
for i in range(num_exp):
    F=[]; G=[]
    dd = d/2
    for _ in range(dd):
        f,g = sample_annulus(new_ring)
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
    f2 = matrix(CC, d, 1, ff)
    g2 = matrix(CC, d, 1, gg)
    FF2 = A*f2
    GG2 = A*g2
    
    s = Is_in_ring(FF2,GG2,d,[a,b])
    print(s)
    if s == True:
        num_suc +=1


# In[65]:


num_suc/num_exp


# In[66]:


num_suc


# In[68]:


num_suc


# In[70]:


num_suc


# In[71]:


num_suc/num_exp


# In[73]:


14673.0/20000


# In[ ]:




