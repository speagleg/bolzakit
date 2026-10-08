# P1 FEM for -Delta on Bolza surface, Klein model octagon with opposite-side pairings
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla, sys
n=int(sys.argv[1]); k_eigs=int(sys.argv[2]) if len(sys.argv)>2 else 30
R=np.arccosh(3+2*np.sqrt(2)); rho=np.tanh(R)
V=[np.array([rho*np.cos(k*np.pi/4), rho*np.sin(k*np.pi/4)]) for k in range(8)]
def tohyp(p): # Klein -> hyperboloid
    r2=p@p; t=1/np.sqrt(1-r2); return np.array([t,t*p[0],t*p[1]])
def tok(X): return X[1:]/X[0]
def hyp_interp(A,B,s): # point at fraction s of hyperbolic arclength from A to B (Klein coords)
    a,b=tohyp(A),tohyp(B); c=a[0]*b[0]-a[1]*b[1]-a[2]*b[2]; L=np.arccosh(c)
    if L<1e-14: return A
    X=(np.sinh((1-s)*L)*a+np.sinh(s*L)*b)/np.sinh(L); return tok(X)
pts=[];idx={}
def add(p):
    key=(round(p[0],12),round(p[1],12))
    if key in idx: return idx[key]
    idx[key]=len(pts); pts.append(p); return idx[key]
tris=[]; O=np.zeros(2)
for k in range(8):
    A=V[k];B=V[(k+1)%8]
    grid={}
    for i in range(n+1):
        # row i: points along hyperbolic geodesic from P=O->A at s_i to Q=O->B at s_i ; use hyperbolic param
        s=i/n
        P=hyp_interp(O,A,s); Q=hyp_interp(O,B,s)
        for j in range(i+1):
            t=j/i if i>0 else 0
            grid[(i,j)]=add(hyp_interp(P,Q,t))
    for i in range(n):
        for j in range(i+1):
            tris.append((grid[(i,j)],grid[(i+1,j)],grid[(i+1,j+1)]))
            if j<i: tris.append((grid[(i,j)],grid[(i+1,j+1)],grid[(i,j+1)]))
pts=np.array(pts); N=len(pts)
# side pairing: boost along direction theta_k by 2d mapping side k+4 -> side k
d=np.arccosh(1+np.sqrt(2))
def boost(th,L):
    c,s=np.cosh(L),np.sinh(L); u=np.array([np.cos(th),np.sin(th)])
    M=np.eye(3); M[0,0]=c; M[0,1:]=s*u; M[1:,0]=s*u; M[1:,1:]=np.eye(2)+(c-1)*np.outer(u,u); return M
parent=list(range(N))
def find(x):
    while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
    return x
from scipy.spatial import cKDTree
tree=cKDTree(pts)
# boundary nodes of side m
def side_nodes(m):
    A=V[m];B=V[(m+1)%8]; return [idx[(round(q[0],12),round(q[1],12))] for q in [hyp_interp(A,B,j/n) for j in range(n+1)]]
for k in range(4):
    th=(2*k+1)*np.pi/8
    M=boost(th,2*d)
    for a in side_nodes(k+4):
        img=tok(M@tohyp(pts[a])); dist,b=tree.query(img)
        assert dist<1e-8, dist
        parent[find(a)]=find(b)
roots=sorted(set(find(i) for i in range(N))); rid={r:i for i,r in enumerate(roots)}
lab=np.array([rid[find(i)] for i in range(N)]); Nf=len(roots)
# assemble
T=np.array(tris); I=[];J=[];Kv=[];Mv=[]
P=pts[T]  # (nt,3,2)
x1,x2,x3=P[:,0],P[:,1],P[:,2]
Dm=np.stack([x2-x1,x3-x1],axis=2) # (nt,2,2) columns
det=np.linalg.det(Dm); area=np.abs(det)/2
Dinv=np.linalg.inv(Dm)
G=np.einsum('tij,kj->tki',np.transpose(Dinv,(0,2,1)),np.array([[-1,-1],[1,0],[0,1]]))  # grads (nt,3,2)
# quadrature: edge midpoints
mids=[(x1+x2)/2,(x2+x3)/2,(x3+x1)/2]
phis=[np.array([.5,.5,0]),np.array([0,.5,.5]),np.array([.5,0,.5])]
Kloc=np.zeros((len(T),3,3)); Mloc=np.zeros((len(T),3,3))
for q,ph in zip(mids,phis):
    r2=np.sum(q*q,1); w=area/3
    A=(1/np.sqrt(1-r2))[:,None,None]*(np.eye(2)[None]-np.einsum('ti,tj->tij',q,q))
    Kloc+=w[:,None,None]*np.einsum('tai,tij,tbj->tab',G,A,G)
    Mloc+=(w*(1-r2)**-1.5)[:,None,None]*np.outer(ph,ph)[None]
L=lab[T]
rows=np.repeat(L,3,axis=1).ravel(); cols=np.tile(L,(1,3)).ravel()
K=sp.csr_matrix((Kloc.ravel(),(rows,cols)),shape=(Nf,Nf)); Mm=sp.csr_matrix((Mloc.ravel(),(rows,cols)),shape=(Nf,Nf))
print("nodes",Nf,"area",Mm.sum(),"4pi",4*np.pi)
vals=sla.eigsh(K,k=k_eigs,M=Mm,sigma=-0.5,which='LM',return_eigenvectors=False)
vals=np.sort(vals); print(np.array2string(vals,precision=6,max_line_width=120))
