"""FixedK profile using exact linear NNLS reduction and bounded scalar rates."""
import math
import numpy as np
from scipy.optimize import minimize_scalar

def profile_minimum(t,y,K):
 e=K/(1+K);mask=np.isfinite(y);b=y[mask]
 def evaluate(lr):
  z=np.exp(-np.exp(lr)*t);w=e*(1-z)
  # total=u+v and initialcis=v, withu,v>=0. This is exactly the
  # fixed-K closed-state constraint, solved by2DNNLS active-set enumeration.
  A=np.column_stack([np.column_stack([w,1-w])[mask],np.column_stack([w+z,1-w-z])[mask]])
  unconstrained=np.linalg.lstsq(A,b,rcond=None)[0]
  candidates=[np.array([max(float(A[:,0]@b/(A[:,0]@A[:,0])),0),0.]),np.array([0.,max(float(A[:,1]@b/(A[:,1]@A[:,1])),0)])]
  if np.all(unconstrained>=0):candidates.append(unconstrained)
  best=min(candidates,key=lambda c:float(np.sum((A@c-b)**2)));sse=float(np.sum((A@best-b)**2));total=float(best.sum())
  return sse,[math.log(total),float(best[1]/total),float(lr)]
 grid=np.linspace(-21,0,61);values=[evaluate(x)[0] for x in grid];index=int(np.argmin(values));brackets={index} if 0<index<len(grid)-1 else set()
 for j in range(1,len(grid)-1):
  if values[j]<min(values[j-1],values[j+1])-1e-12:brackets.add(j)
 candidates=[evaluate(grid[0]),evaluate(grid[-1])]
 for j in brackets:
  opt=minimize_scalar(lambda x:evaluate(x)[0],bounds=(grid[j-1],grid[j+1]),method='bounded',options=dict(xatol=1e-10));candidates.append(evaluate(opt.x))
 return min(candidates,key=lambda pair:pair[0])

