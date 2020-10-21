v0=sol.v
v0Adj=sol.v_adj
v0/=sqrt(v0'*v0)
v0Adj/=conj(v0Adj'*l(sol.params[sol.eigval],1)*v0)
