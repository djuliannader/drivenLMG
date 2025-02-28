module troterization
push!(LOAD_PATH, pwd())
using LinearAlgebra
import diagonalization
export troter


function troter(H0,Jz,J::Int64,nn::Int64,b,om)
 pi=acos(-1)
 T=2*pi/om
 dt=T/nn
#--- building the time independent Hamiltonian ----
 # diagonal matrix elements
 U0=exp(-im*H0*dt/2)
#--- builiding the time dependent Hamiltonian
 ddiag=[1 for i in 0:J]
 Ut=Jz
 U=Diagonal(ddiag)
#---- Troterization ----------------- 
 for i in 1:nn
  facu1=(im*b/om)*(-sin(i*om*dt)+sin((i-1)*om*dt))
  U1=exp(facu1*Ut)
  UM=U0*U1*U0
  U=UM*U
 end
 #println(H1)
 #println(U)
 return U
end




end
