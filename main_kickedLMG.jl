push!(LOAD_PATH, pwd())
using LinearAlgebra
using QuantumOptics
import diagonalization
import troterization
import wigner
import statistics


# ----- Imput parameters  ----
k = 1                      # state of interest
J=50                       # System size
ep=1.0                     # LMG parameter
gx=-0.95                   # LMG parameter
gy= 3*gx                   # LMG parameter
epsilon = 0.001            # strength of the Kick
tau = 2.89134              # period of the kick
NN=100                     # Size of the Grid
name1="wignertest.dat"     # Wigner output file
name2="husimitest.dat"     # Husimi output file
# ------------------------- ----


#---------------- Main body -----------------------
#--------------------------------------------------
#---------  Building Floquet operator -------------#
HH0 = diagonalization.matrixH00(J,ep,gx,gy)
ev0 = eigvals(HH0)
Jz =  diagonalization.matrixJz(J)
Jx =  diagonalization.matrixJx(J)
Kop = Jx + Jz
Floquet = exp(-im*tau*HH0)*exp(-im*epsilon*Kop)
fstates = eigvecs(Floquet)
fev     = eigvals(Floquet)
println("-> Floquet operator obtained")
#---   sorting Floquet states   -------------#
sortf   = statistics.sortingvecf(Floquet,HH0,J)
listvec=[fstates[i,sortf[k]] for i in 1:length(fev)]
println("-> Stationary state ",k,"-th obtained")
# building floquet state for QuantumOptics library
psi = wigner.buildingstate(listvec,J)
# calculating Husimi and Wigner using Quantum optics library
htest = wigner.husimif(psi,NN,name2)
println("-> Husimi function obtained")
wtest = wigner.wignerf(psi,NN,name1)
println("-> Wigner function obtained")
# calculating <F_k|H0|F_k>
psif = listvec
psift = transpose(psif)
fexpval = real(psift*HH0*psif)
#------------------------------------------------------


#
#open("resonances.dat","w") do io
#T=0.1
#tint = 0.01
#Nmax=400
#for ik in 1:Nmax
#   test1 = statistics.expectation(J,ep,gx,gy,epsilon,T,1)
#   println(io,T," ",test1[2])
#   println(ik/Nmax)
#   T=T+tint
#end

#end

println("See file resonances.dat for expectation value vs tau")

# ----  Printing results ------------
println("----- Results ----------")
println("Floquet stationary state: ",k)
println("<E_k|H0|E_k>/J: ",ev0[k]/J)
println("<F_k|H0|F_k>/J: ",fexpval/J)
println("Go to file ",name2," for Husimi function")
println("Go to file ",name1," for Wigner function")
println("-----------------------" )
#-------------------------------------

