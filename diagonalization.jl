module diagonalization
push!(LOAD_PATH, pwd())
using LinearAlgebra
export matrixH0

#J=100
#ep=1
#gx=-4.12
#gy=-12.36



function matrixH0(J,ep,gx,gy)
      v=ep*(gx-gy)/(2.0*(2.0*J-1.0))
      w=ep*(gx+gy)/(2.0*(2.0*J-1.0))
      diag=[ep*(-J+2*(jj))+(w)*(J*(J+1)-(-J+2*jj)^2) for jj in 0:J]  
      HMatrix=Array(Diagonal(diag))
      for i in 0:J-1
        HMatrix[i+1,i+2]=1.0*(v/2.0)*(J*(J+1)-(-J+2*i+2)*((-J+2*i+2)-1))^(1.0/2.0)*(J*(J+1)-((-J+2*i+2)-1)*((-J+2*i+2)-2))^(1.0/2.0)
      end
      for i in 0:J-1
        HMatrix[i+2,i+1]=HMatrix[i+1,i+2]        
      end
      return HMatrix
      end

function matrixH00(J,ep,gx,gy)
      v=ep*(gx-gy)/(2.0*(2.0*J-1.0))
      w=ep*(gx+gy)/(2.0*(2.0*J-1.0))
      diagz=[m for m in -J:J]
      diag0=[0.0 for i in 0:(2*J)]
      diagar = [0 for m in -J:J-1]
      diagab = [(J*(J+1)-(m*(m+1)))^(1/2) for m in -J:J-1]
      sigzop= Diagonal(diagz)
      sigplusop = Tridiagonal(diagar,diag0,diagab)
      sigminusop = transpose(sigplusop)
      Hmatrix = ep*sigzop + 1.0*(v/2)*(sigplusop^2 + sigminusop^2) + 1.0*(w/2)*(sigplusop*sigminusop + sigminusop*sigplusop)
      return Hmatrix
  end

function matrixJz(J)
      diagz=[m for m in -J:J]
      sigzop= Diagonal(diagz)
      return sigzop
end

function matrixJx(J)
     diag0=[0.0 for i in 0:(2*J)]
     diagar = [0 for m in -J:J-1]
     diagab = [(J*(J+1)-(m*(m+1)))^(1/2) for m in -J:J-1]
     sigplusop = Tridiagonal(diagar,diag0,diagab)
     sigminusop = transpose(sigplusop)
     sigmx = (1/2)*(sigplusop + sigminusop)
     return sigmx
end



end