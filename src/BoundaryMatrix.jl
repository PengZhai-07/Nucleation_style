########################################################
#
#	FUNCTION TO BUILD THE 2D BOUNDARIES FOR SEM
#
########################################################

# 			   (R)T
#        5 10 15 20 25 30
# 	     4  9 14 19 24 29
# (B)L   3  8 13 18 23 28   (T)R 
# 	     2  7 12 17 22 27
# 	     1  6 11 16 21 26
# 	           (L)B

## Reflection coefficient: R22 = (Z2-Z1)/(Z2+Z1)
## Transmission coefficient: T21 = 2Z2/(Z2+Z1)
## wave from 2 to 1

function BoundaryMatrix!(NGLL, NelX, NelY, rho1, vs1, rho2, vs2, 
                         dy_deta, dx_dxi, wgll, iglob, side)

	# INPUT: 
	#		wgll = GLL weights (see GetGLL)
	#		NelX, NelY = no. of elements in x or y direction
	#		iglob = global to local index
	#		jac1D = line Jacobian
	#		side = 'L', 'R', 'T', 'B'

	if side == 'B'
		eB = collect(0:NelY-1)*NelX .+ 1
		igll = 1
		jgll = collect(1:NGLL)
        jac1D = dy_deta
        impedance = rho1*vs1          # 9248880: large impedance means absorbing boundary

	elseif side == 'T'
		eB = collect(0:NelY-1)*NelX .+ NelX
		igll = NGLL
		jgll = collect(1:NGLL)
        jac1D = dy_deta
        impedance = rho1*vs1

	elseif side == 'R'
		eB = (NelY-1)*NelX .+ collect(1:NelX)
		igll = collect(1:NGLL)
		jgll = NGLL
        jac1D = dx_dxi
        impedance = rho1*vs1

    # fault line
	else
		eB = collect(1:NelX)
		igll = collect(1:NGLL)
		jgll = 1
        jac1D = dx_dxi
        impedance = 1    # why the impedance of fault line is 1? Actually, the fault line is a reflecting boundary?
	end

	NelB = length(eB)  # number of all elements at one boundary
	ng = NelB*(NGLL-1) .+ 1  # number of all GLL nodes at one boundary
	iB = zeros(Int, ng)
	B = zeros(ng)
	jB = zeros(NGLL, NelB)

	for e = 1:NelB
		ip = (NGLL-1)*(e-1) .+ collect(1:NGLL)
		iB[ip] = iglob[igll, jgll, eB[e]]   # global node index of boundary GLL nodes
		jB[:,e] = ip
		B[ip] .+= jac1D*wgll*impedance
	end

	return B, iB
end
