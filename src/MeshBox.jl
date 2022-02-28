###############################################################################
#		Spectral Element Mesh for Rectangular Box, with internal 
#		Gauss-Legendre-Lobatto (GLL) dub-grids.
#
#	INPUT:	NGLL = no. of GLL nodes
#           Nel = Total no. of elements
#			NELX = no. of elements in x
#			NELY = no. of elements in y
#           FltNglob::Int = NelX*(NGLL - 1) + 1     no. of GLL nodes in x direction
#           dxe = Size of one element along X
#           dye = Size of one element along y

#	OUTPUT:	I = iglob[i,j,e] is the global node index of the 
#							 (i,j)th GLL node internal to the
#							 e-th element.
#			Elements are numbered row by row from  bottom-left
#			to top-right. The table iglob is needed to assemble 
#			global data from local data.

#			x[:] = global x coordinates of all GLL nodes, starting at 0
#			y[:] = global y coordinates of all GLL nodes, starting at 0
###############################################################################


function MeshBox!(NGLL, Nel, NelX, NelY, FltNglob, dxe, dye)

	XGLL = GetGLL(NGLL)[1]
	#  [-1.000000000000000, -0.6546536707079772, 0.000000000000000E+00, 0.6546536707079771, 1.000000000000000]

	iglob = zeros(Int, NGLL, NGLL, Nel)
	nglob = FltNglob*(NelY*(NGLL-1) + 1)           #            nglob= no. of all GLL nodes        

    x::Vector{Float64} = zeros(nglob)
    y::Vector{Float64} = zeros(nglob)

	et = 0
	last_iglob = 0


	# collect: Return an Array of all items in a collection or iterator
	# reashape(iterator,row,column)
	ig = reshape(collect(1:NGLL*NGLL), NGLL, NGLL)        # all GLL in one element 
	igL = reshape(collect(1:NGLL*(NGLL-1)), NGLL-1, NGLL) #  no Left(Bottom) edge
	igB = reshape(collect(1:NGLL*(NGLL-1)), NGLL, NGLL-1) # no Bottom(Left) edge
	igLB = reshape(collect(1:(NGLL-1)*(NGLL-1)), NGLL-1, NGLL-1) # rest of the elements no bottom & left edge

	xgll = repeat(0.5*(1 .+ XGLL), 1, NGLL)    # from [-1,1] to [0,1]
	ygll = dye*xgll'
	xgll = dxe*xgll              # project to local coordinate


	@inbounds for ey = 1:NelY         # number of y elements    row
		@inbounds for ex = 1:NelX     # number of x elements    column

			et = et + 1

			# Redundant nodes at element edges
            
            # NGLL = number of GLL nodes per element

			# ig: index of GLL in each element
            
			if et == 1
				ig = reshape(collect(1:NGLL*NGLL), NGLL, NGLL)     
			else
				if ey == 1       # Bottom Row
					ig[1,:] = iglob[NGLL, :, et-1]   # Left edge
					ig[2:end, :] = last_iglob .+ igL # The rest

				elseif ex == 1	    # Left Column
					ig[:,1] = iglob[:,NGLL,et-NelX]	# Bottom edge
					ig[:,2:end] = last_iglob .+ igB 	# The rest

				else 			# Other Elements
					ig[1,:] = iglob[NGLL, :, et-1]	# Left edge
					ig[:,1] = iglob[:, NGLL, et-NelX]# Bottom edge
					ig[2:end, 2:end] = last_iglob .+ igLB
				end
			end

			iglob[:,:,et] = ig
			last_iglob = ig[NGLL, NGLL]

			# Global coordinates of computational nodes
			# dimension of x and y is 1, index of x/y can be array, but the final result is still a vector seprately !!!
			@inbounds x[ig] .= dxe*(ex-1) .+ xgll
			@inbounds y[ig] .= dye*(ey-1) .+ ygll

		end
	end
	return iglob, x, y
end

