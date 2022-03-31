# calculate the nucleation size of model with fault damage zone

function NucleationSize(P)        #  (Rubin & Ampuero, 2005)
#sigma = 5e7             # unit: Pa   50MPa
sigma = maximum(P[3].Seff)             # unit: Pa   50MPa
#mu = 3e10               # unit: Pa    3GPa
mu = P[2].mu               # unit: Pa    3GPa
#L = 0.008               # unit: m    8mm
L = maximum(P[3].xLf)               # unit: m    8mm
#b = 0.019
b = maximum(P[3].ccb)
#b_a = 0.0041
b_a = -minimum(P[3].cca .- P[3].ccb)

h_hom = 2/pi*mu*L*b/sigma/b_a^2
return h_hom
# halfwidth of fault damage zone
# H = 2000        # unit: m  
#H = P[2].ThickY

end

function CohesiveZoneSize(P)             # (Day, 2005)
    
    #sigma = 5e7             # unit: Pa   50MPa
    sigma = maximum(P[3].Seff)             # unit: Pa   50MPa
    #mu = 3e10               # unit: Pa    30GPa
    mu = P[2].mu               # unit: Pa    3GPa
    #L = 0.008               # unit: m    8mm
    L = maximum(P[3].xLf)               # unit: m    8mm
    #b = 0.019
    b = maximum(P[3].ccb)
    #b_a = 0.0041
    b_a = -minimum(P[3].cca .- P[3].ccb)
    
    CZone = 9*pi/32*mu*L/b/sigma
    return CZone
    # halfwidth of fault damage zone
    # H = 2000        # unit: m  
    #H = P[2].ThickY
    
    end