# calculate the nucleation size of model with fault damage zone
# mu is only is number 
# sigma, L, b, b_a are all vectors


function NucleationSize(P, alpha)        
#sigma = 5e7             # unit: Pa   50MPa
sigma = maximum(P[3].Seff)             # unit: Pa   50MPa
#mu = 3.2e10               # unit: Pa    32 GPa
mu = P[2].mu               # unit: Pa    32 GPa
#L = 0.008               # unit: m    8mm
L = maximum(P[3].xLf)               # unit: m    8mm
#b = 0.019 non-changebale
b = maximum(P[3].ccb)
#b_a = 0.0041
b_a = 0.004

h_hom_host = 2/pi*mu*L*b/sigma/b_a^2        #  (Rubin & Ampuero, 2005)
h_hom_dam  = 2/pi*mu*alpha*L*b/sigma/b_a^2        #  (Rubin & Ampuero, 2005)   smaller!!  downlimit of nucleation size

return h_hom_host, h_hom_dam

end

function CohesiveZoneSize(P, alpha)             # (Day, 2005)
    
    #sigma = 5e7             # unit: Pa   50MPa
    sigma = maximum(P[3].Seff)             # unit: Pa   50MPa
    #mu = 3.2 e10               # unit: Pa    32 GPa
    mu = P[2].mu*alpha             # unit: Pa    32 GPa     
    #L = 0.008               # unit: m    8mm
    L = minimum(P[3].xLf)              # unit: m    8mm
    #b = 0.019
    b = maximum(P[3].ccb)
    #b_a = 0.0041
    
    CZone = 9*pi/32*mu*L/b/sigma

    print("sigma=",sigma,"  mu=",mu,"  L=",L,"  b=",b, "\n")
    return CZone 
    
    end