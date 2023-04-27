using DSP
using Distributions
using IterativeSolvers
using LinearAlgebra
using Noise
using Random
using Plots
using EasyFit


global Noise_lvl = [0.0,1.00]
global Noise_lvcont = LinRange(0.00,1.00,20)
global ev = 1.6*10^(-19)
global kb = 1.38*10^(-23)
global nEng = 50
global hbar = 1.057*10^(-34)

#Properties of the material studied 
global noisenorm = 0.036  #value of the noise calculated from the material's Young's modulus, for Si
global alphaSi = 4.73*10^(-4) #ev/K 
global betaSi = 635 #K 
global theta_D = 645 #Debye temp
global wD = 8.9*10^13
global m = 28.0855 * 1.66*10^(-27)
global Jx = Jy = Jz = 1.4*ev #eV
#other elements 
global alphaGaAs = 4.77 * 10^-4
global betaGaAs = 355

global alphaCdTe = 3.08 * 10^-4
global betaCdTe = 104

global Eg = 1.12*ev #for Si

global a = 5.43*10^(-10) #distance between each site (lattice constant), for Si
#The Noise describes the static disorder of the material, due to defects or the non-crystalline structure for example.
#The noise due to temperature will be directly introduced into the hopping and site energies. 

#for the site energies : we obtain the displacement via the einstein model of lattice vibrations 

function generateH(MaterialName,l)
    #nb of atoms per direction
    if cmp("Si",MaterialName) == 0
        a = 5.43*10^(-10) #angstroms 
        theta_D = 645 #Debye temp
        wD = 13.19*10^(13)
        m = 28.0855 * 1.66*10^(-27)
        noiseT = 1/a * sqrt(kb/(m*wD^2))
        hoppot = [5,5,1.2,1.2,1.2] .* ev
        noiselvlV = 4*10^(-2) 
        noiselvlJ = noiselvlV
    elseif cmp("GaAs",MaterialName)==0 
        a = 5.65*10^(-10) 
        m = 144.6*1.66*10^(-27)
        theta_D = 360
        wD = 4.968*10^(13)
        noiseT = 1/a * sqrt(kb/(m*wD^2))
        hoppot = [4.5,4.9,0.9,1.3,1.7] .* ev 
        noiselvlV = 1.9*10^(-2)
        noiselvlJ = noiselvlV
    elseif cmp("CdTe",MaterialName) == 0
        a = 6.48*10^(-10)
        theta_D = 146
        wD = 2.014*10^(13)
        m = 140.4 * 1.66*10^(-27)
        noiseT = 1/a * sqrt(kb/(m*wD^2))
        hoppot = [2.9,5.1,1.8,2.2,2.6] .* ev
        noiselvlV = 3.4*10^(-3)
        noiselvlJ = noiselvlV
    end 
    H = zeros(Complex{Float64}, 2*l^3, 2*l^3)
    lat = a/2 * [1 1 0; 0 1 1; 1 0 1] 
    pos_atm = Float64[]
    count = 0
    for i in 0:l-1
        for j in 0:l-1
            for k in 0:l-1
                pos = i*lat[1,:] + j*lat[2,:] + k*lat[3,:]
                posp = pos .+ a*0.25
                pos_atm = vcat(pos_atm, pos)
                pos_atm = vcat(pos_atm,posp)
            end
        end
    end
    pos_atm = transpose(reshape(pos_atm,3,2*l^3))
    # pl = scatter(pos_atm[:,1],pos_atm[:,2],pos_atm[:,3]) #-> plots the lattice structure in 3D
    # display(pl)
    for i in 1:2*l^3 #iterating over atoms (defining the "reference atom")
        pos1 = pos_atm[i,:]
        for j in 1:2*l^3 #iterating over all atoms 
            pos2 = pos_atm[j,:]
            if i == j
                if i%2 == 1 
                    H[i, j] += hoppot[1] # on-site energy
                else 
                    H[i, j] += hoppot[2]
                end 
            else 
                dx, dy, dz = abs.(pos1 .- pos2) #implementing boundary conditions 
                dx -= (l-1)*a*round(dx/(l*a))
                dy -= (l-1)*a*round(dy/(l*a))
                dz -= (l-1)*a*round(dz/(l*a))
                dist = norm([dx,dy,dz])
                if dist < 1.01*a*sqrt(3)/4 
                    if i%2==1 && j%2==1
                        H[i, j] += hoppot[3] # A1-A1 hopping
                        count +=1
                    elseif i%2==1 && j%2==0
                        H[i, j] += hoppot[4] # A1-A2 hopping
                        count +=1
                    elseif i%2==0 && j%2==1
                        H[i, j] += hoppot[4] # A2-A1 hopping
                        count +=1
                    elseif i%2==0 && j%2==0
                        H[i, j] += hoppot[5] # A2-A2 hopping
                        count +=1
                    end 
                end 
            end
        end
    end
#print(count, "  ")
return H, noiselvlV, noiselvlJ, (hoppot[1]+hoppot[2])/2,noiseT
end 
   

function WignerNoise!(H,noiselvlV,noiselvlJ,noiseT,T) 
    n = length(H[1,:])
    for i in 1:n
        for j in 1:n
            if H[i,j]!=0 
                if i==j #site energies terms 
                    H[i,j] = (1 + 2*(rand()-0.5)*(noiseT*sqrt(T)))*H[i,j] 
                else #hoppMat terms 
                    H[i,j] = (1 + exp(-hbar*wD/(kb*T))*(rand() - 0.5))*H[i,j]
                end 
            end 
        end 
    end 
return H
end 

function dos3D(nbpoints, energylevels) #updated DOS histogram
    Emin = minimum(energylevels)
    Emax = maximum(energylevels)
    dE = (Emax - Emin) / nbpoints 
    energyrange = LinRange(Emin - 2*dE, Emax, nbpoints)
    dos = zeros(nbpoints)
    for p in 1:length(energyrange)-1
        for Eng in energylevels
            if energyrange[p]<Eng && Eng < energyrange[p+1]
                dos[p] +=1
            end 
        end 
    end
    dos = dos/(dE * length(energylevels))
    return dos
end


function urb_energy(DOS1,DOS2,E1,E2,nbptE) 
    tailstart = 0 
    indtailend = 0
    #finding the index that corresponds to the end of the tail for E2 (mobility edge)
    for i in 1:length(E2) 
        if E2[i] - E1[1] > 0
            indtailend = i 
            break 
        end 
    end 
    DOSmobedge = DOS2[trunc(Int64,indtailend*length(DOS2)/length(E2) +1)]  #nbptE is not * lengthDOS : to be reviewed
    for j in 1:length(DOS2) #need correspondance between DOS and E 
        if DOS2[j] - DOSmobedge/exp(1) < 0
            tailstart = j
            break 
        end 
    end 
    engstart = trunc(Int64,tailstart*length(E2)/length(DOS2) + 1)   
    return E1[1] - E2[1] #E2[indtailend] - E2[engstart]
end 

function urb2(DOS1,DOS2,E1,E2,nbptE)
    return(E2[length(E2)]-E1[length(E1)])
end 

#used to be E2[1] - E1[1] -> gave the good dependency but to a factor of 1/100 ...

function UrbTail(V0,N) 
    #V = zeros(Float64,N)
    #V = wells!(V0,V,N)
    Eurb = Float64[]
    nbptE = 20
    T = LinRange(0.01,1000,3)
    #No noise : 
    #NoisyV = pot3D(V,0.0,T[1]) 
    #H = init_H3D(NoisyV,T[1],0.0) 
    H,noiseV,noiseJ,V0,noiseT = generateH("Si",N) 
    noiselvlV = noiseV * Noise_lvl[1]
    noiselvlJ = noiseJ * Noise_lvl[1]
    WignerNoise!(H,noiselvlV,noiselvlJ,noiseT,T[1])
    E0 = real.((eigvals(H)))
    display(E0)
    DOS0 = dos3D(nbptE,E0)
    x = LinRange(minimum(E0)/ev ,maximum(E0)/ev,nbptE)
    p = plot(x,DOS0,title ="Noise effect on Density of States (3D solid)",xlabel="Energy (eV)",ylabel="DOS",label="Temp Noise (K) \n 0")
    #iterating over the noise levels
    for i in 2:length(T)
        #NoisyV = pot3D(V,0,T[i])
        #H = init_H3D(NoisyV,T[i],0) 
        H,noiseV,noiseJ,V0,noiseT = generateH("Si",N) 
        WignerNoise!(H,noiselvlV,noiselvlJ,noiseT,T[i])
        E = real.((eigvals(H)))
        x = LinRange(minimum(E)/ev ,maximum(E)/ev,nbptE)
        DOS = dos3D(nbptE,E)
        display(E)
        plot!(x,DOS,label=T[i])
    end
    display(p)
end


function UrbvsTemp(V0,N) #returns Urbach energy, for fixed static disorder
    #V = zeros(Float64,N)
    #V = wells!(V0,V,N)
    nbptE = 20 #arbitrarily chosen, might need to be changed
    T = LinRange(0.01,1000,12)
    #iterating over the noise levels
    #NoisyV = pot3D(V,0,T[1])
    #H = init_H3D(NoisyV,T[1],0) 
    H,noiseV,noiseJ,V0,noiseT = generateH("Si",N) 
    noiselvlV = noiseV * Noise_lvl[1]
    noiselvlJ = noiseJ * Noise_lvl[1]
    WignerNoise!(H,noiselvlV,noiselvlJ,noiseT,T[1])
    E0 = real.((eigvals(H)))
    DOS_0 = dos3D(nbptE,E0) 
    Eurb = Float64[0] 
    for j in 2:length(T)
        #NoisyV = pot3D(V,0,T[j])
        #H = init_H3D(NoisyV,T[j],0) 
        H,noiseV,noiseJ,V0,noiseT = generateH("Si",N) 
        noiselvlV = noiseV * Noise_lvl[1]
        noiselvlJ = noiseJ * Noise_lvl[1]
        WignerNoise!(H,noiselvlV,noiselvlJ,noiseT,T[j])
        E = real.((eigvals(H)))
        DOS = dos3D(nbptE,E) 
        push!(Eurb,abs.(urb2(DOS_0,DOS,E0,E,nbptE))/ev)
    end 
    # p=plot(T,Eurb,xlabel="Temperature (K)",ylabel="Urbach Energy (eV)",title="Urbach energy dependence on temperature",label="Algorythm values")
    # plot!(T,alphaSi * T.^2 ./(2*(T .+ betaSi)),label="Varshni model for Si")
    # display(p) 
    return Eurb 
end


V0=5
N=10
#UrbTail(V0,N)
nbrep = 8
T = LinRange(0.001,1000,12)
Eurb = zeros(Float64,length(T))
for i in 1:nbrep 
    Eurb += UrbvsTemp(V0,N)
end 
Eurb ./= nbrep 
p=scatter(T,Eurb,xlabel="Temperature (K)",ylabel="Urbach Energy (eV)",title="Urbach energy dependence on temperature",label="Algorythm values")
plot!(T,alphaSi * T.^2 ./(2*(T .+ betaSi)),label="Varshni model for Si")
display(p) 

#UrbvsTemp(V0,N)

