using Plots 
using LinearAlgebra
using Statistics 
using CurveFit

global ev = 1.6*10^(-19)
global Noise_lvl = [0.00,0.5,1.00]
global Noise_lvcont = LinRange(0,1,10)
global hbar = 1.057*10^(-34)
global kb = 1.38*10^(-23)
#We'll be studying 2 materials, commonly used in solar cells (either regular or thin-films)

function generateH(MaterialName,l)
    #nb of atoms per direction
    if cmp("Si",MaterialName) == 0
        a = 5.43*10^(-10) #angstroms 
        theta_D = 645 #Debye temp
        wD = 8.9*10^(13)
        m = 28.0855 * 1.66*10^(-27)
        noiseT = 1/a * sqrt(kb/(m*wD^2))
        hoppot = [5,5,1.2,1.2,1.2] .* ev
        noiselvlV = 0.04 
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
        noiselvlV = 1.4*10^(-2)
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
    #pl = scatter(pos_atm[:,1],pos_atm[:,2],pos_atm[:,3]) #-> plots the lattice structure in 3D
    #display(pl)
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
                    H[i,j] = (1 + (noiselvlV + noiseT*sqrt(T))*2*(rand()-0.5))*H[i,j] 
                else #hoppMat terms 
                    H[i,j] = (1 + (noiselvlJ + exp(-hbar*4.963*10^(13)/(kb*T)))*2*(rand()-0.5))*H[i,j]
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
    DOSmobedge = DOS2[trunc(Int64,indtailend/nbptE +1)]  
    for j in 1:length(DOS2) #need correspondance between DOS and E 
        if DOS2[j] - DOSmobedge/exp(1) < 0
            tailstart = j
            break 
        end 
    end 
    engstart = trunc(Int64,tailstart/nbptE + 1)   
    return E2[1] - E1[1] #E2[indtailend] - E2[engstart]
end 

function urb2(DOS1,DOS2,E1,E2,nbptE)
    return(E2[length(E2)]-E1[length(E1)])
end 



function UrbTail(l,MaterialName) #only plots DOS 
    Eurb = Float64[]
    nbptE = 20
    T = 300 #Kelvin
    #No noise : 
    H = generateH(MaterialName,l)[1] 
    E = real.((eigvals(H)))
    DOS1 = dos3D(nbptE,E)
    x = LinRange(minimum(E)/ev,maximum(E)/ev,nbptE)
    p = plot(x,DOS1,title ="Noise effect on Density of States (3D solid)",xlabel="Energy (eV)",ylabel="DOS",label="Noise FWHM (eV) \n 0",legend=:best)
    #iterating over the noise levels
    for i in 2:length(Noise_lvl)
        H,noiseV,noiseJ,V0,noiseT = generateH(MaterialName,l) 
        noiselvlV = noiseV * Noise_lvl[i]
        noiselvlJ = noiseJ * Noise_lvl[i]
        WignerNoise!(H,noiselvlV,noiselvlJ,noiseT,T)
        E = real.((eigvals(H)))
        #display(H)
        xp = LinRange(minimum(E)/ev,maximum(E)/ev,nbptE)
        DOS = dos3D(nbptE,E)
        plot!(xp,DOS,label=noiselvlV*V0/ev)
    end
    display(p)
end


function UrbvsNoise(l,MaterialName) #returns Urbach energy,used multiple times to plot an accurate Urbach Energy vs Noise
    Eurb = Float64[0]
    nbptE = 100 #arbitrarily chosen, might need to be changed (higher nbptE -> higher energy resolution)
    T = 300
    #No noise : 
    H = generateH(MaterialName,l)[1]
    E0 = real.((eigvals(H)))
    DOS0 = dos3D(nbptE,E0)
    Eurb = [0.00]
    #iterating over the noise levels
    for i in 2:length(Noise_lvcont)
        H,noiseV,noiseJ,V0,noiseT = generateH(MaterialName,l) 
        noiselvlV = noiseV * Noise_lvcont[i]
        noiselvlJ = noiseJ * Noise_lvcont[i]
        WignerNoise!(H,noiselvlV,noiselvlJ,noiseT,T)
        E = real.((eigvals(H)))
        DOS = dos3D(nbptE,E)
        Urbe = urb2(DOS0,DOS,E0,E,nbptE)/ev
        push!(Eurb,abs.(Urbe))
    end
    return Eurb 
end


# l = 10
# V0 = 5*ev
# MaterialName = "GaAs"
# noisenorm = 0.019
# nbrep = 5 #nb of iterations (to obtain average values of Eurb)
# Eurb = zeros(length(Noise_lvcont))
# for i in 1:nbrep 
#     Eurb = Eurb .+ UrbvsNoise(l,MaterialName) 
# end 
# Eurb = Eurb./nbrep
# x = Noise_lvcont.*V0*(noisenorm + 0.002)/ev
# p = scatter(x,Eurb,title = "Urbach Energy dependance on Noise",xlabel = "FWHM of the Gaussian Noise (eV)",ylabel="Urbach energy (eV)")
# a,b = linear_fit(x,Eurb)
# y = b.*x .+ a 
# plot!(x,y,label="linear fit")
# display(p)

#Visualising DOS 


UrbTail(10,"Si") 





        