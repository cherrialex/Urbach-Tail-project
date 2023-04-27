using LinearAlgebra
using Random
using Plots
using QuadGK
using Polynomials

#Need for more precise Noise values : use Young Modulus for precise materials --> maximum displacement

global Noise_lvdiscr = [0.0,1.0]
global Noise_lvl = LinRange(0,1.00,20)
global kB=8.6173324E-5 #in eV
global J0 = 1 #in eV : value to be changed
global T = 300 #to be changed when we study the effect of temperature  
global V0 = 5 
global ev = 1.6*10^-19
global a = 5.8*10^(-8) #10nm well width, subject to change 

# Calculates number of eigenvalues less than 'sigma' in tridiagonal matrix 



function wells(V0,N) 
    V = zeros(Float64,0)
    for i in 1:N 
        push!(V,V0)
    end  
    return V     
end


#   SiteEnergy - scalar eV; reference for site energy
#   disorder - scalar eV ; amount of Gaussian / normal energetic disorder, for trace of Hamiltonian
#   modelJ(theta) - function, takes degrees, returns eV ; model for the transfer integral 
#   B - scalar ; Thermodynamic (B)eta, used to populate Probability Density Function
#   Z - scalar ; Partition function, weighting for absolute Boltzmann populations
#   U - function; Free energy function, used to generate Bolztmann populations 
#   N - integer ; size of diagonal of Hamiltonian



function modelJ(r,a) 
    return J0*exp.(-r./a)  #a is a characteristic distance :well width for example
end 


function U(r) 
    return V0 * sin(2*r*π/a)^2 
end 

function partitionFct(V0,U) 
    return quadgk(x->exp(V0)*U(x),0,a)[1]  
end 


function randH(SiteEnergy, disorder, modelJ,B,Z,U,N)
    # Random Trace / diagonal elements
        D = SiteEnergy + disorder*randn(N)*V0/25 #the 1/25 (4%) corresponds to the max strain in 1D (arbitrarily chosen)
        positions=Float64[]
        count = 0
        for i in 1:N-1 
            r = 0.0
            while true #do while loop 
                count +=1
                r = a * rand()     #random position in the well (between 0 and a)
                p = exp(-U(r)*B)/Z  #probability by stat mech 
                if p > 1.0 
                    #print("Error \n ")
                elseif p > rand()
                    break   
                end   #rejection sampling of distribution
            end
            push!(positions,r)     
        end
    E = modelJ.(positions,a) 
    #Squared for computational speed
    E = E.^2
    return D,E
end  
    

function sturm(D,E_squared,sigma)  
    t = 0.0  
    countnegatives=0  
    
    t=D[1] - sigma  #first term of sequence calculated differently, to avoid needing E[0], t[0]
    if t < 0.0
        countnegatives = countnegatives + 1
    end
    for i in 2:length(D)
        t = D[i] - sigma - E_squared[i-1]/t   # Sturm sequence, overwriting temporary values...
        if t < 0.0                     # if t<0, we've found another eigenvalue
            countnegatives = countnegatives + 1 
        end
    end

    return countnegatives
end


function getDOSsturm(sigma,V0,Noiselv,modelJ,B,Z,U,N) #we will obtain the DOS between Emin - dE and Emax
    dE = (sigma[length(sigma)]-sigma[1])/length(sigma)
    V = wells(V0,N)
    D,Esqrd= randH(V, Noiselv,modelJ,B,Z,U,N)
    a = 0
    DOS = zeros(Float64,0) #the number of negatives in a given interval corresponds to the number of energy eigenvalues in the interval --> DOS
    for i in 1:length(sigma) 
        totnbnegative = sturm(D,Esqrd,sigma[i])
        push!(DOS,(totnbnegative - a)/dE) #avoid recounting energies from lower sigmas
        a = totnbnegative  
    end 
    return sigma,DOS
end 

function urbachenergy(sigma1,DOS1,sigma2,DOS2) 
    mobility_edge = 0
    lowest_eng = 0
    j = 1
    for i in 1:length(sigma1) 
        if DOS1[i]!=0
            mobility_edge = sigma1[i]
            break
        end 
    end 
    while DOS2[j]==0 
        j+=1
    end 
    lowest_eng = sigma2[j]
    return mobility_edge - lowest_eng #Urbach energy (mobility edge - energy of the least energetic localised state)
end 
    

function DOSplot(N) 
    Emin = 0
    Emax = 10 #random values, will be changed 
    nbE =  50 #large enough to have decent precision (low dE intervals)
    sigma = LinRange(Emin,Emax,nbE)
    B = 1/(kB*T)
    Z = partitionFct(V0,U)  
    Energies0,DOS0 = getDOSsturm(sigma,V0,Noise_lvdiscr[1],modelJ,B,Z,U,N)
    plot(Energies0,DOS0,title="1D Density of States dependance on Noise",xlabel="Energy (eV)",ylabel="DOS",label="Noise FWHM \n 0")
    for i in 2:length(Noise_lvdiscr)
        Energies, DOS = getDOSsturm(sigma,V0,Noise_lvdiscr[i],modelJ,B,Z,U,N)
        
        p = plot!(Energies,DOS,label=Noise_lvdiscr[i]*V0/10,xlims=(2.5,7.5))
        display(p)
    end 
end 


function UrbEvsNoise(N) 
    Emin = 0
    Emax = 10 #random values, will be changed 
    nbE =  50 #large enough to have decent precision (low dE intervals)
    Urb = Float64[0] 
    sigma = LinRange(Emin,Emax,nbE)
    B = 1/(kB*T)
    Z = partitionFct(V0,U)  
    Energies0,DOS0 = getDOSsturm(sigma,V0,Noise_lvl[1],modelJ,B,Z,U,N)
    for i in 2:length(Noise_lvl)
        Energies, DOS = getDOSsturm(sigma,V0,Noise_lvl[i],modelJ,B,Z,U,N)
        push!(Urb,urbachenergy(Energies0,DOS0,Energies,DOS))
        #print("Urbach energy for ",Noise_lvl[i]*V0/25,"eV ",Urb[i-1],"eV \n") 
    end 

    #obtaining the coefficient describing the linear dependence of Urbach Energy in regards to noise 
#     slopemoy = 0
#     for i in 1:length(Urb)-1 
#         slopemoy += (Urb[i+1] - Urb[1])/((Noise_lvl[i+1]-Noise_lvl[1])*V0/25)

#     end 
#     slopemoy /= length(Urb) 
     return(Urb)
end 

#Visualising the DOS
@time DOSplot(1000)



#Obtaining the Urbach Energy dependance on noise 
# nbit = 30
# #slope = 0
# Eurb = zeros(Float64,length(Noise_lvl))
# for i in 1:nbit 
#     #slope += UrbEvsNoise(1000)[1]
#     Eurb = Eurb .+ UrbEvsNoise(1000)
# end 
# Eurb = Eurb ./nbit

# p = plot(Noise_lvl*V0/25,Eurb,title ="Urbach Energy dependence on Noise (1D)",xlabel = "Noise (eV)",ylabel = "Urbach Energy (eV)")
# #print("Average slope value for Eurb vs Noise : ",slope/nbit)
# display(p)