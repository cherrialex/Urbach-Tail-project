using DSP
using Distributions
using IterativeSolvers
using LinearAlgebra
using Noise
using Random
using Plots
using EasyFit


global Noise_lvl = [0.0,1.0]
global ev = 1.6*10^(-19)


function wells!(V0,V,N) 
    for i in 1:2*N 
        if i%2 == 1
            V[i] = 0
        else 
            V[i] = V0*ev
        end
    end  
    return V     
end

function generate_J(NoisyV,V) 
    J = zeros(Float64,trunc(Int64,length(NoisyV)/2))
    p = 1
    for i in 2:2:length(NoisyV) 
        J[p] = ev + V[p] - NoisyV[p]
        p+=1
    end 
return J
end 


function NoiseV(V,Noiselv) 
    NoisyV = zeros(Float64,length(V))
    for i in 1:length(V) 
        if i%2 == 0 
            NoisyV[i] = V[i] + Noiselv*rand()*maximum(V)/50
                end 
    end #Noiselv is between 0 and 1 
return NoisyV
end

function init_H(V,t) 
    Vtrunc = V[2:2:length(V)]
    tboundary = pop!(t)
    H = zeros(Complex{Float64},length(Vtrunc),length(Vtrunc)) 
    for i in 1:length(Vtrunc)
        for j in 1:length(Vtrunc) 
            if i==j 
                H[i,j] = Vtrunc[i] 
            elseif i==j-1 
                H[i,j] = t[i]
            elseif i==j+1
                H[i,j] = t[j] 
            end 
        end 
    end 
    H[length(Vtrunc),1],H[1,length(Vtrunc)] = tboundary,tboundary 
    return H 
end



function getDOSstep(nbptE,E) 
    sort!(E) #g(E) = dN(E)/dE : choose a small dE interval, count the nb of energies and divide 
    DOS = zeros(Float64,trunc(Int64,length(E)*nbptE))
    dE = maximum(E)/(trunc(Int64,nbptE*length(E))) #width of the little energy interval
    j = 1 
    p = 1   
    Eval = minimum(E) - dE
    for i in 1:length(DOS) 
        if Eval <= E[j] && Eval + dE > E[p] 
            while Eval + dE > E[p] && p < length(E)
                p+=1
            end    
        end 
        DOS[i] = (p-j)/dE 
        j = p
        Eval += dE 
    end 
    return DOS
end 



function pos_discont(DOS) #finds the position of the Van Hove singularities
    discpos1 = 1
    discpos2 = 1
    max1 = 0
    max2 = 0
    for i in 1:trunc(Int64,length(DOS)/2) #1st discontinuity (left of the band)
        if DOS[i] > max1
            discpos1 = i 
            max1 = DOS[i]
        end 
    end 
    for j in trunc(Int64,length(DOS)/2):length(DOS) 
        if DOS[j]>max2
            discpos2 = j 
            max2 = DOS[j]
        end 
    end 
    return discpos1,discpos2 
end 

function moy(arr) #mean of an array 
    moy = 0
    for i in 1:length(arr) 
        moy += arr[i] 
    end 
return (moy/length(arr))
end 

function fiturbach(DOS,E,nbptE)
    posleft,posright = pos_discont(DOS)[1],pos_discont(DOS)[2]
    N0l = DOS[posleft] #left side : DOS = N0l*exp(E - Ecl/Eurbl) 
    N0r = DOS[posright] #right side : DOS = N0r*exp(Ecr - E/Eurbr)
    Ecl = E[1 + trunc(Int64,posleft/nbptE)]
    Ecr = E[1 + trunc(Int64,posright/nbptE)]
    El = []
    Er= []

    Eurbl = 0
    Eurbr = 0
    for i in 1:posleft - 1 #obtaining the Urbach energy for the left tail
        if DOS[i] != 0
            push!(El,(Ecl - E[1 + trunc(Int64,i/nbptE)])/log(N0l/DOS[i]))
        end 
    end 
    Eurbl = mean(El) 

    for p in 1:posleft - 1 #fitting the left exponential tail to the DOS
        DOS[p] = N0l * exp((E[1 + trunc(Int64,p/nbptE)] - Ecl)/Eurbl)
    end 


    for j in posright + 1 : length(DOS) #obtaining the Urbach energy for the right tail
        if DOS[j] != 0 
            push!(Er, (E[1 + trunc(Int64,j/nbptE)] -  Ecr)/log(N0r/DOS[j]))
        end 
    end 
    
    Eurbr = moy(Er) 

    for k in posright + 1:length(DOS)-1 #fitting the right tail 
        DOS[k] = N0r * exp((Ecr - E[1 + trunc(Int64,k/nbptE)])/Eurbr)
    end  
    #print(Eurbl/ev ,"   ",Eurbr/ev)
return DOS, Eurbl, Eurbr
end 



function Tailstep(V0,N) 
    V = zeros(Float64,2*N)
    V = wells!(V0,V,N)
    nbptE = 0.6
    E = zeros(Float64,0)
    H = Matrix(1.0I,length(V),length(V))
    Eurbli = 0
    Eurbri = 0
    Eurbr = [0.0] 
    Eurbl = [0.0]
    for i in 1:length(Noise_lvl)
        NoisyV = NoiseV(V,Noise_lvl[i])
        J = generate_J(NoisyV,V) 
        H = init_H(NoisyV,J) 
        E = real.((eigvals(H)))
        DOS = getDOSstep(nbptE,E) 
        x = LinRange(minimum(E)/ev ,maximum(E)/ev,length(DOS))
        xlabel!("Energy (eV)") 
        ylabel!("Density of states") 
        if i==1 
            #plot(x,DOS,title ="Density of states dependance on noise percentage", label = "0%",xlim=(5.26,5.34))
            
        else 
            corrDOS, Eurbli, Eurbri = fiturbach(DOS,E,nbptE)
            push!(Eurbl,Eurbli/ev)
            push!(Eurbr,Eurbri/ev) 
            p = plot(x,DOS, title ="Density of states dependance on noise percentage", label = "2%",xlims=(2.5,5.5))
            # a,b,c = fitexp(x[mbr:length(x)],DOS[mbr:length(x)])
            # xsl = x[mbr:length(x)]
            # plot(xsl,a.*exp.(xsl./b) .+ c)
            display(p)
        end
    end
end 


@time Tailstep(5,1000)



