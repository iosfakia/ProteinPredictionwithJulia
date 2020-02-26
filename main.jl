import Pkg; Pkg.add("DataStructures"); Pkg.add("MAT"); Pkg.add("DynamicalSystems"); Pkg.add("LinearAlgebra");
Pkg.add("Arpack"); Pkg.add("Statistics"); Pkg.add("LightGraphs");
using MAT
using DataStructures
using DynamicalSystems
using LinearAlgebra
using Arpack
using Statistics
using LightGraphs


function cgr(seq)

    len = length(seq);
    x=zeros(Float64, length(seq))
    y=zeros(Float64, length(seq))

    k=counter(seq)
    kati=keys(k)
    count=length(kati.dict)

    angle =2pi/count
    kleidia=[key for key in keys(k)]

    for i=1:count
        if i==1
            global proigoumeno_vector = [0;0]
            global vector=([cos(angle) -sin(angle); sin(angle) cos(angle);]*[1; 1])
            global Dictionarymou=Dict(kleidia[i]=>vector)
        else
            global vector=([cos(angle) -sin(angle); sin(angle) cos(angle);]*(proigoumeno_vector))
            global Dictionarymou=merge(Dictionarymou,Dict(kleidia[i]=>vector))
        end
        proigoumeno_vector= vector
    end

    v=[value for value in values(Dictionarymou)]

    for j= 1:(len-1)
        for i=1:count
            simeio=v[i]
            if seq[j] == kleidia[i]
                if j==1
                    x[j] = 0.5 * (0 + simeio[1])
                    y[j] = 0.5 * (0 + simeio[2])
                end
            x[j+1] = 0.5 * (x[j] + simeio[1])
            y[j+1] = 0.5 * (y[j] + simeio[2])
            end
        end
    end

    pinakas_feature=[]
    append!(pinakas_feature,Recurrence_Analysis(x))
    append!(pinakas_feature,Recurrence_Analysis(y))
    append!(pinakas_feature,hvg(x))
    append!(pinakas_feature,hvg(y))
    append!(pinakas_feature,mean(x))
    append!(pinakas_feature,mean(y))
    return(pinakas_feature)
end

function Recurrence_Analysis(onedim)
    onedim=vec(onedim)
    dime=2
    t=2
    reconstructed=embed(onedim, dime, t)
    println(reconstructed)

    R = RecurrenceMatrix(reconstructed, 0.1; metric = "euclidean");

    apotelesmata_rqa=[recurrencerate(R),determinism(R),dl_average(R),dl_max(R),dl_entropy(R),laminarity(R),trappingtime(R),vl_max(R)]

    return(apotelesmata_rqa)
end

function hvg(y)

    N = length(y)
    A = zeros(N, N)

    for i = 1:N
        for j = i:N-1
            if (y[i] >= y[j+1] && i<N)
                j=j+1
                A[i,j] = 1
                A[j,i] = 1
            else
                j=j+1
                A[i,j] = 1
                A[j,i] = 1
                break
            end
        end
    end
    g = Graph(A)

    cc=global_clustering_coefficient(g)
    vathmos=maximum(degree(g))

    diam=diameter(g)

    energ=eigvals(A)
    energy=abs.(energ)
    total_energy=sum(energy)

    n=nv(g)
    m=ne(g)

    tade=(laplacian_spectrum(g))

    for i=1:length(energy)
        tade[i]= tade[i] - (2m/n)
    end

    laplacian_energy=sum(abs.(tade))

    global paths=([],[])
    global array=[]
    for i=1:length(y)
            paths=(bellman_ford_shortest_paths(g,i))
            append!(array,sum(paths.dists))
    end
    array = array/(length(y)-1)
    avg_short=mean(array)

    hvg_features = [vathmos,diam,total_energy,laplacian_energy,avg_short,cc]
    return(hvg_features)
end

println("Please enter input file directory:\n")

myfile = readline()

global data=[]
open(myfile) do f
    for (i, line) in enumerate(eachline(f))
        global data = vcat(data,line)
    end
end

n=length(data)

for i=1:n
    seq=data[i]

    if i==1
        global pinakas_feature=cgr(seq)
        println(pinakas_feature)
    else
        pinakas_feature=hcat(pinakas_feature,cgr(seq))
    end
end
final_matrix=convert(Array{Float64,2}, pinakas_feature)
println(final_matrix)
