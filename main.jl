using MAT
using DataStructures
using DynamicalSystems
using LinearAlgebra
using Arpack
using Statistics
using LightGraphs
using LIBSVM
using Random

function cgr(seq)

    len = length(seq);
    x=zeros(Float64, length(seq))
    y=zeros(Float64, length(seq))

    k=counter(seq)
    #print(k)
    kati=keys(k)
    count=length(kati.dict)
    #kleidia=collect(kati.dict)
    angle =2pi/count
    #println(angle)
    kleidia=[key for key in keys(k)]
    #println(kleidia)
    #first_point=[cos(angle) -sin(angle); sin(angle) cos(angle);]*[1; 1]
    #println(first_point)

    for i=1:count
        if i==1
            global proigoumeno_vector = [0;0]
            global vector=([cos(angle) -sin(angle); sin(angle) cos(angle);]*[1; 1])
            global Dictionarymou=Dict(kleidia[i]=>vector)
        else
            #print(proigoumeno_vector)
            global vector=([cos(angle) -sin(angle); sin(angle) cos(angle);]*(proigoumeno_vector))
            global Dictionarymou=merge(Dictionarymou,Dict(kleidia[i]=>vector))
        end
        proigoumeno_vector= vector
    end

    #println(vector)
    #println(Dictionarymou)

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
    #rqa=Recurrence_Analysis(kappa)
    #println(kappa)
    flag=false
    #ts = reconstruct(onedim, 1, 10)
    #println(ts)
    #maxlen=length(onedim)
    if isodd(length(onedim))
         push!(onedim, 0)
         flag = true
    end
    #mutInfo = mutualinformation(onedim, ts[:,2])
    #println(mutInfo)
    t = estimate_delay(onedim, "mi_min" , nbins=10)
    println("tau\n",t)
    #kati=estimate_delay(kappa, "ac_min")
    #println(kati)
    if flag
        pop!(onedim)
    end

    if length(onedim) > 400
        n=15
        ds = estimate_dimension(onedim, t, 1:n,"fnn", rtol=10.0, atol=2.0)
    else
        n=10
        ds = estimate_dimension(onedim, t, 1:n,"fnn", rtol=10.0, atol=2.0)
    end
    print("\n",ds)
    global dime = 0
    while dime == 0
        for j=1:n
            if ds[j] == 0
                global dime = j
                break
            end
        end
        for j=1:n-1
            if ds[j] == ds[j+1] && dime == 0
                global dime = j
                break
            end
        end
        for j=1:n-1
            if ds[j+1] > ds[j]  && dime == 0
                global dime = j
                break
            end
        end
        for j=1:n-1
            if (ds[j] - ds[j+1])/100 <= 0.03 && dime == 0
                global dime = j
                break
            end
        end
    end

    println("dimension: ", dime)
    #dim=trunc(Int, D)
    #print(D)
    reconstructed=embed(onedim, dime, t)
    println(reconstructed)
    #println("\n",i)

    R = RecurrenceMatrix(reconstructed, 0.1; metric = "euclidean");

    apotelesmata_rqa=[recurrencerate(R),determinism(R),dl_average(R),dl_max(R),dl_entropy(R),laminarity(R),trappingtime(R),vl_max(R)]
    #println(apotelesmata_rqa)
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

    #adjacency_matrix(A)
    #println(A)

    g = Graph(A)

    cc=global_clustering_coefficient(g)
    vathmos=maximum(degree(g))

    diam=diameter(g)
    #rintln(vathmos,diam)

    energ=eigvals(A)
    energy=abs.(energ)
    total_energy=sum(energy)
    #println(total_energy)
    n=nv(g)
    m=ne(g)
    #println(n,m)

    #laplasianos=laplacian_matrix(g)
    tade=(laplacian_spectrum(g))

    for i=1:length(energy)
        tade[i]= tade[i] - (2m/n)
    end

    laplacian_energy=sum(abs.(tade))
    #println(laplacian_energy)

    global paths=([],[])
    global array=[]
    for i=1:length(y)
            paths=(bellman_ford_shortest_paths(g,i))
            append!(array,sum(paths.dists))
    end
    array = array/(length(y)-1)
    avg_short=mean(array)

    #println(avg_short)

    hvg_features = [vathmos,diam,total_energy,laplacian_energy,avg_short,cc]
    return(hvg_features)
end

file = matopen("C:/Users/johan/Documents/telika_arxeia/data.mat")
data=read(file,"data")
labels=read(file,"labels") # note that this does NOT introduce a variable ``varname`` into scope
close(file)

#seq=""

#seq="CCHHHHHHHHHHHHHHHHHHHHHHHHHCHHCHHCCHHHHHHHHHCCCHHHHHHHHCCEEEEEEEEEEECCHHHHHHEEEEHHHHHHHHHHCC"
#println(length(seq))
#global pinakas_feature=cgr(seq)
n=length(labels)
#println(n)
indexes=collect(381:n)
indexes=indexes[shuffle(1:end), :]
kappa=length(indexes)
nea_labels=labels[indexes]
global tade=0
for i=1:300
    if indexes[i]==572 || indexes[i] == 805 || indexes[i]== 976 || indexes[i]==945 || indexes[i]==1060 ||  indexes[i]==1249 || indexes[i]==1291 || indexes[i]==1351 || indexes[i]==1422 || indexes[i] == 1427
        global tade= tade + 1
        continue
    end
    #seq=data[i]

    seq=data[indexes[i]]
    #println(seq,"\n")
    if i==1
        global pinakas_feature=cgr(seq)
        println(pinakas_feature)
    else
        println(i)
        pinakas_feature=hcat(pinakas_feature,cgr(seq))
    end
end
telikos=convert(Array{Float64,2}, pinakas_feature)
#final_mat=transpose(telikos)
#println(final_mat)

model = svmtrain(telikos[:,1:2:300-tade], nea_labels[1:2:300-tade])
#println(model)

# Test model on the other half of the data.
(predicted_labels, decision_values) = svmpredict(model, telikos[:,2:2:300-tade]);
println(predicted_labels)

#sum(predicted_labels .== labels[1001:1:1601])/length(predicted_labels)
#println( "Accuracy: \n", mean((predicted_labels .== labels[1401:1:1600]))*100)
