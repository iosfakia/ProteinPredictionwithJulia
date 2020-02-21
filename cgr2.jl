using DataStructures
using DynamicalSystems
using MAT

function Recurrence_Analysis(onedim)

    #maxlen=length(onedim)
    #if isodd(maxlen)
    #    maxlen=maxlen-1
    #end

    #t=1
    #Ds = estimate_dimension(onedim, t, 1:5,"fnn")
    #D=maximum(Ds)
    #dim=trunc(Int, D)
    #print("dimension\n",dim,"\n")
    #con_mat=embed(onedim, dim, t)

    #R = RecurrenceMatrix(con_mat, 0.1; metric = "euclidean");

    #apotelesmata_rqa=[recurrencerate(R),determinism(R),dl_average(R),dl_max(R),dl_entropy(R),laminarity(R),trappingtime(R),vl_max(R)]
    #println(apotelesmata_rqa)
    #return(apotelesmata_rqa)
end


file = matopen("C:/Users/johan/Documents/telika_arxeia/data.mat")
data=read(file,"data")
labels=read(file,"labels") # note that this does NOT introduce a variable ``varname`` into scope
close(file)

#seq=("CCHHHHHHHHHHHHHHHHHHHHHHHHCHHCHHCCHHHHHCHHHHCCCHHHHHHHHCCEEEEEE")
#seq=("CCCHHHHHHHHHHHHHHHHHHCCCHHHHHHHHCCCCCEEEEEECCC")
seq =""
for i=225:380
    seq=data[i]
    len = length(seq);
    x=zeros(Float64, 1, len)
    y=zeros(Float64, 1, len)

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

    #println("Timeseries:\n",x)

    #println( "Y is: ", y , "\n")
    flag=false
    onedim=vec(x)
    #rqa=Recurrence_Analysis(kappa)
    #println(kappa)

    #ts = reconstruct(onedim, 1, 10)
    #println(ts)
    #maxlen=length(onedim)
    if isodd(length(onedim))
         push!(onedim, 0)
         flag=true
    end
    #mutInfo = mutualinformation(onedim, ts[:,2])
    #println(mutInfo)
    t = estimate_delay(onedim, "mi_min", nbins=10 )
    println("tau\n",t)
    if flag
        pop!(onedim)
    end
    #kati=estimate_delay(kappa, "ac_min")
    #println(kati)
    println(i)
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
    if dime == 0
        return(i)
    end
    #dim=trunc(Int, D)
    #print(D)
    #x=embed(onedim, dime, t)
    #println(x)
    #println("\n",i)
end
