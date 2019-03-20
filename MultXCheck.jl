using DelimitedFiles
using Plots
using LaTeXStrings

z = [.20,.25,.30,.35,.40,.45,.50,.55,.60,.65,.70,.75,.85]

MultN = readdlm("data/multiplicities_hadron.txt")
# AccN = readdlm("data/.txt")
MultM = readdlm("data/MarcinMult.txt")

pullM = Float64[]
Mult_p = zeros(12,45)
Mult_m = zeros(12,45)

for i in 1:9
    for j in 1:5
        for k in 1:12
            global pull
            if MultN[(i-1)*5*12+(j-1)*12+(k-1)+1,8]>0 && MultM[(i-1)*5*12+(j-1)*12+(k-1)+1,4]>0
                push!(pullM,(MultN[(i-1)*5*12+(j-1)*12+(k-1)+1,8]-MultM[(i-1)*5*12+(j-1)*12+(k-1)+1,4])/MultN[(i-1)*5*12+(j-1)*12+(k-1)+1,8])
                Mult_p[k,(i-1)*5+(j-1)+1] = (MultN[(i-1)*5*12+(j-1)*12+(k-1)+1,8]-MultM[(i-1)*5*12+(j-1)*12+(k-1)+1,4])/MultN[(i-1)*5*12+(j-1)*12+(k-1)+1,8]
            end
            if MultN[(i-1)*5*12+(j-1)*12+(k-1)+1,16]>0 && MultM[(i-1)*5*12+(j-1)*12+(k-1)+1,6]>0
                push!(pullM,(MultN[(i-1)*5*12+(j-1)*5+(k-1)+1,16]-MultM[(i-1)*5*12+(j-1)*12+(k-1)+1,6])/MultN[(i-1)*5*12+(j-1)*12+(k-1)+1,7])
                Mult_m[k,(i-1)*5+(j-1)+1] = (MultN[(i-1)*5*12+(j-1)*12+(k-1)+1,16]-MultM[(i-1)*5*12+(j-1)*12+(k-1)+1,6])/MultN[(i-1)*5*12+(j-1)*12+(k-1)+1,7]
            end
        end
    end
end

histogram(pullM)
savefig("plots/PullXCheck.png")

zmid = zeros(12)
for i in 1:12
    zmid[i] = (z[i+1]+z[i])/2
end

plot(zmid,Mult_p,layout=(5,9),seriestype=:scatter,
                                size = (3200,1800),
                                dpi = 300,
                                ylims = (-1,1.),
                                # ylims = (0.75,1.25),
                                legend=false)
savefig("plots/MultXCplus.png")

plot(zmid,Mult_m,layout=(5,9),seriestype=:scatter,
                                size = (3200,1800),
                                dpi = 300,
                                ylims = (-1.,1.),
                                # ylims = (0.75,1.25),
                                legend=false)
savefig("plots/MultXCminus.png")
