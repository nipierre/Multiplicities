using DelimitedFiles
using DataFrames
using GLM

RC = readdlm("./data/proton_semi_inclusive_RC_backup2.txt")
cqel = readdlm("./data/sigtot_RC_wqel.dat")
cnoqel = readdlm("./data/sigtot_RC_woqel.dat")

QelCorr = zeros(9,6)
z_range = [.1 .225 .275 .325 .375 .425 .475 .525 .575 .625 .675 .725 .8 .925];

for i in 1:9
    for j in 1:6
      if cnoqel[i,2*(j-1)+1]==0
        QelCorr[i,j] = 0
      else
        QelCorr[i,j] = cqel[i,2*(j-1)+1]/cnoqel[i,2*(j-1)+1]
      end
    end
end

for c in 1:2
    for i in 1:9
        for j in 1:6
            for k in 1:14
                 println(QelCorr[i,j])
                 println(RC[(c-1)*9*14*6+(i-1)*14*6+(j-1)*14+k,8])
                 RC[(c-1)*9*14*6+(i-1)*14*6+(j-1)*14+k,8] = RC[(c-1)*9*14*6+(i-1)*14*6+(j-1)*14+k,8]*QelCorr[i,j]
                 println(RC[(c-1)*9*14*6+(i-1)*14*6+(j-1)*14+k,8])
            end
            # x = Float64[]
            # y = Float64[]
            # bin = Int64[]
            # for k in 2:6
            #     if RC[(c-1)*9*14*6+(i-1)*14*6+(j-1)*14+k,8]>0
            #         push!(y,RC[(c-1)*9*14*6+(i-1)*14*6+(j-1)*14+k,8])
            #         push!(x,z_range[k])
            #         push!(bin,k)
            #     end
            # end
            # data = DataFrame(X=x,Y=y)
            # ols = lm(@formula(Y ~ X),data,true  )
            # param = coef(ols)
            # # println(param[1],":",param[2])
            # for k in 1:14
            #     RC[(c-1)*9*14*6+(i-1)*14*6+(j-1)*14+k,8] = z_range[k]*param[2]+param[1]#*QelCorr[i,j]
            # end
        end
    end
end


writedlm("data/newRC.txt", RC)
writedlm("data/elcorr.txt",QelCorr)
