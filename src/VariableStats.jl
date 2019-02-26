function init_VariableStats()
  VarStats = VariableStats(
  zeros(Float64, 1,1), # Histogram
  [0],
  [0.0],
  [0.0],
  [0.0],
  [0.0],
  [0.0],
  [0.0],
  [0.0],
  [0.0],
  [0.0],
  [0.0],
  [0.0],
  [0.0],
  Array{Dict{Tuple{Int,Int,Int}, Float64},1}()
  )
  return VarStats
end

export init_VariableStats

function getPosNegCentroids(CCobj::CCstats, VarStats::VariableStats, QuantileScores::Array{Float32, 3}; threshigh::Float64 = 0.0, threslow::Float64 = 0.0, min_nobs_within_thres::Int = 100)
  ValueWeighting = true
  if size(QuantileScores) != size(CCobj.LabelCube) error("size  of CCobj.LabelCube and QuantileScores not matching") end
  if length(CCobj.affectedLonLatTimeScoreValue) == 0 error("CCobj.affectedLonLatTimeScoreValue is empty") end
  VarStats.High_LONcenters = zeros(Float64, CCobj.numExtremes)
  VarStats.High_LATcenters = zeros(Float64, CCobj.numExtremes)
  VarStats.High_TIMEcenters = zeros(Float64, CCobj.numExtremes)
  sweightsHigh = zeros(Float64, CCobj.numExtremes)
  valweightsHigh = zeros(Float64, CCobj.numExtremes)
  Div = zeros(Int, CCobj.numExtremes)
  VarStats.Low_LONcenters = zeros(Float64, CCobj.numExtremes)
  VarStats.Low_LATcenters = zeros(Float64, CCobj.numExtremes)
  VarStats.Low_TIMEcenters = zeros(Float64, CCobj.numExtremes)
  sweightsLow = zeros(Float64, CCobj.numExtremes)
  valweightsLow = zeros(Float64, CCobj.numExtremes)
  Div2 = zeros(Int, CCobj.numExtremes)
  for i = 1:length(CCobj.affectedLonLatTimeScoreValue)
    for lonlattimeidx = keys(CCobj.affectedLonLatTimeScoreValue[i])
      (lonidx, latidx, timeidx) = lonlattimeidx
      if ValueWeighting
        val = QuantileScores[lonidx, latidx, timeidx]
      else
        val = 1
      end
      if val > threshigh
        innerstCentroids!(sweightsHigh, valweightsHigh, VarStats.High_LATcenters, VarStats.High_LONcenters, VarStats.High_TIMEcenters, CCobj, lonlattimeidx, val, i)
        Div[i] += 1
      elseif val < threslow
        innerstCentroids!(sweightsLow, valweightsLow, VarStats.Low_LATcenters, VarStats.Low_LONcenters, VarStats.Low_TIMEcenters, CCobj, lonlattimeidx, val, i)
        Div2[i] += 1
      end
    end
  end
  broadcast!(/ ,VarStats.High_LATcenters, VarStats.High_LATcenters, sweightsHigh)
  #broadcast!(acosd ,VarStats.High_LATcenters,VarStats.High_LATcenters)
  broadcast!(/, VarStats.High_LONcenters, VarStats.High_LONcenters, valweightsHigh)
  broadcast!(/ ,VarStats.High_TIMEcenters, VarStats.High_TIMEcenters, valweightsHigh)
  broadcast!(/ ,VarStats.Low_LATcenters, VarStats.Low_LATcenters, sweightsLow)
  #broadcast!(acosd ,VarStats.Low_LATcenters,VarStats.Low_LATcenters)
  broadcast!(/, VarStats.Low_LONcenters, VarStats.Low_LONcenters, valweightsLow)
  broadcast!(/ ,VarStats.Low_TIMEcenters, VarStats.Low_TIMEcenters, valweightsLow)
  for j = 1:length(Div)
    if Div[j] < min_nobs_within_thres
      VarStats.High_LATcenters[j] = NaN; VarStats.High_LONcenters[j] = NaN; VarStats.High_TIMEcenters[j] = NaN
    end
    if Div2[j] < min_nobs_within_thres
      VarStats.Low_LATcenters[j] = NaN; VarStats.Low_LONcenters[j] = NaN; VarStats.Low_TIMEcenters[j] = NaN
    end
  end
  return Centroids = (VarStats.High_LONcenters, VarStats.High_LATcenters, VarStats.High_TIMEcenters, VarStats.Low_LONcenters, VarStats.Low_LATcenters, VarStats.Low_TIMEcenters)
end

function spatDist(lon1, lon2, lat1, lat2)
    meanlat = (lat1 .+ lat2) / 2
    spatD = sqrt.((lon1 .- lon2) .* (lon1 .- lon2) + (cosd.(meanlat) .* (lat1 .- lat2)) .* (cosd.(meanlat) .* (lat1 .- lat2))) .* 111.324
    return spatD
end

export spatDist

function getDistances(VarStats::VariableStats)
  VarStats.SpatialDistance = spatDist(VarStats.Low_LONcenters, VarStats.High_LONcenters, VarStats.Low_LATcenters, VarStats.High_LATcenters)
  VarStats.TemporalDistance = (VarStats.High_TIMEcenters-VarStats.Low_TIMEcenters)
  return VarStats.SpatialDistance, VarStats.TemporalDistance # D_day negative: time of negative response > time of positive response, i.e. negative response after positive
end

export getDistances

function getPosNegSums(SumScoresDict::Array{Dict{Tuple{Int64,Int64},Float32},1})
    posres = zeros(eltype(values(SumScoresDict[1])), length(SumScoresDict))
    negres = zeros(eltype(values(SumScoresDict[1])), length(SumScoresDict))
    for lab = 1:length(SumScoresDict)
        if isassigned(SumScoresDict, lab)
            for v = values(SumScoresDict[lab])
              if !isnan(v)
                if v < 0
                    negres[lab] += v
                elseif v < 10^38
                    posres[lab] += v
                end
              end
            end
        end
    end
    bil = posres + negres
    return bil, posres, negres, posres ./ negres
end


function getBalance(VarStats::VariableStats, CCobj::CCstats, inMemVarAnomalies; multiplyByConst::Bool = true)
    SumScoresDict = get_ScoresSumAreaWeighted(CCobj, inMemVarAnomalies);
    VarStats.Balance, VarStats.PosSum, VarStats.NegSum, VarStats.Comp = getPosNegSums(SumScoresDict)
    if multiplyByConst
        kmconst = (111.324 * 1000 * abs(CCobj.Latitudes[1]-CCobj.Latitudes[2])) * (111.324 * 1000 * abs(CCobj.Longitudes[1]-CCobj.Longitudes[2])) * abs(CCobj.Times[2] - CCobj.Times[1])
        VarStats.Balance = VarStats.Balance .* kmconst
        VarStats.PosSum =  VarStats.PosSum .* kmconst
        VarStats.NegSum =  VarStats.NegSum .* kmconst
    end
    return VarStats.Balance, VarStats.NegSum, VarStats.PosSum, VarStats.Comp
end

export getBalance

function getWeightedSumOfValuesPerMask(CCobj::CCstats, landcover, xin)
    if size(landcover) != size(xin)[1:2] || size(CCobj.LabelCube) != size(xin) error("size not matching") end
    xout = Array{Dict{Int32, Float32}, 1}(undef, CCobj.numExtremes)
    concomp = CCobj.LabelCube
    if size(concomp) != size(xin) error("size of input arrays does not match") end
    if size(concomp)[1:2] != size(landcover) error("size of input arrays does not match") end
    for k = 1:size(concomp, 3)
        for j = 1:size(concomp, 2)
            for i = 1:size(concomp, 1)
                if(concomp[i, j, k] > 0)
                    label = concomp[i, j, k]
                    if !isnan(xin[i, j, k])
                        if isassigned(xout, label)
                            if(haskey(xout[label], landcover[i,j]))
                                xout[label][landcover[i,j]] = xout[label][landcover[i,j]] + xin[i, j, k] * CCobj.CosLatitudes[j]
                            else
                                xout[label][landcover[i,j]] = xin[i, j, k] * CCobj.CosLatitudes[j]
                            end
                        else
                            xout[label] = Dict(landcover[i,j] => xin[i, j, k] * CCobj.CosLatitudes[j])
                        end
                    end
                end
            end
        end
    end
    return(xout)
end

export getWeightedSumOfValuesPerMask


function getWeightedAverageOfValuesPerMask(CCobj::CCstats, landcover, xin)
    if size(landcover) != size(xin)[1:2] || size(CCobj.LabelCube) != size(xin) error("size not matching") end
    xout = Array{Dict{Int32, Float32}, 1}(undef, CCobj.numExtremes)
    xdiv = Array{Dict{Int32, Float32}, 1}(undef, CCobj.numExtremes)
    concomp = CCobj.LabelCube
    if size(concomp) != size(xin) error("size of input arrays does not match") end
    if size(concomp)[1:2] != size(landcover) error("size of input arrays does not match") end
    for k = 1:size(concomp, 3)
        for j = 1:size(concomp, 2)
            for i = 1:size(concomp, 1)
                if(concomp[i, j, k] > 0)
                    label = concomp[i, j, k]
                    if !isnan(xin[i, j, k])
                        if isassigned(xout, label)
                            if(haskey(xout[label], landcover[i,j]))
                                xout[label][landcover[i,j]] = xout[label][landcover[i,j]] + xin[i, j, k] * CCobj.CosLatitudes[j]
                                xdiv[label][landcover[i,j]] = xdiv[label][landcover[i,j]] + CCobj.CosLatitudes[j]
                            else
                                xout[label][landcover[i,j]] = xin[i, j, k] * CCobj.CosLatitudes[j]
                                xdiv[label][landcover[i,j]] = CCobj.CosLatitudes[j]
                            end
                        else
                            xout[label] = Dict(landcover[i,j] => xin[i, j, k] * CCobj.CosLatitudes[j])
                            xdiv[label] = Dict(landcover[i,j] => CCobj.CosLatitudes[j])
                        end
                    end
                end
            end
        end
    end
    for l = 1:length(xout)
        if isassigned(xout,l )
            for m = keys(xout[l])
                xout[l][m] = xout[l][m] / xdiv[l][m]
            end
        end
    end
    return(xout)
end

export getWeightedAverageOfValuesPerMask


function get_ScoresSumAreaWeighted(CCobj::CCstats, Scores)
    ScoresSum = Array{Dict{Tuple{Int,Int}, eltype(Scores)}, 1}(undef, CCobj.numExtremes)
    get_ScoresSumAreaWeighted!(ScoresSum, CCobj, Scores)
  return(ScoresSum)
end

get_ScoresSumAreaWeighted

function get_ScoresSumAreaWeighted!(xout::Array{Dict{Tuple{Int64,Int64},tp},1}, CCobj::CCstats, Scores::AbstractArray{tp, 3}) where {tp}
  if CCobj.LabelCube == zeros(Int, 1, 1, 1) error("empty CCobj.LabelCube") end
  #CCobj.affectedLonLatIdxCounts = Array(Dict{Tuple{Int,Int}, Int}, CCobj.numExtremes)
  #xout = CCobj.affectedLonLatIdxCounts
  for k = 1:size(CCobj.LabelCube, 3)
    for j = 1:size(CCobj.LabelCube, 2)
      for i = 1:size(CCobj.LabelCube, 1)
        l = CCobj.LabelCube[i,j,k]
        if l > 0
          if isassigned(xout, l)
            if haskey(xout[l], (i, j))
              xout[l][(i,j)] = xout[l][(i,j)] + Scores[i,j,k] * CCobj.CosLatitudes[j]
            else
              xout[l][(i,j)] = Scores[i,j,k] * CCobj.CosLatitudes[j]
            end
          else
            xout[l] = Dict((i, j) => Scores[i,j,k] * CCobj.CosLatitudes[j]) # lon -> i, lat -> j
          end
        end
      end
    end
  end
  return xout
end

 export get_ScoresSumAreaWeighted!
