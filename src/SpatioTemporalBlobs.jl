module SpatioTemporalBlobs

import StatsBase.median

mutable struct VariableStats
  Histogram::Array{Float64,2}
  ContrastingResponse::Array{Int, 1}
  Low_LATcenters::Array{Float64, 1}
  Low_LONcenters::Array{Float64, 1}
  Low_TIMEcenters::Array{Float64, 1}
  High_LATcenters::Array{Float64, 1}
  High_LONcenters::Array{Float64, 1}
  High_TIMEcenters::Array{Float64, 1}
  SpatialDistance::Array{Float64,1}
  TemporalDistance::Array{Float64,1}
  Balance::Array{Float32,1}
  PosSum::Array{Float32,1}
  NegSum::Array{Float32,1}
  Comp::Array{Float32,1}
  affectedLonLatTimeValue::Array{Dict{Tuple{Int64,Int64,Int64},Int64},1}
end

mutable struct CCstats
  Bools::Array{UInt8, 3}
  LabelCube::Array{Int, 3} # --> LabelCube
  ConnectionFootprint::Array{Bool, 3} # trues(3 ,3,3)
  numExtremes::Int
  minAffectedVoxels::Int
  labels::Array{Int, 1}
  perm::Array{Int,1}
  affectedVoxels::Array{Int, 1} # Volume [counts]
  affectedPixels::Array{Int, 1} # Area [counts]
  affectedTimeSteps::Array{Int, 1} # Duration [counts]
  affectedVolume::Array{Float64, 1} # Volume [days * km^2]
  affectedArea::Array{Float64, 1} # Area [km^2]
  affectedTime::Array{Float64, 1} # Duration [days]
  Tmins::Array{Float64, 1}
  Tmaxs::Array{Float64, 1}
  Tcenters::Array{Float64, 1}
  LATmins::Array{Float64, 1}
  LATmaxs::Array{Float64, 1}
  LATcenters::Array{Float64, 1}
  LONmins::Array{Float64, 1}
  LONmaxs::Array{Float64, 1}
  LONcenters::Array{Float64, 1}
  Latitudes::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
  Longitudes::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
  Times::Array{Float64, 1}
  CosLatitudes::Array{Float64, 1}
  AreaLatitudes::Array{Float64, 1}
  affectedLonLatIdxCounts::Array{Dict{Tuple{Int64,Int64},Int64},1}
  affectedLonLatTimeScoreValue::Array{Dict{Tuple{Int64,Int64,Int64},Float64},1}
  VarStats::Dict{String,VariableStats}
end

include("BasicStats.jl")

include("VariableStats.jl")

include("IO.jl")

end # module
