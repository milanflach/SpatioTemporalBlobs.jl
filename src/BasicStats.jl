using Images

function init_CCobj()
  CCobj = CCstats(zeros(Int, 1,1,1),  # Bools
  zeros(Int, 1,1,1), # Connection Cube
  trues(3,3,3),
  0, # numExtremes
  0, # minAffectedVoxels
  [0], # labels
  [0], # perm
  [0], # affV
  [0], # affA
  [0], # affT
  [0.0], #affV
  [0.0], #affA [km^2]
  [0.0], # affT [days]
  [0.0], # Tmin
  [0.0], # Tmax
  [0.0], # Tcenters
  [0.0], # LATmin
  [0.0], # LATmax
  [0.0], # LATcenters
  [0.0], # LONmin
  [0.0], # LONmax
  [0.0], # LONcenters
  500.0:0.00001:500.0,  # Latitudes
  500.0:0.00001:500.0, # Longitudes
  [0.0], #Times
  [0.0],
  [0.0],
  Array{Dict{Tuple{Int64,Int64},Int64},1}(undef,0),
  Array{Dict{Tuple{Int64,Int64,Int64},Float64},1}(undef,0),
  Dict{String, VariableStats}() # Variable Statistics
  );
  return(CCobj)
end

export init_CCobj


function get_affectedTimes(CCobj::CCstats, TimeStepLength)
  CCobj.affectedTime = zeros(Float64, length(CCobj.Tmins))
  for i = 1:length(CCobj.Tmins)
   CCobj.affectedTime[i] = Float64(CCobj.Tmaxs[i] - CCobj.Tmins[i]) / (3600 * 24 * 1000) + TimeStepLength
  end
  CCobj.affectedTime
end

function getLonLatTimebounds(CCobj::CCstats)
  (minbounds, maxbounds) = getBounds(CCobj.LabelCube, CCobj.numExtremes)

for j = 1:3
  replace!(i -> i > size(CCobj.LabelCube, j) ? size(CCobj.LabelCube, j) : i, view(minbounds,j,:))
  replace!(i -> i < 1 ? 1 : i, view(minbounds,j,:))
  replace!(i -> i > size(CCobj.LabelCube, j) ? size(CCobj.LabelCube, j) : i, view(maxbounds,j,:))
  replace!(i -> i < 1 ? 1 : i, view(maxbounds,j,:))
end

  #println(extrema(minbounds[1,:]))
  CCobj.LONmins = CCobj.Longitudes[minbounds[1,:]]
  CCobj.LATmins = CCobj.Latitudes[minbounds[2,:]]
  CCobj.Tmins = CCobj.Times[minbounds[3,:]]

  CCobj.LONmaxs = CCobj.Longitudes[maxbounds[1,:]]
  CCobj.LATmaxs = CCobj.Latitudes[maxbounds[2,:]]
  CCobj.Tmaxs = CCobj.Times[maxbounds[3,:]]

  CCobj.affectedTimeSteps = maxbounds[3,:] .- minbounds[3,:] .+ 1

  return(minbounds, maxbounds)
end


function countVoxels(LabelCube::AbstractArray{<:Integer, N}, numExtremes::Int = 0) where {N}
  if numExtremes == 0
    numExtremes = maximum(LabelCube)
  end
  lAr=zeros(Int,numExtremes)
  for i=1:length(LabelCube)
    j=LabelCube[i]
    if j>0
      lAr[j]=lAr[j]+1
    end
  end
  return lAr
end

# get latitude index of a floatrange of latitudes
function getLatidx(Latitudes::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}, lat::Float64)
    latidx = 0
    counter = 0
    for i in Latitudes
    counter += 1
        if i == lat
            latidx = counter
        end
    end
    if latidx == 0 error("lat = $lat not in Latitudes = $Latitudes") end
    return latidx
end

# get pixel Area of a given latitude
function pixelArea(CCobj::CCstats, lat::Float64)
    latidx = getLatidx(CCobj.Latitudes, lat)
    km_const = 111.324 #2pi / 360 * 6378.388 (Earth radius in km)
    latgridsize = abs(convert(Float64, CCobj.Latitudes.step))
    longridsize = abs(convert(Float64, CCobj.Longitudes.step))
    dlon = km_const * longridsize
    dlat = km_const * CCobj.CosLatitudes[latidx] * latgridsize
    Area = dlon * dlat
    return Area
end

function pixelArea(CCobj::CCstats, latidx::Int)
    #latidx = getLatidx(CCobj.Latitudes, lat)
    km_const = 111.324 #2pi / 360 * 6378.388 (Earth radius in km)
    latgridsize = abs(convert(Float64, CCobj.Latitudes.step))
    longridsize = abs(convert(Float64, CCobj.Longitudes.step))
    dlon = km_const * longridsize
    dlat = km_const * CCobj.CosLatitudes[latidx] * latgridsize
    Area = dlon * dlat
    return Area
end

function get_affectedLonLatCounter(CCobj::CCstats)
  CCobj.affectedLonLatIdxCounts = Array{Dict{Tuple{Int,Int}, Int},1}(undef, CCobj.numExtremes)
  get_affectedLonLatCounter!(CCobj.affectedLonLatIdxCounts, CCobj)
  return(CCobj.affectedLonLatIdxCounts)
end

export get_affectedLonLatCounter

function get_affectedLonLatCounter!(xout::Array{Dict{Tuple{Int,Int}, Int}, 1}, CCobj::CCstats)
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
              xout[l][(i,j)] += 1
            else
              xout[l][(i,j)] = 1
            end
          else
            xout[l] = Dict((i, j) => 1) # lon -> i, lat -> j
          end
        end
      end
    end
  end
  for i=1:length(xout)
    isassigned(xout,i) || (xout[i]=Dict((1,1)=>1))
  end
  return xout
end

export get_affectedLonLatCounter!

function get_affectedLonLatTimeScoreValues(CCobj::CCstats, Scores::AbstractArray{<:AbstractFloat, 3})
  CCobj.affectedLonLatTimeScoreValue = Array{Dict{Tuple{Int,Int,Int}, Float64},1}(undef, CCobj.numExtremes)
  get_affectedLonLatTimeScoreValues!(CCobj.affectedLonLatTimeScoreValue, CCobj.LabelCube, Scores)
  return(CCobj.affectedLonLatTimeScoreValue)
end

export get_affectedLonLatTimeScoreValues

function get_affectedLonLatTimeScoreValues!(xout::Array{Dict{Tuple{Int,Int,Int}, T},1}, LabelCube::AbstractArray{<:Integer,3}, Scores::AbstractArray{T, 3}; tshift::Int = 0) where {T}
  if LabelCube == zeros(Int, 1, 1, 1) error("empty CCobj.LabelCube") end
  for k = 1:size(LabelCube, 3)
    for j = 1:size(LabelCube, 2)
      for i = 1:size(LabelCube, 1)
        l = LabelCube[i,j,k]
        if l > 0
          if k+tshift <= size(LabelCube, 3) && k+tshift > 0
            if isassigned(xout, l)
              if haskey(xout[l], (i, j, k))
                error("double defined $i $j $k")
              else
                xout[l][(i,j,k)] = Scores[i,j,k+tshift]
              end
            else
              xout[l] = Dict((i, j,k) => Scores[i,j,k+tshift]) # lon -> i, lat -> j, times -> k
            end
          end
        end
      end
    end
  end
  return xout
end

export get_affectedLonLatTimeScoreValues!



function get_affectedArea!(CCobj::CCstats)

  CCobj.affectedLonLatIdxCounts = get_affectedLonLatCounter(CCobj);

  CCobj.affectedPixels = zeros(Int, length(CCobj.affectedLonLatIdxCounts))
  CCobj.affectedArea = zeros(Float64, length(CCobj.affectedLonLatIdxCounts))
  for i = 1:length(CCobj.affectedLonLatIdxCounts)
    CCobj.affectedPixels[i] = length(CCobj.affectedLonLatIdxCounts[i])
    for lonlats = keys(CCobj.affectedLonLatIdxCounts[i])
      (lonidx, latidx) = lonlats
      CCobj.affectedArea[i] = CCobj.affectedArea[i] + CCobj.AreaLatitudes[latidx]
    end
  end
  return (CCobj.affectedArea, CCobj.affectedPixels)
end

function get_affectedAreaAndVolume!(CCobj::CCstats)

  CCobj.affectedLonLatIdxCounts = get_affectedLonLatCounter(CCobj);

  #@show CCobj.affectedLonLatIdxCounts
  timesteplength = median(CCobj.Times[2:end]-CCobj.Times[1:(end-1)])

  CCobj.affectedPixels = zeros(Int, length(CCobj.affectedLonLatIdxCounts))
  CCobj.affectedArea = zeros(Float64, length(CCobj.affectedLonLatIdxCounts))
  CCobj.affectedVolume = zeros(Float64, length(CCobj.affectedLonLatIdxCounts))
  CCobj.affectedVoxels = zeros(Int, length(CCobj.affectedLonLatIdxCounts))
  for i = 1:length(CCobj.affectedLonLatIdxCounts)
    CCobj.affectedPixels[i] = length(CCobj.affectedLonLatIdxCounts[i])
    for lonlats = keys(CCobj.affectedLonLatIdxCounts[i])
      (lonidx, latidx) = lonlats
      #print(" -- lonidx = $lonidx, latidx = $latidx")
      CCobj.affectedArea[i] = CCobj.affectedArea[i] + CCobj.AreaLatitudes[latidx]
      CCobj.affectedVoxels[i] = CCobj.affectedVoxels[i] + CCobj.affectedLonLatIdxCounts[i][lonlats]
      CCobj.affectedVolume[i] = CCobj.affectedVolume[i] + CCobj.AreaLatitudes[latidx] * CCobj.affectedLonLatIdxCounts[i][lonlats] * timesteplength
    end
  end
  return (CCobj.affectedArea, CCobj.affectedPixels, CCobj.affectedVoxels , CCobj.affectedVolume)
end

export get_affectedAreaAndVolume!

# get fill CCobj.AreaLatitudes wit the calculated areas
function getAreaLatitudes(CCobj::CCstats)
    CCobj.AreaLatitudes = zeros(Float64, length(CCobj.Latitudes))
    for i = 1:length(CCobj.Latitudes)
        CCobj.AreaLatitudes[i] = pixelArea(CCobj, i)
    end
    return(CCobj.AreaLatitudes)
end

function countVoxels_weighted(CCobj::CCstats)
  if CCobj.numExtremes == 0
    CCobj.numExtremes = maximum(CCobj.LabelCube)
  end
  lAr=zeros(Float64,CCobj.numExtremes)
  for t=1:size(CCobj.LabelCube, 3)
    for la=1:size(CCobj.LabelCube, 2)
      for lo=1:size(CCobj.LabelCube, 1)
        j=CCobj.LabelCube[lo, la, t]
        if j>0
          if t+1 > length(CCobj.Times)
            timesteplength = Float64(DateTime(year(CCobj.Times[t])+1,1,1)-CCobj.Times[t]) / (24 * 3600 * 1000)
          else
            timesteplength = Float64(CCobj.Times[t+1]-CCobj.Times[t]) / (24 * 3600 * 1000)
          end
          lAr[j] = lAr[j] + (CCobj.AreaLatitudes[la] * timesteplength)
        end
      end
    end
  end
  return lAr
end

function get_affectedVolume!(CCobj::CCstats)
  if CCobj.numExtremes == 0
    CCobj.numExtremes = maximum(CCobj.LabelCube)
  end
  timesteplength = median(CCobj.Times[2:end]-CCobj.Times[1:(end-1)])
  lAr=zeros(Float64,CCobj.numExtremes)
  for t=1:size(CCobj.LabelCube, 3)
    for la=1:size(CCobj.LabelCube, 2)
      for lo=1:size(CCobj.LabelCube, 1)
        j=CCobj.LabelCube[lo, la, t]
        if j>0
          lAr[j] = lAr[j] + (CCobj.AreaLatitudes[la] * timesteplength)
        end
      end
    end
  end
  CCobj.affectedVolume = lAr
  return CCobj.affectedVolume
end

export get_affectedVolume!

function any_loop(sel_labels::AbstractArray{tp, N}, testlabel::tp) where {tp, N}
    @inbounds for j = 1:length(sel_labels)
        if sel_labels[j] == testlabel
            return  true
        end
    end
    return false
end

# loop wrapper around any(sel_labels .!= testlabel)
function any_notloop(sel_labels::AbstractArray{tp, N}, testlabel::tp) where {tp, N}
    @inbounds for j = 1:length(sel_labels)
        if sel_labels[j] != testlabel
            return  true
        end
    end
    return false
end

# Set Bools to zero which are not within the sel_labels of the LabelCube
function clean_bools!(Bools::Array{<:Integer, N}, LabelCube::Array{<:Integer, N}, sel_labels::Array{<:Integer, 1}) where {N}
    for i = 1:length(Bools)
        if Bools[i] != 0
            if !any_loop(sel_labels, LabelCube[i])
                Bools[i] = 0
            end
        end
    end
    return(Bools)
end

# get a connection cube cleaned for the desired minimum of affected voxels (default) or the desired number of extremes
function getLabelCube!(Bools::Array{<:Integer, N}, ConnectionFootprint::Array{<:Integer, N}; mode::String = "Voxels", minAffectedVoxels::Int = 1, numExtremes::Int = 100) where {N}
    #ConnectionFootprint = trues(3,3,3) # consider all 3 * 3 * 3 = 27 connections
    # first run
    CC = label_components(Bools, ConnectionFootprint)
    affVoxels = countVoxels(CC)
    # second run to exclude low connectivity
    if mode == "Voxels"
        sel_labels = find(affVoxels .>= minAffectedVoxels)
    else
        perm = sortperm(affVoxels, rev = true)
        sel_labels = perm[1:numExtremes]
    end
    clean_bools!(Bools, CC, sel_labels)
    CC2 = label_components(Bools, ConnectionFootprint)
    return(Bools, CC2)
end

export getLabelCube!


function getLabelCube!(CCobj::CCstats)
    if CCobj.numExtremes == 0 && CCobj.minAffectedVoxels == 0 error("either numExtremes or minAffectedVoxels has to be != 0") end
    if CCobj.Bools == zeros(Int, 1,1,1) error("you forgot to write your Boolean cube to CCobj.Bools") end
    if all(CCobj.ConnectionFootprint .== 0)  error("you do not consider any connections in CCobj.ConnectionFootprint") end
    #ConnectionFootprint = trues(3,3,3) # consider all 3 * 3 * 3 = 27 connections
    # first run
    print("first run to compute connected components \n")
    @time CCobj.LabelCube = label_components(CCobj.Bools, CCobj.ConnectionFootprint)
    print("compute affected voxels \n")
    @time affVoxels = countVoxels(CCobj.LabelCube)
    print("get permutation \n")
    perm = sortperm(affVoxels, rev = true)
    # second run to exclude low connectivity
    if CCobj.minAffectedVoxels != 0
      print("select only events which are connected to more than $(CCobj.minAffectedVoxels) voxels \n")
        @time sel_labels = findall(affVoxels .>= CCobj.minAffectedVoxels)
    else
      print("select the $(CCobj.numExtremes) Events with the largest number of affected Voxels \n")
        @time begin
        sel_labels = perm[1:CCobj.numExtremes]
        end
    end
    CCobj.perm=perm[1:length(sel_labels)]
    print("remove the other events from Bools \n")
    @time clean_bools!(CCobj.Bools, CCobj.LabelCube, sel_labels)
    print("clean and sort label cube with cleaned bools \n")
    @time CCobj.LabelCube = label_components(CCobj.Bools, CCobj.ConnectionFootprint)
    #cleanAndSort_LabelCube!(CCobj, perm)
    CCobj.numExtremes = maximum(CCobj.LabelCube)
    return(CCobj.LabelCube)
end

export getLabelCube!

function cleanAndSort_LabelCube!(CCobj::CCstats, perm::Array{<:Integer, 1})
  #newlabelCube = similar(CCobj.LabelCube)
  for i = 1:length(CCobj.LabelCube)
    if CCobj.Bools[i] == 0 # Bools are already cleaned --> if zero label cube is also zero
      CCobj.LabelCube[i] = 0
    else
      CCobj.LabelCube[i] = findfirst(isequal(CCobj.LabelCube[i]), perm)
    end
  end
  return CCobj.LabelCube
end


# get minimum and maximum values of each dimension for each label of the LabelCube
function getBounds!(minbounds, maxbounds, LabelCube::Array{<:Integer, N}, numExtremes::Int = 0) where {N}
    if numExtremes == 0
     numExtremes = maximum(LabelCube)
    end
    for k = 1:size(LabelCube, 3)
        for j = 1:size(LabelCube, 2)
            for i = 1:size(LabelCube, 1)
                signat = LabelCube[i,j,k]
                if signat > 0
                    if minbounds[1, signat] > i  minbounds[1, signat] = i end
                    if minbounds[2, signat] > j  minbounds[2, signat] = j end
                    if minbounds[3, signat] > k  minbounds[3, signat] = k end
                    if maxbounds[1, signat] < i  maxbounds[1, signat] = i end
                    if maxbounds[2, signat] < j  maxbounds[2, signat] = j end
                    if maxbounds[3, signat] < k  maxbounds[3, signat] = k end
                end
            end
        end
    end
    return(minbounds, maxbounds)
end

function getBounds(LabelCube::Array{<:Integer, N}, numExtremes::Int = 0) where {N}
    if numExtremes == 0
     numExtremes = maximum(LabelCube)
    end
    minbounds=fill(maximum(size(LabelCube)) + 1,3,numExtremes)
    maxbounds=zeros(Int,3,numExtremes)
    getBounds!(minbounds, maxbounds, LabelCube, numExtremes)
    return(minbounds, maxbounds)
end

function getConCompStats(CCobj::CCstats)
  if CCobj.Times == [0.0] error("CCobj.Times is empty") end
  if minimum(CCobj.Latitudes) == 500.0 error("CCobj.Latitudes is undefined") end
  if minimum(CCobj.Longitudes) == 500.0 error("CCobj.Longitudes is undefined") end
  if CCobj.LabelCube == zeros(Int, 1,1,1) error("CCobj.LabelCube is empty") end

  CCobj.CosLatitudes = cos.(collect(CCobj.Latitudes) / 360.0 * 2pi)
  CCobj.AreaLatitudes = getAreaLatitudes(CCobj)

  CCobj.numExtremes = maximum(CCobj.LabelCube)

  #print("get affected voxels, volume, area, pixels \n")
  get_affectedAreaAndVolume!(CCobj)

  #print("get permutattion according to volume \n")
  CCobj.perm = sortperm(CCobj.affectedVolume, rev = true);
  CCobj.labels = collect(1:CCobj.numExtremes)

  #print("get minima and maxima of lon lat time \n")
  getLonLatTimebounds(CCobj)

  #print("get duration")
  TimeStepLength = median(CCobj.Times[2:end]-CCobj.Times[1:(end-1)])
  get_affectedTimes(CCobj, TimeStepLength)

  return CCobj
end

export getConCompStats


function getCentroids(CCobj::CCstats)
  if length(CCobj.affectedLonLatTimeScoreValue) == 0 error("CCobj.affectedLonLatTimeScoreValue is empty") end
  CCobj.LONcenters = zeros(Float64, CCobj.numExtremes)
  CCobj.LATcenters = zeros(Float64, CCobj.numExtremes)
  CCobj.Tcenters = zeros(Float64, CCobj.numExtremes)
  getCentroids!(CCobj.LONcenters, CCobj.LATcenters, CCobj.Tcenters, CCobj.affectedLonLatTimeScoreValue, CCobj)
  return CCobj.LONcenters, CCobj.LATcenters, CCobj.Tcenters
end

export getCentroids


function getCentroids!(LonCenters::Array{Float64, 1}, LatCenters::Array{Float64, 1}, TimeCenters::Array{Float64, 1}, LonLatTime::Array{Dict{Tuple{Int64,Int64,Int64},Float64},1}, CCobj::CCstats; ValueWeighting::Bool = true)
  sweights=zeros(Float64, CCobj.numExtremes)
  valweights=zeros(Float64, CCobj.numExtremes)
  for i = 1:length(LonLatTime)
    if isassigned(LonLatTime, i)
        for lonlattimeidx = keys(LonLatTime[i])
          if ValueWeighting
            val = LonLatTime[i][lonlattimeidx]
          else
            val = 1
          end
          innerstCentroids!(sweights, valweights, LatCenters, LonCenters, TimeCenters, CCobj, lonlattimeidx, val, i)
        end
    end
  end
  broadcast!(/ ,LatCenters, LatCenters, sweights)
  broadcast!(/, LonCenters, LonCenters, valweights)
  broadcast!(/ ,TimeCenters, TimeCenters, valweights)
  return LonCenters, LatCenters, TimeCenters
end

export getCentroids!

function innerstCentroids!(sweights, valweights, LatCenters, LonCenters, TimeCenters, CCobj, lonlattimeidx, val_orig, i)
  (lonidx, latidx, timeidx) = lonlattimeidx
  # weighting only works with positive values
  val = abs(val_orig)
  if lonidx < 1  lonidx = 1 end
  if latidx < 1  latidx = 1 end
  if timeidx < 1  timeidx = 1 end
  if lonidx > size(CCobj.LabelCube, 1)  lonidx = size(CCobj.LabelCube, 1) end
  if latidx > size(CCobj.LabelCube, 2)  latidx = size(CCobj.LabelCube, 2) end
  if timeidx > size(CCobj.LabelCube, 3)  timeidx = size(CCobj.LabelCube, 3) end
  LatCenters[i] = LatCenters[i] + CCobj.Latitudes[latidx]*CCobj.CosLatitudes[latidx]*val
  sweights[i]  += CCobj.CosLatitudes[latidx]*val
  LonCenters[i] = LonCenters[i] + CCobj.Longitudes[lonidx]*val
  TimeCenters[i] = TimeCenters[i] + CCobj.Times[timeidx]*val
  valweights[i] += val
  return sweights[i], valweights[i], LatCenters[i], LonCenters[i], TimeCenters[i]
end
