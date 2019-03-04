using JLD
using NetCDF

function Array2NetCDF(varname::String, outputname::String, x::AbstractArray{tp, 3}; overwrite::Bool = false) where {tp}
  if overwrite
    if isfile(outputname)
    rm(outputname)
    end
  end
  nccreate(outputname, varname,"dim1_size$(size(x,1))",size(x,1),"dim2_size$(size(x,2))",size(x, 2),"dim3_size$(size(x,3))",size(x,3),t=eltype(x))
  ncwrite(x,outputname, varname)
  ncclose(outputname)
end


function Array2NetCDF(varname::String, outputname::String, x::AbstractArray; overwrite::Bool = false)
  if overwrite
    if isfile(outputname)
      rm(outputname)
    end
  end
  nccreate(outputname,varname,string(varname,"_array__length$(length(x))"),length(x),t=eltype(x))
  ncwrite(reshape(x, length(x)),outputname,varname)
  ncclose(outputname)
end

function Array2NetCDF(varname::String, outputname::String, x::Array{Bool, 3}; overwrite::Bool = false)
  if overwrite
    if isfile(outputname)
      rm(outputname)
    end
  end
  nccreate(outputname, varname,"dim1_size$(size(x,1))",size(x,1),"dim2_size$(size(x,2))",size(x, 2),"dim3_size$(size(x,3))",size(x,3),t=Int32)
  ncwrite(broadcast(Int32, x),outputname,varname)
  ncclose(outputname)
end

function Array2NetCDF(varname::String, outputname::String, x::StepRangeLen; overwrite::Bool = false)
  if overwrite
    if isfile(outputname)
      rm(outputname)
    end
  end
  nccreate(outputname,varname,string(varname,"_array_length$(length(x))"),length(x),t=eltype(x))
  ncwrite(collect(x),outputname,varname)
  ncclose(outputname)
end

function Array2NetCDF(varname::String, outputname::String, x::tp; overwrite::Bool = false) where {tp}
  if overwrite
    if isfile(outputname)
      rm(outputname)
    end
  end
  nccreate(outputname,varname,string(varname,"_array_length$(length(x))"),length(x),t=eltype(x))
  ncwrite([x],outputname,varname)
  ncclose(outputname)
end


function CCobj2NetCDF(outputname, CCobj; overwrite = true, endskip::Int = 3)
  if overwrite
    if isfile(outputname)
      rm(outputname)
    end
  end
  for i = 1:(length(propertynames(CCobj))-endskip)
    varname = string(propertynames(CCobj)[i])
    print("\n write $(varname) ...")
    x = getfield(CCobj, propertynames(CCobj)[i])
    Array2NetCDF(varname, outputname, x, overwrite = false)
  end
  return "-- done -- \n"
end

export CCobj2NetCDF

function NetCDF2CCobj(inputname; endskip::Int = 3)
  CCobj = init_CCobj()
  NetCDF2CCobj!(CCobj, inputname, endskip = endskip)
  return CCobj
end

export NetCDF2CCobj

function NetCDF2CCobjVarStats(inputname, CCobj::CCstats, variable_name; endskip::Int = 1, startskip::Int = 1)
  CCobj.VarStats[variable_name] = init_VariableStats()
  NetCDF2CCobj!(CCobj.VarStats[variable_name], inputname, endskip = endskip, startskip = startskip)
  return CCobj
end

export NetCDF2CCobjVarStats


function NetCDF2CCobj!(CCobj, inputname; endskip::Int = 3, startskip::Int = 0 )
  for i = (1+startskip):(length(propertynames(CCobj))-endskip)
    varname = string(propertynames(CCobj)[i])
    print("\n read $(varname) ...")
    x = ncread(inputname, varname)
    if length(x) == 1 || varname == "Latitudes" || varname == "Longitudes"
        x = broadcast(eltype(getfield(CCobj, propertynames(CCobj)[i])), x)
    else
        x = convert(typeof(getfield(CCobj, propertynames(CCobj)[i])), x)
    end
    if length(x) == 1
      setfield!(CCobj, propertynames(CCobj)[i], x[1])
    elseif varname == "Latitudes"
      setfield!(CCobj, propertynames(CCobj)[i], maximum(x):round(x[2]-x[1], digits=5):minimum(x))
    elseif varname == "Longitudes"
      setfield!(CCobj, propertynames(CCobj)[i], minimum(x):round(x[2]-x[1], digits=5):maximum(x))
    else
      setfield!(CCobj, propertynames(CCobj)[i], x)
    end
  end
  ncclose()
  return CCobj
end

export NetCDF2CCobj!
