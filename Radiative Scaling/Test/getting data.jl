using ZarrDatasets, NCDatasets

# from https://github.com/eurec4a/eurec4a-intake/blob/master/Simulations/DALES/botany_dx100m_nx1536.yaml line 44
url = "https://swift.dkrz.de/v1/dkrz_b42c6bc6-439f-4efb-9491-db3de70ad9ba/botany_1536_radiation/botany_1536_radiation.zarr"
data = ZarrDataset(url; consolidated=true)

#need to somehow turn this into a .nc file i think

#need to actually create the file
NCDataset("dales_radiation_rlut_data_60hours.nc", "c") do nc
    x, y = data["xt"][:], data["yt"][:]
    #and its dimensions (so the 1536x1536 grid and then the 60 hours)
    defDim(nc, "xt", length(x))
    defDim(nc, "yt", length(y))
    defDim(nc, "time", 60)
    #save them all 
    defVar(nc, "xt", x, ("xt",))
    defVar(nc, "yt", y, ("yt",))
    defVar(nc, "time", Float64.(0:59), ("time",); attrib=Dict("units" => "hours since 2020-02-01 01:00:00"))
    #radiative flux
    radiation = defVar(nc, "rlut", Float32, ("xt", "yt", "time"); deflatelevel=1, chunksizes=[256, 256, 1])

    #download the 60 radiation fields 
    for t in 1:60
        radiation[:, :, t] = coalesce.(data["rlut"][:, :, t, 3], NaN32)
    end
end