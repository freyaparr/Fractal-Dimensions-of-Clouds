using ZarrDatasets, NCDatasets

# from https://github.com/eurec4a/eurec4a-intake/blob/master/Simulations/DALES/botany_dx100m_nx1536.yaml line 44
url = "https://swift.dkrz.de/v1/dkrz_b42c6bc6-439f-4efb-9491-db3de70ad9ba/botany_1536_radiation/botany_1536_radiation.zarr"
data = ZarrDataset(url; consolidated=true)

#need to somehow turn this into a .nc file i think


