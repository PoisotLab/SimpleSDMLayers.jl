"""
    SharedSocioeconomicPathway

Enumeration of the four SSPs, which can be listed with
`instances(SharedSocioeconomicPathway)`. These are meant to be used with `CMIP6`
models.
"""
@enum SharedSocioeconomicPathway SSP126 SSP245 SSP370 SSP585

"""
    RepresentativeConcentrationPathway

Enumeration of the four RCPs, which can be listed with
`instances(RepresentativeConcentrationPathway)`. These are meant to be used with
`CMIP5` models.
"""
@enum RepresentativeConcentrationPathway RCP26 RCP45 RCP60 RCP85

"""
    CMIP6

Enumeration of the models from CMIP6, which can be listed with
`instances(CMIP6)`.
"""
@enum CMIP6 BCCCSM2MR CNRMCMS61 CNRMESM21 CanESM5 GFDLESM4 IPSLCM6ALR MIROCES2L MIROC6 MRIESM20 UKESM10LL MPIESM12HR

"""
    CMIP5

Enumeration of the models from CMIP5, which can be listed with
`instances(CMIP5)`.
"""
@enum CMIP5 ACCESS10 BNUESM CCSM4 CESM1BGC CESM1CAM5 CMCCCMS CMCCCM CNRMCM5 CSIROMK360 CanESM2 FGOALSG2 FIOESM GFDLCM3 GFDLESM2G GFDLESM2M GISSE2HCC GISSE2H GISSE2RCC GISSE2R HADGEM2AO HADGEM2CC IPSLCM5ALR IPSLCM5AMR MIROCESMCHEM MIROCESM MIROC5 MPIESMLR MPIESMMR MRICGCM3 MRIESM1 NORESM1M BCCCSM11 INMCM4

# Future paths
_rasterpath(model::CMIP6) = _rasterpath(Val{model})
_rasterpath(model::CMIP5) = _rasterpath(Val{model})
_rasterpath(ssp::SharedSocioeconomicPathway) = _rasterpath(Val{ssp})
_rasterpath(rcp::RepresentativeConcentrationPathway) = _rasterpath(Val{rcp})

# SSP path
_rasterpath(::Type{Val{SSP126}}) = "ssp126"
_rasterpath(::Type{Val{SSP245}}) = "ssp245"
_rasterpath(::Type{Val{SSP370}}) = "ssp370"
_rasterpath(::Type{Val{SSP585}}) = "ssp585"

# RCP path
_rasterpath(::Type{Val{RCP26}}) = "rcp26"
_rasterpath(::Type{Val{RCP45}}) = "rcp45"
_rasterpath(::Type{Val{RCP60}}) = "rcp60"
_rasterpath(::Type{Val{RCP85}}) = "rcp85"

# CMIP6 path
_rasterpath(::Type{Val{BCCCSM2MR}}) = "BCC-CSM2-MR"
_rasterpath(::Type{Val{CNRMCMS61}}) = "CNRM-CM6-1"
_rasterpath(::Type{Val{CNRMESM21}}) = "CNRM-ESM2-1"
_rasterpath(::Type{Val{CanESM5}}) = "CanESM5"
_rasterpath(::Type{Val{GFDLESM4}}) = "GFDL-ESM4"
_rasterpath(::Type{Val{IPSLCM6ALR}}) = "IPSL-CM6A-LR"
_rasterpath(::Type{Val{MIROCES2L}}) = "MIROC-ES2L"
_rasterpath(::Type{Val{MIROC6}}) = "MIROC6"
_rasterpath(::Type{Val{MRIESM20}}) = "MRI-ESM2-0"
_rasterpath(::Type{Val{UKESM10LL}}) = "UKESM1-0-LL"
_rasterpath(::Type{Val{MPIESM12HR}}) = "MPI-ESM1-2-HR"

# CMIP5 path
_rasterpath(::Type{Val{ACCESS10}}) = "ACCESS1-0"
_rasterpath(::Type{Val{BNUESM}}) = "BNU-ESM"
_rasterpath(::Type{Val{CCSM4}}) = "CCSM4"
_rasterpath(::Type{Val{CESM1BGC}}) = "CESM1-BGC"
_rasterpath(::Type{Val{CESM1CAM5}}) = "CESM1-CAM5"
_rasterpath(::Type{Val{CMCCCMS}}) = "CMCC-CMS"
_rasterpath(::Type{Val{CMCCCM}}) = "CMCC-CM"
_rasterpath(::Type{Val{CNRMCM5}}) = "CNRM-CM5"
_rasterpath(::Type{Val{CSIROMK360}}) = "CSIRO-Mk3-6-0"
_rasterpath(::Type{Val{CanESM2}}) = "CanESM2"
_rasterpath(::Type{Val{FGOALSG2}}) = "FGOALS-g2"
_rasterpath(::Type{Val{FIOESM}}) = "FIO-ESM"
_rasterpath(::Type{Val{GFDLCM3}}) = "GFDL-CM3"
_rasterpath(::Type{Val{GFDLESM2G}}) = "GFDL-ESM2G"
_rasterpath(::Type{Val{GFDLESM2M}}) = "GFDL-ESM2M"
_rasterpath(::Type{Val{GISSE2HCC}}) = "GISS-E2-H-CC"
_rasterpath(::Type{Val{GISSE2H}}) = "GISS-E2-H"
_rasterpath(::Type{Val{GISSE2RCC}}) = "GISS-E2-R-CC"
_rasterpath(::Type{Val{GISSE2R}}) = "GISS-E2-R"
_rasterpath(::Type{Val{HADGEM2AO}}) = "HadGEM2-AO"
_rasterpath(::Type{Val{HADGEM2CC}}) = "HadGEM2-CC"
_rasterpath(::Type{Val{IPSLCM5ALR}}) = "IPSL-CM5A-LR"
_rasterpath(::Type{Val{IPSLCM5AMR}}) = "IPSL-CM5A-MR"
_rasterpath(::Type{Val{MIROCESMCHEM}}) = "MIROC-ESM-CHEM"
_rasterpath(::Type{Val{MIROCESM}}) = "MIROC-ESM"
_rasterpath(::Type{Val{MIROC5}}) = "MIROC5"
_rasterpath(::Type{Val{MPIESMLR}}) = "MPI-ESM-LR"
_rasterpath(::Type{Val{MPIESMMR}}) = "MPI-ESM-MR"
_rasterpath(::Type{Val{MRICGCM3}}) = "MRI-CGCM3"
_rasterpath(::Type{Val{MRIESM1}}) = "MRI-ESM1"
_rasterpath(::Type{Val{NORESM1M}}) = "NorESM1-M"
_rasterpath(::Type{Val{BCCCSM11}}) = "bcc-csm1-1"
_rasterpath(::Type{Val{INMCM4}}) = "inmcm4"
