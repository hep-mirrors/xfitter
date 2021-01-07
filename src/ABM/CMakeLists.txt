add_library(ABM SHARED
  alphasvfn.f
  asy-hpcc.f
  asympcoef.f
  CMakeLists.txt
  dishq.f
  dislt.f
  disqpm.f
  f2ccbmsn.f
  gauss.f
  grid.F
  hplog.f
  hqcoef.f
  hqnnlocoef.f
  initgridconst.F
  lpcoef.f
  ome.f
  s2nlo.f
  sf_abkm_wrap.f
  spline.f
  split.f
  WGPLG.f
)
target_link_libraries(ABM xfitter)
install(TARGETS ABM DESTINATION ${DEST_MODULE})
