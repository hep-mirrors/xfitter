add_library(ABM SHARED
  alphasvfn.f
  apsolve.f
  asy-hpcc.f
  asympcoef.f
  CMakeLists.txt
  cq.f
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
  ABM.cc
)
target_link_libraries(ABM xfitter)
install(TARGETS ABM DESTINATION ${DEST_MODULE})
