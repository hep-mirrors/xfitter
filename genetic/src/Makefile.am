AUTOMAKE_OPTIONS = foreign
noinst_LTLIBRARIES = libgenetic.la

if ENABLE_GENETIC
 pkginclude_HEADERS = genetic_search.h mixmax_wrapper.h
 libgenetic_la_SOURCES = mixmax_wrapper.C genetic_search.C
 LDADD = mixmax_r004/libmixmax.a
else 
 pkginclude_HEADERS = genetic_dummy.h
 libgenetic_la_SOURCES = genetic_dummy.C
endif



