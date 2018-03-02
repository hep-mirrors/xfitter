AUTOMAKE_OPTIONS = foreign

noinst_LTLIBRARIES = libgenetic.la
dist_noinst_HEADERS = genetic_search.h mixmax_wrapper.h genetic_dummy.h

if ENABLE_GENETIC
 libgenetic_la_SOURCES = mixmax_wrapper.C genetic_search.C
 LDADD = mixmax_r004/libmixmax.a
else 
 libgenetic_la_SOURCES = genetic_dummy.C
endif



