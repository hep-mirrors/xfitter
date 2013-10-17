AUTOMAKE_OPTIONS = foreign
noinst_LIBRARIES = libgenetic.a

if ENABLE_GENETIC
 pkginclude_HEADERS = genetic_search.h mixmax_wrapper.h
 libgenetic_a_SOURCES = mixmax_wrapper.C genetic_search.C
 LDADD = mixmax_r004/libmixmax.a
else 
 pkginclude_HEADERS = genetic_dummy.h
 libgenetic_a_SOURCES = genetic_dummy.C
endif



