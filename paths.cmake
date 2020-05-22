
#set(CMAKE_BUILD_TYPE DEBUG)

# - if you want to enable debugging symbols in the code please add the following to your paths.cmake:
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
