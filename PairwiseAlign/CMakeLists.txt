cmake_minimum_required(VERSION 3.16)
project(PairwiseAlign C)

set(CMAKE_C_STANDARD 11)

include_directories(src)


add_executable(PairwiseAlign
        src/aliCost.c
        src/aliCost.h
        src/aliDnaAffine.c
        src/aliDnaAffineLongGaps.c
        src/aliDnaLinear.c
        src/aliGetSeq.c
        src/aliGetSeq.h
        src/aliOutSW.h
        src/aliProtAffine.c
        src/mem.c
        src/mem.h
        src/printStrings.c src/aliCalcAE.c src/aliCalcAE.h src/aliOutAE.c src/aliOutAE.h)
#add_executable(aliCalcSW src/aliCalcSW.c src/mem.c src/aliCost.c)
#target_link_libraries(PairwiseAlign m)
#target_link_libraries(aliCalcSW m)
#target_link_libraries(PairwiseAlign src)
#target_link_libraries(aliCalcSW src)
#add_executable(aliOutSW src/aliOutSW.c)

add_library(souf SHARED src/aliGetSeq.c src/aliCost.c src/aliCalcSW.c src/aliOutSW.c)
add_executable(aliDnaLinear src/aliDnaLinear.c src/mem.c src/aliCost.c)

add_library(souf3 SHARED src/aliGetSeq.c src/aliCost.c src/aliCalcAE.c src/aliOutAE.c)
add_executable(aliDnaAffineLongGaps src/aliDnaAffineLongGaps.c src/mem.c src/aliCost.c)

add_library(souf2 SHARED src/aliGetSeq.c src/aliCost.c src/aliCalcAE.c src/aliOutAE.c)
add_executable(aliDnaAffine src/aliDnaAffine.c src/mem.c src/aliCost.c)
add_executable(aliProtAffine src/aliProtAffine.c src/mem.c src/aliCost.c)

target_link_libraries(aliDnaLinear souf)
target_link_libraries(aliDnaLinear m)
target_link_libraries(aliDnaAffine souf2)
target_link_libraries(aliDnaAffine m)
target_link_libraries(aliDnaAffineLongGaps souf3)
target_link_libraries(aliDnaAffineLongGaps m)
target_link_libraries(aliProtAffine souf2)
target_link_libraries(aliProtAffine m)
