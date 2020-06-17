# make file for suvendu's code, calculating 
#-------------------------------------
CC = g++
CFLAGS = -Wall -g 
#-------------------------------------
strain_map.X: p_vec_ten.o main.o 
	${CC} ${CFLAGS} p_vec_ten.o main.o -o strain_map.X

p_vec_ten.o: p_vec_ten.cc p_vec_ten.hh
	${CC} ${CFLAGS} -c p_vec_ten.cc

main.o: main.cpp p_vec_ten.hh p_vec_ten.cc
	${CC} ${CFLAGS} -c main.cpp







