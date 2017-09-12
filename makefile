# makefile for FELAC, implemented in ANSI C
CC = gcc
FLAGS = -Wall -c
NAME = CFLAG
OBJS = a1.c \
gidmsh.c gidres.c gidpre.c ea1b.c ua1a.c ea1a.c startb.c starta.c bet6.c \
beq9g3.c bet3g2.c beq4g2.c aet6.c aeq9g3.c aet3g2.c aeq4g2.c solv.c epgsub.c 
 
#SRC = $(wildcard *.c)
obj = $(patsubst %.c,%.o,$(OBJS))
 
 
a1: $(obj)
	gcc $(obj) -o a1  -L"%felacpath%/lib/" -lfelac -L"%felacpath%/solv/SuperLU/" -lsolv -L"%felacpath%/solv/AztecOO/" -laztec -llapack -lm
$(obj):%.o:%.c
	gcc -std=c99 -c $< -o $@
 
.PHONY: clean
clean:
	del a1 a1.exe *.o *.obj
