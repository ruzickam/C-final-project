# Compile

final: final.c
	gcc -o final final.c -std=c99 -Wall -pedantic -lg2 -lX11 -lm

