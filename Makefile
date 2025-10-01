bin/main: src/*.c
	cc -o bin/main src/*.c -I raylib/src -L raylib/src -lraylib -lGL -lm -lpthread -ldl -lrt -lX11
