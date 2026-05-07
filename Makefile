CC ?= cc
CFLAGS ?= -O3 -std=c11 -Wall -Wextra -Wpedantic
CPPFLAGS ?= -Iinclude
LDFLAGS ?=

SRC := src/main.c src/io.c src/sw_scalar.c src/sw_simd.c
OBJ := $(patsubst src/%.c,build/%.o,$(SRC))
BIN := build/water

.PHONY: all clean test

all: $(BIN)

build:
	mkdir -p build

build/%.o: src/%.c | build
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BIN): $(OBJ)
	$(CC) $(OBJ) $(LDFLAGS) -o $@

clean:
	rm -rf build

test: $(BIN)
	./tests/test_example.sh
