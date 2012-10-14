CC=nasm
LD=ld
SOURCES=foo.asm
OBJ=foo.o
AFLAGS=-f elf32 -g -ggdb -F dwarf
LDFLAGS=-shared -melf_i386
OUTPUT=../bin/libasm.so
all:
	$(CC) $(SOURCES) $(AFLAGS) -o $(OBJ)
	$(LD) $(LDFLAGS) $(OBJ) -o $(OUTPUT)
