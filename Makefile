# The compiler:
CC= gcc

# Flags:
# see https://stackoverflow.com/questions/3375697/useful-gcc-flags-for-c
CFLAGS= @compilerFlags
# C files:
SOURCES=  main.c random.c stat.c adj.c sync.c top.c lists.c nuc.c

# Name of output:
TARGET= net

$(TARGET):
	$(CC) $(SOURCES) $(CFLAGS) -o $(TARGET)

clean:
	rm -f $(TARGET)

valgrind:
	valgrind --leak-check=yes ./$(TARGET)
