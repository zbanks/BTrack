TARGET=btest
SRCS = main.c ../src/BTrack.c ../src/OnsetDetectionFunction.c
OBJS = $(SRCS:%.c=%.o)

CC = gcc
LFLAGS = -lfftw3 -lsndfile -lsamplerate -lm
CFLAGS = -I. -I../src/ -ggdb3 -O0 -Wall -Wextra -Werror -Wno-unused-parameter

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LFLAGS)

.PHONY: clean
clean:
	-rm $(OBJS) $(TARGET)

.PHONY: all
all: $(TARGET)
