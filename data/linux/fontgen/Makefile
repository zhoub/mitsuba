CPPFLAGS = -I ../../../include -I /usr/include/freetype2 -DSINGLE_PRECISION -DSPECTRUM_SAMPLES=3
LDFLAGS = -lfreetype -L../../../dist -lmitsuba-core -lmitsuba-hw -lboost_filesystem -lboost_system -lz
CC = g++

all: bin2c fontgen

bin2c: bin2c.c

fontgen: fontgen.o
	$(CC) $^ $(LDFLAGS) -o $@

clean:
	rm -f fontgen
