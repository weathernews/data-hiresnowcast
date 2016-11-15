CFLAGS = -I/usr/local/include
LDLIBS = -L/usr/local/lib -lpng

all: grib-img grib-tile

grib-img: grib-img.c
grib-tile: grib-tile.c

clean:
	/bin/rm -f grib-img grib-tile
