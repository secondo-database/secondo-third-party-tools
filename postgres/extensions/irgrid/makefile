EXTENSION = irgrid
DATA = irgrid--0.0.1.sql
OBJS = irgrid.o IrregularGrid2D.o
MODULES = irgrid IrregularGrid2D
MODULE_big = irgrid

# postgis
SHLIB_LINK += -L/home/markus/Downloads/postgis-3.0.2/liblwgeom/.libs -llwgeom
SHLIB_LINK += -lgeos_c
SHLIB_LINK += -L/home/markus/Downloads/proj-7.2.0/src/.libs -l:libproj.so

# libpqxx
SHLIB_LINK += -L/usr/lib64 -l:libpqxx-5.0.so

# postgis
PG_LDFLAGS += -L /home/markus/Downloads/postgis-3.0.2/libpgcommon
PG_CPPFLAGS += -I /home/markus/Downloads/postgis-3.0.2/libpgcommon
PG_LDFLAGS += -L /home/markus/Downloads/postgis-3.0.2/liblwgeom
PG_CPPFLAGS += -I /home/markus/Downloads/postgis-3.0.2/liblwgeom

# libpqxx
PG_LDFLAGS += -L /usr/include/pqxx
PG_CPPFLAGS += -I /usr/include/pqxx

PG_CONFIG = pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)

CXXFLAGS = -fPIC -Wall -Werror -g -O0 -lpqxx -lpq
CXX = g++

include $(PGXS)

### libs info	###
### 1) lwgeom_calculate_gbox: L/home/markus/Downloads/postgis-3.0.2/liblwgeom/.libs -llwgeom
### 2) GEOSGeom_destroy: -lgeos_c
### 3) geod_polygon_addpoint: L/home/markus/Downloads/proj-7.2.0/src/.libs -l:libproj.so
