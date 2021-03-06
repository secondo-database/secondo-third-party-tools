# Needed debian packages libpqxx-dev libproj-dev libpq-dev postgresql-server-dev-10

EXTENSION = irgrid
DATA = irgrid--0.0.1.sql
OBJS = irgrid.o IrregularGrid2D.o
MODULE_big = irgrid

PG_CONFIG = pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)

#CXXFLAGS += -fPIC -Wall -Werror -g
#CXXFLAGS += -Wall -Werror -g
SHLIB_LINK += -lpqxx -lpq -llwgeom

include $(PGXS)

