# Building

```
$ apt-get install libpqxx-dev libproj-dev libpq-dev postgresql-server-dev-10
$ make
$ make install
```

# Installation
```
# create extension irgrid;
CREATE EXTENSION

# \dx
                                     List of installed extensions
  Name   | Version |   Schema   |                             Description
---------+---------+------------+---------------------------------------------------------------------
 irgrid  | 0.0.1   | public     | irregular grid
 plpgsql | 1.0     | pg_catalog | PL/pgSQL procedural language
 postgis | 2.4.3   | public     | PostGIS geometry, geography, and raster spatial types and functions
(3 rows)
```

# Uninstall
```
# drop extension irgrid;
DROP EXTENSION

$ make uninstall 
```
