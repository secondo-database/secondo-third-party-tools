CREATE OR REPLACE FUNCTION
createIrgrid(geometry,int,int, text, text) RETURNS int AS 'MODULE_PATHNAME','createIrgrid'
LANGUAGE C STRICT;

CREATE OR REPLACE FUNCTION
getCellnumbers(geometry, text) RETURNS int[] AS 'MODULE_PATHNAME','getCellnumbers'
LANGUAGE C STRICT;
