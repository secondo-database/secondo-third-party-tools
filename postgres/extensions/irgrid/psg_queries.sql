-- login / admin
-- sudo -u postgres psql testdb
-- sudo systemctl start postgresql

select * from pofw where osm_id = '54848858';

-- text representation
select ST_AsEWKT(geog) from pofw where osm_id = '54848858';
select ST_AsText(geog) from pofw where osm_id = '54848858';

-- point (center) represented as text
-- ST_GeomFromEWKT: geography --> geometry
-- ST_Extent: geometry --> BOX
-- ST_Centroid: BOX --> POINT
select ST_AsText(
  ST_Centroid(
    ST_Extent(
      ST_GeomFromEWKT('SRID=4326;MULTIPOLYGON(((7.4695507 51.5008368,7.469708 51.5008497,7.4697333 51.5007225,7.4698341 51.50073,7.4698437 51.5006675,7.4698438 51.5006645,7.4697483 51.5006582,7.4697491 51.5006541,7.4697011 51.5006508,7.4695917 51.5006431,7.4695507 51.5008368)))'))));


-- point (center)
select ST_Centroid(geog, true) from pofw where osm_id = '54848858';
select ST_AsText(geog) from (select ST_Centroid(geog, true) AS geog from pofw where osm_id = '54848858') AS a;

-- point coordinates
select ST_X(geog::geometry), ST_Y(geog::geometry) from (select ST_Centroid(geog, true) AS geog from pofw where osm_id = '54848858') AS a;
-- select ST_X(ST_Centroid(geog, true)::geometry), ST_Y(ST_Centroid(geog, true)::geometry) from (select ST_Envelope(geog::geometry) AS geog from pofw where osm_id = '54848858') as a;

-- table info
SELECT table_name, column_name, data_type
FROM information_schema.columns
WHERE table_name = 'pofw';

-- create irgrid2d
-- select createIrgrid(ST_MakeEnvelope(7.86403009, 8.21026003, 51.35471354 ,51.4913301, 4326),2,3, 'dbname = testdb', 'pofw');
-- select getCellnumbers(ST_MakeEnvelope(8.09280, 8.09281, 51.354718, 51.354719, 4326), 'dbname = testdb');
-- select getCellnumbers(ST_MakeEnvelope(8, 8.08, 51.354718, 51.406, 4326), 'dbname = testdb');
-- select getCellnumbers(ST_MakeEnvelope(7.994002, 8.210258, 51.407456, 52.4911, 4326), 'dbname = testdb');
-- select getCellnumbers(ST_MakeEnvelope(1, 2, 3, 4, 4326), 'dbname = testdb');

 -- types (create, insert, select, delete)
CREATE TYPE grid_cell AS (
    cell_id	integer,
    upper_cell_id	integer,
    points_per_cell	integer,
    c_l	double precision,
    c_r	double precision,
    c_b	double precision,
    c_t	double precision
);

CREATE TABLE irgrid2d (
	ID varchar(16) PRIMARY KEY NOT NULL,
    grid grid_cell[],
    bb_l	double precision,
    bb_r	double precision,
    bb_b	double precision,
    bb_t	double precision,
    no_row	integer,
    no_cell integer
);

-- INSERT INTO irgrid2d(ID, bb_l, bb_r, bb_b, bb_t, no_row, no_cell, grid)
-- VALUES (
--	'g1', 7.86403009, 8.21026003, 51.35471354 ,51.4913301, -1, -1,
--	ARRAY[row(1, -1, 10, 1.0, 2.09098674, 7, 19.51560973)::grid_cell,
--		 row(2, -1, 10, 1.0, 2.09098674, 7, 19.51560973)::grid_cell,
--		 row(3, -1, 10, 1.0, 2.09098674, 7, 22.51560973)::grid_cell]
--);

-- SELECT bb_l, bb_r, bb_b, bb_t FROM irgrid2d WHERE ID='g1';
-- SELECT cell_id, upper_cell_id, c_l, c_r, c_b, c_t, points_per_cell FROM irgrid2d, unnest(grid) WHERE ID='g1';
-- SELECT cell_id, upper_cell_id, c_l, c_r, c_b, c_t, points_per_cell FROM irgrid2d, unnest(grid);
-- SELECT unnest(grid) FROM irgrid2d;

-- DELETE FROM irgrid2d
