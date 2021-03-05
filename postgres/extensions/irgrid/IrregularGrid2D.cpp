#include "IrregularGrid2D.h"

#include <iostream>
#include <map>
#include <iterator>
#include <set>
#include <math.h>

#include <pqxx/pqxx>
using namespace std;
using namespace pqxx;

Cell::Cell() {
  cellId = -1;
  valFrom = -1;
  valTo = -1;
  nbrPoints = -1;
}

void
Cell::setCellId(int cell_id) {
  this->cellId = cell_id;
}

int
Cell::getCellId() {
  return cellId;
}

void
Cell::setValFrom(double val_from) {
  valFrom = val_from;
}

double
Cell::getValFrom() {
  return valFrom;
}

void
Cell::setValTo(double val_to) {
  valTo = val_to;
}

double
Cell::getValTo() {
  return valTo;
}

int
Cell::getNbrOfPoints() {
  return nbrPoints;
}

void
Cell::setNbrOfPoints(int nbr_points) {
  nbrPoints = nbr_points;
}

Cell::~Cell() { }

VCell::VCell() {
}

std::vector<HCell>&
VCell::getRow() {
  return row;
}

VCell::~VCell() { }

HCell::HCell() {
  upper = nullptr;
  upper_cellId = -1;
}

void
HCell::setUpper(HCell* upper_) {
  upper = upper_;
}

void
HCell::setUpperCellId(int ucell_id) {
  this->upper_cellId = ucell_id;
}

HCell*
HCell::getUpper() {
  return upper;
}

HCell::~HCell() { }

template <unsigned dim>
Rectangle<dim>::Rectangle(double min[], double max[]) {
  for( unsigned i = 0; i < dim; i++ ) {
      this->min[i] = min[i];
      this->max[i] = max[i];
  }
}

template <unsigned dim>
bool Rectangle<dim>::Intersects( const Rectangle<dim>& r) const {
  for( unsigned i = 0; i < dim; i++ )
    if( max[i] < r.min[i] || r.max[i] < min[i] )
      return false;

  return true;
}

template <unsigned dim>
Rectangle<dim>
Rectangle<dim>::Intersection( const Rectangle<dim>& r) const {
  double auxmin[dim], auxmax[dim];
  for( unsigned i = 0; i < dim; i++ )
  {
    auxmin[i] = MAX( min[i], r.min[i] );
    auxmax[i] = MIN( max[i], r.max[i] );
  }
  return Rectangle<dim>( auxmin, auxmax );
}

double getAlmostEqualFactor(){
  return 0.00000001;
}

bool AlmostEqualAbsolute( const double &d1, const double &d2,
  const double &epsilon)
{
  double diff = fabs(d1-d2);
  return ( diff < epsilon );
}

template <unsigned dim>
bool Rectangle<dim>::Contains( const Rectangle<dim>& r) const
{
  for( unsigned i = 0; i < dim; i++ ){
    if( min[i] > r.min[i] || max[i] < r.max[i] ){
        if(min[i] > r.min[i] && !::AlmostEqualAbsolute(min[i], r.min[i],
          getAlmostEqualFactor())) return false;
        if(max[i] < r.max[i] && !::AlmostEqualAbsolute(max[i], r.max[i],
          getAlmostEqualFactor())) return false;
    }
  }
  return true;
}

IrregularGrid2D::IrregularGrid2D(double left, double right,
  double bottom, double top,
  int row_count, int cell_count,
  char* conn_str, char* rel_str) {

  double min[2], max[2];
  min[0] = left;
  min[1] = bottom;
  max[0] = right;
  max[1] = top;

  boundingBox = new Rectangle<2>(min, max);

  rowCount = row_count;
  cellCount = cell_count;

  connection_str = conn_str;
  relation_str = rel_str;

  created = false;
}

IrregularGrid2D::~IrregularGrid2D() { }

bool
IrregularGrid2D::generateIrgrid2D() {

  // sort input by y-coordinates
  processInput();

  // create irregular grid 2d by point density
  buildGrid();

  storeIrgrid2d();

  return this->created;
}

void*
initIrgrid2D(double l, double r, double b,  double t,
  int row_count, int cell_count, char* conn_str, char* rel_str)
{
   IrregularGrid2D *irgrid( new IrregularGrid2D(l, r, b, t,
    row_count, cell_count, conn_str, rel_str) );

   return( reinterpret_cast< void* >( irgrid ) );
}

void
destroyIrregularGrid2D( void* irgrid2d )
{
   delete( reinterpret_cast< IrregularGrid2D* >( irgrid2d ) );
}

bool
createIrgrid2D(void *irgrid2d)
{
   bool res = reinterpret_cast< IrregularGrid2D* >( irgrid2d ) ->
    generateIrgrid2D();

   return res;
}

ARInfo
getCellnos(void *irgrid2d, double r_l, double r_r, double r_b, double r_t,
  char *conn_str)
{
     ARInfo res = reinterpret_cast< IrregularGrid2D* >( irgrid2d ) ->
      cellnos(r_l, r_r, r_b, r_t, conn_str);

   return res;
}

Rectangle<2>*
IrregularGrid2D::getBoundingBox() {
  return boundingBox;
}

int
IrregularGrid2D::getRowCount() {
  return this->rowCount;
}

int
IrregularGrid2D::getCellCount() {
  return this->cellCount;
}

void
IrregularGrid2D::storeIrgrid2d() {
  try {
    //std::vector<CellInfo*> cell_info_vect {};

    std::string connection_str_in;
    if (connection_str == nullptr) {
      connection_str_in = "dbname = testdb";
      outTxtPsql("storeIrgrid2d: connection_str is NULL. "\
        "Set to dbname = testdb");
    } else {
      std::string connection_str_in_tmp( connection_str );
      connection_str_in = connection_str_in_tmp;
    }

    std::string sql_in_str = "";
    bool connected = false;
    connection C(connection_str_in);

    std::string sql_del = "DELETE FROM irgrid2d";

    std::string sql_inp1 = "INSERT INTO irgrid2d(ID, bb_l, bb_r, bb_b, bb_t, "\
      "no_row, no_cell, grid) VALUES ('g1', ";
    std::string sql_inp2 = "]);";
    std::string sql_inp_bb = "";
    std::string sql_inp_cv = "";

    if (C.is_open()) {
       connected = true;
    } else {
       outTxtPsql("CONNECTION ERROR to irgrid2d");
    }

    if (connected) {
      // clear irgrid2d
      work W(C);
      W.exec( sql_del );
      //W.commit();

      // bouding box coordinates
        double bb_l = boundingBox->getMinX();
        double bb_r = boundingBox->getMaxX();
        double bb_b = boundingBox->getMinY();
        double bb_t = boundingBox->getMaxY();

        sql_inp_bb = std::to_string(bb_l) + ", " + std::to_string(bb_r) +
         ", " + std::to_string(bb_b) + ", " +  std::to_string(bb_t) +
         ", " + std::to_string(rowCount) +
         ", " +  std::to_string(cellCount) + ", " +
         "ARRAY[";

      std::vector<VCell>* col = &getColumnVector();
      if (col != nullptr) {
        for(size_t colIdx = 0; colIdx < col->size(); colIdx++) {
          VCell* vcell = &col->at(colIdx);
          double c_b = vcell->getValFrom();
          double c_t = vcell->getValTo();

          std::vector<HCell>* row_vect = &vcell->getRow();

          // create cell  grid value statements
          for(size_t cellIdx = 0; cellIdx < row_vect->size(); cellIdx++) {
            HCell* cell = &row_vect->at(cellIdx);
            int cid = cell->getCellId();
            int upper_cid = -1;
            if (cell->getUpper() != nullptr) {
              upper_cid = cell->getUpper()->getCellId();
            }
            int np = cell->getNbrOfPoints();
            double c_l = cell->getValFrom();
            double c_r = cell->getValTo();

            /*CellInfo* ci = new CellInfo(cid, upper_cid,
              np, c_l, c_r, c_b, c_t);*/

            std::string sql_in_str_tmp = "row(" + std::to_string(cid) + ", "
              + std::to_string(upper_cid) + ", " + std::to_string(np) + ", "
              + std::to_string(c_l) + ", " + std::to_string(c_r) + ", "
              + std::to_string(c_b) + ", "
              + std::to_string(c_t) + ")::grid_cell";
              if (!(cellIdx == (row_vect->size())-1 &&
                 colIdx == (col->size())-1)) {
                sql_in_str_tmp += ", ";
              }
              sql_inp_cv += sql_in_str_tmp;
              //cell_info_vect.push_back(ci);
          }
        }
      }
      sql_in_str = sql_inp1 + sql_inp_bb + sql_inp_cv + sql_inp2;
      //outTxtPsql("SQL: " + sql_in_str);
      const char * sql_in = &sql_in_str[0];
      W.exec( sql_in );
      W.commit();

      C.disconnect ();
      this->created=true;
      //return cell_info_vect;
    }
  } catch (const std::exception &e) {
      outTxtPsql((string)e.what());
    }
}

bool
pointComparisonX(RPoint p1, RPoint p2) {
  return p1.x < p2.x;
}

void
IrregularGrid2D::buildGrid() {
  // create grid structure
  double bb_top = boundingBox->getMaxY();
  double bb_bot = boundingBox->getMinY();
  double colOffset = (bb_top - bb_bot) / rowCount;

  double bb_right = boundingBox->getMaxX();
  double bb_left = boundingBox->getMinX();
  double hOffset = (bb_right - bb_left) / cellCount;

  double col_boundary_val = bb_bot;
  int hcell_id = 1;

  for(int c = 0; c < rowCount; c++) {
    VCell vcell = VCell();
    vcell.setValFrom(col_boundary_val);

    col_boundary_val += colOffset;
    vcell.setValTo(col_boundary_val);

    double cell_boundary_val = bb_left;
    for(int r = 0; r < cellCount; r++) {
      HCell hcell = HCell();

      hcell.setValFrom(cell_boundary_val);

      cell_boundary_val += hOffset;
      hcell.setValTo(cell_boundary_val);

      hcell.cellId = hcell_id;
      hcell_id++;

      // will be determined later
      hcell.upper = nullptr;

      vcell.getRow().push_back(hcell);

    }
    columnVector.push_back(vcell);
  }

  // adjust boundaries by point distribution
  int nbrOfPoints = points.size();
  int pointsPerRow = nbrOfPoints / rowCount;
  std::vector<RPoint> tmp_row_points {};

  int pointIdx = 0;
  int point_counter = 0;
  for(int colIdx = 0; colIdx < rowCount; colIdx++) {
    for(size_t dp_idx = pointIdx; dp_idx < points.size(); dp_idx++) {
      RPoint rp = points[dp_idx];
      tmp_row_points.push_back(rp);
      point_counter ++;

      if ((point_counter == pointsPerRow) || (dp_idx == points.size()-1)) {
        // adjust y-boundaries

        if (colIdx > 0) {
          getColumnVector()[colIdx].setValFrom(
            getColumnVector()[colIdx-1].getValTo());
        }

        double new_val_to_y;
        if (dp_idx == points.size()-1) {
          new_val_to_y = bb_top;
        } else {
          new_val_to_y = rp.y;
        }
        getColumnVector()[colIdx].setValTo(new_val_to_y);
        getColumnVector()[colIdx].setNbrOfPoints(point_counter);

        point_counter = 0;

        // adjust x-boundaries
        std::sort(tmp_row_points.begin(), tmp_row_points.end(),
          pointComparisonX);

        std::vector<HCell> c_row =  columnVector[colIdx].getRow();
        int pointsPerCell = tmp_row_points.size() / cellCount;

        int tmpPointIdx = 0;
        for(int h = 0; h < cellCount; h++) {
          std::vector<HCell>* hcel_vec_ptr = &columnVector[colIdx].getRow();

          for(size_t tmp_rp_idx = tmpPointIdx;
              tmp_rp_idx < tmp_row_points.size(); tmp_rp_idx++) {

            RPoint rp_t = tmp_row_points[tmp_rp_idx];
            point_counter ++;

            if((point_counter == pointsPerCell)
              || (tmp_rp_idx == tmp_row_points.size()-1)
              /* || (h == cellCount-1 ) */ ) {

              if (h > 0) {
                hcel_vec_ptr->at(h).setValFrom(
                  hcel_vec_ptr->at(h-1).getValTo());
              }

              double new_val_to_x;
              if ((tmp_rp_idx == tmp_row_points.size()-1)
                  || (h == cellCount-1)) {
                new_val_to_x = bb_right;
              } else {
                 new_val_to_x = rp_t.x;
              }
              hcel_vec_ptr->at(h).setValTo(new_val_to_x);

              tmpPointIdx = tmp_rp_idx + 1;
              hcel_vec_ptr->at(h).setNbrOfPoints(point_counter);

              // Smoothing in case of odd number of points per cell
              if (h+1 < cellCount) {
                pointsPerCell = (tmp_row_points.size()-tmpPointIdx)
                  / (cellCount-1-h);
              }

              point_counter = 0;

              // next cell
              break;
            }
          }
        }

        tmp_row_points.clear();
        pointIdx = dp_idx + 1;

        // Smoothing in case of odd number of points per row
        if (colIdx+1 < rowCount) {
          pointsPerRow = (points.size()-pointIdx)
            / (rowCount-1-colIdx);
        }

        point_counter = 0;

        // one row up
        break;
      }
    }
  }

  // clear aux. vectors
  if (points.size() > 0) {
    points.clear();
  }
  if (tmp_row_points.size() > 0) {
    tmp_row_points.clear();
  }

  // update cell pointer
  if (rowCount > 1 && cellCount > 0) {
    for(int c = 0; c < rowCount-1; c++) {
      std::vector<HCell>* row_lower = &getColumnVector().at(c).getRow();
      std::vector<HCell>* row_upper = &getColumnVector().at(c+1).getRow();

      int pointToCellIdx = 0;
      for(int h = 0; h < cellCount; h++) {
        HCell* lower_cell_ptr = &row_lower->at(h);
        HCell* upper_cell_ptr = &row_upper->at(pointToCellIdx);

        if (lower_cell_ptr->getValFrom() <= upper_cell_ptr->getValTo()) {
          lower_cell_ptr->setUpper(upper_cell_ptr);
        } else {
          HCell* next_upper_cell_ptr;
          do {
            pointToCellIdx ++;
            next_upper_cell_ptr = &row_upper->at(pointToCellIdx);
            lower_cell_ptr->setUpper(next_upper_cell_ptr);
          } while (lower_cell_ptr->getValFrom()
              >= next_upper_cell_ptr->getValTo());
        }
      }
    }
  }
}

bool
pointComparisonY(RPoint p1, RPoint p2) {
  return p1.y < p2.y;
}

RPoint
getRectangleCentre(double c_x, double c_y) {
  RPoint r_c { c_x, c_y };

  return r_c;
}

// check if a point is inside the irgrid2d bounding box
bool
insideBoundingBox(Rectangle<2>* bbox, double p_x, double p_y) {
  double le = bbox->getMinX();
  double ri = bbox->getMaxX();
  double bo = bbox->getMinY();
  double to = bbox->getMaxY();

  if (p_x >= le && p_x <= ri
      && p_y >= bo && p_y <=to) {
    return true;
  }

  return false;
}

void
IrregularGrid2D::processInput() {
  /*connection C("dbname = testdb user = user \
    hostaddr = 127.0.0.1 port = 5432");*/

  // psql connection
  //connection C("dbname = testdb");
  //std::string relation_str_in = "pofw";

  std::string connection_str_in;
  if (connection_str == nullptr) {
    connection_str_in = "dbname = testdb";
    outTxtPsql("processInput: connection_str is NULL. Set to dbname = testdb");
  } else {
    std::string connection_str_in_tmp( connection_str );
    connection_str_in = connection_str_in_tmp;
  }

  std::string relation_str_in;
  if (relation_str == nullptr) {
    relation_str_in = "pofw";
    outTxtPsql("processInput: relation_str is NULL. Set to pofw");
  } else {
    std::string relation_str_in_tmp( relation_str );
    relation_str_in = relation_str_in_tmp;
  }
  connection C(connection_str_in);

  /* std::string sql_p1 = "select ST_X(geog::geometry), ST_Y(geog::geometry) "\
    "from (select ST_Centroid(geog, true) AS geog from "; */

  std::string sql_p1 = "select ST_X(ST_Centroid(geog, true)::geometry), "\
    "ST_Y(ST_Centroid(geog, true)::geometry) from "\
    " (select ST_Envelope(geog::geometry) AS geog from ";

  std::string sql_p2 = ") AS a;";
  std::string sql_con = sql_p1 + relation_str_in + sql_p2;
  const char * sql = &sql_con[0];

  try {
    if (C.is_open()) {
      //outTxtPsql((string)"Connected to :" + C.dbname());

      nontransaction N(C);
      result R( N.exec( sql ));

      std::string in_size = std::to_string(R.size());
      /*outTxtPsql((string)"nbr of input polygons: " + in_size +
        ", relation: " + C.dbname());*/

      double st_x, st_y;
      std::size_t offset = 0;
      for (result::const_iterator c = R.begin(); c != R.end(); ++c) {
        // pont x-coord
        if (!c[0].is_null()) {
          std::string st_x_str = c[0].as<string>();
          st_x = stod(st_x_str, &offset);

          //outTxtPsql("st_x: " + std::to_string(st_x));
        }
        // point y-coord
        if (!c[1].is_null()) {
          std::string st_y_str = c[1].as<string>();
          st_y = stod(st_y_str, &offset);

          //outTxtPsql("st_y: " + std::to_string(st_y));
        }

        if (insideBoundingBox(boundingBox, st_x, st_y)) {
          points.push_back(getRectangleCentre(st_x, st_y));
        }
      }
      C.disconnect ();

      // sort point vector by y-coordinates
      std::sort(points.begin(), points.end(), pointComparisonY);
    } else {
      outTxtPsql("CONNECTION ERROR");
    }
  } catch (const std::exception &e) {
       outTxtPsql((string)e.what());
  }
}

void
IrregularGrid2D::setColumnVector(std::vector<VCell> column_vect) {
  this->columnVector = column_vect;
}

std::vector<VCell>&
IrregularGrid2D::getColumnVector() {
  return this->columnVector;
}

static IrregularGrid2D*
loadGrid(char *grid_id, char *connection_str) {
  //std::string grid_id_in( grid_id );
  IrregularGrid2D *irgrid = nullptr;
  double bb_l, bb_r, bb_b, bb_t;
  double nbr_r, nbr_c;

  int c_id, uc_id;
  double c_l, c_r, c_b, c_t;

  VCell vc;
  VCell* vc_ptr;
  std::vector<VCell> column_vec {};
  // temporary support structures
  std::map<int, int> cellRef;
  std::map<int, HCell*> cellIds;

  std::string connection_str_in;
  if (connection_str == nullptr) {
    connection_str_in = "dbname = testdb";
    IrregularGrid2D::outTxtPsql("loadGrid: connection_str is NULL. "\
      "Set to dbname = testdb");
  } else {
    std::string connection_str_in_tmp( connection_str );
    connection_str_in = connection_str_in_tmp;
  }
  std::string relation_str_in = "irgrid2d";
  connection C(connection_str_in);

  std::string sql_bb = "SELECT bb_l, bb_r, bb_b, bb_t, no_row, no_cell FROM "
    + relation_str_in
    + " WHERE ID='g1';";

  std::string sql_cell = "SELECT cell_id, upper_cell_id, c_l, c_r, c_b, c_t, "\
    "points_per_cell FROM "
    + relation_str_in
    + ", unnest(grid)"
    + " WHERE ID='g1';";

  const char * sql = &sql_bb[0];
  const char * sql_c = &sql_cell[0];

  try {
    if (C.is_open()) {
      nontransaction N(C);
      result R( N.exec( sql ));

      std::size_t offset = 0;
      for (result::const_iterator c = R.begin(); c != R.end(); ++c) {
        // bb_l
        if (!c[0].is_null()) {
          std::string bb_l_str = c[0].as<string>();
          bb_l = stod(bb_l_str, &offset);
        }
        // bb_r
        if (!c[1].is_null()) {
          std::string bb_r_str = c[1].as<string>();
          bb_r = stod(bb_r_str, &offset);
        }
        // bb_b
        if (!c[2].is_null()) {
          std::string bb_b_str = c[2].as<string>();
          bb_b = stod(bb_b_str, &offset);
        }
        //bb_t
        if (!c[3].is_null()) {
          std::string bb_t_str = c[3].as<string>();
          bb_t = stod(bb_t_str, &offset);
        }
        // nbr_r
        if (!c[4].is_null()) {
          nbr_r = c[4].as<int>();
        }
        //nbr_c
        if (!c[5].is_null()) {
          nbr_c = c[5].as<int>();
        }
      }

      // columnVector
      int row_cnt = 0;
      result RC( N.exec( sql_c ));
      for (result::const_iterator c = RC.begin(); c != RC.end(); ++c) {
        row_cnt ++;

        // cell_id
        if (!c[0].is_null()) {
          c_id = c[0].as<int>();
        }
        // upper_cell_id
        if (!c[1].is_null()) {
          uc_id = c[1].as<int>();
        }
        //points_per_cell
        // c_l
        if (!c[2].is_null()) {
          std::string c_l_str = c[2].as<string>();
          c_l = stod(c_l_str, &offset);
        }
        //c_r
        if (!c[3].is_null()) {
          std::string c_r_str = c[3].as<string>();
          c_r = stod(c_r_str, &offset);
        }
        // c_b
        if (!c[4].is_null()) {
          std::string c_b_str = c[4].as<string>();
          c_b = stod(c_b_str, &offset);
        }
        //c_t
        if (!c[5].is_null()) {
          std::string c_t_str = c[5].as<string>();
          c_t = stod(c_t_str, &offset);
        }

        if (row_cnt == 1) {
          vc = VCell();
          vc.setValFrom(c_b);
          vc.setValTo(c_t);
          column_vec.push_back(vc);
          vc_ptr = &(column_vec.back());
        }

        HCell* hc = new HCell();
        hc->setValFrom(c_l);
        hc->setValTo(c_r);
        hc->setCellId(c_id);
        hc->setUpperCellId(uc_id);
        // will be determined later
        hc->setUpper(nullptr);
        vc_ptr->getRow().push_back(*hc);

        // update support structures
        if (uc_id != -1) {
          cellRef.insert(std::make_pair(c_id, uc_id));
        }
        cellIds.insert(std::make_pair(c_id, hc));

        // next vcell
        if (row_cnt == nbr_c) {
          row_cnt = 0;
        }
      }

      if (nbr_r > 0 && nbr_c > 0) {
        // update pointer
        if (nbr_r > 0 && nbr_c > 0) {
          for(int colIdx = 0; colIdx < nbr_r-1; colIdx++) {
            VCell* vcell = &column_vec.at(colIdx);
            std::vector<HCell>* row_vect = &vcell->getRow();
            for(int cIdx = 0; cIdx < nbr_c; cIdx++) {
              HCell* hcell = &(*row_vect).at(cIdx);
              if(cellRef.find(hcell->getCellId()) != cellRef.end()) {
                int cell_ref = cellRef.at(hcell->getCellId());
                if(cellIds.find(cell_ref) != cellIds.end()) {
                  hcell->setUpper(cellIds.at(cell_ref));
                }
              }
            }
          }
        }
      }

      irgrid = ( new IrregularGrid2D(bb_l, bb_r, bb_b, bb_t,
        nbr_r, nbr_c, nullptr, nullptr) );

      irgrid->setColumnVector(column_vec);

      C.disconnect ();
    }
  }  catch (const std::exception &e) {
    IrregularGrid2D::outTxtPsql((string)e.what());
  }

  return irgrid;
}

template <class C>
bool
InCell(C cell, double val) {
  return (val >= cell.getValFrom()
    && val < cell.getValTo());
}

template <class C>
bool
GtCell(C cell, double val) {
  return (cell.getValFrom() >= val);
}

template <class C>
int
CellBS(const std::vector<C>* c_vec, int start, int end, const double val) {
  if (start > end) {
    return -1;
  }

  const int mid = start + ((end - start) / 2);

  if (InCell(c_vec->at(mid), val)) {
    return mid;
  } else if (GtCell(c_vec->at(mid), val)) {
    return CellBS(c_vec, start, mid-1, val);
  }

  return CellBS(c_vec, mid+1, end, val);
}

ARInfo
IrregularGrid2D::cellnos(double r_l, double r_r, double r_b, double r_t,
  char *conn_str) {

  double min[2], max[2];
  min[0] = r_l;
  min[1] = r_b;
  max[0] = r_r;
  max[1] = r_t;

  Rectangle<2> *search_window_ptr = new Rectangle<2>(min, max);

  IrregularGrid2D *input_irgrid2d_ptr = loadGrid(nullptr, conn_str);

  std::set<int> cell_ids;
  if (input_irgrid2d_ptr != nullptr && search_window_ptr != nullptr) {
    bool skipSearch = false;

    Rectangle<2> *b_box = input_irgrid2d_ptr->getBoundingBox();
    if (!search_window_ptr->Intersects(*b_box)) {
      cell_ids.insert(0);
      skipSearch = true;
    }

    // 'truncate' search window in case of partial cutting
    if (!(*b_box).Contains(*search_window_ptr)) {
      search_window_ptr = new Rectangle<2>(
        search_window_ptr->Intersection(*b_box));

      cell_ids.insert(0);
      outTxtPsql("INSERT (PC): 0");
    }

    if (!skipSearch) {
      std::vector<VCell>* col = &input_irgrid2d_ptr->getColumnVector();

      double le = search_window_ptr->getMinX();
      double ri = search_window_ptr->getMaxX();
      double bo = search_window_ptr->getMinY();
      double to = search_window_ptr->getMaxY();

      int pos_bo = CellBS(col, 0, col->size(), bo);
      if (pos_bo != -1) {
        VCell vCell = col->at(pos_bo);
        std::vector<HCell>* row = &vCell.getRow();

        int pos_le = CellBS(row, 0, row->size(), le);
        if (pos_le != -1) {

          // collect ids
          unsigned int cellIdx = pos_le;
          while (cellIdx < row->size()) {
            HCell i = row->at(cellIdx);
            cell_ids.insert(i.getCellId());
            //outTxtPsql("INSERT: " + to_string(i.getCellId()));

            if ((ri >= i.getValFrom() && ri < i.getValTo())
              || (cellIdx == row->size()-1  && ri >= i.getValFrom()))  {
                HCell fi = row->at(pos_le);
                if (to >= vCell.getValTo() && fi.getUpper() != nullptr) {
                  vCell = col->at(++pos_bo);
                  row = &vCell.getRow();

                  HCell * u = fi.getUpper();
                  int nbr_cpr = input_irgrid2d_ptr->getCellCount();
                  int cid_pos = (u->getCellId()) % nbr_cpr;
                  pos_le = cid_pos == 0 ? nbr_cpr-1 : cid_pos-1;

                  cellIdx = pos_le-1;
              } else if (to < vCell.getValTo() || fi.getUpper() == nullptr) {
                break;
              }
            }
            cellIdx++;
          }
        }
      }
    }
  }

  int* arr = nullptr;
  if (cell_ids.size() > 0) {
    arr = new int[cell_ids.size()];
    set<int>::iterator it;
    int cnt = 0;
    for (it = cell_ids.begin(); it != cell_ids.end(); ++it) {
      arr[cnt] = *it;
      cnt++;
    }
  }
  ARInfo arInfo {arr, cell_ids.size()};

  return arInfo;
}

int main() {
    return 0;
}
