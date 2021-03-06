#ifndef IRREGULARGRID2D_H_
#define IRREGULARGRID2D_H_

#ifdef __cplusplus
extern "C" {
#endif
#include "postgres.h"
#ifdef __cplusplus
}
#endif

#ifndef MAX
#define MAX( a, b ) ((a) > (b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN( a, b ) ((a) < (b) ? (a) : (b))
#endif

//Represents a point of a rectangle corner
typedef struct {
  double x;
  double y;
} RPoint;

// struct for passing getCellnos result to postgres
typedef struct {
    int* res;
    long unsigned int res_size;
} ARInfo;


#ifdef __cplusplus
#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>

template <unsigned dim>
class Rectangle;

class Cell {
  public:
    Cell();
    ~Cell();

    int getCellId();
    void setValFrom(double valFrom);
    double getValFrom();
    void setValTo(double valTo);
    double getValTo();
    int getNbrOfPoints();
    void setNbrOfPoints(int nbrPoints);

    void setCellId(int cell_id);

    friend class IrregularGrid2D;

  private:
    int cellId;
    double valFrom;
    double valTo;
    int nbrPoints;
};

class HCell : public Cell {
  public:
    HCell();
    ~HCell();

    void setUpperCellId(int upper_cell_id);
    void setUpper(HCell* upper_);

    HCell* getUpper();

    friend class IrregularGrid2D;

  private:
    HCell* upper;
    int upper_cellId;
};

class VCell : public Cell {
  public:
    VCell();
    ~VCell();

    std::vector<HCell> &getRow();

  private:
    // row at a specific level
    std::vector<HCell> row {};
};

template <unsigned dim>
class Rectangle {
  public:
    Rectangle(double min[], double max[]);
    ~Rectangle();

  double getMinX() {
    return min[0];
  }

  double getMaxX() {
    return max[0];
  }

  double getMinY() {
    return dim < 2 ? 0 : min[1];
  }

  double getMaxY() {
    return dim < 2 ? 0 : max[1];
  }

  bool Intersects( const Rectangle<dim>& r ) const;
  Rectangle<dim> Intersection( const Rectangle<dim>& b) const;
  bool Contains( const Rectangle<dim>& r) const;

  private:
    double min[dim];
    double max[dim];
};
// Represents a row cell as rectangle
/* struct CellInfo {
  int cellId;
  int upperCellId;
  Rectangle<2>* cell;
  int statNbrOfPoints;

  CellInfo (int c_id, int uc_id, int nbr_points,
    double left, double right, double bottom, double top) {
    cellId = c_id;
    upperCellId = uc_id;
    statNbrOfPoints = nbr_points;

    double min[2], max[2];
    min[0] = left;
    min[1] = bottom;
    max[0] = right;
    max[1] = top;

    cell = new Rectangle<2>(min, max);
  }
}; */

class IrregularGrid2D {
  public:
    IrregularGrid2D(double l, double r, double b, double t,
      int row_count, int cell_count, char *conn_str, char *rel_str);
    ~IrregularGrid2D();

    Rectangle<2>* getBoundingBox();
    int getRowCount();
    int getCellCount();

    void* initIrgrid2D(double l, double r, double b, double t,
      int row_count, int cell_count, char *conn_str, char *rel_str);
    void destroyIrregularGrid2D( void *irgrid2d );

    bool createIrgrid2D( void *irgrid2d );
    bool generateIrgrid2D();

    ARInfo getCellnos(void *irgrid2d, double r_l, double r_r, double r_b,
     double r_t, char *conn_str);
    ARInfo cellnos(double r_l, double r_r, double r_b, double r_t,
     char *conn_str);

    // print output in psql terminal
    static void outTxtPsql(std::string txt_str) {
      if (!txt_str.empty()) {
        txt_str += "\n";
        elog(INFO, "Info: %s", txt_str.c_str());
      }
    }

    void setColumnVector(std::vector<VCell> column_vect);
    std::vector<VCell> &getColumnVector();

    void storeIrgrid2d();
    //static IrregularGrid2D* loadGrid(char *grid_id, char *conn_str);

  private:
    // points sorted by y-coordinates
  std::vector<RPoint> points{};
  // irregular grid bounding box
  Rectangle<2> * boundingBox;
  // number of rows and cells per row
  int rowCount, cellCount;
  // column (y-axis) vector with access to the section vectors (x-axes)
  std::vector<VCell> columnVector {};

  bool created;

  char * connection_str;
  char * relation_str;

  // void createIrgrid2D(Stream<Rectangle<2>> rStream);
  void processInput();
  void buildGrid();
};
#endif

#ifdef __cplusplus
extern "C" {
#endif
void* initIrgrid2D(double l, double r, double b, double t,
  int row_count, int cell_count, char *conn_str, char *rel_str);

void destroyIrregularGrid2D( void *irgrid2d );

bool createIrgrid2D( void *irgrid2d );

ARInfo getCellnos(void *irgrid2d, double r_l, double r_r, double r_b,
 double r_t, char *conn_str);
#ifdef __cplusplus
}; 
#endif

#endif /* IRREGULARGRID2D_H_ */
