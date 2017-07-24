///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
//
// MCG4139 Final Project
//
// Study of the effect on aerodynamic drag of the angle of
// inclination of a cab-mounted deflector from a tractor-
// trailer assembly
//
// Project done by Jonathan Charbonneau
// Student number 7199186
//
// Submitted to Professor James McDonald
//
// April 28th, 2017
// University of Ottawa
// Faculty of Engineering
// Department of Mechanical Engineering
//
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

#include<iostream>
#include<iomanip>
#include<string>
#include<sstream>
#include<fstream>
#include<functional>
#include<vector>
#include<limits>
#include<utility>
// This code used the "eigen3" linear-algebra library.
// It is available at eigen.tuxfamily.org/
#include<Eigen/Core>
#include<Eigen/Dense>
#include<stdexcept>
#include<chrono>
#include<thread>
#include"splitTime.h"


///////////////////////////////////////////////////////////////////////
//                    Some useful types
///////////////////////////////////////////////////////////////////////
struct Cell {
  int node0;
  int node1;
  int node2;
};

using Node2D = Eigen::Vector2d;
enum class location {truck, freestream, road, other};

struct Edge {
  int node0;
  int node1;
  int cell_l;
  int cell_r;
  location locate;
  Node2D fluid_n_hat;
  double l;
};

double max_x = -1e8;
double max_y = -1e8;
double min_x = 1e8;
double min_y = 1e8;

///////////////////////////////////////////////////////////////////////
//                    Simple Geometry stuff
///////////////////////////////////////////////////////////////////////
auto n_hat_and_length(const Node2D& n0, const Node2D& n1) {
  auto length = (n1-n0).norm();
  Node2D n_hat;
  n_hat.x() =  (n1.y()-n0.y())/length;
  n_hat.y() = -(n1.x()-n0.x())/length;
  return std::make_pair(n_hat, length);
}

///////////////////////////////////////////////////////////////////////
//                    Read gmsh
///////////////////////////////////////////////////////////////////////
// TODO: Put this function in a separate file
auto read_gmsh_file(const std::string& filename) {
 
  std::ifstream fin(filename);
  if(!fin) {
    throw std::runtime_error("Cannot open file: " + filename + ".");
  }

  ///////////////////////////////////////////////////////
  // Lambda function to consume lines that are expected
  // but should be ignored
  auto expect_line = [filename] (std::ifstream& fin, const std::string& expected) {
    std::string s;
    do {
      std::getline(fin,s);
    } while (s.empty());

    if(s != expected) {
      throw std::runtime_error("Error reading file: " + filename + ".\n" +
                               "Expected \"" + expected + "\", but got \"" + s +
                               "\".");
    }
  };


  std::cout << "Reading gmsh file: " << filename << '\n';

  ///////////////////////////////////////////////////////
  // Read file
  expect_line(fin, "$MeshFormat");
  expect_line(fin, "2.2 0 8");
  expect_line(fin, "$EndMeshFormat");
  expect_line(fin, "$Nodes");

  int number_of_nodes;
  
  fin >> number_of_nodes;

  std::vector<Node2D> nodes(number_of_nodes);

  // Loop through all nodes
  for(int i = 0; i < number_of_nodes; ++i) {
    int    dummy_index;
    double dummy_z_coordinate;
    fin >> dummy_index
        >> nodes[i].x()
        >> nodes[i].y()
        >> dummy_z_coordinate;
    // Check for boundaries
    if(nodes[i].x() > max_x){
      max_x = nodes[i].x();
    }
    if(nodes[i].x() < min_x){
      min_x = nodes[i].x();
    }
    if(nodes[i].y() > max_y){
      max_y = nodes[i].y();
    }
    if(nodes[i].y() < min_y){
      min_y = nodes[i].y();
    }
    
    if(dummy_index != i+1) {
      throw std::runtime_error("Error with node index.");
    }
    if(fabs(dummy_z_coordinate) > 1.0e-12) {
      throw std::runtime_error("Error, node has z component.");
    }
  }

  expect_line(fin, "$EndNodes");
  expect_line(fin, "$Elements");

  int number_of_elements; //not all will be cells
  fin >> number_of_elements;

  std::vector<Cell> cells;
  for(int i = 0; i < number_of_elements; ++i) {
    std::string s;
    do {
      std::getline(fin,s);
    } while (s.empty());

    std::istringstream ss(s);

    int element_num;
    ss >> element_num;
    if(element_num - 1 != i) {
      throw std::runtime_error("Error reading element number.");
    }

    int element_type;
    ss >> element_type;

    if(element_type == 2) {
      //triangular cell

      int dummy;
      Cell c;

      ss >> dummy >> dummy >> dummy >> c.node0 >> c.node1 >> c.node2;

      c.node0 -= 1;
      c.node1 -= 1;
      c.node2 -= 1;

      cells.push_back(c);

    }
  }

  expect_line(fin, "$EndElements");

  std::cout << "done.\n";
  std::cout << "\nFound boundaries: \n\tmax_x = " << max_x
	    << "\n\tmin_x = " << min_x
	    << "\n\tmax_y = " << max_y
	    << "\n\tmin_y = " << min_y << "\n\n";

  return std::make_pair(nodes, cells);
  
}

///////////////////////////////////////////////////////////////////////
//                      Compute Edges
///////////////////////////////////////////////////////////////////////
auto compute_edges(const std::vector<Node2D>& nodes, const std::vector<Cell>& cells,
		   const double max_x, const double max_y, const double min_x, const double min_y){

  std::cout << "Computing Edges.\n";
  double target = 0.01;
  double tolerance = 1e-6;

  std::vector<Edge> edges;

  for(int i = 0; i < cells.size(); ++i) {

    const auto& cell = cells[i];

    auto e = std::find_if(edges.begin(), edges.end(), [&cell] (const Edge& edge) {
        if((cell.node0 == edge.node0) && (cell.node1 == edge.node1)) return true;
        if((cell.node0 == edge.node1) && (cell.node1 == edge.node0)) return true;
        return false;
      });

    if(e != edges.end()) { //it was found
      e->cell_r = i;
    } else { // new edge
      edges.push_back({cell.node0, cell.node1, i, -1, location::other});
    }

    e = std::find_if(edges.begin(), edges.end(), [&cell] (const Edge& edge) {
        if((cell.node1 == edge.node0) && (cell.node2 == edge.node1)) return true;
        if((cell.node1 == edge.node1) && (cell.node2 == edge.node0)) return true;
        return false;
      });

    if(e != edges.end()) { //it was found
      e->cell_r = i;
    } else { // new edge
      edges.push_back({cell.node1, cell.node2, i, -1, location::other});
    }

    e = std::find_if(edges.begin(), edges.end(), [&cell] (const Edge& edge) {
        if((cell.node2 == edge.node0) && (cell.node0 == edge.node1)) return true;
        if((cell.node2 == edge.node1) && (cell.node0 == edge.node0)) return true;
        return false;
      });

    if(e != edges.end()) { //it was found
      e->cell_r = i;
    } else { // new edge
      edges.push_back({cell.node2, cell.node0, i, -1, location::other});
    }

    if(static_cast<double>(i+1)/static_cast<double>(cells.size()) >= target) {
      std::cout << i+1 << "/" << cells.size() << " cells processed.\n";
      std::cout.flush();
      target += 0.01;
    }

  }

  std::cout << "All cells processed.\n";

  //make sure all edges are propoerly aligned
  // (that the unit normal points from cell_l to cell_r)
    
  for(auto& edge: edges) {

    using std::swap;

    const auto edge_centroid = 0.5*(nodes[edge.node0]+nodes[edge.node1]);

    const auto& cell_left = cells[edge.cell_l];
    const auto left_cell_centroid = ( nodes[cell_left.node0]
                                     +nodes[cell_left.node1]
                                     +nodes[cell_left.node2])/3.0;

    const auto vec1 = left_cell_centroid - edge_centroid;

    auto n_hat_l = n_hat_and_length(nodes[edge.node0],nodes[edge.node1]);
    const auto& n_hat = n_hat_l.first;

    if(n_hat.dot(vec1) > 0.0) swap(edge.cell_l, edge.cell_r);
    
    
    //////////////////////////////////////////////////////////
    // Locate Edges in the control volume

    const Node2D& n0 = nodes[edge.node0];
    const Node2D& n1 = nodes[edge.node1];

    if((edge.cell_l == -1) || (edge.cell_r == -1)){
      if((n1.x() > (max_x - tolerance) || (n1.x() < (min_x + tolerance)))){
	// if it is the front or the back, set freestream
	edge.locate = location::freestream;
      }
      else if(n1.y() > (max_y - tolerance)){
	// if it is the top edge, set freestream
	edge.locate = location::freestream;
      }
      else if(n1.y() < (min_y + tolerance)){
	//if it is the bottom edge, reflect values
	edge.locate = location::road;
      }
      else{
	// if it's not a boundary, it's the truck
	edge.locate = location::truck;	
      }
    }
    else{
      edge.locate = location::other;
    }

  if(edge.locate == location::truck){


	const auto edge_centroid = 0.5*(nodes[edge.node0]+nodes[edge.node1]);

	if(edge.cell_r == -1){
	  const auto& cell_left = cells[edge.cell_l];
	  const auto left_cell_centroid = ( nodes[cell_left.node0]
					    +nodes[cell_left.node1]
					    +nodes[cell_left.node2])/3.0;

	  const auto vec2 = left_cell_centroid - edge_centroid;
	  std::tie(edge.fluid_n_hat,edge.l) = n_hat_and_length(nodes[edge.node0], nodes[edge.node1]);
	  
	  
	  if(edge.fluid_n_hat.dot(vec2) > 0.0) edge.fluid_n_hat.x() *= -1.0;
	} else{
	  const auto& cell_right = cells[edge.cell_r];
	  const auto right_cell_centroid = ( nodes[cell_right.node0]
					    +nodes[cell_right.node1]
					    +nodes[cell_right.node2])/3.0;

	  const auto vec3 = right_cell_centroid - edge_centroid;

	  std::tie(edge.fluid_n_hat,edge.l) = n_hat_and_length(nodes[edge.node0], nodes[edge.node1]);
	  
	  if(edge.fluid_n_hat.dot(vec3) > 0.0) edge.fluid_n_hat.x() *= -1.0;
	  
	}

  }// if statement
  }
  return edges;

}



///////////////////////////////////////////////////////////////////////
//                    Local Lax-Friedrichs
///////////////////////////////////////////////////////////////////////
template<typename PDE_type, typename Vec_type>
auto Local_Lax_Friedrichs(const Vec_type& Ul, const Vec_type& Ur) {
  
  const auto Fl = PDE_type::Fx(Ul);
  const auto Fr = PDE_type::Fx(Ur);
  const auto max_lambda = std::max(PDE_type::max_lambda_x(Ul),PDE_type::max_lambda_x(Ur));
  typename PDE_type::Vector_type F = 0.5*(Fl+Fr-max_lambda*(Ur-Ul));
  return F;
}

///////////////////////////////////////////////////////////////////////
//                      Euler Equations
///////////////////////////////////////////////////////////////////////

class EulerEquations {
public:
 
  using Vector_type = Eigen::Vector4d;
  static constexpr int number_of_unknowns = 4;

  ////////////////////////////////////
  // density of the air
  template<typename Vec_in>
  static double rho(const Vec_in& U){
    return U[0];
  }

  ////////////////////////////////////
  // ux
  template<typename Vec_in>
  static double ux(const Vec_in& U){
    return U[1]/U[0];
  }

  ////////////////////////////////////
  // uy
  template<typename Vec_in>
  static double uy(const Vec_in& U){
    return U[2]/U[0];
  }

  ////////////////////////////////////
  // pressure
  template<typename Vec_in>
  static double p(const Vec_in& U){
    return (U[3] - 0.5*(U[1]*U[1]+U[2]*U[2])/U[0])*(k-1);
  }

  /////////////////////////////////////
  // Free-stream conditions
  double fs_rho = 1.225;
  double fs_ux  = 105.0/3.6;
  double fs_p   = 101325.0;
  Eigen::Vector4d U_fs = {fs_rho, fs_rho*fs_ux, 0.0,fs_p/(k-1)+0.5*fs_rho*fs_ux*fs_ux};
   
  ////////////////////////////////////
  //  Flux x
  template<typename Vec_in>
  static Vector_type Fx(const Vec_in& U) {
    Vector_type F;
    auto u_x = ux(U);
    auto u_y = uy(U);
    auto pressure = p(U);
    auto dens = rho(U);
    F[0] = U[1];
    F[1] = dens*u_x*u_x + pressure;
    F[2] = dens*u_x*u_y;
    F[3] = u_x*(0.5*dens*(u_x*u_x + u_y*u_y) + k*pressure/(k-1));
    return F;
  }

  ////////////////////////////////////
  //  Rotate
  template<typename Vec_in>
  static Vector_type rotate(const Vec_in& U, Node2D n_hat) {
    Vector_type U_rot = U;
    U_rot[1] =  U[1]*n_hat.x() + U[2]*n_hat.y();
    U_rot[2] = -U[1]*n_hat.y() + U[2]*n_hat.x();
    return U_rot;
  }

  ////////////////////////////////////
  //  reflect_x
  template<typename Vec_in>
  static Vector_type reflect_x(const Vec_in& U) {
    Vector_type U_ref = U;
    U_ref[1] *= -1.0;
    return U_ref;
  }

  ////////////////////////////////////
  //  max_lambda_x
  template<typename Vec_in>
  static double max_lambda_x(const Vec_in& U) {
    return fabs(U[1]/U[0]) + sqrt(k*p(U)/U[0]);
  }

private:
  static constexpr double k = 1.4;
};


///////////////////////////////////////////////////////////////////////
//             Unstructured Finite-Volume Scheme
///////////////////////////////////////////////////////////////////////
template<typename PDE_type>
class Unstructured_FV_Solver {
public:

  using SolutionVector_type = typename PDE_type::Vector_type;

  Unstructured_FV_Solver() = default;
  Unstructured_FV_Solver(const Unstructured_FV_Solver&) = default;
  Unstructured_FV_Solver(Unstructured_FV_Solver&&) = default;
  Unstructured_FV_Solver& operator=(const Unstructured_FV_Solver&) = default;
  Unstructured_FV_Solver& operator=(Unstructured_FV_Solver&&) = default;

  EulerEquations ee;

  Unstructured_FV_Solver(const std::string& mesh_filename) {

    find_number_of_cores();
    
    std::tie(nodes,cells) = read_gmsh_file(mesh_filename);
    edges = compute_edges(nodes,cells,max_x,max_y,min_x,min_y);

    Global_U.resize(number_of_cells() * PDE_type::number_of_unknowns);
    Global_dUdt.resize(number_of_cells() * PDE_type::number_of_unknowns);
  
    for(int i = 0; i < number_of_cells(); ++i) {
      U(i) = ee.U_fs; // Initialize computation with free-stream conditions
    }
    compute_areas();
    compute_frontal_area();
    time = 0.0;
    iteration = 0;
  }

  auto U(int i) {
    return Global_U.segment<PDE_type::number_of_unknowns>(i*PDE_type::number_of_unknowns);
  }

  auto dUdt(int i) {
    return Global_dUdt.segment<PDE_type::number_of_unknowns>(i*PDE_type::number_of_unknowns);
  }

  auto number_of_cells() {return cells.size();}
  auto number_of_nodes() {return nodes.size();}
  auto number_of_edges() {return edges.size();}

  void write_to_vtk(const std::string& filename) {

    std::ofstream fout(filename);
    if(!fout) {
      throw std::runtime_error("error opening vtk file.");
    }

    fout << "# vtk DataFile Version 2.0\n"
         << "Unstructured Solver\n"
         << "ASCII\n"
         << "DATASET UNSTRUCTURED_GRID\n"
         << "POINTS " << number_of_nodes() << " double\n";

    for(const auto& node : nodes) {
      fout << node.x() << " " << node.y() << " 0\n";
    }

    fout << "\nCELLS " << number_of_cells() << " " << 4*number_of_cells() << '\n';
    for(const auto& cell : cells) {
      fout << "3 " << cell.node0 << " " << cell.node1 << " " << cell.node2 << '\n';
    }

    fout << "\nCELL_TYPES " << number_of_cells() << '\n';
    for(int i = 0; i < number_of_cells(); ++i) {
      fout << "5\n";
    }

    fout << "\nCELL_DATA " << number_of_cells() << '\n'
	 << "SCALARS rho double 1\n"
	 << "LOOKUP_TABLE default\n";
    for(int i = 0; i < number_of_cells(); ++i){
      fout << U(i)[0] << '\n';
    }

    fout << "\nVECTORS u double\n";
    for(int i = 0; i < number_of_cells(); ++i){
      fout << U(i)[1]/U(i)[0] << " " << U(i)[2]/U(i)[0] << " 0.0\n";
    }

    fout << "\nSCALARS p double 1\n"
	 << "LOOkUP_TABLE default\n";
    for(int i = 0; i< number_of_cells(); ++i){
      fout << ee.p(U(i)) << '\n';
    }
    
  }
  
  void time_march_to_time(double final_time, double CFL, int i) {

    const double tolerance = 1.0e-12;      
    while(time < final_time - tolerance) {

      Global_dUdt.fill(0.0);
      double dt = std::numeric_limits<double>::max();

      for(const auto& edge : edges) {

        double l;
        Node2D n_hat;
	std::tie(n_hat,l) = n_hat_and_length(nodes[edge.node0], nodes[edge.node1]);
	

        typename PDE_type::Vector_type Ul;
        typename PDE_type::Vector_type Ur;
        typename PDE_type::Vector_type Ul_rot;
        typename PDE_type::Vector_type Ur_rot;

        if(edge.cell_l != -1) {
          Ul = U(edge.cell_l);
          Ul_rot = PDE_type::rotate(Ul, n_hat);
        }

        if(edge.cell_r != -1) {
          Ur = U(edge.cell_r);
          Ur_rot = PDE_type::rotate(Ur, n_hat);
        }
	
        if(edge.cell_l == -1) {
	  
	  switch(edge.locate){
	  case location::road: Ul_rot = PDE_type::reflect_x(Ur_rot); break;
	  case location::freestream: Ul_rot = PDE_type::rotate(ee.U_fs, n_hat); break;
	  case location::truck: Ul_rot = PDE_type::reflect_x(Ur_rot); break;
	  case location::other: throw std::runtime_error("Problem with edges!"); break;
	  default: throw std::runtime_error ("Error with boundary checks!"); 	
	  } // switch 
   	  
        }// if statement
	  
        if(edge.cell_r == -1) {
	  
	  switch(edge.locate){
	  case location::road: Ur_rot = PDE_type::reflect_x(Ul_rot); break;
	  case location::freestream: Ur_rot = PDE_type::rotate(ee.U_fs, n_hat); break;
	  case location::truck: Ur_rot = PDE_type::reflect_x(Ul_rot); break;
	  case location::other: throw std::runtime_error("Problem with edges!"); break;
	  default: throw std::runtime_error("Error with boundary checks!");
	  } // switch
	  
        }// if statement
		
        auto F_rot = Local_Lax_Friedrichs<PDE_type>(Ul_rot,Ur_rot);
	
        n_hat.y() *= -1.0;
        auto F = PDE_type::rotate(F_rot,n_hat);        

        if(edge.cell_l != -1) {
          dUdt(edge.cell_l) -= F*l/areas[edge.cell_l];
          dt = std::min(dt, CFL*(areas[edge.cell_l]/l)/PDE_type::max_lambda_x(Ul));
        }

        if(edge.cell_r != -1) {
          dUdt(edge.cell_r) += F*l/areas[edge.cell_r];
          dt = std::min(dt, CFL*(areas[edge.cell_r]/l)/PDE_type::max_lambda_x(Ur));
        }
		
      } //end loop over edges

      Global_U += dt*Global_dUdt;
      time += dt;
      ++ iteration;

      std::cout << "Time = " << time
                << "    dt = " << dt
                << "    residual = " << Global_dUdt.norm() << '\n';
      
    }
    postProcessing(i);
    
  }// end time march to time

  void postProcessing(int i){
    
    double pressure = 0;
    double force = 0;
    double cd = 0;
    
    for(const auto& edge : edges){
      if(edge.locate == location::truck){

	if(edge.cell_r == -1){
	  pressure = ee.p(U(edge.cell_l));
	} else{
	  pressure = ee.p(U(edge.cell_r));
	}
	
	force += 2.0 * pressure * edge.fluid_n_hat.x() * edge.l;
      }// if statement
    }// for loop

    cd = fabs(force)/(0.5 * ee.fs_rho * ee.fs_ux * ee.fs_ux * frontal_area*2.0);
    
    std::ofstream fout("force_output.txt", std::ofstream::app);
    if(!fout){
      throw std::runtime_error("Could not write force_output.txt");
    }
    
    fout  << i << "\t" << force << '\t' << "Cd = " << cd << '\n';
    
  }
  


  void make_movie(double final_time, double CFL, int num_frames, const std::string& filename_base) {
    const auto dt = final_time/static_cast<double>(num_frames-1);

    auto get_name = [&filename_base] (int i) {
      std::stringstream ss;
      ss << filename_base << std::setfill('0') << std::setw(5) << i << ".vtk";
      return ss.str();
    };

    // The reason for post-processing is to make sure that there is zero drag
    // when the initial conditions are set.
    postProcessing(0);
    std::cout << "Writing initial conditions to file " << get_name(0) <<'\n'; 
    write_to_vtk(get_name(0));

    // Start a timer to calculate the total computation time
    std::chrono::time_point<std::chrono::system_clock> startTime, stopTime;
    startTime = std::chrono::system_clock::now();

    // Meat and potatoes of the code
    for(int i = 1; i < num_frames; ++i) {
      std::cout << "Time Marching to time =" << dt*static_cast<double>(i) << '\n';
      time_march_to_time(dt*static_cast<double>(i), CFL, i);
      std::cout << "Writing output to file " << get_name(i) <<'\n';
      write_to_vtk(get_name(i));
      }

    // Stop the timer and display the time
    stopTime = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = stopTime - startTime;
    double timeArray[4] = {0};
    splitTime(elapsed_seconds.count(),timeArray);
    std::cout << "Total computation time = "
	      << floor(elapsed_seconds.count())
	      << " seconds or "
	      << timeArray [0] << " days " << timeArray [1] << " hours, "
	      << timeArray [2] << " minutes and " << timeArray [3] << " seconds.\n";

  }

  

private:
  std::vector<Node2D>    nodes;
  std::vector<Cell>      cells;
  std::vector<Edge>      edges;
  std::vector<double>    areas;
  Eigen::VectorXd     Global_U;
  Eigen::VectorXd  Global_dUdt;
  double                  time;
  double             iteration;
  double          frontal_area;
  unsigned int               n;
  
  void compute_areas() {
       
    areas.resize(number_of_cells());
    for(int i = 0; i < number_of_cells(); ++i) {
      const Cell& cell = cells[i];
      const Node2D& n0 = nodes[cell.node0];
      const Node2D& n1 = nodes[cell.node1];
      const Node2D& n2 = nodes[cell.node2];
      areas[i] = 0.5*fabs( n0.x()*n1.y()-n0.y()*n1.x()
                          +n1.x()*n2.y()-n1.y()*n2.x()
                          +n2.x()*n0.y()-n2.y()*n0.x());
    }
  }

  void compute_frontal_area(){
    for(auto& edge: edges) {
      if (edge.locate == location::truck){
	double l;
        Node2D n_hat;
	std::tie(n_hat,l) = n_hat_and_length(nodes[edge.node0], nodes[edge.node1]);
	if (n_hat.x() < 0){
	  frontal_area += fabs(n_hat.x() * l);
	}
      }
    }
    std::cout << "Frontal area = " << frontal_area << '\n';
  }

  void find_number_of_cores(){
  n = std::thread::hardware_concurrency();
  std::cout << n << " concurrent threads are supported.\n";
  }

};


///////////////////////////////////////////////////////////////////////
//                         Main
///////////////////////////////////////////////////////////////////////
int main() {

  std::cout << " --------------------------------------------------------\n"
               " |                                                      |\n"
               " |                  MCG 4139 - Project                  |\n"
               " |                                                      |\n"
               " --------------------------------------------------------\n\n";

  double final_time = 6.0;
  double CFL = 0.75;
  int num_frames = 401;

  // Get mesh filename
  std::string mesh_filename;
  std::cout << "Please enter the filename for analysis: \n";
  getline(std::cin,mesh_filename);

  // Actual program
  auto solver = Unstructured_FV_Solver<EulerEquations>(mesh_filename);
  solver.make_movie(final_time, CFL, num_frames,"movie");
      
  return 0;
}


  
