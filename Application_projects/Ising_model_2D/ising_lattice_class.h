/*******************************************************************************
*
* Lattice class definition
*
*******************************************************************************/

//--- Preprocessor directives --------------------------------------------------

#ifndef ISING_LATTICE_CLASS_H
#define ISING_LATTICE_CLASS_H

#include <iostream>
#include <fstream>

#include <string>
#include <vector>
#include <random>
#include <cmath>

using namespace std;

// Define the PRNG
#define SEED 42
random_device device;
mt19937 generator(device());
//mt19937 generator(SEED);

// Define precise constants
constexpr double pi = 3.14159265358979323846;

//--- Contents -----------------------------------------------------------------

class lattice {

private:
    int tot_lenght_;
    vector<int> lattice_;
    vector<vector<int>> nearest_neighbors_;

public:
    const int side_lenght, geometry_flag, initial_flag;

    lattice(const int &SIDE, const int &G_FLAG, const int &I_FLAG):
       /* Class lattice constructor **********************************

       Inputs:

       SIDE = number of sites per side.

       G_FLAG =  geometry_flag for the topology.
                 1 for 1D lattice with PBC.
                 2 for 2D square lattice with PBC.
                 else gives "Error: not implemented yet".

       I_FLAG = initial_flag for the initial condition.
                0 for initialising cold lattice.
                else initialise hot lattice.

        *************************************************************/
        side_lenght(SIDE),
        geometry_flag(G_FLAG),
        initial_flag(I_FLAG)
        {// CONSTRUCTION BEGIN

        double random_number;
        vector<int> nearest_list;

        // Defining topology and nearest neighbors
        if (geometry_flag == 1) {
            tot_lenght_ = side_lenght;
            lattice_.reserve(tot_lenght_);
            nearest_neighbors_.reserve(tot_lenght_);
            nearest_list.reserve(2);

            // Fill nearest_neighbors_ with nearest_list for each site
            for (int i = 0; i < tot_lenght_; i++){
                // Fill nearest_list for site i
                nearest_list.clear();
                nearest_list.push_back(i - 1);
                nearest_list.push_back(i + 1);

                // Implementing PBC
                if (i == 0){
                    // Correct nearest_list first element
                    nearest_list[0] = tot_lenght_-1;
                } else if (i == tot_lenght_-1){
                    // Correct nearest_list last element
                    nearest_list[1] = 0;
                }

                // Push nearest_list into the nearest_neighbors_ vector
                nearest_neighbors_.push_back(nearest_list);
            }
            nearest_list.clear();
        } else if (geometry_flag == 2) {
            tot_lenght_ = side_lenght * side_lenght;
            lattice_.reserve(tot_lenght_);
            nearest_neighbors_.reserve(tot_lenght_);
            nearest_list.reserve(4);

            // Fill nearest_neighbors_ with nearest_list for each site
            for (int i = 0; i < tot_lenght_; i++){
                // Fill nearest_list for site i
                nearest_list.clear();
                nearest_list.push_back(i - 1);
                nearest_list.push_back(i + 1);
                nearest_list.push_back(i - side_lenght);
                nearest_list.push_back(i + side_lenght);

                // Implementing PBC
                if (i % side_lenght  == 0){
                    // Correct nearest_list first culomn elements
                    nearest_list[0] += side_lenght;
                } else if ((i + 1) % side_lenght  == 0){
                    // Correct nearest_list last culomn elements
                    nearest_list[1] -= side_lenght;
                }
                if (i < side_lenght){
                    // Correct nearest_list first row elements
                    nearest_list[2] += tot_lenght_;
                } else if (i >= tot_lenght_ - side_lenght){
                    // Correct nearest_list last row elements
                    nearest_list[3] -= tot_lenght_;
                }

                // Push nearest_list into the nearest_neighbors_ vector
                nearest_neighbors_.push_back(nearest_list);
            }
            nearest_list.clear();
        } else {
            cerr << "Error: geometry not implemented." << endl;
            exit(1);
        }



        // Initialising the lattice configuration
        if (initial_flag == 0) {
            for (int i = 0; i < tot_lenght_; i++) lattice_.push_back(1);
        } else {
            for (int i = 0; i < tot_lenght_; i++){
                random_number = rand_double();
                if (random_number < 0.5){
                    lattice_.push_back(-1);
                } else {
                    lattice_.push_back(1);
                }
            }
        }
    }// CONSTRUCTION END

    int rand_int(){
        /* Generate a random index for the lattice */

        uniform_int_distribution<long int> randomint(0, tot_lenght_ - 1);
        return randomint(generator);
    }

    double rand_double(){
        /* Generate a random double */

        uniform_real_distribution<double> prng(0.0, 1.0);
        return prng(generator);
    }

    double energy(double extfield){
        /* Compute the energy of the present configuration */

        double sum, ener = 0.;
        for(int i = 0; i < tot_lenght_; i++) {
           sum = 0.;
           for(auto nn : nearest_neighbors_[i]) sum += lattice_[nn];
           ener += -0.5 * lattice_[i] * sum - extfield * lattice_[i];
        }
        ener = ener / tot_lenght_;
        return ener;
    }

    double magnetization(){
        /* Compute the magnetization of the present configuration */

        double sum = 0;
        for(int i = 0; i < tot_lenght_; i++) sum += lattice_[i];
        sum = sum / tot_lenght_;
        return sum;
    }

    void show_configuration(){
        /* Print the lattice configuration */

        int value;
        cout << "Lattice configuration: " << endl;
        if (geometry_flag == 1){
            for (int i = 0; i < tot_lenght_; i++){
                value = lattice_[i];
                if (value == -1) value = 0;
                cout << value << " ";
            }
            cout << endl << endl;
        } else if (geometry_flag == 2){
            for (int i = 0; i < tot_lenght_; i++){
                value = lattice_[i];
                if (value == -1) value = 0;
                cout << value << " ";
                if ((i + 1) % side_lenght == 0){
                    cout << endl;
                }
            }
            cout << endl;
        }
    }

    void save_configuration(string file_state = "ising_state.dat"){
        /* Save the current configuration */

        ofstream file;
        file.open(file_state);
        for (int i = 0; i < tot_lenght_; i++){
            file << lattice_[i] << endl;
        }
        file.close();
    }

    void load_configuration(string file_state = "ising_state.dat"){
        /* Load configuration from path */

        int i = 0;
        ifstream file(file_state);
        cout << "Loading initial configuration..." << endl;
        if (file.is_open()) {
            while((file >> lattice_[i]) && (i < tot_lenght_)) i++;
            file.close();
        } else {
            cerr << "Error: unable to open the file." << endl;
            exit(1);
        }

        if (i != tot_lenght_){
            cerr << "Error: different number of sites." << endl;
            exit(1);
        }
    }

    void show_nearest_index(const int &index){
        /* Print the nearest neighbors of the site index */

        int size;
        size = nearest_neighbors_[index].size();
        for (int i = 0; i < size; i++){
            cout << nearest_neighbors_[index][i] << " ";
        }
    }

    void show_nearest_neighbors(){
        /* Print the list of all nearest neighbors to each lattice site */

        for (int i = 0; i < tot_lenght_; i++){
            cout << "Nearest neighbors to " << i << " are: ";
            show_nearest_index(i);
            cout << endl;
        }
        cout << endl;
    }

    void update(double beta, double extfield){
        /* Update the state of the lattice with a MC step */

        int ind;
        double random_number, force;

        // Metropolis-Hastings algorithm per Ising
        for(int it = 0; it < tot_lenght_; it++){
           force = 0.;
           ind = rand_int();
           random_number = rand_double();

           // MC step for the site with index ind
           for(auto nn : nearest_neighbors_[ind]) force += lattice_[nn];

           force = beta * (force + extfield);

           force = exp(-2.0 * force * lattice_[ind]);

           // accept or reject step
           if (random_number < force) lattice_[ind] = -lattice_[ind];
        }
    }

};

#endif
