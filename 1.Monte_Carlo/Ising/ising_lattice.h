/**
* Lattice class definition
*/

//----Preprocessor directives--------------------------------------
#ifndef ISING_LATTICE_H
#define ISING_LATTICE_H

#include <iostream>
#include <fstream>

#include <vector>
#include <random>
#include <cmath>

using namespace std;

//----Contents-----------------------------------------------------

class lattice {
private:
    int tot_lenght_;
    vector<int> lattice_;
    vector<vector<int>> nearest_neighbors_;

public:
    const int side_lenght, temperature_flag, geometry_flag;
    lattice(const int &L = 5, const int &G_FLAG = 0, const int &T_FLAG = 0):
        side_lenght(L),
        temperature_flag(T_FLAG),
        geometry_flag(G_FLAG)
    {
        /* Initialising the lattice */
        vector<int> nearest_list;
        if (geometry_flag == 0) {
            /* Initialising 1D lattice and nearest neighbors */
            tot_lenght_ = side_lenght;
            for (int i = 0; i < tot_lenght_; i++){
                nearest_list.clear();
                nearest_list.push_back(i - 1);
                nearest_list.push_back(i + 1);
                if (i == 0){
                  nearest_list[0] = tot_lenght_-1;
                } else if (i == tot_lenght_-1){
                  nearest_list[1] = 0;
                }
                nearest_neighbors_.push_back(nearest_list);
            }

        } else if (geometry_flag == 1) {
            /* Initialising 2D square lattice and nearest neighbors */
            tot_lenght_ = side_lenght * side_lenght;
        } else {
            /* Not implemented */
        }

        lattice_.reserve(tot_lenght_);

        if (temperature_flag == 0) {
            /* Initialising cold lattice */
            for (int i = 0; i < tot_lenght_; i++){
                lattice_.push_back(0.);
            }

        } else if (temperature_flag == 1) {
            /* Initialising hot lattice */
            for (int i = 0; i < tot_lenght_; i++){
                lattice_.push_back(1.);
            }

        } else {
            /* Not implemented */
        }

    }

    void show_configuration(){
        cout << "Lattice configuration: " << endl;
        for (int i = 0; i < tot_lenght_; i++){
            cout << lattice_[i] << " ";
        }
        cout << endl;
    }

    void show_nearest_index(const int &index){
        int size;
        size = nearest_neighbors_[index].size();

        for (int i = 0; i < size; i++){
            cout << nearest_neighbors_[index][i] << " ";
        }
    }

    void show_nearest_neighbors(){
        for (int i = 0; i < tot_lenght_; i++){
            cout << "The nearest neighbors to element " << i << " are: ";
            show_nearest_index(i);
            cout << endl;
        }
    }

};

#endif
