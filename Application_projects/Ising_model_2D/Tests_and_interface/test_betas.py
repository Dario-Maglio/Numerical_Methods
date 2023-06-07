"""*****************************************************************************
*
* Test fit intervals
*
*****************************************************************************"""

import os

import numpy as np

#--- Contents ------------------------------------------------------------------

def load_data():
    """ Load data produced by analysis.cpp """

    data = {}
    for side in [20, 30]:
        # define data file path
        filename = f"side_{side}_data.dat"
        file_path = os.path.join("Data_analysis", filename)
        file_path = os.path.join("..", file_path)
        print("Loading " + file_path)
        # load data from each side file
        if os.path.isfile(file_path):
            data[side] = np.loadtxt(file_path, unpack='True')

    return data

#--- Main ----------------------------------------------------------------------

if __name__ == '__main__':

    data = load_data()
    print("Loading complete! \n")

    x, _, _, _, _, _, _, y, y_err = data[20]
    x, y, y_err = zip(*sorted(zip(x, y, y_err)))

    for pair in enumerate(x):
        print(pair)
