import source
import time
import numpy as np
import pymaster as nmt

def start(repository_name):

        print ' ---- START ----'
        start = time.time()

        NSIDE = 1024

        # MASK INPUT
        mask_name = repository_name + '/mock_cats/Y1redLSSmaskv1.04096ring.dat'
        source.make_mask(repository_name, mask_name, NSIDE)

        #HALOGEN MOCK CATALOGS
        nmocks = 1800

        # REDSHIFT BINNING
        photoz_min = 0.6
        photoz_max = 1.0
        nshells = 4

        # BINNING ELLS'S INPUTS
        lmin = 20
        lmax = 360
        dl = 20

        # Create binning scheme
        #bpws, ells, weights = binning(repository_name, lmin, lmax, dl)

        # Or reading binning scheme
        bin_name = repository_name + '/bin.dat'
        bpws, ells, weights = np.loadtxt(bin_name, unpack=True)

	Bin = nmt.NmtBin(NSIDE, ells=ells, bpws=bpws, weights=weights)
        #Bin = nmt.NmtBin(NSIDE, nlb=15, lmax=400)

        #OVERDENSITY MAPS FROM MOCKS
        #source.make_compute_all_serie(repository_name, NSIDE, Bin, photoz_min, photoz_max, nshells, nmocks)
        source.make_compute_all_parallel(repository_name, NSIDE, Bin, photoz_min, photoz_max, nshells, nmocks)

        # MAKING GHAPH PLOTS
        source.make_plots(repository_name, Bin, nshells, nmocks)

        # COMPUTING COVARIANCE MATRIX
        #source.make_covmat(repository_name, Bin, nshells, nmocks)


        end = time.time()
        dt = end - start
        print 'This job took %.1f seconds' %dt

        print ' ---- END ---'

        return


#if __name__ == '__main__':
#       start()

#Change for your directory!
repository_name = '/home/drc01/sobreira/andluiz/halogen'
start(repository_name)
