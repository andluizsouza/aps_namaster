''' Wrapper to aps_measurements_namaster '''

# basic libraries
import os
import sys
import time
import mpi4py.MPI as MPI

# scientific libraries
import healpy as hp
import numpy as np
import shlex, subprocess
import matplotlib as plt
plt.use('Agg')
import skymapper as skm
import pymaster as nmt
import matplotlib.colors as colors


def plot_mask(mask, NSIDE, figname):

	sep = 15
	
	fig, ax, proj = skm.plotHealpix(mask, NSIDE, sep=sep, cb_label='Fracgood pixel')
        fig.savefig(figname)

	return

def plot_over(i, RA, DEC, NSIDE, figname):

	# returns a list of counts (in units of 1/arcmin^2)
	# get count in healpix cells (as vertices), restrict to non-empty cells
	bc, ra, dec, vertices = skm.getCountAtLocations(RA, DEC, nside=NSIDE, return_vertices=True)

	bc_mean = np.sum(bc)/len(bc)
	bc = bc/bc_mean - 1

        # setup figure
        import matplotlib.cm as cm
        cmap = cm.YlGnBu
        fig = plt.figure(figsize=(10,8))
        ax = fig.add_subplot(111, aspect='equal')

        # setup map: define AEA map optimal for given RA/Dec
        proj = skm.createConicMap(ax, ra, dec, proj_class=skm.AlbersEqualAreaProjection)
        # add lines and labels for meridians/parallels
        sep = 15
        meridians = np.arange(-90, 90+sep, sep)
        parallels = np.arange(0, 360+sep, sep)
        skm.setMeridianPatches(ax, proj, meridians, linestyle='-', lw=0.5, alpha=0.3, zorder=2)
        skm.setParallelPatches(ax, proj, parallels, linestyle='-', lw=0.5, alpha=0.3, zorder=2)
        skm.setMeridianLabels(ax, proj, meridians, loc="left", fmt=skm.pmDegFormatter)
        skm.setParallelLabels(ax, proj, parallels, loc="bottom")
	plt.title('Galaxy overdensity map - shell ' + str(i))

        # add vertices as polygons
        vmin, vmax = np.percentile(bc,[10,90])
        poly = skm.addPolygons(vertices, proj, ax, color=bc, vmin=vmin, vmax=vmax, cmap=cmap, zorder=3, rasterized=True)

        # add colorbar
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.0)
        cb = fig.colorbar(poly, cax=cax)
        cb.set_label('$\delta_g$ [overdensity]')
        cb.solids.set_edgecolor("face")

        # show (and save) ...
        fig.tight_layout()
        fig.savefig(figname)

	return

def make_mask(dir_name, mask_name, NSIDE):

	print 'Reading original mask'

        index_old, fracdet_old = np.loadtxt(mask_name, unpack=True)

        Npix_old = hp.nside2npix(4096)
        
	mask_old = np.zeros(Npix_old)
	#mask_old = np.ones(Npix_old)*hp.UNSEEN

        for i in range(len(index_old)):
                mask_old[int(index_old[i])] = fracdet_old[i]

        # Downgrade old mask to get the fracdets of new mask
        print 'Writing new mask'
        mask_new = hp.ud_grade(mask_old, NSIDE)
        #dir_name = os.getcwd()
	output_name = dir_name + '/outputs/mask.fits'
        hp.write_map(output_name, mask_new)

        Npix = hp.nside2npix(NSIDE)
        fsky = np.sum(mask_new)/Npix
        print "Fraction of survey sky:", fsky

	print 'Plotting mask'
        output_name = dir_name + '/outputs/mask.png'
        plot_mask(mask_new, NSIDE, output_name)

	return

def make_compute_all_serie(dir_name, NSIDE, Bin, photoz_min, photoz_max, nshells, nmocks):


        
	input_name = dir_name + '/outputs/mask.fits'
	#mask = hp.read_map('/home/drc01/sobreira/andluiz/halogen/outputs/mask.fits', verbose=False)
	mask = hp.read_map(input_name, verbose=False)


	# ******** Apodize mask ********
        # The following function calls create apodized versions of the raw mask
	# with an apodization scale of 1.0 degrees using three different methods
        # mask_apod = nmt.mask_apodization(mask, 1.0, apotype="Smooth")


	# Computing coupling matrix
	
	print 'Computing coupling matrix'
	wsp = coupling_matrix(NSIDE, mask, Bin)

	ell = Bin.get_effective_ells()
        np.savetxt(dir_name + '/outputs/cldata/ell.dat', ell.T)

	for i in range(nmocks):

		print ' --- MOCK '+ str(i+1) + ' --- ' 

		if i < 10:
			catalog_name = dir_name + '/mock_cats/HALOGENlamps_V2.0.0_DNF_mock000' + str(i) + '_masked.dat'
		elif i >= 10 and i < 100:
                        catalog_name = dir_name + '/mock_cats/HALOGENlamps_V2.0.0_DNF_mock00' + str(i) + '_masked.dat'
		elif i >= 100 and i < 1000:
                        catalog_name = dir_name + '/mock_cats/HALOGENlamps_V2.0.0_DNF_mock0' + str(i) + '_masked.dat'
		else:
                        catalog_name = dir_name + '/mock_cats/HALOGENlamps_V2.0.0_DNF_mock' + str(i) + '_masked.dat'

		#Getting overdensity maps for i-th mock
		overmaps = make_overmap_single(catalog_name, mask, NSIDE, photoz_min, photoz_max, nshells, i)

		print 'Computing namaster'
		for j in range(nshells):

	                dmap = overmaps[j,:]
			f0 = nmt.NmtField(mask, [dmap])

	                cl_decoupled = nmt.compute_full_master(f0, f0, Bin, cl_noise=None, cl_guess=None, workspace=wsp)[0]
		        output_name = dir_name + '/outputs/cldata/cl_mock' + str(i) + '_shell'  + str(j) + '.dat '
			np.savetxt(output_name, cl_decoupled)			

	return



def make_compute_all_parallel(dir_name, NSIDE, Bin, photoz_min, photoz_max, nshells, nmocks):

        input_name = dir_name + '/outputs/mask.fits'
        #mask = hp.read_map('/home/drc01/sobreira/andluiz/halogen/outputs/mask.fits', verbose=False)
        mask = hp.read_map(input_name, verbose=False)


        # ******** Apodize mask ********
        # The following function calls create apodized versions of the raw mask
        # with an apodization scale of 1.0 degrees using three different methods
        # mask_apod = nmt.mask_apodization(mask, 1.0, apotype="Smooth")
  

        # Computing coupling matrix
        print 'Computing coupling matrix'
        wsp = coupling_matrix(NSIDE, mask, Bin)


	# Creating array with all the measures
	nell = Bin.get_n_bands()
	cl_all = np.zeros(shape=(nmocks,nshells,nell))

	# Parallel version

  	# mpi4py has the notion of a "communicator" - a collection of processors all operating together, usually on the same program.
	# Each processor in the communicator is identified by a number, its rank,  We'll use that number to split the tasks
	comm = MPI.COMM_WORLD

	# find out which number processor this particular instance is, and how many there are in total
	rank = comm.Get_rank()
	size = comm.Get_size()

	# the enumerate function gives us a number i in addition to the task. 
	# (In this specific case i is the same as task! But that's not true usually)

	task_list = np.arange(nmocks)


	for i, task in enumerate(task_list):
		# This is how we split up the jobs.
		# The % sign is a modulus, and the "continue" means "skip the rest of this bit and go to the next time through the loop"
		# If we had e.g. 4 processors, this would mean that
		# processor zero did tasks 0, 4, 8, 12, 16, ...
		# and processor one did tasks 1, 5, 9, 13, 17, ... and do on.

		if i%size != rank: continue
		# print "Task number %d (%d) being done by processor %d of %d" % (i, task, rank, size)
		make_compute_single_parallel(dir_name, mask, wsp, NSIDE, photoz_min, photoz_max, nshells, Bin, task)

	
	ell = Bin.get_effective_ells()
        np.savetxt(dir_name + '/outputs/cldata/ell.dat', ell.T)

	return


def make_compute_single_parallel(dir_name, mask, wsp, NSIDE, photoz_min, photoz_max, nshells, Bin, i):

	print ' --- MOCK '+ str(i+1) + ' --- ' 

        if i < 10:
        	catalog_name = dir_name + '/mock_cats/HALOGENlamps_V2.0.0_DNF_mock000' + str(i) + '_masked.dat'
        elif i >= 10 and i < 100:
        	catalog_name = dir_name + '/mock_cats/HALOGENlamps_V2.0.0_DNF_mock00' + str(i) + '_masked.dat'
	elif i >= 100 and i < 1000:
        	catalog_name = dir_name + '/mock_cats/HALOGENlamps_V2.0.0_DNF_mock0' + str(i) + '_masked.dat'
	else:
        	catalog_name = dir_name + '/mock_cats/HALOGENlamps_V2.0.0_DNF_mock' + str(i) + '_masked.dat'

	#Getting overdensity maps for mock i
        overmaps = make_overmap_single(catalog_name, mask, NSIDE, photoz_min, photoz_max, nshells, i)

	#Setting cl of this mock
	nell = Bin.get_n_bands()
	cl_mock = np.zeros(shape=(nshells, nell))

        print 'Computing namaster'
        for j in range(nshells):

		#Overdensity map for mock i and redshift shell j
		dmap = overmaps[j,:]


		"""
		Computes the full MASTER estimate of the power spectrum of two fields (f1 and f2). This is equivalent to successively calling:

		    - :func:`pymaster.NmtWorkspace.compute_coupling_matrix`
                    - :func:`pymaster.deprojection_bias`
   		    - :func:`pymaster.compute_coupled_cell`
	  	    - :func:`pymaster.NmtWorkspace.decouple_cell`

		    :param NmtField f1,f2: fields to correlate
		    :param NmtBin b: binning scheme defining output bandpower
		    :param cl_noise: noise bias (i.e. angular power spectrum of masked noise realizations) (optional).
		    :param cl_guess: set of power spectra corresponding to a best-guess of the true power spectra of f1 and f2. 
			Needed only to compute the contaminant cleaning bias (optional).
		    :param NmtWorkspace workspace: object containing the mode-coupling matrix associated with an incomplete sky coverage. 
			If provided, the function will skip the computation of the mode-coupling matrix and use the given information.

		    :return: set of decoupled bandpowers

		"""

	        f0 = nmt.NmtField(mask, [dmap])
		cl_decoupled = nmt.compute_full_master(f0, f0, Bin, cl_noise=None, cl_guess=None, workspace=wsp)

		cl_mock[j,:] = cl_decoupled[0]

	#Saving measures on this mock i
        #dir_name = os.getcwd()
	#output_name = '/home/drc01/sobreira/andluiz/halogen/outputs/cldata/cl_mock' + str(i) + '.dat '
	output_name = dir_name + '/outputs/cldata/cl_mock' + str(i) + '.dat'
	np.savetxt(output_name, cl_mock.T)

        return


def make_overmap_single(catalog_name, mask, NSIDE, photoz_min, photoz_max, nshells, mock):

	Npix = hp.nside2npix(NSIDE)
	fsky = np.sum(mask)/Npix

        overmaps = np.zeros(shape=(nshells, Npix))

	print 'Reading mock'

	RA, DEC, ZPHOTO = np.loadtxt(catalog_name, unpack=True, usecols=(0,1,2))

        #Below are the correct ways to go from RA and DEC to Healpix's theta and phi
        theta = (180-(DEC+90))*(np.pi/180)
        phi = RA*(np.pi/180)

	#Redshift infos
        #print 'Min redshift:', np.amin(ZPHOTO)
	#print 'Mean redshift:', np.mean(ZPHOTO)
        #print 'Max redshfit:', np.amax(ZPHOTO)

	dz = (photoz_max - photoz_min)/nshells
 
        print 'Making overdensity maps'

	for i in range(nshells):

                #Start or restart counting in the shell
                ra = []
                dec = []

                photoz_1 = photoz_min + i*dz
                photoz_2 = photoz_min + (i+1)*dz
                #print ' -- Redshift shell', i+1, '--'
                #print photoz_1, '< photo-z <', photoz_2

                for j in range(len(RA)):
                        if (ZPHOTO[j]>photoz_1 and ZPHOTO[j]<photoz_2):
               		        ra.append(RA[j])
                                dec.append(DEC[j])
                               

                raarray = np.array(ra,'f')
                decarray = np.array(dec,'f')

                #print 'Plotting overdensity map'
	        #dir_name = os.getcwd()
                #figname = '/home/drc01/sobreira/andluiz/halogen/outputs/over_' + str(i) + '.png'
                #figname = dir_name + '/outputs/over_' + str(i) + '.png'
                #plot_over(i+1, raarray, decarray, NSIDE, figname)

                #Below are the correct ways to go from RA and DEC to Healpix's theta and phi (see Healpix documentation for those definitions)
                thetaarray = (180-(decarray+90))*(np.pi/180)
                phiarray = raarray*(np.pi/180)

                #Command below reads in thetas and phis from previous step into a list of Healpix pixel numbers. 
                #'healpixarray' below is a list of pixel numbers that is as long as the number of galaxies in the data
                healpixarray = hp.ang2pix(NSIDE, thetaarray, phiarray)

                #print 'Counting galaxies in each pixel'
                #index of occupied pixels and galaxies counting in each pixel
                map_index, map_count = np.unique(healpixarray, return_counts=True)

                #print 'Making overdensity map'

                #total number of galaxies
                Ngal = np.sum(map_count)
                #print 'Total number of galaxies:', Ngal
                
                #mean density per pixel
                mean_count_pix = Ngal/fsky/Npix

                #initializing overdensity		
		over = -np.ones(Npix)
		
		#inside of mask and there is counting in pixel
		for j in range(len(map_index)):
			pix = int(map_index[j])
			if mask[pix] != 0:
				over[pix] = map_count[j]/mean_count_pix/mask[pix] - 1.0

		#name_over = 'outputs/overmaps/over_mock' + str(mock) + '_shell' + str(i) + '.fits'
                #print 'Writing file over.fits'
                #hp.write_map(name_over, over)                        

		#Getting overdensity maps together
		overmaps[i,:] = over


        return overmaps


'''

https://namaster.readthedocs.io/en/latest/sample_workspaces.html

#This script showcases the use of NmtWorkspace objects to speed up
the computation of power spectra for many pairs of fields with the same masks.

'''

def binning(dir_name, lmin, lmax, dl):

        print 'Creating binning ells'

        nbands = int((lmax - lmin)/dl)


        for i in range(nbands+1):
                if i == 0:
                        ell = np.arange(lmin - int(dl/2) + i*dl, lmin - int(dl/2) + (i+1)*dl + 1, 1)
                        bandpowers = i*np.ones(dl+1)
                        weights = np.ones(dl+1)
                else:
                        ell = np.concatenate((ell, np.arange(lmin - int(dl/2) + i*dl, lmin - int(dl/2) + (i+1)*dl + 1, 1)))
                        bandpowers = np.concatenate((bandpowers, i*np.ones(dl+1)))
                        weights = np.concatenate((weights, np.ones(dl+1)))


      
	bpws = np.array(bandpowers, dtype='int32')
	ls = np.array(ell, dtype='int32')
	w = np.array(weights, dtype='float32')/(dl+1)


	#name_output = '/home/drc01/sobreira/andluiz/halogen/outputs/bin.dat'
        #dir_name = os.getcwd()
	name_output = dir_name + '/bin.dat'
	output = np.column_stack((bpws,ls,w))
	np.savetxt(name_output, output)

        return bpws, ls, w

def coupling_matrix(NSIDE, mask, Bin):


	# Initialize spin-0 field
	field = np.ones_like(mask)
	f0 = nmt.NmtField(mask, [field])

	# We then generate an NmtWorkspace object that we use to compute and store the mode coupling matrix. 
	#Note that this matrix depends ONLY on the masks of the two fields to correlate, but not on the maps themselves.
	w=nmt.NmtWorkspace()
	w.compute_coupling_matrix(f0, f0, Bin)

	return w;

	#The function defined below will compute the power spectrum between two NmtFields f_a and f_b, 
	#using the coupling matrix stored in the NmtWorkspace wsp.
	#Note that the most expensive operations in the MASTER algorithm are the computation of the coupling matrix and the deprojection bias. 
	#Since these two objects are precomputed, this function should be pretty fast!

def compute_namaster(mask, field, wsp, Bin):

	#Compute the power spectrum (a la anafast) of the masked fields
	#You can use n_iter=0 here to speed up the computation, but the default value of 3 is recommended in general.
	f0 = nmt.NmtField(mask, [field])
	cl_coupled = nmt.compute_coupled_cell(f0, f0)
	
	#Decouple power spectrum into bandpowers inverting the coupling matrix
	cl_decoupled = wsp.decouple_cell(cl_coupled, cl_bias=None, cl_noise=None)

	ell = Bin.get_effective_ells()
	cl = cl_decoupled[0]

	output = np.zeros(shape=(len(ell),2))
	output[:,0] = ell
	output[:,1] = cl

	return output



def make_plots(dir_name, Bin, nshells, nmocks):


        ell = Bin.get_effective_ells()
	nell = len(ell)

        #dir_name = os.getcwd()


	cl_mean = np.zeros(shape=(nshells, nell))


	for i in range(nmocks):
	        output_name = dir_name + '/outputs/cldata/cl_mock' + str(i) + '.dat'
		cl_mean += np.loadtxt(output_name).T

	cl_mean = cl_mean/nmocks

	for i in range(nshells):

		print 'Plotting Cls shell ' + str(i+1)

	
                plt.pyplot.figure(figsize = (6,5), dpi = 1024)
                plt.pyplot.semilogy(ell, cl_mean[i,:], 'r.', label = "NaMaster data")
                plt.pyplot.title("Angular power spectrum - shell " + str(i+1))
                plt.pyplot.xlabel("$\ell$")
                plt.pyplot.ylabel("$C_\ell$")
                plt.pyplot.legend(loc = "upper right", shadow = True)
                plt.pyplot.grid(color = 'k', linestyle = '--', linewidth = 0.2)
                plt.pyplot.xlim((0.8*np.amin(ell), 1.2*np.amax(ell)))
		#plt.pyplot.ylim((0.8*np.amin(cl_mean[i,:]), 1.2*np.amax(cl_mean[i,:])))
                plt.pyplot.savefig(dir_name +  '/outputs/cl_' + str(i) + '.png')



	return

def make_covmat(dir_name, Bin, nshells, nmocks):


        ell = Bin.get_effective_ells()
        nell = len(ell)

	#Computing average Cl's at photo-z bin on all the mocks
        cl_mean = np.zeros(shape=(nshells, nell))

        #dir_name = os.getcwd()

        for k in range(nmocks):
                input_name = dir_name + '/outputs/cldata/cl_mock' + str(k) + '.dat'
                cl_mean += np.loadtxt(input_name).T

        cl_mean = cl_mean/nmocks
	#cl_mean[i, :] is the average Cl's at photo-z i-bin


	#Making the full covariance matrix with all the photo-z bins

	print 'Computing and plotting...'

	for shell_i in range(nshells):

		for shell_j in range(shell_i, nshells):

			print 'Covariance ' + str(shell_i) + '-' + str(shell_j) 

			#Covariance matriz between photo-z shell_i and shell_j
			cov_partial = np.zeros(shape=(nell,nell))

			#Reading and computing covariance on each mock
			for k in range(nmocks):				

				input_name = dir_name + '/outputs/cldata/cl_mock' + str(k) + '.dat'
		                cl_mock = np.loadtxt(input_name).T

				for ell_i in range(nell):
					for ell_j in range(nell):
						cov_partial[ell_i][ell_j] += (cl_mock[shell_i][ell_i] - cl_mean[shell_i][ell_i])*(cl_mock[shell_j][ell_j] - cl_mean[shell_j][ell_j])


			cov_partial = cov_partial/(nmocks-1)
	                output_name = dir_name +  '/outputs/covmatrix/cov_' + str(shell_i) + '_' + str(shell_j) + '.dat'
			np.savetxt(output_name, cov_partial)

			#Plotting matrix

			v_min, v_max = np.amin(cov_partial), np.amax(cov_partial)			
			lin_thresh = v_max*5E-3

			fig, ax = plt.pyplot.subplots(1, 1)
			pcm = ax.pcolormesh(ell, ell, cov_partial, norm=colors.SymLogNorm(linthresh=lin_thresh, linscale=1.0, vmin=v_min, vmax=v_max),  cmap='YlGnBu')
			#pcm = ax.pcolormesh(ell, ell, cov_partial, norm=colors.LogNorm(),  cmap='YlGnBu')
			fig.colorbar(pcm, ax=ax, extend='max', orientation='vertical')
			plt.pyplot.title("Covariance matrix - bins (" + str(shell_i+1) + ', ' + str(shell_j+1) + ')')
	                plt.pyplot.xlabel("$\ell$")
        	        plt.pyplot.ylabel("$\ell$")
			output_name = dir_name +  '/outputs/covmatrix/cov_' + str(shell_i) + '_' + str(shell_j) + '.png'
			fig.savefig(output_name)

	return 

