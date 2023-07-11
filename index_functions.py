import bagpipes as pipes
import numpy as np
import matplotlib.pyplot as plt

def d4000(spec,obj,plot):

	#Convert flux to frequency units as per Bruzual 83
	spec[:,1] *=  1e10 * spec[:,0]**2/3e+8

	if spec.shape[1] == 3:
		spec[:,2] *=  1e10 * spec[:,0]**2/3e+8 #Convert errors to frequency units also

	mask1 = (spec[:,0] > 3850) & (spec[:,0] < 3950)
	mask2 = (spec[:,0] > 4000) & (spec[:,0] < 4100)

	left = spec[mask1, :]
	right = spec[mask2, :]

	m_left = np.mean(left[:,1])
	m_right = np.mean(right[:,1])

	ratio = m_right/m_left

	#Plot ----------------------------------------------------------------------------
	if plot == True:
		#plot
		plt.title(obj)
		s = spec[spec[:,0]>3800]
		s = s[s[:,0]<4150]
		plt.plot(s[:,0],s[:,1])

		plt.axvspan(3850, 3950, color='limegreen', alpha=0.2)
		plt.axvspan(4000, 4100, color='limegreen', alpha=0.2)

		#plt.savefig(path + obj + '_D4000.png', format = 'png')
		plt.show()
	#-------------------------------------------------------------------------------------

	#Calculate errors if they exist
	if spec.shape[1] == 3:

		#Continuum (add in quadrature then scale like the mean)
		left_err = np.sqrt(np.sum(left[:,2]**2))/float(left.shape[0])
		right_err = np.sqrt(np.sum(left[:,2]**2))/float(left.shape[0])

		#ratio (add fractional errors in quadrature then scale by the absolute value)
		ratio_err = np.abs(ratio) * np.sqrt( (left_err/m_left)**2 + (right_err/m_right)**2 )
	
		return (round(ratio,5), round(ratio_err,5))

	else:
		return round(ratio,5)

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def hdelta(spec,obj,plot):

	#Continuum
	mask1 = (spec[:,0] > 4041.60) & (spec[:,0] < 4079.75)
	mask2 = (spec[:,0] > 4128.50) & (spec[:,0] < 4161.00)

	left = spec[mask1, :]
	right = spec[mask2, :]

	m_left = np.mean(left[:,1])
	m_right = np.mean(right[:,1])

	cont_fluxes = np.array([m_left, m_right]) #Put mean continuum fluxes into array

	#Feature
	maskf = (spec[:,0] > 4083.50) & (spec[:,0] < 4122.25)
	feature = spec[maskf, :]
	feature_flux = np.mean(feature[:,1])

	#EW
	feature_width = 4122.25 - 4083.50 #feature[-1,0] - feature[0,0]
	cont_flux = np.mean(cont_fluxes) #making continuum be stright line
	area = feature_width*(cont_flux - feature_flux)
	ew = area/cont_flux

	#Plot ----------------------------------------------------------------------------
	if plot == True:
		#plot
		plt.title(obj)
		s = spec[spec[:,0]>4000]
		s = s[s[:,0]<4200]
		plt.plot(s[:,0],s[:,1])

		plt.axvspan(4041.6, 4079.75, color='limegreen', alpha=0.2)
		plt.axvspan(4128.5, 4161.0, color='limegreen', alpha=0.2)
		plt.axvspan(4083.50, 4122.25, color='red', alpha=0.2)

		#plt.savefig(path + obj + '_H-delta.png', format = 'png')
		plt.show()
	#-------------------------------------------------------------------------------------

	#Calculate errors if they exist
	if spec.shape[1] == 3:

		#Continuum error (add in quadrature, scale like the mean and then same again for final err)
		left_err = np.sqrt(np.sum(left[:,2]**2))/float(left.shape[0])
		right_err = np.sqrt(np.sum(left[:,2]**2))/float(right.shape[0])
		
		cont_err = 0.5 * np.sqrt(left_err**2 + right_err**2)
		
		#Feature flux error
		feature_err = np.sqrt(np.sum(feature[:,2]**2))/float(feature.shape[0])

		#Area error (add fractional errors in quadrature then scale by feature width)
		area_err = feature_width * np.sqrt( cont_err**2 + feature_err**2 )

		#EW error (add fractional errors in quadrature)
		ew_err = np.abs(ew) * np.sqrt( (area_err/area)**2 + (cont_err/cont_flux)**2 )
	
		return (ew, ew_err)

	else:
		return ew

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def halpha(spec,obj,plot):
	#Continuum
	mask1 = (spec[:,0] > 6510) & (spec[:,0] < 6530)
	mask2 = (spec[:,0] > 6600) & (spec[:,0] < 6640)

	left = spec[mask1, :]
	right = spec[mask2, :]

	m_left = np.mean(left[:,1])
	m_right = np.mean(right[:,1])

	cont_fluxes = np.array([m_left, m_right]) #Put mean continuum fluxes into array

	#Feature
	maskf = (spec[:,0] > 6531.99) & (spec[:,0] < 6594.0)
	feature = spec[maskf, :]
	feature_flux = np.mean(feature[:,1])

	#EW
	feature_width = 6594.0 - 6531.99 #feature[-1,0] - feature[0,0]
	cont_flux = np.mean(cont_fluxes) #making continuum be stright line
	area = feature_width*(cont_flux - feature_flux)
	ew = area/cont_flux

	#Plot ----------------------------------------------------------------------------
	if plot == True:
		#plot
		plt.title(obj)
		s = spec[spec[:,0]>6500]
		s = s[s[:,0]<6700]
		plt.plot(s[:,0],s[:,1])

		plt.axvspan(6510, 6530, color='limegreen', alpha=0.2)
		plt.axvspan(6600, 6640, color='limegreen', alpha=0.2)
		plt.axvspan(6531.99, 6594, color='red', alpha=0.2)

		#plt.savefig(path + obj + '_H-alpha.png', format = 'png')
		plt.show()
	#-------------------------------------------------------------------------------------

	#Calculate errors if they exist
	if spec.shape[1] == 3:

		#Continuum error (add in quadrature, scale like the mean and then same again for final err)
		left_err = np.sqrt(np.sum(left[:,2]**2))/float(left.shape[0])
		right_err = np.sqrt(np.sum(left[:,2]**2))/float(right.shape[0])
		
		cont_err = 0.5 * np.sqrt(left_err**2 + right_err**2)
		
		#Feature flux error
		feature_err = np.sqrt(np.sum(feature[:,2]**2))/float(feature.shape[0])

		#Area error (add fractional errors in quadrature then scale by feature width)
		area_err = feature_width * np.sqrt( cont_err**2 + feature_err**2 )

		#EW error (add fractional errors in quadrature)
		ew_err = np.abs(ew) * np.sqrt( (area_err/area)**2 + (cont_err/cont_flux)**2 )
	
		return (ew, ew_err)

	else:
		return ew
