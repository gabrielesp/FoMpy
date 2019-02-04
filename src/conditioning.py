	

# # class interpolator(metaclass=ABCMeta):
# # 	''
# # 	def __init__(self, **kwargs):
# # 		self.dataset = []
# # 		# self.n_sims = None

# # class polynomial_interpolation(interpolator):
# # 	''
# # 	def __init__(self,**kwargs):
# # 		interpolator.__init__(self, **kwargs)

# # class spline_interpolation(interpolator):
# # 	''
# # 	def __init__(self,**kwargs):
# # 		interpolator.__init__(self, **kwargs)


# 	def normalize(self,data,norm):

# 	       norm_fl= float(norm)
# 	       normalized_data = data/norm_fl
# 	       return normalized_data

# 	def smooth_data(self,data,arg1, arg2):

# 	       print("under construction")
# 	       return smoothed_data

# def interpol(x,y,n,s):
# 	from scipy.interpolate import CubicSpline

# # ‘natural’: The second derivative at curve ends are zero.
# #  Assuming a 1D y, bc_type=((2, 0.0), (2, 0.0)) is the same condition.
# 	#print(x)
# 	#print(y)
# 	cs_natural = CubicSpline(x, y, bc_type='natural')
# 	# print(cs_natural.derivative())
# 	# # cs_natural_2 = CubicSpline(x, y, bc_type=((1, 0.0), (1, 0.0)),bc_type=((2, 0.0), (2, 0.0)))

# 	xvals = np.linspace(x[0], x[-1], n)
# 	# fig,ax1 = plt.subplots()
# 	# ax2 = ax1.twinx()
# 	# ax2.scatter(x,y, color='r',s=10)
# 	# ax2.plot(xvals, cs_natural(xvals), label="natural")

# 	# ax2.set_ylim(np.min(cs_natural(xvals)),np.max(cs_natural(xvals)))

# 	# ax1.set_yscale('log')	

# 	# ax1.scatter(x,y, color='r',s=10)
# 	# ax1.plot(xvals, cs_natural(xvals), label="natural_log")
# 	# ax1.set_ylim(np.min(cs_natural(xvals)),np.max(cs_natural(xvals)))

# 	# ax1.legend()
# 	# # ax2.legend()
# 	# plt.show()

# # #splrep creates a k degree of the spline fit
# # 	if s == 0:
# # 		spline = interpolate.splrep(x, y, k=5)
# # 	else :
# # 		spline = interpolate.splrep(x, y, s)

# # 	#splev evaluates the spline fit into the array stored in xvals (n-long)
# # 	xvals = np.linspace(x[0], x[-1], n)
# # 	yinterp = interpolate.splev(xvals, spline, der=0)

# 	return (xvals, cs_natural(xvals))