

def find_closest(A, target):
    #A must be sorted
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A)-1)
    left = A[idx-1]
    right = A[idx]
    idx -= target - left < right - target
    return idx

def get_diff(arr, order = 1, diff = 'central'):
    x = arr[:,0]
    y = arr[:,1]

    dx = np.gradient(x)
    dy = np.gradient(y)
    d1 = dy/dx

    d2y = np.gradient(d1)
    d2 = d2y/dx

    d3y = np.gradient(d2y)
    d3 = d3y/dx

    if order == 1:
        gm_c = np.column_stack((x[2:], d1[2:]))

        if diff == 'central':
            return(gm_c)

    if order == 2:
        gm2_c = np.column_stack((x[2:], d2[2:]))

        if diff == 'central':
            return(gm2_c)

    if order == 3:
        gm3_c = np.column_stack((x[2:], d3[2:]))
        if diff == 'central':
            return(gm3_c)

def checkPath(cwd):
	try:
		os.makedirs(cwd)
	except OSError:
		if not os.path.isdir(cwd):
			raise
